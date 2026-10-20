!------------------------------------------------------------------------------
! Module: G_DUMP
! Purpose: GPU reductions for global heat-release and gas-mass output data.
!------------------------------------------------------------------------------
MODULE G_DUMP
    USE PRECISION_PARAMETERS, ONLY: EB,FB
    USE GLOBAL_CONSTANTS
    USE MESH_POINTERS
    USE OUTPUT_DATA
    USE ISOSMOKE, ONLY: SMOKE3D_TO_FILE
    USE DEVICE_VARIABLES, ONLY: DEVICE_TYPE,SUBDEVICE_TYPE,DEVICE,N_DEVC
    USE G_CONS, ONLY: GPU_DEBUG_H2D_D2H,GPUON,GPU_IO_ON
    USE G_MESH_POINTERS
    USE CUDAFOR
    USE G_ERRO, ONLY: GPU_WAIT
    IMPLICIT NONE

    CONTAINS

        ! The GPU path is deliberately limited to the three automatic
        ! reaction Smoke3D fields.  Custom SM3D definitions and cut-cell
        ! geometry continue through the unchanged CPU implementation.
        LOGICAL FUNCTION GPU_SMOKE3D_SUPPORTED(NM)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: NM

            GPU_SMOKE3D_SUPPORTED = .FALSE.
            IF (.NOT.GPUON .OR. .NOT.GPU_IO_ON) RETURN
            IF (.NOT.ALLOCATED(SMOKE3D_FILE)) RETURN
            IF (N_SMOKE3D/=3 .OR. CC_IBM) RETURN
            IF (.NOT.ALLOCATED(GPU_MESHES)) RETURN
            IF (NM<LBOUND(GPU_MESHES,1) .OR. NM>UBOUND(GPU_MESHES,1)) RETURN
            IF (.NOT.ALLOCATED(GPU_MESHES(NM)%DUMP)) RETURN
            IF (.NOT.ALLOCATED(GPU_MESHES(NM)%DUMP%SMOKE_PACK_D)) RETURN

            ! READ_SM3D creates these automatic fields in this order.
            IF (SMOKE3D_FILE(1)%QUANTITY_INDEX/=1) RETURN
            IF (SMOKE3D_FILE(1)%Y_INDEX<=0 .AND. SMOKE3D_FILE(1)%Z_INDEX<=0) RETURN
            IF (SMOKE3D_FILE(2)%QUANTITY_INDEX/=11) RETURN
            IF (SMOKE3D_FILE(3)%QUANTITY_INDEX/=64) RETURN
            GPU_SMOKE3D_SUPPORTED = .TRUE.
        END FUNCTION GPU_SMOKE3D_SUPPORTED

        SUBROUTINE DUMP_SMOKE3D_GPU(T,DT,NM)
            IMPLICIT NONE
            REAL(EB), INTENT(IN) :: T,DT
            INTEGER, INTENT(IN) :: NM
            INTEGER :: N,NVALS,MODE,Y_INDEX,Z_INDEX
            REAL(FB) :: STIME
            REAL(EB) :: FR_C
            TYPE(SMOKE3D_TYPE), POINTER :: S3

            IF (.NOT.GPU_SMOKE3D_SUPPORTED(NM)) RETURN
            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)
            NVALS = (IBAR+1)*(JBAR+1)*(KBAR+1)
            STIME = REAL(T_BEGIN + (T-T_BEGIN)*TIME_SHRINK_FACTOR,FB)
            FR_C = RTE_SOURCE_CORRECTION_FACTOR**0.25_EB

            DO N=1,N_SMOKE3D
                S3 => SMOKE3D_FILE(N)
                SELECT CASE(N)
                   CASE(1); MODE = 1
                   CASE(2); MODE = 2
                   CASE(3); MODE = 3
                   CASE DEFAULT; CYCLE
                END SELECT
                Y_INDEX = S3%Y_INDEX
                Z_INDEX = S3%Z_INDEX
                CALL SMOKE3D_KERNEL<<<DIM3((NVALS+GPU_DUMP_REDUCTION_BLOCK-1)/GPU_DUMP_REDUCTION_BLOCK,1,1), &
                                      GPU_DUMP_REDUCTION_BLOCK,0,GM%STRM>>>(GM%DUMP%SMOKE_PACK_D,GM%RHOD,GM%TMPD,GM%QD, &
                    GM%CHI_RD,GM%ZZD,SPECIES_Z2YD,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,N_TRACKED_SPECIES,N_SPECIES, &
                    MODE,Y_INDEX,Z_INDEX,TMPM,QR_CLIP,FR_C,RTE_SOURCE_CORRECTION)
                CALL GPU_WAIT(GM%STRM)
                GM%DUMP%SMOKE_PACK_H = GM%DUMP%SMOKE_PACK_D
                CALL SMOKE3D_TO_FILE(NM,STIME,N,GM%DUMP%SMOKE_PACK_H,IBP1,JBP1,KBP1)
            ENDDO
        END SUBROUTINE DUMP_SMOKE3D_GPU

        INTEGER FUNCTION GPU_DEVC_REDUCTION_KIND(DV)
            IMPLICIT NONE
            TYPE(DEVICE_TYPE), INTENT(IN) :: DV
            INTEGER :: IND,AXIS

            ! Map each DEVC quantity to its GPU reduction.
            GPU_DEVC_REDUCTION_KIND = 0
            IF (DV%N_QUANTITY/=1) RETURN
            IF (.NOT.ALLOCATED(DV%QUANTITY_INDEX)) RETURN

            IND = DV%QUANTITY_INDEX(1)

            SELECT CASE(IND)
               CASE(91:93)
                  IF (DV%SPATIAL_STATISTIC/='AREA INTEGRAL') RETURN
                  IF (DV%Y_INDEX>0 .OR. DV%Z_INDEX>0) RETURN
                  AXIS = ABS(DV%IOR_ASSUMED)
                  IF (AXIS==0) AXIS = ABS(DV%IOR)
                  IF (AXIS<1 .OR. AXIS>3) RETURN
                  GPU_DEVC_REDUCTION_KIND = 1 ! MASS FLUX X/Y/Z
               CASE(111:113)
                  IF (DV%SPATIAL_STATISTIC/='AREA INTEGRAL') RETURN
                  IF (DV%Y_INDEX>0 .OR. DV%Z_INDEX>0) RETURN
                  AXIS = ABS(DV%IOR_ASSUMED)
                  IF (AXIS==0) AXIS = ABS(DV%IOR)
                  IF (AXIS<1 .OR. AXIS>3) RETURN
                  IF (SIM_MODE/=DNS_MODE) GPU_DEVC_REDUCTION_KIND = 2 ! ENTHALPY FLUX X/Y/Z
            END SELECT
        END FUNCTION GPU_DEVC_REDUCTION_KIND

        LOGICAL FUNCTION GPU_DEVC_REDUCTION_SUPPORTED(DV)
            IMPLICIT NONE
            TYPE(DEVICE_TYPE), INTENT(IN) :: DV

            GPU_DEVC_REDUCTION_SUPPORTED = GPU_DEVC_REDUCTION_KIND(DV)>0
        END FUNCTION GPU_DEVC_REDUCTION_SUPPORTED

        SUBROUTINE UPDATE_DEVC_REDUCTIONS_GPU(NM)
            IMPLICIT NONE

            INTEGER, INTENT(IN) :: NM
            INTEGER :: N,IND,AREA_AXIS,FLUX_AXIS,N_ACTIVE,GRID_BLOCKS,REDUCTION_KIND
            TYPE(DEVICE_TYPE), POINTER :: DV
            TYPE(SUBDEVICE_TYPE), POINTER :: SDV

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)
            N_ACTIVE = 0

            DO N=1,N_DEVC
                DV => DEVICE(N)
                REDUCTION_KIND = GPU_DEVC_REDUCTION_KIND(DV)
                IF (REDUCTION_KIND==0) CYCLE
                IF (DV%SUBDEVICE_INDEX(NM)==0) CYCLE
                SDV => DV%SUBDEVICE(DV%SUBDEVICE_INDEX(NM))
                IF (SDV%MESH/=NM) CYCLE
                SDV%GPU_VALUE_READY = .FALSE.

                AREA_AXIS = ABS(DV%IOR_ASSUMED)
                IF (AREA_AXIS==0) AREA_AXIS = ABS(DV%IOR)
                IND = DV%QUANTITY_INDEX(1)
                SELECT CASE(IND)
                   CASE(91,111); FLUX_AXIS = 1
                   CASE(92,112); FLUX_AXIS = 2
                   CASE(93,113); FLUX_AXIS = 3
                END SELECT

                IF (SDV%I1<0 .OR. SDV%J1<0 .OR. SDV%K1<0) CYCLE
                IF (SDV%I2<SDV%I1 .OR. SDV%J2<SDV%J1 .OR. SDV%K2<SDV%K1) CYCLE
                GRID_BLOCKS = MAX(1,((SDV%I2-SDV%I1+1)*(SDV%J2-SDV%J1+1)*(SDV%K2-SDV%K1+1) &
                                    + GPU_DUMP_REDUCTION_BLOCK-1)/GPU_DUMP_REDUCTION_BLOCK)

                SELECT CASE(REDUCTION_KIND)
                   CASE(1,2)
                      CALL DEVC_FLUX_REDUCTION_KERNEL<<<DIM3(GRID_BLOCKS,1,1),GPU_DUMP_REDUCTION_BLOCK,0,GM%STRM>>>( &
                          GM%DUMP%DEVC_PARTIAL_D,GM%UD,GM%VD,GM%WD,GM%RHOD,GM%TMPD,GM%MUD,GM%ZZD,H_SENS_ZD, &
                          GM%DXD,GM%DYD,GM%DZD,GM%RDXND,GM%RDYND,GM%RDZND,GM%RRND, &
                          SDV%I1,SDV%I2,SDV%J1,SDV%J2,SDV%K1,SDV%K2,AREA_AXIS,FLUX_AXIS,IND, &
                          IBP1,JBP1,KBP1,N_TRACKED_SPECIES,I_MAX_TEMP,SIM_MODE,DNS_MODE,TMPA,CPOPR, &
                          DV%QUANTITY_RANGE(1),DV%QUANTITY_RANGE(2))
                      CALL FINALIZE_DEVC_REDUCTION_KERNEL<<<DIM3(1,1,1),GPU_DUMP_REDUCTION_BLOCK,0,GM%STRM>>>( &
                          GM%DUMP%DEVC_PARTIAL_D,GM%DUMP%DEVC_RESULT_D,N,GRID_BLOCKS)
                      N_ACTIVE = N_ACTIVE + 1
                END SELECT
            ENDDO

            IF (N_ACTIVE==0) RETURN
            CALL GPU_WAIT(GM%STRM)
            GM%DUMP%DEVC_RESULT_H = GM%DUMP%DEVC_RESULT_D

            DO N=1,N_DEVC
                DV => DEVICE(N)
                IF (.NOT.GPU_DEVC_REDUCTION_SUPPORTED(DV)) CYCLE
                IF (DV%SUBDEVICE_INDEX(NM)==0) CYCLE
                SDV => DV%SUBDEVICE(DV%SUBDEVICE_INDEX(NM))
                IF (SDV%MESH/=NM) CYCLE
                IF (SDV%I1<0 .OR. SDV%J1<0 .OR. SDV%K1<0) CYCLE
                IF (SDV%I2<SDV%I1 .OR. SDV%J2<SDV%J1 .OR. SDV%K2<SDV%K1) CYCLE
                SDV%VALUE_1 = GM%DUMP%DEVC_RESULT_H(N)
                SDV%VALUE_2 = 0._EB
                SDV%VALUE_3 = 0._EB
                SDV%VALUE_4 = 0._EB
                SDV%GPU_VALUE_READY = .TRUE.
            ENDDO
        END SUBROUTINE UPDATE_DEVC_REDUCTIONS_GPU

        SUBROUTINE UPDATE_GLOBAL_OUTPUT_GPU(T,DT,NM)
            IMPLICIT NONE

            REAL(EB), INTENT(IN) :: T,DT
            INTEGER, INTENT(IN) :: NM
            INTEGER :: N_WALL_TOTAL,N_WALL_OUTPUTS,N_MASS_OUTPUTS
            TYPE(DIM3) :: CELL_GRID,WALL_GRID
            REAL(EB) :: ENTHALPY_SUM_OLD

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)

            IF (GPU_DEBUG_H2D_D2H) CALL SYNC_GLOBAL_OUTPUT_INPUTS_H2D(NM)

            N_WALL_TOTAL = N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS
            N_WALL_OUTPUTS = 3 + N_TRACKED_SPECIES
            N_MASS_OUTPUTS = 1 + N_SPECIES + N_TRACKED_SPECIES
            CELL_GRID = DIM3(GM%DUMP%CELL_BLOCKS,1,1)
            WALL_GRID = DIM3(GM%DUMP%WALL_BLOCKS,N_WALL_OUTPUTS,1)

            ENTHALPY_SUM_OLD = ENTHALPY_SUM(NM)

            CALL HRR_VOLUME_REDUCTION_KERNEL<<<CELL_GRID,GPU_DUMP_REDUCTION_BLOCK,0,GM%STRM>>>( &
                GM%DUMP%HRR_PARTIAL_D,GM%QD,GM%QRD,GM%D_PBAR_DTD,GM%D_PBAR_DT_SD,GM%PRESSURE_ZONED, &
                GM%RHOD,GM%TMPD,GM%ZZD,H_SENS_ZD,GM%DXD,GM%DYD,GM%DZD,GM%CELL_INDEXD,GM%CELL_SOLID, &
                GM%INTERPOLATED_MESHD,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,N_TRACKED_SPECIES,I_MAX_TEMP, &
                GM%DUMP%CELL_BLOCKS,MERGE(1,0,NM>1))
            CALL FINALIZE_OUTPUT_REDUCTION_KERNEL<<<DIM3(4,1,1),GPU_DUMP_REDUCTION_BLOCK,0,GM%STRM>>>( &
                GM%DUMP%HRR_PARTIAL_D,GM%DUMP%HRR_RESULT_D,GM%DUMP%CELL_BLOCKS,4)

            CALL HRR_WALL_REDUCTION_KERNEL<<<WALL_GRID,GPU_DUMP_REDUCTION_BLOCK,0,GM%STRM>>>( &
                GM%DUMP%WALL_PARTIAL_D,GM%WALL_TOPO_AOS,GM%UD,GM%VD,GM%WD,GM%ZZD,GM%WALL_TMP_F, &
                GM%WALL_RHO_F,GM%WALL_Q_CON_F,GM%WALL_ZZ_F,GM%WALL_RHO_D_F,GM%WALL_RDN, &
                GM%WALL_M_DOT_G_PP_ADJUST,GM%WALL_AREA,H_SENS_ZD,GM%INTERPOLATED_MESHD, &
                N_WALL_TOTAL,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,N_TRACKED_SPECIES,I_MAX_TEMP, &
                GM%DUMP%WALL_BLOCKS,N_WALL_OUTPUTS,SOLID_BOUNDARY,OPEN_BOUNDARY,MERGE(1,0,NM>1))
            CALL FINALIZE_OUTPUT_REDUCTION_KERNEL<<<DIM3(N_WALL_OUTPUTS,1,1),GPU_DUMP_REDUCTION_BLOCK,0,GM%STRM>>>( &
                GM%DUMP%WALL_PARTIAL_D,GM%DUMP%WALL_RESULT_D,GM%DUMP%WALL_BLOCKS,N_WALL_OUTPUTS)

            IF (MASS_FILE) THEN
                CELL_GRID = DIM3(GM%DUMP%CELL_BLOCKS,N_MASS_OUTPUTS,1)
                CALL MASS_VOLUME_REDUCTION_KERNEL<<<CELL_GRID,GPU_DUMP_REDUCTION_BLOCK,0,GM%STRM>>>( &
                    GM%DUMP%MASS_PARTIAL_D,GM%RHOD,GM%ZZD,SPECIES_Z2YD,GM%DXD,GM%DYD,GM%DZD, &
                    GM%CELL_INDEXD,GM%CELL_SOLID,GM%INTERPOLATED_MESHD,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1, &
                    N_SPECIES,N_TRACKED_SPECIES,GM%DUMP%CELL_BLOCKS,N_MASS_OUTPUTS,MERGE(1,0,NM>1))
                CALL FINALIZE_OUTPUT_REDUCTION_KERNEL<<<DIM3(N_MASS_OUTPUTS,1,1), &
                    GPU_DUMP_REDUCTION_BLOCK,0,GM%STRM>>>(GM%DUMP%MASS_PARTIAL_D,GM%DUMP%MASS_RESULT_D, &
                    GM%DUMP%CELL_BLOCKS,N_MASS_OUTPUTS)
            ENDIF

            CALL GPU_WAIT(GM%STRM)

            GM%DUMP%HRR_RESULT_H = GM%DUMP%HRR_RESULT_D
            GM%DUMP%WALL_RESULT_H = GM%DUMP%WALL_RESULT_D
            IF (MASS_FILE) GM%DUMP%MASS_RESULT_H = GM%DUMP%MASS_RESULT_D

            Q_DOT(1) = Q_DOT(1) + GM%DUMP%HRR_RESULT_H(1)
            Q_DOT(3) = Q_DOT(3) + GM%DUMP%HRR_RESULT_H(2)
            Q_DOT(7) = Q_DOT(7) + GM%DUMP%HRR_RESULT_H(3)
            ENTHALPY_SUM(NM) = GM%DUMP%HRR_RESULT_H(4)
            IF (ICYC>0) THEN
                Q_DOT(9) = Q_DOT(9) + (ENTHALPY_SUM(NM)-ENTHALPY_SUM_OLD)/DT
            ELSE
                Q_DOT(9) = 0._EB
            ENDIF

            Q_DOT(4:6) = Q_DOT(4:6) + GM%DUMP%WALL_RESULT_H(1:3)
            M_DOT(1:N_TRACKED_SPECIES) = M_DOT(1:N_TRACKED_SPECIES) + &
                                         GM%DUMP%WALL_RESULT_H(4:3+N_TRACKED_SPECIES)
            IF (MASS_FILE) MASS_DT(0:N_SPECIES+N_TRACKED_SPECIES) = &
                MASS_DT(0:N_SPECIES+N_TRACKED_SPECIES) + DT*GM%DUMP%MASS_RESULT_H(1:N_MASS_OUTPUTS)

            IF (NM==UPPER_MESH_INDEX) THEN
                Q_DOT_SUM = Q_DOT_SUM + DT*Q_DOT
                M_DOT_SUM = M_DOT_SUM + DT*M_DOT
            ENDIF

        END SUBROUTINE UPDATE_GLOBAL_OUTPUT_GPU


        SUBROUTINE SYNC_GLOBAL_OUTPUT_INPUTS_H2D(NM)
            IMPLICIT NONE

            INTEGER, INTENT(IN) :: NM
            INTEGER :: IW,N,N_WALL_TOTAL
            TYPE(WALL_TYPE), POINTER :: WC
            TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
            REAL(EB), ALLOCATABLE :: WALL_TMP_F_H(:),WALL_RHO_F_H(:),WALL_Q_CON_F_H(:)
            REAL(EB), ALLOCATABLE :: WALL_ZZ_F_H(:,:),WALL_RHO_D_F_H(:,:),WALL_M_DOT_ADJUST_H(:,:)

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)

            GM%UD = U
            GM%VD = V
            GM%WD = W
            GM%RHOD = RHO
            GM%TMPD = TMP
            GM%ZZD(:,:,:,1:N_TRACKED_SPECIES) = ZZ(:,:,:,1:N_TRACKED_SPECIES)
            GM%QD = Q
            GM%QRD = QR
            GM%D_PBAR_DTD = D_PBAR_DT
            GM%D_PBAR_DT_SD = D_PBAR_DT_S

            N_WALL_TOTAL = N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS
            IF (N_WALL_TOTAL<=0) RETURN

            ALLOCATE(WALL_TMP_F_H(GM%WALL_LB:GM%WALL_UB),WALL_RHO_F_H(GM%WALL_LB:GM%WALL_UB), &
                     WALL_Q_CON_F_H(GM%WALL_LB:GM%WALL_UB))
            ALLOCATE(WALL_ZZ_F_H(GM%WALL_LB:GM%WALL_UB,1:N_TRACKED_SPECIES), &
                     WALL_RHO_D_F_H(GM%WALL_LB:GM%WALL_UB,1:N_TRACKED_SPECIES), &
                     WALL_M_DOT_ADJUST_H(GM%WALL_LB:GM%WALL_UB,1:N_TRACKED_SPECIES))

            WALL_TMP_F_H = 0._EB
            WALL_RHO_F_H = 0._EB
            WALL_Q_CON_F_H = 0._EB
            WALL_ZZ_F_H = 0._EB
            WALL_RHO_D_F_H = 0._EB
            WALL_M_DOT_ADJUST_H = 0._EB
            IF (N_TRACKED_SPECIES>=1) WALL_ZZ_F_H(:,1) = 1._EB

            DO IW=1,N_WALL_TOTAL
                WC => WALL(IW)
                IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY .OR. WC%B1_INDEX<=0) CYCLE
                B1 => BOUNDARY_PROP1(WC%B1_INDEX)
                WALL_TMP_F_H(IW) = B1%TMP_F
                WALL_RHO_F_H(IW) = B1%RHO_F
                WALL_Q_CON_F_H(IW) = B1%Q_CON_F
                IF (ALLOCATED(B1%ZZ_F)) WALL_ZZ_F_H(IW,1:N_TRACKED_SPECIES) = B1%ZZ_F(1:N_TRACKED_SPECIES)
                IF (ALLOCATED(B1%RHO_D_F)) WALL_RHO_D_F_H(IW,1:N_TRACKED_SPECIES) = B1%RHO_D_F(1:N_TRACKED_SPECIES)
                IF (ALLOCATED(B1%M_DOT_G_PP_ADJUST)) WALL_M_DOT_ADJUST_H(IW,1:N_TRACKED_SPECIES) = &
                    B1%M_DOT_G_PP_ADJUST(1:N_TRACKED_SPECIES)
            ENDDO

            GM%WALL_TMP_F = WALL_TMP_F_H
            GM%WALL_RHO_F = WALL_RHO_F_H
            GM%WALL_Q_CON_F = WALL_Q_CON_F_H
            DO N=1,N_TRACKED_SPECIES
                GM%WALL_ZZ_F(:,N) = WALL_ZZ_F_H(:,N)
                GM%WALL_RHO_D_F(:,N) = WALL_RHO_D_F_H(:,N)
                GM%WALL_M_DOT_G_PP_ADJUST(:,N) = WALL_M_DOT_ADJUST_H(:,N)
            ENDDO

            DEALLOCATE(WALL_TMP_F_H,WALL_RHO_F_H,WALL_Q_CON_F_H,WALL_ZZ_F_H,WALL_RHO_D_F_H, &
                       WALL_M_DOT_ADJUST_H)

        END SUBROUTINE SYNC_GLOBAL_OUTPUT_INPUTS_H2D


        ATTRIBUTES(GLOBAL) SUBROUTINE HRR_VOLUME_REDUCTION_KERNEL(PARTIAL,Q,QR,D_PBAR_DT,D_PBAR_DT_S, &
            PRESSURE_ZONE,RHO,TMP,ZZ,H_SENS_Z,DX,DY,DZ,CELL_INDEX,CELL_SOLID,INTERPOLATED_MESH, &
            IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,N_TRACKED_SPECIES,I_MAX_TEMP,N_BLOCKS,SKIP_INTERPOLATED)
            IMPLICIT NONE

            INTEGER, VALUE :: IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,N_TRACKED_SPECIES,I_MAX_TEMP,N_BLOCKS
            INTEGER, VALUE :: SKIP_INTERPOLATED
            REAL(EB), DEVICE :: PARTIAL(N_BLOCKS,4)
            REAL(EB), DEVICE :: Q(0:IBP1,0:JBP1,0:*),QR(0:IBP1,0:JBP1,0:*)
            REAL(EB), DEVICE :: RHO(0:IBP1,0:JBP1,0:*),TMP(0:IBP1,0:JBP1,0:*)
            REAL(EB), DEVICE :: ZZ(0:IBP1,0:JBP1,0:KBP1,1:*)
            REAL(EB), DEVICE :: D_PBAR_DT(0:*),D_PBAR_DT_S(0:*),H_SENS_Z(0:I_MAX_TEMP,1:*)
            REAL(EB), DEVICE :: DX(0:*),DY(0:*),DZ(0:*)
            INTEGER, DEVICE :: PRESSURE_ZONE(0:IBP1,0:JBP1,0:*),CELL_INDEX(0:IBP1,0:JBP1,0:*)
            INTEGER, DEVICE :: INTERPOLATED_MESH(1:IBAR,1:JBAR,1:*)
            LOGICAL, DEVICE :: CELL_SOLID(0:*)
            REAL(EB), SHARED :: SQ(GPU_DUMP_REDUCTION_BLOCK),SQR(GPU_DUMP_REDUCTION_BLOCK)
            REAL(EB), SHARED :: SP(GPU_DUMP_REDUCTION_BLOCK),SH(GPU_DUMP_REDUCTION_BLOCK)
            INTEGER :: TID,GID,I,J,K,REM,OFFSET,N,IC,IZ,ITMP
            REAL(EB) :: VC,Q_LOCAL,QR_LOCAL,P_LOCAL,H_LOCAL,H_S,TMP_LOCAL,FRAC

            TID = THREADIDX%X
            GID = (BLOCKIDX%X-1)*BLOCKDIM%X + TID
            Q_LOCAL = 0._EB
            QR_LOCAL = 0._EB
            P_LOCAL = 0._EB
            H_LOCAL = 0._EB

            IF (GID<=IBAR*JBAR*KBAR) THEN
                K = (GID-1)/(IBAR*JBAR) + 1
                REM = MOD(GID-1,IBAR*JBAR)
                J = REM/IBAR + 1
                I = MOD(REM,IBAR) + 1
                IC = CELL_INDEX(I,J,K)
                IF (.NOT.CELL_SOLID(IC) .AND. &
                    (SKIP_INTERPOLATED==0 .OR. INTERPOLATED_MESH(I,J,K)==0)) THEN
                    VC = DX(I)*DY(J)*DZ(K)
                    Q_LOCAL = Q(I,J,K)*VC
                    QR_LOCAL = QR(I,J,K)*VC
                    IZ = PRESSURE_ZONE(I,J,K)
                    P_LOCAL = 0.5_EB*(D_PBAR_DT_S(IZ)+D_PBAR_DT(IZ))*VC
                    TMP_LOCAL = TMP(I,J,K)
                    ITMP = MAX(0,MIN(I_MAX_TEMP-1,INT(TMP_LOCAL)))
                    FRAC = TMP_LOCAL-REAL(ITMP,EB)
                    H_S = 0._EB
                    DO N=1,N_TRACKED_SPECIES
                        H_S = H_S + ZZ(I,J,K,N)*(H_SENS_Z(ITMP,N) + &
                              FRAC*(H_SENS_Z(ITMP+1,N)-H_SENS_Z(ITMP,N)))
                    ENDDO
                    H_LOCAL = RHO(I,J,K)*H_S*VC
                ENDIF
            ENDIF

            SQ(TID) = Q_LOCAL
            SQR(TID) = QR_LOCAL
            SP(TID) = P_LOCAL
            SH(TID) = H_LOCAL
            CALL SYNCTHREADS()
            OFFSET = BLOCKDIM%X/2
            DO WHILE (OFFSET>=1)
                IF (TID<=OFFSET) THEN
                    SQ(TID) = SQ(TID) + SQ(TID+OFFSET)
                    SQR(TID) = SQR(TID) + SQR(TID+OFFSET)
                    SP(TID) = SP(TID) + SP(TID+OFFSET)
                    SH(TID) = SH(TID) + SH(TID+OFFSET)
                ENDIF
                CALL SYNCTHREADS()
                OFFSET = OFFSET/2
            ENDDO
            IF (TID==1) THEN
                PARTIAL(BLOCKIDX%X,1) = SQ(1)
                PARTIAL(BLOCKIDX%X,2) = SQR(1)
                PARTIAL(BLOCKIDX%X,3) = SP(1)
                PARTIAL(BLOCKIDX%X,4) = SH(1)
            ENDIF

        END SUBROUTINE HRR_VOLUME_REDUCTION_KERNEL


        ATTRIBUTES(GLOBAL) SUBROUTINE HRR_WALL_REDUCTION_KERNEL(PARTIAL,WALL_TOPO,U,V,W,ZZ,WALL_TMP_F, &
            WALL_RHO_F,WALL_Q_CON_F,WALL_ZZ_F,WALL_RHO_D_F,WALL_RDN,WALL_M_DOT,WALL_AREA,H_SENS_Z, &
            INTERPOLATED_MESH,N_WALL_TOTAL,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,N_TRACKED_SPECIES,I_MAX_TEMP, &
            N_BLOCKS,N_OUTPUTS,SOLID_BOUNDARY,OPEN_BOUNDARY,SKIP_INTERPOLATED)
            IMPLICIT NONE

            INTEGER, VALUE :: N_WALL_TOTAL,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,N_TRACKED_SPECIES,I_MAX_TEMP
            INTEGER, VALUE :: N_BLOCKS,N_OUTPUTS,SOLID_BOUNDARY,OPEN_BOUNDARY,SKIP_INTERPOLATED
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            REAL(EB), DEVICE :: PARTIAL(N_BLOCKS,N_OUTPUTS)
            REAL(EB), DEVICE :: U(0:IBP1,0:JBP1,0:*),V(0:IBP1,0:JBP1,0:*),W(0:IBP1,0:JBP1,0:*)
            REAL(EB), DEVICE :: ZZ(0:IBP1,0:JBP1,0:KBP1,1:*)
            REAL(EB), DEVICE :: WALL_TMP_F(0:*),WALL_RHO_F(0:*),WALL_Q_CON_F(0:*),WALL_RDN(0:*),WALL_AREA(0:*)
            REAL(EB), DEVICE :: WALL_ZZ_F(0:N_WALL_TOTAL,1:*),WALL_RHO_D_F(0:N_WALL_TOTAL,1:*)
            REAL(EB), DEVICE :: WALL_M_DOT(0:N_WALL_TOTAL,1:*)
            REAL(EB), DEVICE :: H_SENS_Z(0:I_MAX_TEMP,1:*)
            INTEGER, DEVICE :: INTERPOLATED_MESH(1:IBAR,1:JBAR,1:*)
            REAL(EB), SHARED :: SUM_SHARED(GPU_DUMP_REDUCTION_BLOCK)
            INTEGER :: TID,GID,IW,IIG,JJG,KKG,IOR,BTYPE,OUTPUT_INDEX,N,ITMP,OFFSET
            REAL(EB) :: VALUE,AREA,TMP_LOCAL,FRAC,H_S,H_S_ALPHA,H_S_J_ALPHA,U_N

            TID = THREADIDX%X
            GID = (BLOCKIDX%X-1)*BLOCKDIM%X + TID
            OUTPUT_INDEX = BLOCKIDX%Y
            VALUE = 0._EB

            IF (GID<=N_WALL_TOTAL) THEN
                IW = GID
                BTYPE = WALL_TOPO(IW)%BOUNDARY_TYPE
                IIG = WALL_TOPO(IW)%IIG
                JJG = WALL_TOPO(IW)%JJG
                KKG = WALL_TOPO(IW)%KKG
                IF ((SKIP_INTERPOLATED==0 .OR. INTERPOLATED_MESH(IIG,JJG,KKG)==0) .AND. &
                    (BTYPE==SOLID_BOUNDARY .OR. BTYPE==OPEN_BOUNDARY)) THEN
                    AREA = WALL_AREA(IW)
                    SELECT CASE(OUTPUT_INDEX)
                        CASE(1)
                            IOR = WALL_TOPO(IW)%IOR
                            SELECT CASE(IOR)
                                CASE( 1); U_N = -U(IIG-1,JJG,KKG)
                                CASE(-1); U_N =  U(IIG,JJG,KKG)
                                CASE( 2); U_N = -V(IIG,JJG-1,KKG)
                                CASE(-2); U_N =  V(IIG,JJG,KKG)
                                CASE( 3); U_N = -W(IIG,JJG,KKG-1)
                                CASE(-3); U_N =  W(IIG,JJG,KKG)
                            END SELECT
                            TMP_LOCAL = WALL_TMP_F(IW)
                            ITMP = MAX(0,MIN(I_MAX_TEMP-1,INT(TMP_LOCAL)))
                            FRAC = TMP_LOCAL-REAL(ITMP,EB)
                            H_S = 0._EB
                            DO N=1,N_TRACKED_SPECIES
                                H_S = H_S + WALL_ZZ_F(IW,N)*(H_SENS_Z(ITMP,N) + &
                                      FRAC*(H_SENS_Z(ITMP+1,N)-H_SENS_Z(ITMP,N)))
                            ENDDO
                            VALUE = -U_N*WALL_RHO_F(IW)*H_S*AREA
                        CASE(2)
                            VALUE = -WALL_Q_CON_F(IW)*AREA
                        CASE(3)
                            IF (N_TRACKED_SPECIES>1) THEN
                                TMP_LOCAL = WALL_TMP_F(IW)
                                ITMP = MAX(0,MIN(I_MAX_TEMP-1,INT(TMP_LOCAL)))
                                FRAC = TMP_LOCAL-REAL(ITMP,EB)
                                H_S_J_ALPHA = 0._EB
                                DO N=1,N_TRACKED_SPECIES
                                    H_S_ALPHA = H_SENS_Z(ITMP,N) + &
                                        FRAC*(H_SENS_Z(ITMP+1,N)-H_SENS_Z(ITMP,N))
                                    H_S_J_ALPHA = H_S_J_ALPHA + 2._EB*H_S_ALPHA*WALL_RHO_D_F(IW,N)* &
                                        (ZZ(IIG,JJG,KKG,N)-WALL_ZZ_F(IW,N))*WALL_RDN(IW)
                                ENDDO
                                VALUE = -H_S_J_ALPHA*AREA
                            ENDIF
                        CASE DEFAULT
                            IF (BTYPE==SOLID_BOUNDARY) VALUE = WALL_M_DOT(IW,OUTPUT_INDEX-3)*AREA
                    END SELECT
                ENDIF
            ENDIF

            SUM_SHARED(TID) = VALUE
            CALL SYNCTHREADS()
            OFFSET = BLOCKDIM%X/2
            DO WHILE (OFFSET>=1)
                IF (TID<=OFFSET) SUM_SHARED(TID) = SUM_SHARED(TID) + SUM_SHARED(TID+OFFSET)
                CALL SYNCTHREADS()
                OFFSET = OFFSET/2
            ENDDO
            IF (TID==1) PARTIAL(BLOCKIDX%X,OUTPUT_INDEX) = SUM_SHARED(1)

        END SUBROUTINE HRR_WALL_REDUCTION_KERNEL


        ATTRIBUTES(GLOBAL) SUBROUTINE MASS_VOLUME_REDUCTION_KERNEL(PARTIAL,RHO,ZZ,Z2Y,DX,DY,DZ,CELL_INDEX, &
            CELL_SOLID,INTERPOLATED_MESH,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,N_SPECIES,N_TRACKED_SPECIES, &
            N_BLOCKS,N_OUTPUTS,SKIP_INTERPOLATED)
            IMPLICIT NONE

            INTEGER, VALUE :: IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,N_SPECIES,N_TRACKED_SPECIES,N_BLOCKS,N_OUTPUTS
            INTEGER, VALUE :: SKIP_INTERPOLATED
            REAL(EB), DEVICE :: PARTIAL(N_BLOCKS,N_OUTPUTS),RHO(0:IBP1,0:JBP1,0:*)
            REAL(EB), DEVICE :: ZZ(0:IBP1,0:JBP1,0:KBP1,1:*),Z2Y(1:N_SPECIES,1:*)
            REAL(EB), DEVICE :: DX(0:*),DY(0:*),DZ(0:*)
            INTEGER, DEVICE :: CELL_INDEX(0:IBP1,0:JBP1,0:*),INTERPOLATED_MESH(1:IBAR,1:JBAR,1:*)
            LOGICAL, DEVICE :: CELL_SOLID(0:*)
            REAL(EB), SHARED :: SUM_SHARED(GPU_DUMP_REDUCTION_BLOCK)
            INTEGER :: TID,GID,I,J,K,REM,IC,OUTPUT_INDEX,N,OFFSET
            REAL(EB) :: VALUE,VC,Y_VALUE

            TID = THREADIDX%X
            GID = (BLOCKIDX%X-1)*BLOCKDIM%X + TID
            OUTPUT_INDEX = BLOCKIDX%Y
            VALUE = 0._EB

            IF (GID<=IBAR*JBAR*KBAR) THEN
                K = (GID-1)/(IBAR*JBAR) + 1
                REM = MOD(GID-1,IBAR*JBAR)
                J = REM/IBAR + 1
                I = MOD(REM,IBAR) + 1
                IC = CELL_INDEX(I,J,K)
                IF (.NOT.CELL_SOLID(IC) .AND. &
                    (SKIP_INTERPOLATED==0 .OR. INTERPOLATED_MESH(I,J,K)==0)) THEN
                    VC = DX(I)*DY(J)*DZ(K)
                    IF (OUTPUT_INDEX==1) THEN
                        VALUE = RHO(I,J,K)*VC
                    ELSEIF (OUTPUT_INDEX<=1+N_SPECIES) THEN
                        Y_VALUE = 0._EB
                        DO N=1,N_TRACKED_SPECIES
                            Y_VALUE = Y_VALUE + Z2Y(OUTPUT_INDEX-1,N)*ZZ(I,J,K,N)
                        ENDDO
                        Y_VALUE = MIN(1._EB,MAX(0._EB,Y_VALUE))
                        VALUE = RHO(I,J,K)*Y_VALUE*VC
                    ELSE
                        N = OUTPUT_INDEX-1-N_SPECIES
                        VALUE = RHO(I,J,K)*ZZ(I,J,K,N)*VC
                    ENDIF
                ENDIF
            ENDIF

            SUM_SHARED(TID) = VALUE
            CALL SYNCTHREADS()
            OFFSET = BLOCKDIM%X/2
            DO WHILE (OFFSET>=1)
                IF (TID<=OFFSET) SUM_SHARED(TID) = SUM_SHARED(TID) + SUM_SHARED(TID+OFFSET)
                CALL SYNCTHREADS()
                OFFSET = OFFSET/2
            ENDDO
            IF (TID==1) PARTIAL(BLOCKIDX%X,OUTPUT_INDEX) = SUM_SHARED(1)

        END SUBROUTINE MASS_VOLUME_REDUCTION_KERNEL


        ATTRIBUTES(GLOBAL) SUBROUTINE FINALIZE_OUTPUT_REDUCTION_KERNEL(PARTIAL,RESULT,N_BLOCKS,N_OUTPUTS)
            IMPLICIT NONE

            INTEGER, VALUE :: N_BLOCKS,N_OUTPUTS
            REAL(EB), DEVICE :: PARTIAL(N_BLOCKS,N_OUTPUTS),RESULT(N_OUTPUTS)
            REAL(EB), SHARED :: SUM_SHARED(GPU_DUMP_REDUCTION_BLOCK)
            INTEGER :: TID,OUTPUT_INDEX,INDEX,OFFSET
            REAL(EB) :: VALUE

            TID = THREADIDX%X
            OUTPUT_INDEX = BLOCKIDX%X
            VALUE = 0._EB
            INDEX = TID
            DO WHILE (INDEX<=N_BLOCKS)
                VALUE = VALUE + PARTIAL(INDEX,OUTPUT_INDEX)
                INDEX = INDEX + BLOCKDIM%X
            ENDDO
            SUM_SHARED(TID) = VALUE
            CALL SYNCTHREADS()
            OFFSET = BLOCKDIM%X/2
            DO WHILE (OFFSET>=1)
                IF (TID<=OFFSET) SUM_SHARED(TID) = SUM_SHARED(TID) + SUM_SHARED(TID+OFFSET)
                CALL SYNCTHREADS()
                OFFSET = OFFSET/2
            ENDDO
            IF (TID==1) RESULT(OUTPUT_INDEX) = SUM_SHARED(1)

        END SUBROUTINE FINALIZE_OUTPUT_REDUCTION_KERNEL


        ATTRIBUTES(GLOBAL) SUBROUTINE DEVC_FLUX_REDUCTION_KERNEL(PARTIAL,U,V,W,RHO,TMP,MU,ZZ,H_SENS_Z, &
            DX,DY,DZ,RDXN,RDYN,RDZN,RRN,I1,I2,J1,J2,K1,K2,AREA_AXIS,FLUX_AXIS,IND, &
            IBP1,JBP1,KBP1,N_TRACKED_SPECIES,I_MAX_TEMP,SIM_MODE_LOCAL,DNS_MODE_LOCAL,TMPA_IN,CPOPR_IN, &
            QUANTITY_MIN,QUANTITY_MAX)
            IMPLICIT NONE

            INTEGER, VALUE :: I1,I2,J1,J2,K1,K2,AREA_AXIS,FLUX_AXIS,IND
            INTEGER, VALUE :: IBP1,JBP1,KBP1,N_TRACKED_SPECIES,I_MAX_TEMP,SIM_MODE_LOCAL,DNS_MODE_LOCAL
            REAL(EB), VALUE :: TMPA_IN,CPOPR_IN,QUANTITY_MIN,QUANTITY_MAX
            REAL(EB), DEVICE :: PARTIAL(*)
            REAL(EB), DEVICE :: U(0:IBP1,0:JBP1,0:KBP1),V(0:IBP1,0:JBP1,0:KBP1),W(0:IBP1,0:JBP1,0:KBP1)
            REAL(EB), DEVICE :: RHO(0:IBP1,0:JBP1,0:KBP1),TMP(0:IBP1,0:JBP1,0:KBP1),MU(0:IBP1,0:JBP1,0:KBP1)
            REAL(EB), DEVICE :: ZZ(0:IBP1,0:JBP1,0:KBP1,1:*),H_SENS_Z(0:I_MAX_TEMP,1:*)
            REAL(EB), DEVICE :: DX(0:*),DY(0:*),DZ(0:*),RDXN(0:*),RDYN(0:*),RDZN(0:*),RRN(0:*)
            REAL(EB), SHARED :: SUM_SHARED(GPU_DUMP_REDUCTION_BLOCK)
            INTEGER :: TID,GID,NI,NJ,NK,NITEMS,REM,I,J,K,IP,JP,KP,OFFSET,N,ITMP,ITMPA
            REAL(EB) :: VALUE,AREA,VEL,R_DN,TMP_TC,FRAC,FRAC_A,H_G,H_G_SUM,H0,H1,ZAVG

            TID = THREADIDX%X
            GID = (BLOCKIDX%X-1)*BLOCKDIM%X + TID
            NI = I2-I1+1
            NJ = J2-J1+1
            NK = K2-K1+1
            NITEMS = NI*NJ*NK
            VALUE = 0._EB

            IF (GID<=NITEMS) THEN
                K = K1 + (GID-1)/(NI*NJ)
                REM = MOD(GID-1,NI*NJ)
                J = J1 + REM/NI
                I = I1 + MOD(REM,NI)

                IP = I
                JP = J
                KP = K
                SELECT CASE(FLUX_AXIS)
                   CASE(1)
                      IP = I+1
                      VEL = U(I,J,K)
                      R_DN = RDXN(I)
                   CASE(2)
                      JP = J+1
                      VEL = V(I,J,K)
                      R_DN = RDYN(J)
                   CASE(3)
                      KP = K+1
                      VEL = W(I,J,K)
                      R_DN = RDZN(K)
                END SELECT

                IF (IND<=93) THEN
                   VALUE = 0.5_EB*(RHO(I,J,K)+RHO(IP,JP,KP))*VEL
                ELSE
                   TMP_TC = 0.5_EB*(TMP(I,J,K)+TMP(IP,JP,KP))
                   ITMP = MAX(0,MIN(I_MAX_TEMP-1,INT(TMP_TC)))
                   FRAC = TMP_TC-REAL(ITMP,EB)
                   ITMPA = MAX(0,MIN(I_MAX_TEMP-1,INT(TMPA_IN)))
                   FRAC_A = TMPA_IN-REAL(ITMPA,EB)
                   H_G_SUM = 0._EB
                   H_G = 0._EB
                   DO N=1,N_TRACKED_SPECIES
                      ZAVG = 0.5_EB*(ZZ(I,J,K,N)+ZZ(IP,JP,KP,N))
                      H0 = H_SENS_Z(ITMP,N)
                      H1 = H_SENS_Z(ITMP+1,N)
                      H_G_SUM = H_G_SUM + ZAVG*(H0+FRAC*(H1-H0))
                      H0 = H_SENS_Z(ITMPA,N)
                      H1 = H_SENS_Z(ITMPA+1,N)
                      H_G = H_G + ZAVG*(H0+FRAC_A*(H1-H0))
                   ENDDO
                   VALUE = VEL*0.5_EB*(RHO(I,J,K)+RHO(IP,JP,KP))*(H_G_SUM-H_G)
                   VALUE = VALUE - MU(I,J,K)*CPOPR_IN*(TMP(IP,JP,KP)-TMP(I,J,K))*R_DN
                   VALUE = VALUE*0.001_EB
                ENDIF

                SELECT CASE(AREA_AXIS)
                   CASE(1)
                      AREA = DY(J)*DZ(K)
                      IF (ABS(RRN(I))>0._EB) AREA = AREA/RRN(I)
                   CASE(2)
                      AREA = DX(I)*DZ(K)
                   CASE(3)
                      AREA = DX(I)*DY(J)
                      IF (ABS(RRN(I))>0._EB) AREA = AREA/RRN(I)
                   CASE DEFAULT
                      AREA = 0._EB
                END SELECT

                IF (VALUE<QUANTITY_MIN .OR. VALUE>QUANTITY_MAX) VALUE = 0._EB
                VALUE = VALUE*AREA
             ENDIF

            SUM_SHARED(TID) = VALUE
            CALL SYNCTHREADS()
            OFFSET = BLOCKDIM%X/2
            DO WHILE (OFFSET>=1)
                IF (TID<=OFFSET) SUM_SHARED(TID) = SUM_SHARED(TID) + SUM_SHARED(TID+OFFSET)
                CALL SYNCTHREADS()
                OFFSET = OFFSET/2
            ENDDO
            IF (TID==1) PARTIAL(BLOCKIDX%X) = SUM_SHARED(1)

        END SUBROUTINE DEVC_FLUX_REDUCTION_KERNEL


        ATTRIBUTES(GLOBAL) SUBROUTINE FINALIZE_DEVC_REDUCTION_KERNEL(PARTIAL,RESULT,RESULT_INDEX,N_BLOCKS)
            IMPLICIT NONE

            INTEGER, VALUE :: RESULT_INDEX,N_BLOCKS
            REAL(EB), DEVICE :: PARTIAL(*),RESULT(*)
            REAL(EB), SHARED :: SUM_SHARED(GPU_DUMP_REDUCTION_BLOCK)
            INTEGER :: TID,INDEX,OFFSET
            REAL(EB) :: VALUE

            TID = THREADIDX%X
            VALUE = 0._EB
            INDEX = TID
            DO WHILE (INDEX<=N_BLOCKS)
                VALUE = VALUE + PARTIAL(INDEX)
                INDEX = INDEX + BLOCKDIM%X
            ENDDO
            SUM_SHARED(TID) = VALUE
            CALL SYNCTHREADS()
            OFFSET = BLOCKDIM%X/2
            DO WHILE (OFFSET>=1)
                IF (TID<=OFFSET) SUM_SHARED(TID) = SUM_SHARED(TID) + SUM_SHARED(TID+OFFSET)
                CALL SYNCTHREADS()
                OFFSET = OFFSET/2
            ENDDO
            IF (TID==1) RESULT(RESULT_INDEX) = SUM_SHARED(1)

        END SUBROUTINE FINALIZE_DEVC_REDUCTION_KERNEL

        ATTRIBUTES(GLOBAL) SUBROUTINE SMOKE3D_KERNEL(PACK,RHO,TMP,Q,CHI_R,ZZ,Z2Y,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1, &
            N_TRACKED_SPECIES,N_SPECIES,MODE,Y_INDEX,Z_INDEX,TMPM_IN,QR_CLIP_IN,FR_C,RTE_CORRECTION)
            IMPLICIT NONE

            INTEGER, VALUE :: IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,N_TRACKED_SPECIES,N_SPECIES
            INTEGER, VALUE :: MODE,Y_INDEX,Z_INDEX
            LOGICAL, VALUE :: RTE_CORRECTION
            REAL(EB), VALUE :: TMPM_IN,QR_CLIP_IN,FR_C
            REAL(FB), DEVICE :: PACK(*)
            REAL(EB), DEVICE :: RHO(0:IBP1,0:JBP1,0:KBP1),TMP(0:IBP1,0:JBP1,0:KBP1)
            REAL(EB), DEVICE :: Q(0:IBP1,0:JBP1,0:KBP1),CHI_R(0:IBP1,0:JBP1,0:KBP1)
            REAL(EB), DEVICE :: ZZ(0:IBP1,0:JBP1,0:KBP1,1:*),Z2Y(1:N_SPECIES,1:*)
            INTEGER :: GID,I,J,K,C,CYCLE_N,II,JJ,KK,REM
            REAL(EB) :: VALUE,CELL_VALUE,Y_VALUE

            GID = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
            IF (GID>(IBP1*JBP1*KBP1)) RETURN

            REM = GID-1
            K = REM/(IBP1*JBP1)
            REM = MOD(REM,IBP1*JBP1)
            J = REM/IBP1
            I = MOD(REM,IBP1)
            VALUE = 0._EB

            ! Each output thread owns one node and averages its eight adjacent
            ! cell-center values, matching the CPU DUMP_SMOKE3D interpolation.
            DO C=0,7
                II = I + MOD(C,2)
                JJ = J + MOD(C/2,2)
                KK = K + C/4
                SELECT CASE(MODE)
                   CASE(1) ! soot density = rho * reconstructed soot mass fraction
                      Y_VALUE = 0._EB
                      IF (Y_INDEX>0) THEN
                         DO CYCLE_N=1,N_TRACKED_SPECIES
                            Y_VALUE = Y_VALUE + Z2Y(Y_INDEX,CYCLE_N)*ZZ(II,JJ,KK,CYCLE_N)
                         ENDDO
                      ELSEIF (Z_INDEX>0) THEN
                         Y_VALUE = ZZ(II,JJ,KK,Z_INDEX)
                      ENDIF
                      CELL_VALUE = RHO(II,JJ,KK)*Y_VALUE
                   CASE(2) ! HRRPUV is written in kW/m3
                      CELL_VALUE = Q(II,JJ,KK)*0.001_EB
                   CASE(3) ! effective flame temperature, in the CPU units
                      CELL_VALUE = TMP(II,JJ,KK)-TMPM_IN
                      IF (RTE_CORRECTION .AND. CHI_R(II,JJ,KK)*Q(II,JJ,KK)>QR_CLIP_IN) &
                          CELL_VALUE = (CELL_VALUE+TMPM_IN)*FR_C-TMPM_IN
                   CASE DEFAULT
                      CELL_VALUE = 0._EB
                END SELECT
                VALUE = VALUE + CELL_VALUE
            ENDDO
            PACK(GID) = REAL(0.125_EB*VALUE,FB)
        END SUBROUTINE SMOKE3D_KERNEL

END MODULE G_DUMP
