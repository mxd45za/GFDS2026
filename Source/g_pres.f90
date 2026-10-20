!------------------------------------------------------------------------------
! Module: G_PRES
! Purpose: GPU pressure assembly, baroclinic correction, and velocity-error evaluation.
!------------------------------------------------------------------------------
MODULE G_PRES
    USE G_ERRO, ONLY: GPU_WAIT
    USE G_MESH_POINTERS
    USE G_CONS
    USE G_POIS
    USE G_FUNC
    USE CUDAFOR
    CONTAINS

        !--------------------------------------------------------------------------
        ! Pressure-equation host orchestration
        !--------------------------------------------------------------------------
        SUBROUTINE BAROCLINIC_CORRECTION_GPU(T,NM)

            REAL(EB), INTENT(IN) :: T
            INTEGER, INTENT(IN) :: NM
            REAL(EB), POINTER, DIMENSION(:,:,:) :: RHOP,HP,P,RRHO
            REAL(EB), POINTER, DEVICE,DIMENSION(:,:,:) :: RHOPD,HPD,PD,RRHOD
            INTEGER  :: I,J,K
            REAL(EB) :: T_NOW


            INTEGER :: BLKD_1D , GRDD_1D
            TYPE(DIM3) :: BLKD , GRDD , BLKD_2D , GRDD_2D_X , GRDD_2D_Y , GRDD_2D_Z
            INTEGER , PARAMETER :: BS = 8

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)
            BLKD_1D = GPU_BLKD_1D
            GRDD_1D = GM%GRDD_EXTERNAL_WALL_1D

            BLKD_2D = DIM3(8,8,1)
            GRDD_2D_X = dim3((JBAR+2 +  BLKD_2D%X - 1) / BLKD_2D%X,(KBAR+2 +  BLKD_2D%Y - 1) / BLKD_2D%Y ,1)
            GRDD_2D_Y = dim3((IBAR+2 +  BLKD_2D%X - 1) / BLKD_2D%X,(KBAR+2 +  BLKD_2D%Y - 1) / BLKD_2D%Y ,1)
            GRDD_2D_Z = dim3((IBAR+2 +  BLKD_2D%X - 1) / BLKD_2D%X,(JBAR+2 +  BLKD_2D%Y - 1) / BLKD_2D%Y ,1)

            BLKD = GM%BLKD_3D_BP1
            GRDD = GM%GRDD_3D_BP1


            IF (SOLID_PHASE_ONLY .OR. FREEZE_VELOCITY) RETURN

            T_NOW = CURRENT_TIME()

            ! If the baroclinic torque term has been added to the momentum equation RHS, subtract it off.

            P    => WORK1 ! p=rho*(H-K)
            RRHO => WORK2 ! reciprocal of rho

            PD => GM%WORK1D
            RRHOD => GM%WORK2D

            IF (PREDICTOR) THEN
                RHOP=>RHO
                HP => H

                RHOPD => GM%RHOD
                HPD => GM%HD

                IF (GPU_DEBUG_H2D_D2H) THEN
                    GM%RHOD = RHO
                    GM%HD   = H
                ENDIF

            ELSE
                RHOP=>RHOS
                HP => HS

                RHOPD => GM%RHOSD
                HPD => GM%HSD

                IF (GPU_DEBUG_H2D_D2H) THEN
                    GM%RHOSD = RHOS
                    GM%HSD   = HS
                ENDIF

            ENDIF

            IF (GPU_DEBUG_H2D_D2H) THEN
                GM%FVXD  = FVX
                GM%FVYD  = FVY
                GM%FVZD  = FVZ
                GM%KRESD = KRES
            ENDIF

            IF (BAROCLINIC_TERMS_ATTACHED) THEN
                IF (GPU_DEBUG_H2D_D2H) THEN
                    GM%FVX_BD = FVX_B
                    GM%FVY_BD = FVY_B
                    GM%FVZ_BD = FVZ_B
                ENDIF
            ENDIF


            IF (BAROCLINIC_TERMS_ATTACHED) THEN
                CALL REMOVE_FB_TO_FA<<<GRDD,BLKD,0,GM%STRM>>>(IBP1,JBP1,KBP1,GM%FVXD,GM%FVYD,GM%FVZD,GM%FVX_BD,GM%FVY_BD,GM%FVZ_BD, PD , RHOPD , GM%KRESD , RRHOD , HPD)
            ENDIF
            CALL COMPUTE_PD_RRHOD<<<GRDD,BLKD,0,GM%STRM>>>(IBP1,JBP1,KBP1,GM%FVXD,GM%FVYD,GM%FVZD,GM%FVX_BD,GM%FVY_BD,GM%FVZ_BD, PD , RHOPD , GM%KRESD , RRHOD ,HPD )

            CALL ADD_FB_TO_F_XYZ<<<GRDD,BLKD,0,GM%STRM>>>(PD,RHOPD,RRHOD,GM%FVXD,GM%FVYD,GM%FVZD, &
                GM%FVX_BD,GM%FVY_BD,GM%FVZ_BD,GM%RDXND,GM%RDYND,GM%RDZND,IBAR,JBAR,KBAR,TWO_D)

            IF (GPU_DEBUG_H2D_D2H) THEN
                CALL GPU_WAIT(GM%STRM)
                FVX_B = GM%FVX_BD
                FVY_B = GM%FVY_BD
                FVZ_B = GM%FVZ_BD
                FVX   = GM%FVXD
                FVY   = GM%FVYD
                FVZ   = GM%FVZD
            ENDIF

            T_USED(4) = T_USED(4) + CURRENT_TIME() - T_NOW

            BAROCLINIC_TERMS_ATTACHED = .TRUE.

        END SUBROUTINE BAROCLINIC_CORRECTION_GPU

        SUBROUTINE NO_FLUX_GPU(DT,NM)
            REAL(EB), INTENT(IN) :: DT
            INTEGER, INTENT(IN) :: NM
            INTEGER :: BLKD_1D,GRDD_1D,ISTAT,IW,NOM
            INTEGER, PARAMETER :: BS_1D = 64
            REAL(EB) :: RFODT,T_NOW
            REAL(EB), ALLOCATABLE :: WALL_U_NORMAL_H_TMP(:), WALL_U_NORMAL_S_H_TMP(:)
            REAL(EB), DEVICE, POINTER, CONTIGUOUS :: OM_HPD(:,:,:)
            TYPE(OMESH_TYPE), POINTER :: OM
            TYPE(GPU_OMESH_TYPE), POINTER :: GOM_LOCAL

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)
            IF (SOLID_PHASE_ONLY .OR. FREEZE_VELOCITY) RETURN

            T_NOW = CURRENT_TIME()
            RFODT = RELAXATION_FACTOR/DT
            BLKD_1D = GPU_BLKD_1D
            GRDD_1D = GM%GRDD_TOTAL_WALL_1D

            IF (GPU_DEBUG_H2D_D2H) THEN
                GM%UD   = U
                GM%VD   = V
                GM%WD   = W
                GM%USD  = US
                GM%VSD  = VS
                GM%WSD  = WS
                GM%FVXD = FVX
                GM%FVYD = FVY
                GM%FVZD = FVZ
                IF (PREDICTOR) THEN
                    GM%HD = H
                ELSE
                    GM%HSD = HS
                ENDIF
                DO NOM=1,NMESHES
                    IF (.NOT.ALLOCATED(GM%OMESH)) EXIT
                    GOM_LOCAL => GM%OMESH(NOM)
                    IF (.NOT.GOM_LOCAL%ACTIVE) CYCLE
                    OM => MESHES(NM)%OMESH(NOM)
                    IF (PREDICTOR) THEN
                        GOM_LOCAL%H_D = OM%H
                    ELSE
                        GOM_LOCAL%HS_D = OM%HS
                    ENDIF
                ENDDO

                ALLOCATE(WALL_U_NORMAL_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)),STAT=ISTAT)
                ALLOCATE(WALL_U_NORMAL_S_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)),STAT=ISTAT)
                WALL_U_NORMAL_H_TMP   = 0._EB
                WALL_U_NORMAL_S_H_TMP = 0._EB
                DO IW=LBOUND(WALL,1),UBOUND(WALL,1)
                    IF (WALL(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY .OR. WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY) CYCLE
                    IF (IW>N_EXTERNAL_WALL_CELLS .AND. WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE
                    IF (IW<=N_EXTERNAL_WALL_CELLS) THEN
                        IF (EXTERNAL_WALL(IW)%NOM==0 .AND. WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY .AND. &
                            WALL(IW)%BOUNDARY_TYPE/=NULL_BOUNDARY) CYCLE
                    ELSE
                        IF (WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY .AND. WALL(IW)%BOUNDARY_TYPE/=NULL_BOUNDARY) CYCLE
                    ENDIF
                    WALL_U_NORMAL_H_TMP(IW)   = BOUNDARY_PROP1(WALL(IW)%B1_INDEX)%U_NORMAL
                    WALL_U_NORMAL_S_H_TMP(IW) = BOUNDARY_PROP1(WALL(IW)%B1_INDEX)%U_NORMAL_S
                ENDDO
                GM%WALL_U_NORMAL   = WALL_U_NORMAL_H_TMP
                GM%WALL_U_NORMAL_S = WALL_U_NORMAL_S_H_TMP
                DEALLOCATE(WALL_U_NORMAL_H_TMP,WALL_U_NORMAL_S_H_TMP,STAT=ISTAT)
            ENDIF

            ! Fill pressure-potential ghost cells from the destination-owned OMESH
            ! mirrors prepared by GPU_MESH_EXCHANGE(5).
            DO NOM=1,NMESHES
                IF (.NOT.ALLOCATED(GM%OMESH)) EXIT
                GOM_LOCAL => GM%OMESH(NOM)
                IF (.NOT.GOM_LOCAL%ACTIVE) CYCLE
                IF (PREDICTOR) THEN
                    OM_HPD => GOM_LOCAL%H_D
                ELSE
                    OM_HPD => GOM_LOCAL%HS_D
                ENDIF
                IF (PREDICTOR) THEN
                    CALL NO_FLUX_FILL_HP_GHOST_GPU<<<GM%GRDD_EXTERNAL_WALL_1D,GPU_BLKD_1D,0,GM%STRM>>>(GM%HD,OM_HPD,GM%EXTERNAL_WALL_NOMD, &
                        GM%WALL_TOPO_AOS,GM%WALL_IIO,GM%WALL_JJO,GM%WALL_KKO,GM%WALL_IIO_MAX,GM%WALL_JJO_MAX,GM%WALL_KKO_MAX, &
                        N_EXTERNAL_WALL_CELLS,NOM,LBOUND(OM_HPD,1),LBOUND(OM_HPD,2),LBOUND(OM_HPD,3), &
                        IBP1,JBP1,SIZE(OM_HPD,1),SIZE(OM_HPD,2))
                ELSE
                    CALL NO_FLUX_FILL_HP_GHOST_GPU<<<GM%GRDD_EXTERNAL_WALL_1D,GPU_BLKD_1D,0,GM%STRM>>>(GM%HSD,OM_HPD,GM%EXTERNAL_WALL_NOMD, &
                        GM%WALL_TOPO_AOS,GM%WALL_IIO,GM%WALL_JJO,GM%WALL_KKO,GM%WALL_IIO_MAX,GM%WALL_JJO_MAX,GM%WALL_KKO_MAX, &
                        N_EXTERNAL_WALL_CELLS,NOM,LBOUND(OM_HPD,1),LBOUND(OM_HPD,2),LBOUND(OM_HPD,3), &
                        IBP1,JBP1,SIZE(OM_HPD,1),SIZE(OM_HPD,2))
                ENDIF
            ENDDO

            IF (PREDICTOR) THEN
                CALL NO_FLUX_SOLID_SOLID_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1,0,GM%STRM>>>( &
                    GM%FVXD,GM%FVYD,GM%FVZD,GM%UD,GM%VD,GM%WD,GM%USD,GM%VSD,GM%WSD,GM%HD,GM%CELL_INDEXD,GM%CELL_SOLID,GM%CELL_EXTERIOR, &
                    GM%RDXND,GM%RDYND,GM%RDZND,IBAR,JBAR,KBAR,RFODT,1)
                CALL NO_FLUX_WALL_GPU<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>(GM%FVXD, GM%FVYD, GM%FVZD, GM%UD, GM%VD, GM%WD, GM%USD, GM%VSD, GM%WSD, GM%HD, GM%WALL_TOPO_AOS, &
                    GM%EXTERNAL_WALL_NOMD, GM%WALL_U_NORMAL, GM%WALL_U_NORMAL_S, GM%RDXND, GM%RDYND, GM%RDZND, N_EXTERNAL_WALL_CELLS, &
                    N_INTERNAL_WALL_CELLS,RFODT,PRES_FLAG,1,IBP1,JBP1)
            ELSE
                CALL NO_FLUX_SOLID_SOLID_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1,0,GM%STRM>>>( &
                    GM%FVXD,GM%FVYD,GM%FVZD,GM%UD,GM%VD,GM%WD,GM%USD,GM%VSD,GM%WSD,GM%HSD,GM%CELL_INDEXD,GM%CELL_SOLID,GM%CELL_EXTERIOR, &
                    GM%RDXND,GM%RDYND,GM%RDZND,IBAR,JBAR,KBAR,RFODT,0)
                CALL NO_FLUX_WALL_GPU<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>(GM%FVXD, GM%FVYD, GM%FVZD, GM%UD, GM%VD, GM%WD, GM%USD, GM%VSD, GM%WSD, GM%HSD, GM%WALL_TOPO_AOS, &
                    GM%EXTERNAL_WALL_NOMD, GM%WALL_U_NORMAL, GM%WALL_U_NORMAL_S, GM%RDXND, GM%RDYND, GM%RDZND, N_EXTERNAL_WALL_CELLS, &
                    N_INTERNAL_WALL_CELLS,RFODT,PRES_FLAG,0,IBP1,JBP1)
            ENDIF

            IF (GPU_DEBUG_H2D_D2H) THEN
                CALL GPU_WAIT(GM%STRM)
                FVX = GM%FVXD
                FVY = GM%FVYD
                FVZ = GM%FVZD
                IF (PREDICTOR) THEN
                    H = GM%HD
                ELSE
                    HS = GM%HSD
                ENDIF
            ENDIF

            T_USED(4)=T_USED(4)+CURRENT_TIME()-T_NOW

        END SUBROUTINE NO_FLUX_GPU

        SUBROUTINE PRESSURE_SOLVER_COMPUTE_RHS_GPU(T,DT,NM)

            USE MESH_POINTERS
            USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
            USE GLOBAL_CONSTANTS

            INTEGER, INTENT(IN) :: NM
            REAL(EB), INTENT(IN) :: T,DT
            REAL(EB), POINTER,DEVICE , DIMENSION(:,:,:) :: UUD,VVD,WWD,HPD,RHOPD

            TYPE(DIM3) :: BLKD , GRDD
            INTEGER , PARAMETER :: BS = 8
            INTEGER :: ISTAT,IW,BLKD_1D,GRDD_1D,NOM,N
            REAL(EB) :: TNOW,H0_PRESSURE
            REAL(EB), ALLOCATABLE, DIMENSION(:) :: WALL_DUNDT_H_TMP,WALL_WORK1_H_TMP
            REAL(EB), DEVICE, POINTER, CONTIGUOUS :: DX_OTHERD(:),DY_OTHERD(:),DZ_OTHERD(:)

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)
            BLKD = GM%BLKD_3D_BP1
            GRDD = GM%GRDD_3D_BP1

            IF (SOLID_PHASE_ONLY) RETURN
            IF (FREEZE_VELOCITY)  RETURN

            TNOW=CURRENT_TIME()
            BLKD_1D = 256
            GRDD_1D = (N_EXTERNAL_WALL_CELLS + BLKD_1D - 1) / BLKD_1D
            H0_PRESSURE = 0.5_EB*(U0**2 + V0**2 + W0**2)

            IF (GPU_DEBUG_H2D_D2H .OR. .NOT. WALL_BC_ON .OR. .NOT. DIVERGENCE_PART_1_ON) THEN
                ALLOCATE(WALL_DUNDT_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)),STAT=ISTAT)
                WALL_DUNDT_H_TMP = 0._EB
                DO IW=1,N_EXTERNAL_WALL_CELLS
                    WALL_DUNDT_H_TMP(IW) = EXTERNAL_WALL(IW)%DUNDT
                ENDDO
                GM%WALL_DUNDTD = WALL_DUNDT_H_TMP
                DEALLOCATE(WALL_DUNDT_H_TMP,STAT=ISTAT)
            ENDIF

            IF (GPU_DEBUG_H2D_D2H) THEN
                GM%FVXD  = FVX
                GM%FVYD  = FVY
                GM%FVZD  = FVZ
                GM%DDDTD = DDDT
                GM%KRESD = KRES
                DO N=LBOUND(BXS,2),UBOUND(BXS,2)
                    GM%BXSD(:,N) = BXS(:,N)
                    GM%BXFD(:,N) = BXF(:,N)
                ENDDO
                DO N=LBOUND(BYS,2),UBOUND(BYS,2)
                    GM%BYSD(:,N) = BYS(:,N)
                    GM%BYFD(:,N) = BYF(:,N)
                ENDDO
                DO N=LBOUND(BZS,2),UBOUND(BZS,2)
                    GM%BZSD(:,N) = BZS(:,N)
                    GM%BZFD(:,N) = BZF(:,N)
                ENDDO

                ALLOCATE(WALL_WORK1_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)),STAT=ISTAT)
                WALL_WORK1_H_TMP = 0._EB
                DO IW=LBOUND(WALL,1),UBOUND(WALL,1)
                    IF (IW >= LBOUND(WALL_WORK1,1) .AND. IW <= UBOUND(WALL_WORK1,1)) WALL_WORK1_H_TMP(IW) = WALL_WORK1(IW)
                ENDDO
                GM%WALL_WORK1D = WALL_WORK1_H_TMP
                DEALLOCATE(WALL_WORK1_H_TMP,STAT=ISTAT)
            ENDIF

            IF (PRESSURE_ITERATIONS==1) THEN
                GM%WALL_WORK1D = 0._EB
                IF (GPU_DEBUG_H2D_D2H) WALL_WORK1 = 0._EB
            ENDIF

            GM%PRHSD = 0._EB
            IF (PREDICTOR) THEN
                IF (GPU_DEBUG_H2D_D2H) THEN
                    GM%UD   = U
                    GM%VD   = V
                    GM%WD   = W
                    GM%HD   = H
                    GM%RHOD = RHO
                ENDIF
            ELSE
                IF (GPU_DEBUG_H2D_D2H) THEN
                    GM%USD   = US
                    GM%VSD   = VS
                    GM%WSD   = WS
                    GM%HSD   = HS
                    GM%RHOSD = RHOS
                ENDIF
            ENDIF

            IF (PREDICTOR) THEN

                UUD => GM%UD
                VVD => GM%VD
                WWD => GM%WD
                HPD => GM%HD
                RHOPD => GM%RHOD
            ELSE

                UUD => GM%USD
                VVD => GM%VSD
                WWD => GM%WSD
                HPD => GM%HSD
                RHOPD => GM%RHOSD
            ENDIF

            CALL COMPUTE_PRESSURE_BCS_GPU<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>(GM%BXSD, GM%BXFD, GM%BYSD, GM%BYFD, GM%BZSD, GM%BZFD, GM%WALL_TOPO_AOS, &
                GM%WALL_PRESSURE_BC_TYPE, GM%WALL_IIO, GM%WALL_JJO, GM%WALL_KKO, GM%WALL_DUNDTD, GM%WALL_WORK1D, GM%WALL_P_EXTERNAL_OVER_RHOD, &
                GM%FVXD, GM%FVYD, GM%FVZD, UUD, VVD, WWD, HPD, GM%KRESD, GM%HXD, GM%HYD, GM%HZD, GM%DXD, GM%DYD, GM%DZD, N_EXTERNAL_WALL_CELLS, IBAR, &
                JBAR, KBAR, IBP1, JBP1, KBP1, H0_PRESSURE)

            DO NOM=1,NMESHES
                IF (.NOT.ALLOCATED(GM%OMESH)) EXIT
                IF (.NOT.GM%OMESH(NOM)%ACTIVE) CYCLE
                DX_OTHERD => GPU_MESHES(NOM)%DXD
                DY_OTHERD => GPU_MESHES(NOM)%DYD
                DZ_OTHERD => GPU_MESHES(NOM)%DZD
                CALL COMPUTE_PRESSURE_BCS_INTERPOLATED_GPU<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>(GM%BXSD,GM%BXFD,GM%BYSD,GM%BYFD,GM%BZSD,GM%BZFD, &
                    GM%WALL_TOPO_AOS,GM%WALL_PRESSURE_BC_TYPE,GM%EXTERNAL_WALL_NOMD,GM%WALL_IIO,GM%WALL_JJO,GM%WALL_KKO, &
                    GM%WALL_WORK1D,HPD, &
                    GM%DXD,GM%DYD,GM%DZD,DX_OTHERD,DY_OTHERD,DZ_OTHERD,N_EXTERNAL_WALL_CELLS,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,NOM, &
                    GM%WALL_LB,SIZE(GM%BXSD,1),LBOUND(GM%BXSD,1),LBOUND(GM%BXSD,2), &
                    SIZE(GM%BYSD,1),LBOUND(GM%BYSD,1),LBOUND(GM%BYSD,2), &
                    SIZE(GM%BZSD,1),LBOUND(GM%BZSD,1),LBOUND(GM%BZSD,2))
            ENDDO

            CALL COMPUTE_PRHS_GPU<<<GRDD,BLKD,0,GM%STRM>>>(GM%PRHSD,GM%FVXD,GM%FVYD,GM%FVZD,GM%DDDTD,GM%RDXD,GM%RDYD,GM%RDZD,IBAR,JBAR,KBAR)

            IF (GPU_DEBUG_H2D_D2H) THEN
                CALL GPU_WAIT(GM%STRM)
                CALL COPY_PRESSURE_RHS_GPU(NM,.FALSE.)
                DO N=LBOUND(BXS,2),UBOUND(BXS,2)
                    BXS(:,N) = GM%BXSD(:,N)
                    BXF(:,N) = GM%BXFD(:,N)
                ENDDO
                DO N=LBOUND(BYS,2),UBOUND(BYS,2)
                    BYS(:,N) = GM%BYSD(:,N)
                    BYF(:,N) = GM%BYFD(:,N)
                ENDDO
                DO N=LBOUND(BZS,2),UBOUND(BZS,2)
                    BZS(:,N) = GM%BZSD(:,N)
                    BZF(:,N) = GM%BZFD(:,N)
                ENDDO
            ENDIF

            T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW
        END SUBROUTINE PRESSURE_SOLVER_COMPUTE_RHS_GPU

        SUBROUTINE PRESSURE_SOLVER_FFT_INSTEAD(NM)

            INTEGER, INTENT(IN) :: NM
            REAL(EB),POINTER ,DIMENSION(:,:,:):: HP
            REAL(EB),POINTER ,DEVICE,DIMENSION(:,:,:):: HPD
            REAL(EB) :: DX0,DY0,DZ0
            INTEGER :: N

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)
            IF (PREDICTOR) THEN
                HP => H
                HPD => GM%HD
            ELSE
                HP => HS
                HPD => GM%HSD
            ENDIF

            DX0 = DX(1)
            DY0 = DY(1)
            DZ0 = DZ(1)

            IF (GPU_DEBUG_H2D_D2H) THEN
                CALL COPY_PRESSURE_RHS_GPU(NM,.TRUE.)
                DO N=LBOUND(BXS,2),UBOUND(BXS,2)
                    GM%BXSD(:,N) = BXS(:,N)
                    GM%BXFD(:,N) = BXF(:,N)
                ENDDO
                DO N=LBOUND(BYS,2),UBOUND(BYS,2)
                    GM%BYSD(:,N) = BYS(:,N)
                    GM%BYFD(:,N) = BYF(:,N)
                ENDDO
                DO N=LBOUND(BZS,2),UBOUND(BZS,2)
                    GM%BZSD(:,N) = BZS(:,N)
                    GM%BZFD(:,N) = BZF(:,N)
                ENDDO
                HPD = HP
            ENDIF

            CALL PRESSURE_SOLVER_GPU_DISPATCH(NM,HPD,GM%PRHSD,GM%BXSD,GM%BXFD,GM%BYSD,GM%BYFD,GM%BZSD,GM%BZFD,DX0,DY0,DZ0,LBC,MBC,NBC)

            IF (GPU_DEBUG_H2D_D2H) THEN
                CALL GPU_WAIT(GM%STRM)
                HP = HPD
            ENDIF
        END SUBROUTINE PRESSURE_SOLVER_FFT_INSTEAD

        SUBROUTINE PRESSURE_SOLVER_GPU_DISPATCH(NM,HPD,PRHSD,BXSD,BXFD,BYSD,BYFD,BZSD,BZFD,DX0,DY0,DZ0,LBC,MBC,NBC)
            INTEGER, INTENT(IN) :: NM
            REAL(EB), DEVICE, INTENT(INOUT), DIMENSION(0:IBP1,0:JBP1,0:KBP1) :: HPD
            REAL(EB), DEVICE, INTENT(IN), DIMENSION(1:IBP1,1:JBP1,1:KBP1) :: PRHSD
            REAL(EB), DEVICE, INTENT(IN), DIMENSION(:,:) :: BXSD,BXFD,BYSD,BYFD,BZSD,BZFD
            REAL(EB), INTENT(IN) :: DX0,DY0,DZ0
            INTEGER, INTENT(IN) :: LBC,MBC,NBC

            IF (ALLOCATED(GPU_MESHES(NM)%POIS%STRETCH)) THEN
                CALL STRETCHED_POISSON_SOLVER_GPU(NM,HPD,PRHSD,BXSD,BXFD,BYSD,BYFD,BZSD,BZFD)
            ELSE
                CALL FISHPAK_SOLVER_GPU_SAVE1_FUSED_YZ_GLOBAL_FT( &
                    NM,HPD,PRHSD,BXSD,BXFD,BYSD,BYFD,BZSD,BZFD,DX0,DY0,DZ0,LBC,MBC,NBC,SAVE1)
            ENDIF
        END SUBROUTINE PRESSURE_SOLVER_GPU_DISPATCH


        SUBROUTINE PRESSURE_SOLVER_CHECK_RESIDUALS_GPU(NM)

            INTEGER, INTENT(IN) :: NM
            INTEGER , PARAMETER :: BS = 8
            INTEGER , PARAMETER :: BS_1D = 256
            TYPE(DIM3) :: BLKD , GRDD
            REAL(EB),DEVICE,POINTER,DIMENSION(:,:,:)::HPD,RHOPD,PD
            INTEGER :: NBLK

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)
            IF (.NOT.ITERATE_BAROCLINIC_TERM) RETURN
            IF (SOLID_PHASE_ONLY) RETURN
            IF (FREEZE_VELOCITY) RETURN

            IF (PREDICTOR) THEN
                HPD => GM%HD
                RHOPD => GM%RHOD
            ELSE
                HPD => GM%HSD
                RHOPD => GM%RHOSD
            ENDIF

            IF (GPU_DEBUG_H2D_D2H) THEN
                GM%RDXD  = RDX
                GM%RDYD  = RDY
                GM%RDZD  = RDZ
                GM%RDXND = RDXN
                GM%RDYND = RDYN
                GM%RDZND = RDZN
                GM%FVXD  = FVX
                GM%FVYD  = FVY
                GM%FVZD  = FVZ
                GM%FVX_BD = FVX_B
                GM%FVY_BD = FVY_B
                GM%FVZ_BD = FVZ_B
                GM%DDDTD = DDDT
                GM%KRESD = KRES
                GM%HD    = H
                GM%HSD   = HS
                GM%RHOD  = RHO
                GM%RHOSD = RHOS
            ENDIF

            PD => GM%WORK7D
            BLKD = GM%BLKD_3D_BP1
            GRDD = GM%GRDD_3D_BP1
            NBLK = (IBAR*JBAR*KBAR + BS_1D - 1) / BS_1D

            IF (.NOT.ALLOCATED(GM%RESIDUAL_BLOCKMAX_D)) THEN
                ALLOCATE(GM%RESIDUAL_BLOCKMAX_D(NBLK),GM%RESIDUAL_BLOCKMAX_H(NBLK))
                ALLOCATE(GM%RESIDUAL_BLOCKI_D(NBLK),GM%RESIDUAL_BLOCKI_H(NBLK))
                ALLOCATE(GM%RESIDUAL_BLOCKJ_D(NBLK),GM%RESIDUAL_BLOCKJ_H(NBLK))
                ALLOCATE(GM%RESIDUAL_BLOCKK_D(NBLK),GM%RESIDUAL_BLOCKK_H(NBLK))
            ELSEIF (SIZE(GM%RESIDUAL_BLOCKMAX_D) < NBLK) THEN
                DEALLOCATE(GM%RESIDUAL_BLOCKMAX_D,GM%RESIDUAL_BLOCKMAX_H, &
                           GM%RESIDUAL_BLOCKI_D,GM%RESIDUAL_BLOCKI_H, &
                           GM%RESIDUAL_BLOCKJ_D,GM%RESIDUAL_BLOCKJ_H, &
                           GM%RESIDUAL_BLOCKK_D,GM%RESIDUAL_BLOCKK_H)
                ALLOCATE(GM%RESIDUAL_BLOCKMAX_D(NBLK),GM%RESIDUAL_BLOCKMAX_H(NBLK))
                ALLOCATE(GM%RESIDUAL_BLOCKI_D(NBLK),GM%RESIDUAL_BLOCKI_H(NBLK))
                ALLOCATE(GM%RESIDUAL_BLOCKJ_D(NBLK),GM%RESIDUAL_BLOCKJ_H(NBLK))
                ALLOCATE(GM%RESIDUAL_BLOCKK_D(NBLK),GM%RESIDUAL_BLOCKK_H(NBLK))
            ENDIF

            CALL P_CALCULATE_KERNEL<<<GRDD,BLKD,0,GM%STRM>>>(PD,RHOPD,HPD,GM%KRESD, IBP1,JBP1,KBP1)
            CALL RESIDUALS_KERNEL<<<NBLK,BS_1D,0,GM%STRM>>>( &
                GM%WORK8D,GM%FVXD,GM%FVX_BD,GM%FVYD,GM%FVY_BD,GM%FVZD,GM%FVZ_BD,GM%DDDTD,PD,GM%KRESD,GM%WORK10D,GM%WORK11D,GM%WORK12D, &
                GM%RDXND,GM%RDYND,GM%RDZND,GM%RDXD,GM%RDYD,GM%RDZD,GM%RD,GM%RRND,IBAR,JBAR,KBAR)
            CALL RESIDUALS_BLOCKMAX_KERNEL<<<NBLK,BS_1D,0,GM%STRM>>>( &
                GM%RESIDUAL_BLOCKMAX_D, GM%RESIDUAL_BLOCKI_D, GM%RESIDUAL_BLOCKJ_D, GM%RESIDUAL_BLOCKK_D, &
                GM%WORK8D,IBAR,JBAR,KBAR)

            IF (.NOT.ALLOCATED(GM%FINAL_BLOCK_RESULT_D)) ALLOCATE(GM%FINAL_BLOCK_RESULT_D(4))

            CALL FINAL_BLOCKMAX_REDUCE_KERNEL<<<1,BS_1D,0,GM%STRM>>>( &
                GM%FINAL_BLOCK_RESULT_D, &
                GM%RESIDUAL_BLOCKMAX_D,GM%RESIDUAL_BLOCKI_D,GM%RESIDUAL_BLOCKJ_D,GM%RESIDUAL_BLOCKK_D,NBLK)

            CALL GPU_WAIT(GM%STRM)
            GM%FINAL_BLOCK_RESULT_H = GM%FINAL_BLOCK_RESULT_D
            PRESSURE_ERROR_MAX(NM) = GM%FINAL_BLOCK_RESULT_H(1)
            PRESSURE_ERROR_MAX_LOC(1:3,NM) = NINT(GM%FINAL_BLOCK_RESULT_H(2:4))

        END SUBROUTINE PRESSURE_SOLVER_CHECK_RESIDUALS_GPU

        SUBROUTINE COMPUTE_VELOCITY_ERROR_GPU(DT,NM)

            REAL(EB), INTENT(IN) :: DT
            INTEGER, INTENT(IN) :: NM
            INTEGER, PARAMETER :: BS_1D = 64
            INTEGER :: BLKD_1D,GRDD_1D,ISTAT,NBLK,IW,NOM,AXIS
            REAL(EB) :: ITERATIVE_FACTOR,TNOW
            REAL(EB), ALLOCATABLE, DIMENSION(:) :: WALL_WORK1_GPU,WALL_VEL_ERR_H_TMP,WALL_U_NORMAL_H_TMP, &
                                                    WALL_U_NORMAL_S_H_TMP,WALL_RDN_H_TMP
            REAL(EB), DEVICE, ALLOCATABLE, DIMENSION(:) :: WALL_VEL_ERRD
            REAL(EB), DEVICE, POINTER, CONTIGUOUS :: OM_VELD(:,:,:),OM_VEL_SD(:,:,:),OM_FLUXD(:,:,:),OM_HPD(:,:,:)
            REAL(EB), DEVICE, POINTER, CONTIGUOUS :: RDN_OTHERD(:)
            REAL(EB), DEVICE, POINTER, CONTIGUOUS :: HPD(:,:,:)
            TYPE(OMESH_TYPE), POINTER :: OM
            TYPE(GPU_OMESH_TYPE), POINTER :: GOM_LOCAL

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)
            IF (SOLID_PHASE_ONLY) RETURN
            IF (FREEZE_VELOCITY)  RETURN

            TNOW = CURRENT_TIME()

            IF (PREDICTOR) THEN
                ITERATIVE_FACTOR = 0.25_EB
            ELSE
                ITERATIVE_FACTOR = 0.50_EB
            ENDIF

            BLKD_1D = GPU_BLKD_1D
            GRDD_1D = GM%GRDD_TOTAL_WALL_1D

            IF (GPU_DEBUG_H2D_D2H) THEN
                GM%UD   = U
                GM%VD   = V
                GM%WD   = W
                GM%USD  = US
                GM%VSD  = VS
                GM%WSD  = WS
                GM%FVXD = FVX
                GM%FVYD = FVY
                GM%FVZD = FVZ
                GM%HD   = H
                GM%HSD  = HS
                GM%RDXND = RDXN
                GM%RDYND = RDYN
                GM%RDZND = RDZN

                ALLOCATE(WALL_U_NORMAL_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)),STAT=ISTAT)
                ALLOCATE(WALL_U_NORMAL_S_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)),STAT=ISTAT)
                ALLOCATE(WALL_RDN_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)),STAT=ISTAT)
                WALL_U_NORMAL_H_TMP   = 0._EB
                WALL_U_NORMAL_S_H_TMP = 0._EB
                WALL_RDN_H_TMP        = 1._EB
                DO IW=LBOUND(WALL,1),UBOUND(WALL,1)
                    IF (WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE
                    WALL_U_NORMAL_H_TMP(IW)   = BOUNDARY_PROP1(WALL(IW)%B1_INDEX)%U_NORMAL
                    WALL_U_NORMAL_S_H_TMP(IW) = BOUNDARY_PROP1(WALL(IW)%B1_INDEX)%U_NORMAL_S
                    WALL_RDN_H_TMP(IW)        = BOUNDARY_PROP1(WALL(IW)%B1_INDEX)%RDN
                ENDDO
                GM%WALL_U_NORMAL   = WALL_U_NORMAL_H_TMP
                GM%WALL_U_NORMAL_S = WALL_U_NORMAL_S_H_TMP
                GM%WALL_RDN        = WALL_RDN_H_TMP
                DEALLOCATE(WALL_U_NORMAL_H_TMP,WALL_U_NORMAL_S_H_TMP,WALL_RDN_H_TMP,STAT=ISTAT)

                DO NOM=1,NMESHES
                    IF (.NOT.ALLOCATED(GM%OMESH)) EXIT
                    GOM_LOCAL => GM%OMESH(NOM)
                    IF (.NOT.GOM_LOCAL%ACTIVE) CYCLE
                    OM => MESHES(NM)%OMESH(NOM)
                    GOM_LOCAL%U_D = OM%U
                    GOM_LOCAL%V_D = OM%V
                    GOM_LOCAL%W_D = OM%W
                    GOM_LOCAL%US_D = OM%US
                    GOM_LOCAL%VS_D = OM%VS
                    GOM_LOCAL%WS_D = OM%WS
                    GOM_LOCAL%FVX_D = OM%FVX
                    GOM_LOCAL%FVY_D = OM%FVY
                    GOM_LOCAL%FVZ_D = OM%FVZ
                    GOM_LOCAL%H_D = OM%H
                    GOM_LOCAL%HS_D = OM%HS
                ENDDO
            ENDIF

            ALLOCATE(WALL_VEL_ERRD(LBOUND(WALL,1):UBOUND(WALL,1)))
            GM%WALL_WORK1D = 0._EB
            WALL_VEL_ERRD = 0._EB

            IF (PREDICTOR) THEN
                HPD => GM%HD
            ELSE
                HPD => GM%HSD
            ENDIF

            CALL COMPUTE_VELOCITY_ERROR_LOCAL_GPU<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>(GM%WALL_WORK1D,WALL_VEL_ERRD,GM%WALL_TOPO_AOS, &
                GM%WALL_U_NORMAL,GM%WALL_U_NORMAL_S,GM%WALL_RDN,GM%EXTERNAL_WALL_NOMD,GM%EXTERNAL_WALL_AREA_RATIOD,GM%UD,GM%VD,GM%WD,GM%USD,GM%VSD,GM%WSD, &
                GM%FVXD,GM%FVYD,GM%FVZD,HPD,GM%RDXND,GM%RDYND,GM%RDZND,DT,ITERATIVE_FACTOR,N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS, &
                PRES_FLAG,INTERPOLATED_BOUNDARY,MERGE(1,0,PREDICTOR),IBP1,JBP1)

            DO NOM=1,NMESHES
                IF (.NOT.ALLOCATED(GM%OMESH)) EXIT
                GOM_LOCAL => GM%OMESH(NOM)
                IF (.NOT.GOM_LOCAL%ACTIVE) CYCLE
                IF (PREDICTOR) THEN
                    OM_HPD => GOM_LOCAL%H_D
                ELSE
                    OM_HPD => GOM_LOCAL%HS_D
                ENDIF
                DO AXIS=1,3
                    SELECT CASE(AXIS)
                        CASE(1)
                            OM_VELD => GOM_LOCAL%U_D
                            OM_VEL_SD => GOM_LOCAL%US_D
                            OM_FLUXD => GOM_LOCAL%FVX_D
                            RDN_OTHERD => GPU_MESHES(NOM)%RDXND
                        CASE(2)
                            OM_VELD => GOM_LOCAL%V_D
                            OM_VEL_SD => GOM_LOCAL%VS_D
                            OM_FLUXD => GOM_LOCAL%FVY_D
                            RDN_OTHERD => GPU_MESHES(NOM)%RDYND
                        CASE(3)
                            OM_VELD => GOM_LOCAL%W_D
                            OM_VEL_SD => GOM_LOCAL%WS_D
                            OM_FLUXD => GOM_LOCAL%FVZ_D
                            RDN_OTHERD => GPU_MESHES(NOM)%RDZND
                    END SELECT
                    CALL COMPUTE_VELOCITY_ERROR_INTERPOLATED_GPU<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>(GM%WALL_WORK1D,WALL_VEL_ERRD, &
                        GM%WALL_TOPO_AOS,GM%WALL_RDN,OM_VELD,OM_VEL_SD,OM_FLUXD,OM_HPD,RDN_OTHERD,GM%EXTERNAL_WALL_NOMD,GM%WALL_IIO, &
                        GM%WALL_JJO,GM%WALL_KKO,GM%WALL_IIO_MAX,GM%WALL_JJO_MAX,GM%WALL_KKO_MAX,GM%EXTERNAL_WALL_AREA_RATIOD,DT, &
                        ITERATIVE_FACTOR,N_EXTERNAL_WALL_CELLS,INTERPOLATED_BOUNDARY,MERGE(1,0,PREDICTOR),NOM,AXIS, &
                        SIZE(OM_VELD,1),SIZE(OM_VELD,2),LBOUND(OM_VELD,1),LBOUND(OM_VELD,2),LBOUND(OM_VELD,3), &
                        LBOUND(RDN_OTHERD,1),GM%WALL_LB)
                ENDDO
            ENDDO

            IF (GPU_DEBUG_H2D_D2H) THEN
                CALL GPU_WAIT(GM%STRM)
                ALLOCATE(WALL_WORK1_GPU(LBOUND(WALL,1):UBOUND(WALL,1)))
                ALLOCATE(WALL_VEL_ERR_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)))
                WALL_WORK1_GPU = GM%WALL_WORK1D
                WALL_VEL_ERR_H_TMP = WALL_VEL_ERRD
                WALL_WORK1 = WALL_WORK1_GPU(LBOUND(WALL_WORK1,1):UBOUND(WALL_WORK1,1))
                DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
                    IF (WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY .AND. &
                        WALL(IW)%BOUNDARY_TYPE/=INTERPOLATED_BOUNDARY) CYCLE
                    IF (WALL(IW)%CUT_FACE_INDEX>0) CYCLE
                    IF (WALL(IW)%BOUNDARY_TYPE==INTERPOLATED_BOUNDARY .AND. EXTERNAL_WALL(IW)%AREA_RATIO<0.9_EB) CYCLE
                    BOUNDARY_PROP1(WALL(IW)%B1_INDEX)%VEL_ERR_NEW = WALL_VEL_ERR_H_TMP(IW)
                ENDDO
                DEALLOCATE(WALL_WORK1_GPU,WALL_VEL_ERR_H_TMP)
            ENDIF
            VELOCITY_ERROR_MAX(NM) = 0._EB
            VELOCITY_ERROR_MAX_LOC(:,NM) = 0

            NBLK = GRDD_1D
            IF (.NOT.ALLOCATED(GM%VERR_BLOCKMAX_D)) THEN
                ALLOCATE(GM%VERR_BLOCKMAX_D(NBLK),GM%VERR_BLOCKMAX_H(NBLK))
                ALLOCATE(GM%VERR_BLOCKI_D(NBLK),GM%VERR_BLOCKI_H(NBLK))
                ALLOCATE(GM%VERR_BLOCKJ_D(NBLK),GM%VERR_BLOCKJ_H(NBLK))
                ALLOCATE(GM%VERR_BLOCKK_D(NBLK),GM%VERR_BLOCKK_H(NBLK))
            ELSEIF (SIZE(GM%VERR_BLOCKMAX_D) < NBLK) THEN
                DEALLOCATE(GM%VERR_BLOCKMAX_D,GM%VERR_BLOCKMAX_H, &
                           GM%VERR_BLOCKI_D,GM%VERR_BLOCKI_H, &
                           GM%VERR_BLOCKJ_D,GM%VERR_BLOCKJ_H, &
                           GM%VERR_BLOCKK_D,GM%VERR_BLOCKK_H)
                ALLOCATE(GM%VERR_BLOCKMAX_D(NBLK),GM%VERR_BLOCKMAX_H(NBLK))
                ALLOCATE(GM%VERR_BLOCKI_D(NBLK),GM%VERR_BLOCKI_H(NBLK))
                ALLOCATE(GM%VERR_BLOCKJ_D(NBLK),GM%VERR_BLOCKJ_H(NBLK))
                ALLOCATE(GM%VERR_BLOCKK_D(NBLK),GM%VERR_BLOCKK_H(NBLK))
            ENDIF

            CALL VELOCITY_ERROR_SCAN_BLOCKMAX_KERNEL<<<NBLK,BLKD_1D,0,GM%STRM>>>(GM%VERR_BLOCKMAX_D, GM%VERR_BLOCKI_D, GM%VERR_BLOCKJ_D, &
                GM%VERR_BLOCKK_D, WALL_VEL_ERRD, GM%WALL_TOPO_AOS, N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS, &
                INTERPOLATED_BOUNDARY,GM%EXTERNAL_WALL_AREA_RATIOD)

            IF (.NOT.ALLOCATED(GM%FINAL_BLOCK_RESULT_D)) ALLOCATE(GM%FINAL_BLOCK_RESULT_D(4))

            CALL FINAL_BLOCKMAX_REDUCE_KERNEL<<<1,256,0,GM%STRM>>>( &
                GM%FINAL_BLOCK_RESULT_D, &
                GM%VERR_BLOCKMAX_D,GM%VERR_BLOCKI_D,GM%VERR_BLOCKJ_D,GM%VERR_BLOCKK_D,NBLK)

            CALL GPU_WAIT(GM%STRM)
            GM%FINAL_BLOCK_RESULT_H = GM%FINAL_BLOCK_RESULT_D
            VELOCITY_ERROR_MAX(NM) = GM%FINAL_BLOCK_RESULT_H(1)
            VELOCITY_ERROR_MAX_LOC(1:3,NM) = NINT(GM%FINAL_BLOCK_RESULT_H(2:4))

            DEALLOCATE(WALL_VEL_ERRD)

            T_USED(5)=T_USED(5)+CURRENT_TIME()-TNOW

        END SUBROUTINE COMPUTE_VELOCITY_ERROR_GPU

        !--------------------------------------------------------------------------
        ! Pressure, no-flux, residual, and velocity-error kernels
        !--------------------------------------------------------------------------
        ATTRIBUTES(GLOBAL) SUBROUTINE COMPUTE_VELOCITY_ERROR_LOCAL_GPU(WALL_WORK1D,WALL_VEL_ERRD,WALL_TOPO,WALL_U_NORMAL, &
            WALL_U_NORMAL_S,WALL_RDN,EXTERNAL_WALL_NOM,EXTERNAL_WALL_AREA_RATIO,U,V,W,US,VS,WS,FVX,FVY,FVZ,HP,RDXN,RDYN, &
            RDZN,DT,ITERATIVE_FACTOR,N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,PRES_FLAG_LOCAL, &
            INTERPOLATED_BOUNDARY_VALUE,PREDICTOR_FLAG,IBP1,JBP1)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)

            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,PRES_FLAG_LOCAL,INTERPOLATED_BOUNDARY_VALUE
            INTEGER, VALUE :: PREDICTOR_FLAG,IBP1,JBP1
            REAL(EB), VALUE :: DT,ITERATIVE_FACTOR
            REAL(EB), DEVICE, DIMENSION(0:*) :: WALL_WORK1D,WALL_VEL_ERRD,WALL_U_NORMAL,WALL_U_NORMAL_S,WALL_RDN
            INTEGER, DEVICE, DIMENSION(1:*) :: EXTERNAL_WALL_NOM
            REAL(EB), DEVICE, DIMENSION(1:*) :: EXTERNAL_WALL_AREA_RATIO
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: U,V,W,US,VS,WS,FVX,FVY,FVZ,HP
            REAL(EB), DEVICE, DIMENSION(0:*) :: RDXN,RDYN,RDZN
            INTEGER :: IW,II,JJ,KK,IOR,BTYPE
            REAL(EB) :: UN_NEW,UN_NEW_OTHER,VELOCITY_ERROR,DHFCT,SGN_IOR

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS) RETURN
            BTYPE = WALL_TOPO(IW)%BOUNDARY_TYPE
            IF (BTYPE /= SOLID_BOUNDARY .AND. BTYPE /= INTERPOLATED_BOUNDARY_VALUE) RETURN
            IF (WALL_TOPO(IW)%CUT_FACE_INDEX > 0) RETURN
            IF (BTYPE == INTERPOLATED_BOUNDARY_VALUE) THEN
                IF (IW > N_EXTERNAL_WALL_CELLS) RETURN
                IF (EXTERNAL_WALL_NOM(IW) <= 0) RETURN
                IF (EXTERNAL_WALL_AREA_RATIO(IW) < 0.9_EB) RETURN
            ENDIF

            II  = WALL_TOPO(IW)%II
            JJ  = WALL_TOPO(IW)%JJ
            KK  = WALL_TOPO(IW)%KK
            IOR = WALL_TOPO(IW)%IOR
            SGN_IOR = SIGN(1._EB,REAL(IOR,EB))
            DHFCT = 1._EB
            IF (BTYPE == SOLID_BOUNDARY) THEN
                SELECT CASE(PRES_FLAG_LOCAL)
                    CASE(UGLMAT_FLAG,ULMAT_FLAG)
                        DHFCT = 0._EB
                    CASE(GLMAT_FLAG)
                        IF (IW <= N_EXTERNAL_WALL_CELLS) DHFCT = 0._EB
                END SELECT
            ENDIF

            IF (PREDICTOR_FLAG == 1) THEN
                SELECT CASE(IOR)
                    CASE( 1)
                        UN_NEW = U(II,JJ,KK)   - DT*(FVX(II,JJ,KK)   + RDXN(II)  *(HP(II+1,JJ,KK)-HP(II,JJ,KK))*DHFCT)
                    CASE(-1)
                        UN_NEW = U(II-1,JJ,KK) - DT*(FVX(II-1,JJ,KK) + RDXN(II-1)*(HP(II,JJ,KK)-HP(II-1,JJ,KK))*DHFCT)
                    CASE( 2)
                        UN_NEW = V(II,JJ,KK)   - DT*(FVY(II,JJ,KK)   + RDYN(JJ)  *(HP(II,JJ+1,KK)-HP(II,JJ,KK))*DHFCT)
                    CASE(-2)
                        UN_NEW = V(II,JJ-1,KK) - DT*(FVY(II,JJ-1,KK) + RDYN(JJ-1)*(HP(II,JJ,KK)-HP(II,JJ-1,KK))*DHFCT)
                    CASE( 3)
                        UN_NEW = W(II,JJ,KK)   - DT*(FVZ(II,JJ,KK)   + RDZN(KK)  *(HP(II,JJ,KK+1)-HP(II,JJ,KK))*DHFCT)
                    CASE(-3)
                        UN_NEW = W(II,JJ,KK-1) - DT*(FVZ(II,JJ,KK-1) + RDZN(KK-1)*(HP(II,JJ,KK)-HP(II,JJ,KK-1))*DHFCT)
                END SELECT
            ELSE
                SELECT CASE(IOR)
                    CASE( 1)
                        UN_NEW = 0.5_EB*(U(II,JJ,KK)+US(II,JJ,KK)     -DT*(FVX(II,JJ,KK)  +RDXN(II)  *(HP(II+1,JJ,KK)-HP(II,JJ,KK))*DHFCT))
                    CASE(-1)
                        UN_NEW = 0.5_EB*(U(II-1,JJ,KK)+US(II-1,JJ,KK) -DT*(FVX(II-1,JJ,KK)+RDXN(II-1)*(HP(II,JJ,KK)-HP(II-1,JJ,KK))*DHFCT))
                    CASE( 2)
                        UN_NEW = 0.5_EB*(V(II,JJ,KK)+VS(II,JJ,KK)     -DT*(FVY(II,JJ,KK)  +RDYN(JJ)  *(HP(II,JJ+1,KK)-HP(II,JJ,KK))*DHFCT))
                    CASE(-2)
                        UN_NEW = 0.5_EB*(V(II,JJ-1,KK)+VS(II,JJ-1,KK) -DT*(FVY(II,JJ-1,KK)+RDYN(JJ-1)*(HP(II,JJ,KK)-HP(II,JJ-1,KK))*DHFCT))
                    CASE( 3)
                        UN_NEW = 0.5_EB*(W(II,JJ,KK)+WS(II,JJ,KK)     -DT*(FVZ(II,JJ,KK)  +RDZN(KK)  *(HP(II,JJ,KK+1)-HP(II,JJ,KK))*DHFCT))
                    CASE(-3)
                        UN_NEW = 0.5_EB*(W(II,JJ,KK-1)+WS(II,JJ,KK-1) -DT*(FVZ(II,JJ,KK-1)+RDZN(KK-1)*(HP(II,JJ,KK)-HP(II,JJ,KK-1))*DHFCT))
                END SELECT
            ENDIF

            IF (BTYPE == INTERPOLATED_BOUNDARY_VALUE) THEN
                WALL_VEL_ERRD(IW) = UN_NEW
                RETURN
            ENDIF

            IF (PREDICTOR_FLAG == 1) THEN
                UN_NEW_OTHER = -SGN_IOR*WALL_U_NORMAL_S(IW)
            ELSE
                UN_NEW_OTHER = -SGN_IOR*WALL_U_NORMAL(IW)
            ENDIF

            VELOCITY_ERROR = UN_NEW - UN_NEW_OTHER
            WALL_VEL_ERRD(IW) = VELOCITY_ERROR
            WALL_WORK1D(IW) = -SGN_IOR*ITERATIVE_FACTOR*VELOCITY_ERROR/(WALL_RDN(IW)*DT)

        END SUBROUTINE COMPUTE_VELOCITY_ERROR_LOCAL_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE COMPUTE_VELOCITY_ERROR_INTERPOLATED_GPU(WALL_WORK1D,WALL_VEL_ERRD,WALL_TOPO,WALL_RDN, &
            OM_VEL,OM_VEL_S,OM_FLUX,OM_HP,RDN_OTHER,EXTERNAL_WALL_NOM,WALL_IIO,WALL_JJO,WALL_KKO,WALL_IIO_MAX,WALL_JJO_MAX, &
            WALL_KKO_MAX,EXTERNAL_WALL_AREA_RATIO,DT,ITERATIVE_FACTOR,N_EXTERNAL_WALL_CELLS,INTERPOLATED_BOUNDARY_VALUE, &
            PREDICTOR_FLAG,NOM_FILTER,AXIS,OM_NI,OM_NJ,OM_LB_I,OM_LB_J,OM_LB_K,RDN_LB,WALL_LB)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(*)

            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,INTERPOLATED_BOUNDARY_VALUE,PREDICTOR_FLAG,NOM_FILTER,AXIS
            INTEGER, VALUE :: OM_NI,OM_NJ,OM_LB_I,OM_LB_J,OM_LB_K,RDN_LB,WALL_LB
            REAL(EB), VALUE :: DT,ITERATIVE_FACTOR
            REAL(EB), DEVICE, DIMENSION(*) :: WALL_WORK1D,WALL_VEL_ERRD,WALL_RDN
            REAL(EB), DEVICE, DIMENSION(*) :: OM_VEL,OM_VEL_S,OM_FLUX,OM_HP,RDN_OTHER
            INTEGER, DEVICE, DIMENSION(*) :: EXTERNAL_WALL_NOM,WALL_IIO,WALL_JJO,WALL_KKO
            INTEGER, DEVICE, DIMENSION(*) :: WALL_IIO_MAX,WALL_JJO_MAX,WALL_KKO_MAX
            REAL(EB), DEVICE, DIMENSION(*) :: EXTERNAL_WALL_AREA_RATIO
            INTEGER :: IW,IW_WALL,IOR,IIO,JJO,KKO,IIO1,IIO2,JJO1,JJO2,KKO1,KKO2,N_INT_CELLS
            INTEGER :: IDX_C,IDX_V,STRIDE,RDN_COORD,RDN_INDEX
            REAL(EB) :: UN_NEW,UN_NEW_OTHER,VELOCITY_ERROR,DU_DT,SGN_IOR

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS) RETURN
            IW_WALL = IW - WALL_LB + 1
            IF (WALL_TOPO(IW_WALL)%BOUNDARY_TYPE /= INTERPOLATED_BOUNDARY_VALUE) RETURN
            IF (WALL_TOPO(IW_WALL)%CUT_FACE_INDEX > 0) RETURN
            IF (EXTERNAL_WALL_NOM(IW) /= NOM_FILTER) RETURN
            IF (EXTERNAL_WALL_AREA_RATIO(IW) < 0.9_EB) RETURN

            IOR = WALL_TOPO(IW_WALL)%IOR
            IF (ABS(IOR) /= AXIS) RETURN
            SGN_IOR = SIGN(1._EB,REAL(IOR,EB))
            SELECT CASE(AXIS)
                CASE(1)
                    STRIDE = 1
                CASE(2)
                    STRIDE = OM_NI
                CASE(3)
                    STRIDE = OM_NI*OM_NJ
            END SELECT

            IIO1 = WALL_IIO(IW)
            JJO1 = WALL_JJO(IW)
            KKO1 = WALL_KKO(IW)
            IIO2 = WALL_IIO_MAX(IW)
            JJO2 = WALL_JJO_MAX(IW)
            KKO2 = WALL_KKO_MAX(IW)
            N_INT_CELLS = (IIO2-IIO1+1)*(JJO2-JJO1+1)*(KKO2-KKO1+1)
            IF (N_INT_CELLS <= 0) RETURN

            UN_NEW_OTHER = 0._EB
            DO KKO=KKO1,KKO2
                DO JJO=JJO1,JJO2
                    DO IIO=IIO1,IIO2
                        IDX_C = 1 + (IIO-OM_LB_I) + OM_NI*((JJO-OM_LB_J) + OM_NJ*(KKO-OM_LB_K))
                        SELECT CASE(AXIS)
                            CASE(1)
                                RDN_COORD = IIO
                            CASE(2)
                                RDN_COORD = JJO
                            CASE(3)
                                RDN_COORD = KKO
                        END SELECT
                        IF (IOR > 0) THEN
                            IDX_V = IDX_C
                            RDN_INDEX = RDN_COORD - RDN_LB + 1
                            DU_DT = -OM_FLUX(IDX_V) - RDN_OTHER(RDN_INDEX)*(OM_HP(IDX_C+STRIDE)-OM_HP(IDX_C))
                        ELSE
                            IDX_V = IDX_C - STRIDE
                            RDN_INDEX = RDN_COORD - 1 - RDN_LB + 1
                            DU_DT = -OM_FLUX(IDX_V) - RDN_OTHER(RDN_INDEX)*(OM_HP(IDX_C)-OM_HP(IDX_C-STRIDE))
                        ENDIF
                        IF (PREDICTOR_FLAG == 1) THEN
                            UN_NEW_OTHER = UN_NEW_OTHER + OM_VEL(IDX_V) + DT*DU_DT
                        ELSE
                            UN_NEW_OTHER = UN_NEW_OTHER + 0.5_EB*(OM_VEL(IDX_V)+OM_VEL_S(IDX_V)+DT*DU_DT)
                        ENDIF
                    ENDDO
                ENDDO
            ENDDO
            UN_NEW_OTHER = UN_NEW_OTHER/REAL(N_INT_CELLS,EB)

            UN_NEW = WALL_VEL_ERRD(IW_WALL)
            VELOCITY_ERROR = UN_NEW - UN_NEW_OTHER
            WALL_VEL_ERRD(IW_WALL) = VELOCITY_ERROR
            WALL_WORK1D(IW_WALL) = -SGN_IOR*ITERATIVE_FACTOR*VELOCITY_ERROR/(WALL_RDN(IW_WALL)*DT)

        END SUBROUTINE COMPUTE_VELOCITY_ERROR_INTERPOLATED_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE VELOCITY_ERROR_SCAN_BLOCKMAX_KERNEL(BLOCKMAX, BLOCKI, BLOCKJ, BLOCKK, WALL_VEL_ERRD, &
            WALL_TOPO,N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,INTERPOLATED_BOUNDARY_VALUE,EXTERNAL_WALL_AREA_RATIO)
            IMPLICIT NONE
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)

            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,INTERPOLATED_BOUNDARY_VALUE
            REAL(EB), DEVICE, DIMENSION(*) :: BLOCKMAX
            INTEGER, DEVICE, DIMENSION(*) :: BLOCKI,BLOCKJ,BLOCKK
            REAL(EB), DEVICE, DIMENSION(0:*) :: WALL_VEL_ERRD
            REAL(EB), DEVICE, DIMENSION(1:*) :: EXTERNAL_WALL_AREA_RATIO

            REAL(EB), SHARED :: SMAX(64)
            INTEGER, SHARED :: SI(64), SJ(64), SK(64)

            INTEGER :: TID,IW,NWALL,II,JJ,KK,IOR
            LOGICAL :: ACTIVE_WALL
            REAL(EB) :: VAL

            TID = THREADIDX%X
            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + TID
            NWALL = N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS

            VAL = -1._EB
            II = 0
            JJ = 0
            KK = 0
            ACTIVE_WALL = .FALSE.

            IF (IW <= NWALL) THEN
                IF (WALL_TOPO(IW)%CUT_FACE_INDEX == 0) THEN
                    IF (WALL_TOPO(IW)%BOUNDARY_TYPE == SOLID_BOUNDARY) ACTIVE_WALL = .TRUE.
                    IF (WALL_TOPO(IW)%BOUNDARY_TYPE == INTERPOLATED_BOUNDARY_VALUE .AND. IW <= N_EXTERNAL_WALL_CELLS) THEN
                        IF (EXTERNAL_WALL_AREA_RATIO(IW) >= 0.9_EB) ACTIVE_WALL = .TRUE.
                    ENDIF
                ENDIF
                IF (ACTIVE_WALL) THEN
                    VAL = ABS(WALL_VEL_ERRD(IW))
                    II = WALL_TOPO(IW)%II
                    JJ = WALL_TOPO(IW)%JJ
                    KK = WALL_TOPO(IW)%KK
                    IOR = WALL_TOPO(IW)%IOR
                    SELECT CASE(IOR)
                        CASE(-1) ; II = II - 1
                        CASE(-2) ; JJ = JJ - 1
                        CASE(-3) ; KK = KK - 1
                    END SELECT
                ENDIF
            ENDIF

            SMAX(TID) = VAL
            SI(TID) = II
            SJ(TID) = JJ
            SK(TID) = KK
            CALL SYNCTHREADS()

            IF (TID <= 32) THEN
                IF (SMAX(TID+32) > SMAX(TID)) THEN
                    SMAX(TID) = SMAX(TID+32)
                    SI(TID) = SI(TID+32)
                    SJ(TID) = SJ(TID+32)
                    SK(TID) = SK(TID+32)
                ENDIF
            ENDIF
            CALL SYNCTHREADS()

            IF (TID <= 16) THEN
                IF (SMAX(TID+16) > SMAX(TID)) THEN
                    SMAX(TID) = SMAX(TID+16)
                    SI(TID) = SI(TID+16)
                    SJ(TID) = SJ(TID+16)
                    SK(TID) = SK(TID+16)
                ENDIF
            ENDIF
            CALL SYNCTHREADS()

            IF (TID <= 8) THEN
                IF (SMAX(TID+8) > SMAX(TID)) THEN
                    SMAX(TID) = SMAX(TID+8)
                    SI(TID) = SI(TID+8)
                    SJ(TID) = SJ(TID+8)
                    SK(TID) = SK(TID+8)
                ENDIF
            ENDIF
            CALL SYNCTHREADS()

            IF (TID <= 4) THEN
                IF (SMAX(TID+4) > SMAX(TID)) THEN
                    SMAX(TID) = SMAX(TID+4)
                    SI(TID) = SI(TID+4)
                    SJ(TID) = SJ(TID+4)
                    SK(TID) = SK(TID+4)
                ENDIF
            ENDIF
            CALL SYNCTHREADS()

            IF (TID <= 2) THEN
                IF (SMAX(TID+2) > SMAX(TID)) THEN
                    SMAX(TID) = SMAX(TID+2)
                    SI(TID) = SI(TID+2)
                    SJ(TID) = SJ(TID+2)
                    SK(TID) = SK(TID+2)
                ENDIF
            ENDIF
            CALL SYNCTHREADS()

            IF (TID == 1) THEN
                IF (SMAX(2) > SMAX(1)) THEN
                    SMAX(1) = SMAX(2)
                    SI(1) = SI(2)
                    SJ(1) = SJ(2)
                    SK(1) = SK(2)
                ENDIF
                BLOCKMAX(BLOCKIDX%X) = SMAX(1)
                BLOCKI(BLOCKIDX%X) = SI(1)
                BLOCKJ(BLOCKIDX%X) = SJ(1)
                BLOCKK(BLOCKIDX%X) = SK(1)
            ENDIF

        END SUBROUTINE VELOCITY_ERROR_SCAN_BLOCKMAX_KERNEL

        ATTRIBUTES(GLOBAL) SUBROUTINE NO_FLUX_FILL_HP_GHOST_GPU(HP, OM_HP, EXTERNAL_WALL_NOM, WALL_TOPO, WALL_IIO, WALL_JJO, &
            WALL_KKO,WALL_IIO_MAX,WALL_JJO_MAX,WALL_KKO_MAX,N_EXTERNAL_WALL_CELLS,NOM_FILTER,OM_LB_I,OM_LB_J,OM_LB_K, &
            IBP1,JBP1,OM_NI,OM_NJ)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,NOM_FILTER,OM_LB_I,OM_LB_J,OM_LB_K,IBP1,JBP1,OM_NI,OM_NJ
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: HP
            REAL(EB), DEVICE, DIMENSION(0:OM_NI-1,0:OM_NJ-1,0:*) :: OM_HP
            INTEGER, DEVICE, DIMENSION(1:*) :: EXTERNAL_WALL_NOM,WALL_IIO,WALL_JJO,WALL_KKO
            INTEGER, DEVICE, DIMENSION(1:*) :: WALL_IIO_MAX,WALL_JJO_MAX,WALL_KKO_MAX
            INTEGER :: IW,II,JJ,KK,IIO,JJO,KKO,IIO1,IIO2,JJO1,JJO2,KKO1,KKO2,N_INT_CELLS
            REAL(EB) :: H_OTHER

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS) RETURN
            IF (EXTERNAL_WALL_NOM(IW) /= NOM_FILTER) RETURN

            II  = WALL_TOPO(IW)%II
            JJ  = WALL_TOPO(IW)%JJ
            KK  = WALL_TOPO(IW)%KK
            IIO1 = WALL_IIO(IW)
            JJO1 = WALL_JJO(IW)
            KKO1 = WALL_KKO(IW)
            IIO2 = WALL_IIO_MAX(IW)
            JJO2 = WALL_JJO_MAX(IW)
            KKO2 = WALL_KKO_MAX(IW)

            H_OTHER = 0._EB
            DO KKO=KKO1,KKO2
                DO JJO=JJO1,JJO2
                    DO IIO=IIO1,IIO2
                        H_OTHER = H_OTHER + OM_HP(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)
                    ENDDO
                ENDDO
            ENDDO
            N_INT_CELLS = (IIO2-IIO1+1)*(JJO2-JJO1+1)*(KKO2-KKO1+1)
            IF (N_INT_CELLS <= 0) RETURN
            HP(II,JJ,KK) = H_OTHER/REAL(N_INT_CELLS,EB)

        END SUBROUTINE NO_FLUX_FILL_HP_GHOST_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE NO_FLUX_SOLID_SOLID_GPU(FVX,FVY,FVZ,U,V,W,US,VS,WS,HP, &
                                                               CELL_INDEXD,CELL_SOLID,CELL_EXTERIOR, &
                                                               RDXN,RDYN,RDZN,IBAR,JBAR,KBAR,RFODT,PREDICTOR_FLAG)
            INTEGER, VALUE :: IBAR, JBAR, KBAR, PREDICTOR_FLAG
            REAL(EB), VALUE :: RFODT
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: FVX,FVY,FVZ,U,V,W,US,VS,WS,HP
            INTEGER, DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: CELL_INDEXD
            LOGICAL, DEVICE, DIMENSION(0:*) :: CELL_SOLID,CELL_EXTERIOR
            REAL(EB), DEVICE, DIMENSION(0:*) :: RDXN,RDYN,RDZN
            INTEGER :: I,J,K,IC1,IC2
            LOGICAL :: SOLID1,SOLID2
            REAL(EB) :: DUUDT,DVVDT,DWWDT

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBAR) RETURN
            IF (J < 0 .OR. J > JBAR) RETURN
            IF (K < 0 .OR. K > KBAR) RETURN

            ! Static OBST hydro-only path: drive velocity faces embedded between two solid cells to zero.
            IF (I >= 0 .AND. I <= IBAR .AND. J >= 1 .AND. J <= JBAR .AND. K >= 1 .AND. K <= KBAR) THEN
                IC1 = CELL_INDEXD(I,J,K)
                IC2 = CELL_INDEXD(I+1,J,K)
                SOLID1 = CELL_SOLID(IC1) .AND. .NOT. CELL_EXTERIOR(IC1)
                SOLID2 = CELL_SOLID(IC2) .AND. .NOT. CELL_EXTERIOR(IC2)
                IF (SOLID1 .AND. SOLID2) THEN
                    IF (PREDICTOR_FLAG == 1) THEN
                        DUUDT = -RFODT*U(I,J,K)
                    ELSE
                        DUUDT = -RFODT*(U(I,J,K)+US(I,J,K))
                    ENDIF
                    FVX(I,J,K) = -RDXN(I)*(HP(I+1,J,K)-HP(I,J,K)) - DUUDT
                ENDIF
            ENDIF

            IF (I >= 1 .AND. I <= IBAR .AND. J >= 0 .AND. J <= JBAR .AND. K >= 1 .AND. K <= KBAR) THEN
                IC1 = CELL_INDEXD(I,J,K)
                IC2 = CELL_INDEXD(I,J+1,K)
                SOLID1 = CELL_SOLID(IC1) .AND. .NOT. CELL_EXTERIOR(IC1)
                SOLID2 = CELL_SOLID(IC2) .AND. .NOT. CELL_EXTERIOR(IC2)
                IF (SOLID1 .AND. SOLID2) THEN
                    IF (PREDICTOR_FLAG == 1) THEN
                        DVVDT = -RFODT*V(I,J,K)
                    ELSE
                        DVVDT = -RFODT*(V(I,J,K)+VS(I,J,K))
                    ENDIF
                    FVY(I,J,K) = -RDYN(J)*(HP(I,J+1,K)-HP(I,J,K)) - DVVDT
                ENDIF
            ENDIF

            IF (I >= 1 .AND. I <= IBAR .AND. J >= 1 .AND. J <= JBAR .AND. K >= 0 .AND. K <= KBAR) THEN
                IC1 = CELL_INDEXD(I,J,K)
                IC2 = CELL_INDEXD(I,J,K+1)
                SOLID1 = CELL_SOLID(IC1) .AND. .NOT. CELL_EXTERIOR(IC1)
                SOLID2 = CELL_SOLID(IC2) .AND. .NOT. CELL_EXTERIOR(IC2)
                IF (SOLID1 .AND. SOLID2) THEN
                    IF (PREDICTOR_FLAG == 1) THEN
                        DWWDT = -RFODT*W(I,J,K)
                    ELSE
                        DWWDT = -RFODT*(W(I,J,K)+WS(I,J,K))
                    ENDIF
                    FVZ(I,J,K) = -RDZN(K)*(HP(I,J,K+1)-HP(I,J,K)) - DWWDT
                ENDIF
            ENDIF

        END SUBROUTINE NO_FLUX_SOLID_SOLID_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE NO_FLUX_WALL_GPU(FVX, FVY, FVZ, U, V, W, US, VS, WS, HP, WALL_TOPO, EXTERNAL_WALL_NOM, &
            WALL_U_NORMAL, WALL_U_NORMAL_S, RDXN, RDYN, RDZN, N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS, RFODT, &
            PRES_FLAG_LOCAL,PREDICTOR_FLAG,IBP1,JBP1)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,PRES_FLAG_LOCAL,PREDICTOR_FLAG,IBP1,JBP1
            REAL(EB), VALUE :: RFODT
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: FVX,FVY,FVZ,U,V,W,US,VS,WS,HP
            REAL(EB), DEVICE, DIMENSION(0:*) :: WALL_U_NORMAL,WALL_U_NORMAL_S,RDXN,RDYN,RDZN
            INTEGER , DEVICE, DIMENSION(1:*) :: EXTERNAL_WALL_NOM
            INTEGER :: IW,II,JJ,KK,IOR,BTYPE,NOM
            REAL(EB) :: UN,DUUDT,DVVDT,DWWDT,DHFCT,SGN_IOR

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS) RETURN

            BTYPE = WALL_TOPO(IW)%BOUNDARY_TYPE
            IF (BTYPE == OPEN_BOUNDARY .OR. BTYPE == INTERPOLATED_BOUNDARY) RETURN
            IF (IW <= N_EXTERNAL_WALL_CELLS) THEN
                NOM = EXTERNAL_WALL_NOM(IW)
            ELSE
                NOM = 0
            ENDIF
            IF (IW > N_EXTERNAL_WALL_CELLS .AND. BTYPE == NULL_BOUNDARY .AND. NOM == 0) RETURN

            II  = WALL_TOPO(IW)%II
            JJ  = WALL_TOPO(IW)%JJ
            KK  = WALL_TOPO(IW)%KK
            IOR = WALL_TOPO(IW)%IOR

            DHFCT = 1._EB
            SELECT CASE(PRES_FLAG_LOCAL)
                CASE(UGLMAT_FLAG,ULMAT_FLAG)
                    DHFCT = 0._EB
                CASE(GLMAT_FLAG)
                    IF (IW <= N_EXTERNAL_WALL_CELLS) DHFCT = 0._EB
            END SELECT

            IF (NOM /= 0 .OR. BTYPE == SOLID_BOUNDARY .OR. BTYPE == NULL_BOUNDARY) THEN
                SGN_IOR = SIGN(1._EB,REAL(IOR,EB))
                IF (PREDICTOR_FLAG == 1) THEN
                    UN = -SGN_IOR*WALL_U_NORMAL_S(IW)
                ELSE
                    UN = -SGN_IOR*WALL_U_NORMAL(IW)
                ENDIF

                SELECT CASE(IOR)
                    CASE( 1)
                        IF (PREDICTOR_FLAG == 1) THEN
                            DUUDT = RFODT*(UN-U(II,JJ,KK))
                        ELSE
                            DUUDT = 2._EB*RFODT*(UN-0.5_EB*(U(II,JJ,KK)+US(II,JJ,KK)))
                        ENDIF
                        FVX(II,JJ,KK) = -RDXN(II)*(HP(II+1,JJ,KK)-HP(II,JJ,KK))*DHFCT - DUUDT
                    CASE(-1)
                        IF (PREDICTOR_FLAG == 1) THEN
                            DUUDT = RFODT*(UN-U(II-1,JJ,KK))
                        ELSE
                            DUUDT = 2._EB*RFODT*(UN-0.5_EB*(U(II-1,JJ,KK)+US(II-1,JJ,KK)))
                        ENDIF
                        FVX(II-1,JJ,KK) = -RDXN(II-1)*(HP(II,JJ,KK)-HP(II-1,JJ,KK))*DHFCT - DUUDT
                    CASE( 2)
                        IF (PREDICTOR_FLAG == 1) THEN
                            DVVDT = RFODT*(UN-V(II,JJ,KK))
                        ELSE
                            DVVDT = 2._EB*RFODT*(UN-0.5_EB*(V(II,JJ,KK)+VS(II,JJ,KK)))
                        ENDIF
                        FVY(II,JJ,KK) = -RDYN(JJ)*(HP(II,JJ+1,KK)-HP(II,JJ,KK))*DHFCT - DVVDT
                    CASE(-2)
                        IF (PREDICTOR_FLAG == 1) THEN
                            DVVDT = RFODT*(UN-V(II,JJ-1,KK))
                        ELSE
                            DVVDT = 2._EB*RFODT*(UN-0.5_EB*(V(II,JJ-1,KK)+VS(II,JJ-1,KK)))
                        ENDIF
                        FVY(II,JJ-1,KK) = -RDYN(JJ-1)*(HP(II,JJ,KK)-HP(II,JJ-1,KK))*DHFCT - DVVDT
                    CASE( 3)
                        IF (PREDICTOR_FLAG == 1) THEN
                            DWWDT = RFODT*(UN-W(II,JJ,KK))
                        ELSE
                            DWWDT = 2._EB*RFODT*(UN-0.5_EB*(W(II,JJ,KK)+WS(II,JJ,KK)))
                        ENDIF
                        FVZ(II,JJ,KK) = -RDZN(KK)*(HP(II,JJ,KK+1)-HP(II,JJ,KK))*DHFCT - DWWDT
                    CASE(-3)
                        IF (PREDICTOR_FLAG == 1) THEN
                            DWWDT = RFODT*(UN-W(II,JJ,KK-1))
                        ELSE
                            DWWDT = 2._EB*RFODT*(UN-0.5_EB*(W(II,JJ,KK-1)+WS(II,JJ,KK-1)))
                        ENDIF
                        FVZ(II,JJ,KK-1) = -RDZN(KK-1)*(HP(II,JJ,KK)-HP(II,JJ,KK-1))*DHFCT - DWWDT
                END SELECT
            ENDIF

            IF (BTYPE == MIRROR_BOUNDARY) THEN
                SELECT CASE(IOR)
                    CASE( 1) ; FVX(II  ,JJ  ,KK  ) = 0._EB
                    CASE(-1) ; FVX(II-1,JJ  ,KK  ) = 0._EB
                    CASE( 2) ; FVY(II  ,JJ  ,KK  ) = 0._EB
                    CASE(-2) ; FVY(II  ,JJ-1,KK  ) = 0._EB
                    CASE( 3) ; FVZ(II  ,JJ  ,KK  ) = 0._EB
                    CASE(-3) ; FVZ(II  ,JJ  ,KK-1) = 0._EB
                END SELECT
            ENDIF
        END SUBROUTINE NO_FLUX_WALL_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE REMOVE_FB_TO_FA(IBP1,JBP1,KBP1,FVXD,FVYD,FVZD,FVX_BD,FVY_BD,FVZ_BD, PD , RHOPD , KRESD , RRHOD ,HPD )
            INTEGER , VALUE :: IBP1,JBP1,KBP1
            REAL(EB) , DEVICE , DIMENSION(0:IBP1,0:JBP1,0:*) :: FVXD,FVYD,FVZD,FVX_BD,FVY_BD,FVZ_BD , PD , RHOPD , KRESD , RRHOD , HPD

            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBP1) RETURN
            IF (J < 0 .OR. J > JBP1) RETURN
            IF (K < 0 .OR. K > KBP1) RETURN
            FVXD(I,J,K) = FVXD(I,J,K) - FVX_BD(I,J,K)
            FVYD(I,J,K) = FVYD(I,J,K) - FVY_BD(I,J,K)
            FVZD(I,J,K) = FVZD(I,J,K) - FVZ_BD(I,J,K)

            PD(I,J,K) = RHOPD(I,J,K)*(HPD(I,J,K)-KRESD(I,J,K))
            RRHOD(I,J,K) = 1._EB/RHOPD(I,J,K)
        END SUBROUTINE REMOVE_FB_TO_FA

        ATTRIBUTES(GLOBAL) SUBROUTINE COMPUTE_PD_RRHOD(IBP1,JBP1,KBP1,FVXD,FVYD,FVZD,FVX_BD,FVY_BD,FVZ_BD, PD , RHOPD , KRESD , RRHOD ,HPD )
            INTEGER , VALUE :: IBP1,JBP1,KBP1
            REAL(EB) , DEVICE , DIMENSION(0:IBP1,0:JBP1,0:*) :: FVXD,FVYD,FVZD,FVX_BD,FVY_BD,FVZ_BD , PD , RHOPD , KRESD , RRHOD , HPD

            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBP1) RETURN
            IF (J < 0 .OR. J > JBP1) RETURN
            IF (K < 0 .OR. K > KBP1) RETURN

            PD(I,J,K) = RHOPD(I,J,K)*(HPD(I,J,K)-KRESD(I,J,K))
            RRHOD(I,J,K) = 1._EB/RHOPD(I,J,K)
        END SUBROUTINE COMPUTE_PD_RRHOD

        ATTRIBUTES(GLOBAL) SUBROUTINE ADD_FB_TO_F_XYZ(PD,RHOPD,RRHOD,FVXD,FVYD,FVZD,FVX_BD,FVY_BD,FVZ_BD, &
            RDXND,RDYND,RDZND,IBAR,JBAR,KBAR,TWO_D_LOCAL)
            INTEGER , VALUE :: IBAR,JBAR,KBAR
            LOGICAL, VALUE :: TWO_D_LOCAL
            REAL(EB) , DEVICE , DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: PD,RHOPD,RRHOD
            REAL(EB) , DEVICE , DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: FVXD,FVYD,FVZD,FVX_BD,FVY_BD,FVZ_BD
            REAL(EB) , DEVICE , DIMENSION(0:*) :: RDXND,RDYND,RDZND
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I>=0 .AND. I<=IBAR .AND. J>=1 .AND. J<=JBAR .AND. K>=1 .AND. K<=KBAR) THEN
                FVX_BD(I,J,K) = -(PD(I,J,K)*RHOPD(I+1,J,K)+PD(I+1,J,K)*RHOPD(I,J,K))* &
                    (RRHOD(I+1,J,K)-RRHOD(I,J,K))*RDXND(I)/(RHOPD(I+1,J,K)+RHOPD(I,J,K))
                FVXD(I,J,K) = FVXD(I,J,K) + FVX_BD(I,J,K)
            ENDIF

            IF (.NOT.TWO_D_LOCAL .AND. I>=1 .AND. I<=IBAR .AND. J>=0 .AND. J<=JBAR .AND. K>=1 .AND. K<=KBAR) THEN
                FVY_BD(I,J,K) = -(PD(I,J,K)*RHOPD(I,J+1,K)+PD(I,J+1,K)*RHOPD(I,J,K))* &
                    (RRHOD(I,J+1,K)-RRHOD(I,J,K))*RDYND(J)/(RHOPD(I,J+1,K)+RHOPD(I,J,K))
                FVYD(I,J,K) = FVYD(I,J,K) + FVY_BD(I,J,K)
            ENDIF

            IF (I>=1 .AND. I<=IBAR .AND. J>=1 .AND. J<=JBAR .AND. K>=0 .AND. K<=KBAR) THEN
                FVZ_BD(I,J,K) = -(PD(I,J,K)*RHOPD(I,J,K+1)+PD(I,J,K+1)*RHOPD(I,J,K))* &
                    (RRHOD(I,J,K+1)-RRHOD(I,J,K))*RDZND(K)/(RHOPD(I,J,K+1)+RHOPD(I,J,K))
                FVZD(I,J,K) = FVZD(I,J,K) + FVZ_BD(I,J,K)
            ENDIF
        END SUBROUTINE ADD_FB_TO_F_XYZ

        ATTRIBUTES(GLOBAL) SUBROUTINE COMPUTE_PRESSURE_BCS_GPU(BXS_D, BXF_D, BYS_D, BYF_D, BZS_D, BZF_D, WALL_TOPO, &
            WALL_PRESSURE_BC_TYPE, WALL_IIO, WALL_JJO, WALL_KKO, WALL_DUNDT, WALL_WORK1, WALL_P_EXTERNAL_OVER_RHO, FVX, &
            FVY, FVZ, UU, VV, WW, HP, KRES, HX, HY, HZ, DX, DY, DZ, N_EXTERNAL_WALL_CELLS, IBAR, JBAR, KBAR, IBP1, JBP1, &
            KBP1, H0_PRESSURE)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)

            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1
            REAL(EB), VALUE :: H0_PRESSURE
            REAL(EB), DEVICE, DIMENSION(1:MERGE(JBP1,1,JBAR>1),1:*) :: BXS_D,BXF_D
            REAL(EB), DEVICE, DIMENSION(1:IBP1,1:*) :: BYS_D,BYF_D,BZS_D,BZF_D
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: FVX,FVY,FVZ,UU,VV,WW,HP,KRES
            REAL(EB), DEVICE, DIMENSION(0:*) :: WALL_DUNDT,WALL_WORK1,WALL_P_EXTERNAL_OVER_RHO,HX,HY,HZ,DX,DY,DZ
            INTEGER, DEVICE, DIMENSION(0:*) :: WALL_PRESSURE_BC_TYPE
            INTEGER, DEVICE, DIMENSION(1:*) :: WALL_IIO,WALL_JJO,WALL_KKO

            INTEGER :: IW,I,J,K,IOR,BTYPE,PBC,IIO,JJO,KKO
            REAL(EB) :: P_EXTERNAL_OVER_RHO,DX_OTHER,DY_OTHER,DZ_OTHER

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS) RETURN

            I = WALL_TOPO(IW)%II
            J = WALL_TOPO(IW)%JJ
            K = WALL_TOPO(IW)%KK
            IOR = WALL_TOPO(IW)%IOR
            BTYPE = WALL_TOPO(IW)%BOUNDARY_TYPE
            PBC = WALL_PRESSURE_BC_TYPE(IW)

            IF (PBC == NEUMANN) THEN
                SELECT CASE(IOR)
                    CASE( 1) ; BXS_D(J,K) = HX(0)   *(-FVX(0,J,K)    + WALL_DUNDT(IW))
                    CASE(-1) ; BXF_D(J,K) = HX(IBP1)*(-FVX(IBAR,J,K) - WALL_DUNDT(IW))
                    CASE( 2) ; BYS_D(I,K) = HY(0)   *(-FVY(I,0,K)    + WALL_DUNDT(IW))
                    CASE(-2) ; BYF_D(I,K) = HY(JBP1)*(-FVY(I,JBAR,K) - WALL_DUNDT(IW))
                    CASE( 3) ; BZS_D(I,J) = HZ(0)   *(-FVZ(I,J,0)    + WALL_DUNDT(IW))
                    CASE(-3) ; BZF_D(I,J) = HZ(KBP1)*(-FVZ(I,J,KBAR) - WALL_DUNDT(IW))
                END SELECT
                RETURN
            ENDIF

            IF (PBC /= DIRICHLET) RETURN

            IF (BTYPE == INTERPOLATED_BOUNDARY) RETURN

            IF (BTYPE == OPEN_BOUNDARY) THEN
                P_EXTERNAL_OVER_RHO = WALL_P_EXTERNAL_OVER_RHO(IW)
                SELECT CASE(IOR)
                    CASE( 1)
                        IF (UU(0,J,K)<0._EB) THEN
                            BXS_D(J,K) = P_EXTERNAL_OVER_RHO + KRES(1,J,K)
                        ELSE
                            BXS_D(J,K) = P_EXTERNAL_OVER_RHO + H0_PRESSURE
                        ENDIF
                    CASE(-1)
                        IF (UU(IBAR,J,K)>0._EB) THEN
                            BXF_D(J,K) = P_EXTERNAL_OVER_RHO + KRES(IBAR,J,K)
                        ELSE
                            BXF_D(J,K) = P_EXTERNAL_OVER_RHO + H0_PRESSURE
                        ENDIF
                    CASE( 2)
                        IF (VV(I,0,K)<0._EB) THEN
                            BYS_D(I,K) = P_EXTERNAL_OVER_RHO + KRES(I,1,K)
                        ELSE
                            BYS_D(I,K) = P_EXTERNAL_OVER_RHO + H0_PRESSURE
                        ENDIF
                    CASE(-2)
                        IF (VV(I,JBAR,K)>0._EB) THEN
                            BYF_D(I,K) = P_EXTERNAL_OVER_RHO + KRES(I,JBAR,K)
                        ELSE
                            BYF_D(I,K) = P_EXTERNAL_OVER_RHO + H0_PRESSURE
                        ENDIF
                    CASE( 3)
                        IF (WW(I,J,0)<0._EB) THEN
                            BZS_D(I,J) = P_EXTERNAL_OVER_RHO + KRES(I,J,1)
                        ELSE
                            BZS_D(I,J) = P_EXTERNAL_OVER_RHO + H0_PRESSURE
                        ENDIF
                    CASE(-3)
                        IF (WW(I,J,KBAR)>0._EB) THEN
                            BZF_D(I,J) = P_EXTERNAL_OVER_RHO + KRES(I,J,KBAR)
                        ELSE
                            BZF_D(I,J) = P_EXTERNAL_OVER_RHO + H0_PRESSURE
                        ENDIF
                END SELECT
                RETURN
            ENDIF

            SELECT CASE(IOR)
                CASE( 1) ; BXS_D(J,K) = 0.5_EB*(HP(0,J,K)   + HP(1,J,K))    + WALL_WORK1(IW)
                CASE(-1) ; BXF_D(J,K) = 0.5_EB*(HP(IBAR,J,K)+ HP(IBP1,J,K)) + WALL_WORK1(IW)
                CASE( 2) ; BYS_D(I,K) = 0.5_EB*(HP(I,0,K)   + HP(I,1,K))    + WALL_WORK1(IW)
                CASE(-2) ; BYF_D(I,K) = 0.5_EB*(HP(I,JBAR,K)+ HP(I,JBP1,K)) + WALL_WORK1(IW)
                CASE( 3) ; BZS_D(I,J) = 0.5_EB*(HP(I,J,0)   + HP(I,J,1))    + WALL_WORK1(IW)
                CASE(-3) ; BZF_D(I,J) = 0.5_EB*(HP(I,J,KBAR)+ HP(I,J,KBP1)) + WALL_WORK1(IW)
            END SELECT

        END SUBROUTINE COMPUTE_PRESSURE_BCS_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE COMPUTE_PRESSURE_BCS_INTERPOLATED_GPU(BXS_D,BXF_D,BYS_D,BYF_D,BZS_D,BZF_D,WALL_TOPO, &
            WALL_PRESSURE_BC_TYPE,EXTERNAL_WALL_NOM,WALL_IIO,WALL_JJO,WALL_KKO,WALL_WORK1,HP,DX,DY,DZ, &
            DX_OTHER,DY_OTHER,DZ_OTHER,N_EXTERNAL_WALL_CELLS,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,NOM_FILTER,WALL_LB, &
            BX_N1,BX_LB1,BX_LB2,BY_N1,BY_LB1,BY_LB2,BZ_N1,BZ_LB1,BZ_LB2)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE, DIMENSION(*) :: WALL_TOPO

            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,NOM_FILTER,WALL_LB
            INTEGER, VALUE :: BX_N1,BX_LB1,BX_LB2,BY_N1,BY_LB1,BY_LB2,BZ_N1,BZ_LB1,BZ_LB2
            REAL(EB), DEVICE, DIMENSION(*) :: BXS_D,BXF_D,BYS_D,BYF_D,BZS_D,BZF_D
            REAL(EB), DEVICE :: HP(0:IBP1,0:JBP1,0:KBP1)
            REAL(EB), DEVICE, DIMENSION(*) :: WALL_WORK1,DX,DY,DZ,DX_OTHER,DY_OTHER,DZ_OTHER
            INTEGER, DEVICE, DIMENSION(*) :: WALL_PRESSURE_BC_TYPE
            INTEGER, DEVICE, DIMENSION(*) :: EXTERNAL_WALL_NOM,WALL_IIO,WALL_JJO,WALL_KKO

            INTEGER :: IW,IW_IDX,I,J,K,IOR,IIO,JJO,KKO,B_IDX
            REAL(EB) :: D_OTHER

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS) RETURN
            IF (EXTERNAL_WALL_NOM(IW) /= NOM_FILTER) RETURN
            IW_IDX = 1 + IW - WALL_LB
            IF (WALL_TOPO(IW_IDX)%BOUNDARY_TYPE /= INTERPOLATED_BOUNDARY) RETURN
            IF (WALL_PRESSURE_BC_TYPE(IW_IDX) /= DIRICHLET) RETURN

            I = WALL_TOPO(IW_IDX)%II
            J = WALL_TOPO(IW_IDX)%JJ
            K = WALL_TOPO(IW_IDX)%KK
            IOR = WALL_TOPO(IW_IDX)%IOR
            IIO = WALL_IIO(IW)
            JJO = WALL_JJO(IW)
            KKO = WALL_KKO(IW)

            SELECT CASE(IOR)
                CASE( 1)
                    D_OTHER = DX_OTHER(IIO+1)
                    B_IDX = 1 + (J-BX_LB1) + (K-BX_LB2)*BX_N1
                    BXS_D(B_IDX) = (D_OTHER*HP(1,J,K)+DX(2)*HP(0,J,K))/(DX(2)+D_OTHER)+WALL_WORK1(IW_IDX)
                CASE(-1)
                    D_OTHER = DX_OTHER(IIO+1)
                    B_IDX = 1 + (J-BX_LB1) + (K-BX_LB2)*BX_N1
                    BXF_D(B_IDX) = (D_OTHER*HP(IBAR,J,K)+DX(IBAR+1)*HP(IBP1,J,K))/(DX(IBAR+1)+D_OTHER)+WALL_WORK1(IW_IDX)
                CASE( 2)
                    D_OTHER = DY_OTHER(JJO+1)
                    B_IDX = 1 + (I-BY_LB1) + (K-BY_LB2)*BY_N1
                    BYS_D(B_IDX) = (D_OTHER*HP(I,1,K)+DY(2)*HP(I,0,K))/(DY(2)+D_OTHER)+WALL_WORK1(IW_IDX)
                CASE(-2)
                    D_OTHER = DY_OTHER(JJO+1)
                    B_IDX = 1 + (I-BY_LB1) + (K-BY_LB2)*BY_N1
                    BYF_D(B_IDX) = (D_OTHER*HP(I,JBAR,K)+DY(JBAR+1)*HP(I,JBP1,K))/(DY(JBAR+1)+D_OTHER)+WALL_WORK1(IW_IDX)
                CASE( 3)
                    D_OTHER = DZ_OTHER(KKO+1)
                    B_IDX = 1 + (I-BZ_LB1) + (J-BZ_LB2)*BZ_N1
                    BZS_D(B_IDX) = (D_OTHER*HP(I,J,1)+DZ(2)*HP(I,J,0))/(DZ(2)+D_OTHER)+WALL_WORK1(IW_IDX)
                CASE(-3)
                    D_OTHER = DZ_OTHER(KKO+1)
                    B_IDX = 1 + (I-BZ_LB1) + (J-BZ_LB2)*BZ_N1
                    BZF_D(B_IDX) = (D_OTHER*HP(I,J,KBAR)+DZ(KBAR+1)*HP(I,J,KBP1))/(DZ(KBAR+1)+D_OTHER)+WALL_WORK1(IW_IDX)
            END SELECT

        END SUBROUTINE COMPUTE_PRESSURE_BCS_INTERPOLATED_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE COMPUTE_PRHS_GPU(PRHSD,FVXD,FVYD,FVZD,DDDTD,RDXD,RDYD,RDZD,IBAR,JBAR,KBAR)
            INTEGER , VALUE :: IBAR,JBAR,KBAR
            REAL(EB) , DEVICE , DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: FVXD, FVYD , FVZD , DDDTD
            REAL(EB) , DEVICE , DIMENSION(1:IBAR+1,1:JBAR+1,1:*) :: PRHSD
            REAL(EB) , DEVICE , DIMENSION(0:*) ::  RDXD , RDYD , RDZD
            REAL(EB) :: TRM1 , TRM2 , TRM3 , TRM4
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            TRM1 = (FVXD(I-1,J,K)-FVXD(I,J,K))*RDXD(I)
            TRM2 = (FVYD(I,J-1,K)-FVYD(I,J,K))*RDYD(J)
            TRM3 = (FVZD(I,J,K-1)-FVZD(I,J,K))*RDZD(K)
            TRM4 = -DDDTD(I,J,K)
            PRHSD(I,J,K) = TRM1 + TRM2 + TRM3 + TRM4
        END SUBROUTINE COMPUTE_PRHS_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE P_CALCULATE_KERNEL(P,RHOP,HP,KRES, IBP1,JBP1,KBP1)
            INTEGER , VALUE :: IBP1,JBP1,KBP1
            REAL(EB) , DEVICE , DIMENSION(0:IBP1,0:JBP1,0:*) :: P,RHOP,HP,KRES

            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBP1) RETURN
            IF (J < 0 .OR. J > JBP1) RETURN
            IF (K < 0 .OR. K > KBP1) RETURN

            P(I,J,K) = RHOP(I,J,K)*(HP(I,J,K)-KRES(I,J,K))


        END SUBROUTINE P_CALCULATE_KERNEL

        ATTRIBUTES(GLOBAL) SUBROUTINE RESIDUALS_KERNEL(RESIDUAL,FVX,FVX_B,FVY,FVY_B,FVZ,FVZ_B,DDDT,P,KRES,RRHO_X,RRHO_Y,RRHO_Z, &
                                                       RDXN,RDYN,RDZN,RDX,RDY,RDZ,R,RRN,IBAR,JBAR,KBAR)
            IMPLICIT NONE

            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: RESIDUAL,FVX,FVX_B,FVY,FVY_B,FVZ,FVZ_B,DDDT,P,KRES,RRHO_X,RRHO_Y,RRHO_Z
            REAL(EB), DEVICE, DIMENSION(0:*) :: RDXN,RDYN,RDZN,RDX,RDY,RDZ,R,RRN

            INTEGER :: GID, NCELL, TMP, REM, I, J, K
            REAL(EB) :: RHSS, LHSS

            GID = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            NCELL = IBAR*JBAR*KBAR
            IF (GID > NCELL) RETURN

            TMP = GID - 1
            K = TMP / (IBAR*JBAR) + 1
            REM = MOD(TMP, IBAR*JBAR)
            J = REM / IBAR + 1
            I = MOD(REM, IBAR) + 1

            RHSS = (R(I-1)*(FVX(I-1,J,K)-FVX_B(I-1,J,K)) - R(I)*(FVX(I,J,K)-FVX_B(I,J,K))) * RDX(I)*RRN(I) &
                 + (((FVY(I,J-1,K)-FVY_B(I,J-1,K)) - (FVY(I,J,K)-FVY_B(I,J,K))) * RDY(J)) &
                 + (((FVZ(I,J,K-1)-FVZ_B(I,J,K-1)) - (FVZ(I,J,K)-FVZ_B(I,J,K))) * RDZ(K)) &
                 - DDDT(I,J,K)

            LHSS = ((P(I+1,J,K)-P(I,J,K))*RDXN(I)*R(I)       *RRHO_X(I,J,K) - &
                    (P(I,J,K)-P(I-1,J,K))*RDXN(I-1)*R(I-1)   *RRHO_X(I-1,J,K))*RDX(I)*RRN(I) &
                 + ((P(I,J+1,K)-P(I,J,K))*RDYN(J)   *RRHO_Y(I,J,K) - &
                    (P(I,J,K)-P(I,J-1,K))*RDYN(J-1) *RRHO_Y(I,J-1,K))*RDY(J) &
                 + ((P(I,J,K+1)-P(I,J,K))*RDZN(K)   *RRHO_Z(I,J,K) - &
                    (P(I,J,K)-P(I,J,K-1))*RDZN(K-1) *RRHO_Z(I,J,K-1))*RDZ(K) &
                 + ((KRES(I+1,J,K)-KRES(I,J,K))*RDXN(I)*R(I) - (KRES(I,J,K)-KRES(I-1,J,K))*RDXN(I-1)*R(I-1))*RDX(I)*RRN(I) &
                 + ((KRES(I,J+1,K)-KRES(I,J,K))*RDYN(J)   - (KRES(I,J,K)-KRES(I,J-1,K))*RDYN(J-1))*RDY(J) &
                 + ((KRES(I,J,K+1)-KRES(I,J,K))*RDZN(K)   - (KRES(I,J,K)-KRES(I,J,K-1))*RDZN(K-1))*RDZ(K)

            RESIDUAL(I,J,K) = ABS(RHSS-LHSS)
        END SUBROUTINE RESIDUALS_KERNEL


        ATTRIBUTES(GLOBAL) SUBROUTINE RESIDUALS_BLOCKMAX_KERNEL(BLOCKMAX,BLOCKI,BLOCKJ,BLOCKK,RESIDUAL,IBAR,JBAR,KBAR)
            IMPLICIT NONE

            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), DEVICE, DIMENSION(*) :: BLOCKMAX
            INTEGER, DEVICE, DIMENSION(*) :: BLOCKI, BLOCKJ, BLOCKK
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: RESIDUAL

            REAL(EB), SHARED :: SMAX(256)
            INTEGER, SHARED :: SI(256), SJ(256), SK(256)

            INTEGER :: TID, GID, NCELL, TMP, REM, I, J, K, OFFSET
            REAL(EB) :: VAL

            TID = THREADIDX%X
            GID = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            NCELL = IBAR*JBAR*KBAR

            VAL = -1._EB
            I = 0
            J = 0
            K = 0
            IF (GID <= NCELL) THEN
                TMP = GID - 1
                K = TMP / (IBAR*JBAR) + 1
                REM = MOD(TMP, IBAR*JBAR)
                J = REM / IBAR + 1
                I = MOD(REM, IBAR) + 1
                VAL = RESIDUAL(I,J,K)
            ENDIF

            SMAX(TID) = VAL
            SI(TID) = I
            SJ(TID) = J
            SK(TID) = K
            CALL SYNCTHREADS()

            OFFSET = BLOCKDIM%X/2
            DO WHILE (OFFSET >= 1)
                IF (TID <= OFFSET) THEN
                    IF (SMAX(TID+OFFSET) > SMAX(TID)) THEN
                        SMAX(TID) = SMAX(TID+OFFSET)
                        SI(TID) = SI(TID+OFFSET)
                        SJ(TID) = SJ(TID+OFFSET)
                        SK(TID) = SK(TID+OFFSET)
                    ENDIF
                ENDIF
                CALL SYNCTHREADS()
                OFFSET = OFFSET/2
            ENDDO

            IF (TID == 1) THEN
                BLOCKMAX(BLOCKIDX%X) = SMAX(1)
                BLOCKI(BLOCKIDX%X) = SI(1)
                BLOCKJ(BLOCKIDX%X) = SJ(1)
                BLOCKK(BLOCKIDX%X) = SK(1)
            ENDIF

        END SUBROUTINE RESIDUALS_BLOCKMAX_KERNEL

END MODULE G_PRES
