!------------------------------------------------------------------------------
! Module: G_RADI
! Purpose: GPU DOM radiation setup, transport sweep, and radiation-source assembly.
!------------------------------------------------------------------------------
MODULE G_RADI
    USE G_MESH_POINTERS
    USE G_CONS, ONLY: GPU_DEBUG_H2D_D2H
    USE GLOBAL_CONSTANTS
    USE MESH_POINTERS
    USE COMP_FUNCTIONS
    USE G_FUNC, ONLY: EVALUATE_RAMP_DEVICE
    USE G_ERRO, ONLY: GPU_WAIT_ALL,GPU_MPI_IRECV,GPU_MPI_ISEND,GPU_MPI_WAITALL
    USE RADCONS
    IMPLICIT NONE

    CONTAINS

        !--------------------------------------------------------------------------
        ! Radiation-state lifecycle and initialization
        !--------------------------------------------------------------------------
        SUBROUTINE RELEASE_GPU_RADIATION(NM)
            INTEGER, INTENT(IN) :: NM

            GM => GPU_MESHES(NM)
            IF (.NOT.ALLOCATED(GM%RADI)) RETURN
            DEALLOCATE(GM%RADI)
        END SUBROUTINE RELEASE_GPU_RADIATION

        SUBROUTINE GPU_INIT_VALUES_RADIATION(NM)
            INTEGER, INTENT(IN) :: NM
            INTEGER :: IW,NOM
            INTEGER :: SID,N_SLICE,I,J,K,JSTART,JEND,IDX,NCELL
            INTEGER, ALLOCATABLE :: WALL_NIC_MIN_H(:),WALL_NIC_MAX_H(:),IJK_SLICE_MAP_H(:,:)
            TYPE(OMESH_TYPE), POINTER :: OM

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)
            IF (ALLOCATED(GM%RADI)) DEALLOCATE(GM%RADI)
            ALLOCATE(GM%RADI)

            ALLOCATE(GM%RADI%OMESH(1:NMESHES))
            DO NOM=1,NMESHES
                OM => MESHES(NM)%OMESH(NOM)
                IF (ALLOCATED(OM%IL_S)) THEN
                    ALLOCATE(GM%RADI%OMESH(NOM)%IL_S_D(LBOUND(OM%IL_S,1):UBOUND(OM%IL_S,1), &
                        LBOUND(OM%IL_S,2):UBOUND(OM%IL_S,2),LBOUND(OM%IL_S,3):UBOUND(OM%IL_S,3)))
                    GM%RADI%OMESH(NOM)%IL_S_D = OM%IL_S
                ENDIF
                IF (ALLOCATED(OM%IL_R)) THEN
                    ALLOCATE(GM%RADI%OMESH(NOM)%IL_R_D(LBOUND(OM%IL_R,1):UBOUND(OM%IL_R,1), &
                        LBOUND(OM%IL_R,2):UBOUND(OM%IL_R,2),LBOUND(OM%IL_R,3):UBOUND(OM%IL_R,3)))
                    GM%RADI%OMESH(NOM)%IL_R_D = OM%IL_R
                ENDIF
            ENDDO

            GM%RADI%NSLICES = IBAR + JBAR + KBAR - 2
            NCELL = IBAR*JBAR*KBAR
            ALLOCATE(GM%RADI%SLICE_OFFSET_H(1:GM%RADI%NSLICES+1))
            ALLOCATE(IJK_SLICE_MAP_H(1:NCELL,1:3))

            IDX = 1
            DO SID=1,GM%RADI%NSLICES
                GM%RADI%SLICE_OFFSET_H(SID) = IDX
                N_SLICE = SID + 2
                DO K=1,KBAR
                    JSTART = MAX(1,N_SLICE-K-IBAR)
                    JEND = MIN(JBAR,N_SLICE-K-1)
                    IF (JSTART > JEND) CYCLE
                    DO J=JSTART,JEND
                        I = N_SLICE-J-K
                        IJK_SLICE_MAP_H(IDX,1) = I
                        IJK_SLICE_MAP_H(IDX,2) = J
                        IJK_SLICE_MAP_H(IDX,3) = K
                        IDX = IDX + 1
                    ENDDO
                ENDDO
            ENDDO
            GM%RADI%SLICE_OFFSET_H(GM%RADI%NSLICES+1) = IDX

            ALLOCATE(GM%RADI%IJK_SLICE_MAP_D(1:NCELL,1:3))
            ALLOCATE(GM%RADI%KAPPA_GAS_D(0:IBP1,0:JBP1,0:KBP1))
            ALLOCATE(GM%RADI%UII_D(0:IBP1,0:JBP1,0:KBP1))
            ALLOCATE(GM%RADI%UIID_D(0:IBP1,0:JBP1,0:KBP1,1:UIIDIM))
            ALLOCATE(GM%RADI%QR_W_D(0:IBP1,0:JBP1,0:KBP1))
            ALLOCATE(GM%RADI%EMISSION_D(0:IBP1,0:JBP1,0:KBP1))
            ALLOCATE(GM%RADI%ABSORPTION_D(0:IBP1,0:JBP1,0:KBP1))
            ALLOCATE(GM%RADI%RAD_Q_SUM_D(1),GM%RADI%KFST4_SUM_D(1))
            GM%RADI%IJK_SLICE_MAP_D = IJK_SLICE_MAP_H
            DEALLOCATE(IJK_SLICE_MAP_H)
            GM%RADI%KAPPA_GAS_D = 0._EB
            GM%RADI%UII_D = 0._EB
            GM%RADI%UIID_D = UIID
            GM%RADI%QR_W_D = 0._EB
            IF (N_LP_ARRAY_INDICES>0 .AND. ASSOCIATED(QR_W)) GM%RADI%QR_W_D = QR_W
            GM%RADI%EMISSION_D = 0._EB
            GM%RADI%ABSORPTION_D = 0._EB
            IF (STORE_RADIATION_TERMS) THEN
                GM%RADI%EMISSION_D = RADIATION_EMISSION
                GM%RADI%ABSORPTION_D = RADIATION_ABSORPTION
            ENDIF
            GM%RADI%RAD_Q_SUM_D = 0._EB
            GM%RADI%KFST4_SUM_D = 0._EB

            ALLOCATE(GM%RADI%WALL_BR_INDEX_D(GM%WALL_LB:GM%WALL_UB))
            ALLOCATE(GM%RADI%WALL_NIC_MIN_D(1:MAX(1,N_EXTERNAL_WALL_CELLS)))
            ALLOCATE(GM%RADI%WALL_NIC_MAX_D(1:MAX(1,N_EXTERNAL_WALL_CELLS)))
            ALLOCATE(GM%RADI%WALL_WORK1_D(LBOUND(WALL_WORK1,1):UBOUND(WALL_WORK1,1)))
            ALLOCATE(GM%RADI%WALL_WORK2_D(LBOUND(WALL_WORK2,1):UBOUND(WALL_WORK2,1)))
            ALLOCATE(GM%RADI%CELL_ILW_D(LBOUND(CELL_ILW,1):UBOUND(CELL_ILW,1),LBOUND(CELL_ILW,2):UBOUND(CELL_ILW,2)))
            ALLOCATE(GM%RADI%BOUNDARY_RADIA_ILW_D(1:MAX(1,N_BOUNDARY_RADIA_DIM),1:NUMBER_SPECTRAL_BANDS, &
                                                   1:NUMBER_RADIATION_ANGLES))
            ALLOCATE(GM%RADI%WALL_BR_INDEX_H(GM%WALL_LB:GM%WALL_UB))
            ALLOCATE(GM%RADI%WALL_Q_RAD_IN_H(GM%WALL_LB:GM%WALL_UB))
            ALLOCATE(GM%RADI%WALL_Q_RAD_OUT_H(GM%WALL_LB:GM%WALL_UB))
            ALLOCATE(WALL_NIC_MIN_H(1:MAX(1,N_EXTERNAL_WALL_CELLS)))
            ALLOCATE(WALL_NIC_MAX_H(1:MAX(1,N_EXTERNAL_WALL_CELLS)))
            ALLOCATE(GM%RADI%BOUNDARY_RADIA_ILW_H(1:MAX(1,N_BOUNDARY_RADIA_DIM),1:NUMBER_SPECTRAL_BANDS, &
                                                   1:NUMBER_RADIATION_ANGLES))
            GM%RADI%WALL_BR_INDEX_H = 0
            GM%RADI%WALL_Q_RAD_IN_H = 0._EB
            GM%RADI%WALL_Q_RAD_OUT_H = 0._EB
            WALL_NIC_MIN_H = 0
            WALL_NIC_MAX_H = -1
            DO IW=1,GM%WALL_UB
                GM%RADI%WALL_BR_INDEX_H(IW) = WALL(IW)%BR_INDEX
                IF (WALL(IW)%B1_INDEX>0) THEN
                    GM%RADI%WALL_Q_RAD_IN_H(IW) = BOUNDARY_PROP1(WALL(IW)%B1_INDEX)%Q_RAD_IN
                    GM%RADI%WALL_Q_RAD_OUT_H(IW) = BOUNDARY_PROP1(WALL(IW)%B1_INDEX)%Q_RAD_OUT
                ENDIF
            ENDDO
            DO IW=1,N_EXTERNAL_WALL_CELLS
                WALL_NIC_MIN_H(IW) = EXTERNAL_WALL(IW)%NIC_MIN
                WALL_NIC_MAX_H(IW) = EXTERNAL_WALL(IW)%NIC_MAX
            ENDDO
            CALL RADI_PACK_BOUNDARY_RADIA_ILW_HOST
            GM%RADI%WALL_BR_INDEX_D = GM%RADI%WALL_BR_INDEX_H
            GM%RADI%WALL_NIC_MIN_D = WALL_NIC_MIN_H
            GM%RADI%WALL_NIC_MAX_D = WALL_NIC_MAX_H
            GM%WALL_Q_RAD_IN = GM%RADI%WALL_Q_RAD_IN_H
            GM%WALL_Q_RAD_OUT = GM%RADI%WALL_Q_RAD_OUT_H
            GM%RADI%WALL_WORK1_D = WALL_WORK1
            GM%RADI%WALL_WORK2_D = WALL_WORK2
            IF (.NOT.ALLOCATED(DLN_D)) THEN
                ALLOCATE(DLN_D(-3:3,1:NUMBER_RADIATION_ANGLES))
                ALLOCATE(DLX_D(1:NUMBER_RADIATION_ANGLES),DLY_D(1:NUMBER_RADIATION_ANGLES))
                ALLOCATE(DLZ_D(1:NUMBER_RADIATION_ANGLES),RSA_D(1:NUMBER_RADIATION_ANGLES))
                ALLOCATE(DLM_D(1:NUMBER_RADIATION_ANGLES,1:3))
            ENDIF
            DLN_D = DLN
            DLX_D = DLX
            DLY_D = DLY
            DLZ_D = DLZ
            RSA_D = RSA
            DLM_D = DLM
            GM%RADI%CELL_ILW_D = CELL_ILW
            GM%RADI%BOUNDARY_RADIA_ILW_D = GM%RADI%BOUNDARY_RADIA_ILW_H
            DEALLOCATE(WALL_NIC_MIN_H,WALL_NIC_MAX_H)

            IF (N_RADCAL_ARRAY_SIZE>0 .AND. ALLOCATED(Z2RADCAL_SPECIES) .AND. ALLOCATED(RADCAL_SPECIES2KAPPA) .AND. &
                .NOT.ALLOCATED(RADI_RADCAL_SPECIES_INDEX_D)) THEN
                ALLOCATE(RADI_RADCAL_SPECIES_INDEX_D(1:N_RADCAL_ARRAY_SIZE))
                ALLOCATE(RADI_Z2RADCAL_SPECIES_D(1:N_RADCAL_ARRAY_SIZE,1:N_TRACKED_SPECIES))
                ALLOCATE(RADI_RADCAL_SPECIES2KAPPA_D(1:N_RADCAL_ARRAY_SIZE,0:N_KAPPA_Y,0:N_KAPPA_T,1:NUMBER_SPECTRAL_BANDS))
                RADI_RADCAL_SPECIES_INDEX_D = RADCAL_SPECIES_INDEX(1:N_RADCAL_ARRAY_SIZE)
                RADI_Z2RADCAL_SPECIES_D = Z2RADCAL_SPECIES
                RADI_RADCAL_SPECIES2KAPPA_D = RADCAL_SPECIES2KAPPA
            ENDIF

        END SUBROUTINE GPU_INIT_VALUES_RADIATION

        !--------------------------------------------------------------------------
        ! Radiation H2D/D2H synchronization helpers
        !--------------------------------------------------------------------------
        SUBROUTINE RADI_SYNC_WALL_RADIATION_H2D
            INTEGER :: IW

            GM%RADI%WALL_Q_RAD_IN_H  = 0._EB
            GM%RADI%WALL_Q_RAD_OUT_H = 0._EB
            DO IW=1,GM%WALL_UB
                IF (WALL(IW)%B1_INDEX<=0) CYCLE
                GM%RADI%WALL_Q_RAD_IN_H(IW)  = BOUNDARY_PROP1(WALL(IW)%B1_INDEX)%Q_RAD_IN
                GM%RADI%WALL_Q_RAD_OUT_H(IW) = BOUNDARY_PROP1(WALL(IW)%B1_INDEX)%Q_RAD_OUT
            ENDDO
            GM%WALL_Q_RAD_IN  = GM%RADI%WALL_Q_RAD_IN_H
            GM%WALL_Q_RAD_OUT = GM%RADI%WALL_Q_RAD_OUT_H

        END SUBROUTINE RADI_SYNC_WALL_RADIATION_H2D

        SUBROUTINE RADI_SYNC_WALL_RADIATION_D2H
            INTEGER :: IW

            GM%RADI%WALL_Q_RAD_IN_H  = GM%WALL_Q_RAD_IN
            GM%RADI%WALL_Q_RAD_OUT_H = GM%WALL_Q_RAD_OUT
            DO IW=1,GM%WALL_UB
                IF (WALL(IW)%B1_INDEX<=0) CYCLE
                BOUNDARY_PROP1(WALL(IW)%B1_INDEX)%Q_RAD_IN  = GM%RADI%WALL_Q_RAD_IN_H(IW)
                BOUNDARY_PROP1(WALL(IW)%B1_INDEX)%Q_RAD_OUT = GM%RADI%WALL_Q_RAD_OUT_H(IW)
            ENDDO

        END SUBROUTINE RADI_SYNC_WALL_RADIATION_D2H

        SUBROUTINE RADI_PACK_BOUNDARY_RADIA_ILW_HOST
            INTEGER :: IBR,IB,N,NLIM

            IF (.NOT.ALLOCATED(GM%RADI%BOUNDARY_RADIA_ILW_H)) RETURN

            GM%RADI%BOUNDARY_RADIA_ILW_H = 0._EB
            DO IBR=1,N_BOUNDARY_RADIA_DIM
                IF (.NOT.ALLOCATED(BOUNDARY_RADIA(IBR)%BAND)) CYCLE
                DO IB=1,NUMBER_SPECTRAL_BANDS
                    IF (.NOT.ALLOCATED(BOUNDARY_RADIA(IBR)%BAND(IB)%ILW)) CYCLE
                    NLIM = MIN(NUMBER_RADIATION_ANGLES,SIZE(BOUNDARY_RADIA(IBR)%BAND(IB)%ILW))
                    DO N=1,NLIM
                        GM%RADI%BOUNDARY_RADIA_ILW_H(IBR,IB,N) = BOUNDARY_RADIA(IBR)%BAND(IB)%ILW(N)
                    ENDDO
                ENDDO
            ENDDO

        END SUBROUTINE RADI_PACK_BOUNDARY_RADIA_ILW_HOST

        SUBROUTINE RADI_SYNC_BOUNDARY_RADIA_ILW_H2D

            CALL RADI_PACK_BOUNDARY_RADIA_ILW_HOST
            GM%RADI%BOUNDARY_RADIA_ILW_D = GM%RADI%BOUNDARY_RADIA_ILW_H

        END SUBROUTINE RADI_SYNC_BOUNDARY_RADIA_ILW_H2D

        SUBROUTINE RADI_SYNC_BOUNDARY_RADIA_ILW_D2H
            INTEGER :: IBR,IB,N,NLIM

            IF (.NOT.ALLOCATED(GM%RADI%BOUNDARY_RADIA_ILW_H)) RETURN

            GM%RADI%BOUNDARY_RADIA_ILW_H = GM%RADI%BOUNDARY_RADIA_ILW_D
            DO IBR=1,N_BOUNDARY_RADIA_DIM
                IF (.NOT.ALLOCATED(BOUNDARY_RADIA(IBR)%BAND)) CYCLE
                DO IB=1,NUMBER_SPECTRAL_BANDS
                    IF (.NOT.ALLOCATED(BOUNDARY_RADIA(IBR)%BAND(IB)%ILW)) CYCLE
                    NLIM = MIN(NUMBER_RADIATION_ANGLES,SIZE(BOUNDARY_RADIA(IBR)%BAND(IB)%ILW))
                    DO N=1,NLIM
                        BOUNDARY_RADIA(IBR)%BAND(IB)%ILW(N) = GM%RADI%BOUNDARY_RADIA_ILW_H(IBR,IB,N)
                    ENDDO
                ENDDO
            ENDDO

        END SUBROUTINE RADI_SYNC_BOUNDARY_RADIA_ILW_D2H

        !--------------------------------------------------------------------------
        ! Inter-mesh radiation-intensity exchange
        !--------------------------------------------------------------------------
        SUBROUTINE GPU_RADIATION_EXCHANGE
            INTEGER :: NSRC,NDST,N_REQ,N_REQ_MAX
            TYPE(MPI_REQUEST), ALLOCATABLE :: REQUEST(:)

            N_REQ_MAX = 0
            IF (N_MPI_PROCESSES>1) THEN
                DO NDST=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                    IF (.NOT.ALLOCATED(GPU_MESHES(NDST)%RADI)) CYCLE
                    IF (.NOT.ALLOCATED(GPU_MESHES(NDST)%RADI%OMESH)) CYCLE
                    DO NSRC=1,NMESHES
                        IF (PROCESS(NSRC)==MY_RANK) CYCLE
                        IF (ALLOCATED(GPU_MESHES(NDST)%RADI%OMESH(NSRC)%IL_R_D)) N_REQ_MAX = N_REQ_MAX + 1
                    ENDDO
                ENDDO
                DO NSRC=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                    IF (.NOT.ALLOCATED(GPU_MESHES(NSRC)%RADI)) CYCLE
                    IF (.NOT.ALLOCATED(GPU_MESHES(NSRC)%RADI%OMESH)) CYCLE
                    DO NDST=1,NMESHES
                        IF (PROCESS(NDST)==MY_RANK) CYCLE
                        IF (ALLOCATED(GPU_MESHES(NSRC)%RADI%OMESH(NDST)%IL_S_D)) N_REQ_MAX = N_REQ_MAX + 1
                    ENDDO
                ENDDO
            ENDIF
            IF (N_REQ_MAX>0) ALLOCATE(REQUEST(N_REQ_MAX))

            ! Post remote receives before any rank sends.
            N_REQ = 0
            DO NDST=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                IF (.NOT.ALLOCATED(GPU_MESHES(NDST)%RADI)) CYCLE
                IF (.NOT.ALLOCATED(GPU_MESHES(NDST)%RADI%OMESH)) CYCLE
                DO NSRC=1,NMESHES
                    IF (PROCESS(NSRC)==MY_RANK) CYCLE
                    IF (.NOT.ALLOCATED(GPU_MESHES(NDST)%RADI%OMESH(NSRC)%IL_R_D)) CYCLE
                    N_REQ = N_REQ + 1
                    CALL GPU_MPI_IRECV(GPU_MESHES(NDST)%RADI%OMESH(NSRC)%IL_R_D, &
                        SIZE(GPU_MESHES(NDST)%RADI%OMESH(NSRC)%IL_R_D), &
                        PROCESS(NSRC),NSRC,REQUEST(N_REQ))
                ENDDO
            ENDDO

            ! Same-rank neighbors retain the direct device-to-device path.
            DO NSRC=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                IF (.NOT.ALLOCATED(GPU_MESHES(NSRC)%RADI)) CYCLE
                IF (.NOT.ALLOCATED(GPU_MESHES(NSRC)%RADI%OMESH)) CYCLE
                DO NDST=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                    IF (PROCESS(NDST)/=MY_RANK) CYCLE
                    IF (.NOT.ALLOCATED(GPU_MESHES(NSRC)%RADI%OMESH(NDST)%IL_S_D)) CYCLE
                    IF (.NOT.ALLOCATED(GPU_MESHES(NDST)%RADI)) CYCLE
                    IF (.NOT.ALLOCATED(GPU_MESHES(NDST)%RADI%OMESH)) CYCLE
                    IF (.NOT.ALLOCATED(GPU_MESHES(NDST)%RADI%OMESH(NSRC)%IL_R_D)) CYCLE

                    GPU_MESHES(NDST)%RADI%OMESH(NSRC)%IL_R_D = GPU_MESHES(NSRC)%RADI%OMESH(NDST)%IL_S_D
                ENDDO
            ENDDO

            ! The radiation kernels have finished populating IL_S_D before MPI reads it.
            CALL GPU_WAIT_ALL

            DO NSRC=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                IF (.NOT.ALLOCATED(GPU_MESHES(NSRC)%RADI)) CYCLE
                IF (.NOT.ALLOCATED(GPU_MESHES(NSRC)%RADI%OMESH)) CYCLE
                DO NDST=1,NMESHES
                    IF (PROCESS(NDST)==MY_RANK) CYCLE
                    IF (.NOT.ALLOCATED(GPU_MESHES(NSRC)%RADI%OMESH(NDST)%IL_S_D)) CYCLE
                    N_REQ = N_REQ + 1
                    CALL GPU_MPI_ISEND(GPU_MESHES(NSRC)%RADI%OMESH(NDST)%IL_S_D, &
                        SIZE(GPU_MESHES(NSRC)%RADI%OMESH(NDST)%IL_S_D), &
                        PROCESS(NDST),NSRC,REQUEST(N_REQ))
                ENDDO
            ENDDO

            IF (N_REQ>0) THEN
                CALL GPU_MPI_WAITALL(REQUEST(1:N_REQ))
            ENDIF

            CALL GPU_WAIT_ALL

            IF (GPU_DEBUG_H2D_D2H) THEN
                ! Expose local send and receive mirrors to the following CPU stages.
                DO NSRC=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                    IF (.NOT.ALLOCATED(GPU_MESHES(NSRC)%RADI)) CYCLE
                    DO NDST=1,NMESHES
                        IF (.NOT.ALLOCATED(GPU_MESHES(NSRC)%RADI%OMESH(NDST)%IL_S_D)) CYCLE
                        MESHES(NSRC)%OMESH(NDST)%IL_S = GPU_MESHES(NSRC)%RADI%OMESH(NDST)%IL_S_D
                    ENDDO
                ENDDO
                DO NDST=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                    IF (.NOT.ALLOCATED(GPU_MESHES(NDST)%RADI)) CYCLE
                    DO NSRC=1,NMESHES
                        IF (.NOT.ALLOCATED(GPU_MESHES(NDST)%RADI%OMESH(NSRC)%IL_R_D)) CYCLE
                        MESHES(NDST)%OMESH(NSRC)%IL_R = GPU_MESHES(NDST)%RADI%OMESH(NSRC)%IL_R_D
                    ENDDO
                ENDDO
            ENDIF

            IF (ALLOCATED(REQUEST)) DEALLOCATE(REQUEST)

        END SUBROUTINE GPU_RADIATION_EXCHANGE


        LOGICAL FUNCTION RADI_NEEDS_UII_D2H()
            INTEGER :: N

            RADI_NEEDS_UII_D2H = .FALSE.
            DO N=1,N_SLCF
                IF (SLICE(N)%INDEX==18) THEN ! INTEGRATED INTENSITY
                    RADI_NEEDS_UII_D2H = .TRUE.
                    RETURN
                ENDIF
            ENDDO

        END FUNCTION RADI_NEEDS_UII_D2H

        !--------------------------------------------------------------------------
        ! Radiation-stage host orchestration and strict DOM sweep dispatch
        !--------------------------------------------------------------------------
        SUBROUTINE COMPUTE_RADIATION_GPU(T,NM,RAD_ITER)
            REAL(EB), INTENT(IN) :: T
            INTEGER, INTENT(IN) :: NM,RAD_ITER
            REAL(EB) :: TNOW

            TNOW = CURRENT_TIME()

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)

            IF (RADIATION) THEN
                RADIATION_COMPLETED(NM) = .FALSE.
                GM%RADI%T = T
                GM%RADI%RAD_ITER = RAD_ITER
                CALL RADIATION_FVM_GPU
            ELSE
                RADIATION_COMPLETED(NM) = .TRUE.
                IF (GPU_DEBUG_H2D_D2H) THEN
                    GM%QD = Q
                    GM%CHI_RD = CHI_R
                ENDIF
                IF (INIT_HRRPUV .AND. N_REACTIONS==0) THEN
                    GM%QD = 0._EB
                    GM%CHI_RD = 0._EB
                ENDIF
                IF (N_REACTIONS>0) THEN
                    CALL RADIATION_OFF_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1>>>(GM%QRD,GM%CHI_RD,GM%QD,IBP1,JBP1,KBP1)
                ENDIF
                IF (INIT_HRRPUV) CALL RADIATION_APPLY_INIT_HRRPUV_GPU(T,0)
                IF (GPU_DEBUG_H2D_D2H) THEN
                    Q = GM%QD
                    QR = GM%QRD
                ENDIF
            ENDIF

            T_USED(9) = T_USED(9) + CURRENT_TIME() - TNOW

        END SUBROUTINE COMPUTE_RADIATION_GPU

        SUBROUTINE RADIATION_FVM_GPU
            USE MIEV
            USE MATH_FUNCTIONS, ONLY : INTERPOLATE1D, EVALUATE_RAMP
            USE TRAN, ONLY : GET_IJK
            USE COMPLEX_GEOMETRY, ONLY : CC_CGSC,CC_SOLID
            USE PHYSICAL_FUNCTIONS, ONLY : GET_VOLUME_FRACTION, GET_MASS_FRACTION
            IMPLICIT NONE
            REAL(EB) :: RAP, AX, AXU, AXD, AY, AYU, AYD, AZ, AZU, AZD, VC, RU, R_DOWN, RP, AFD, &
                ILXU, ILYU, ILZU, QVAL, BBF, BBFA, NCSDROP, RSA_RAT,EFLUX,SOOT_MASS_FRACTION, &
                AIU_SUM,A_SUM,VOL,VC1,AY1,AZ1,DLO,COS_DLO,AILFU, &
                RAD_Q_SUM_PARTIAL,KFST4_SUM_PARTIAL,ALPHA_CC

            INTEGER  :: N,NN,IIG,JJG,KKG,I,J,K,IW,ICF,II,JJ,KK,IOR,IC,IWUP,IWDOWN, &
                ISTART, IEND, ISTEP, JSTART, JEND, JSTEP, &
                KSTART, KEND, KSTEP, NSTART, NEND, NSTEP, &
                I_UIID, N_UPDATES, IBND, NOM, ARRAY_INDEX,NRA, &
                IMIN, JMIN, KMIN, IMAX, JMAX, KMAX, N_SLICE, M_IJK, IJK, LL, &
                SID, GRDD_1D
            INTEGER  :: IADD,IFACE,INDCF
            REAL(EB) :: XID,YJD,ZKD,KAPPA_PART_SINGLE,DLF,DLA(3),TSI,TMP_EXTERIOR,TEMP_ORIENTATION(3)
            REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_GET
            INTEGER :: IID,JJD,KKD,IP
            LOGICAL :: UPDATE_INTENSITY
            REAL(EB), POINTER, DIMENSION(:,:,:) :: KFST4_GAS,IL,UIIOLD,KAPPA_PART,KFST4_PART,EXTCOE,SCAEFF,SCAEFF_G,IL_UP
            REAL(EB), POINTER, DIMENSION(:)     :: OUTRAD_W,INRAD_W,OUTRAD_F,INRAD_F,IL_F
            REAL(EB), DEVICE, POINTER, DIMENSION(:,:,:) :: KFST4_GAS_D,IL_D,UIIOLD_D,EXTCOE_D,KAPPA_PART_D,SCAEFF_D,SCAEFF_G_D, &
                KFST4_PART_D,IL_UP_D
            REAL(EB), DEVICE, POINTER, DIMENSION(:)     :: OUTRAD_W_D,INRAD_W_D
            TYPE(SURFACE_TYPE), POINTER :: SF
            TYPE(VENTS_TYPE), POINTER :: VT
            TYPE(RAD_FILE_TYPE), POINTER :: RF
            TYPE(LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
            TYPE(LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
            TYPE(CC_CUTFACE_TYPE), POINTER :: CF
            TYPE(CFACE_TYPE), POINTER :: CFA
            TYPE(BOUNDARY_RADIA_TYPE), POINTER :: BR,BR_UP,BR_DOWN
            TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
            TYPE(INITIALIZATION_TYPE), POINTER :: IN
            CHARACTER(20) :: FORMT

            ! WSGG host work variables.
            REAL(EB) :: X_H2O, X_CO2, MOL_RAT,PARTIAL_P,R_MIXTURE,TOTAL_P
            REAL(EB), ALLOCATABLE, DIMENSION(:) :: Z_ARRAY

            ALLOCATE(Z_ARRAY(N_TRACKED_SPECIES))

            KFST4_GAS  => WORK1
            IL         => WORK2
            UIIOLD     => WORK3
            EXTCOE     => WORK4
            KAPPA_PART => WORK5
            SCAEFF     => WORK6
            KFST4_PART => WORK7
            IL_UP      => WORK8
            SCAEFF_G   => WORK9
            OUTRAD_W   => WALL_WORK1
            INRAD_W    => WALL_WORK2
            KFST4_GAS_D => GM%WORK1D
            IL_D         => GM%WORK2D
            UIIOLD_D     => GM%WORK3D
            EXTCOE_D     => GM%WORK4D
            KAPPA_PART_D => GM%WORK5D
            SCAEFF_D     => GM%WORK6D
            KFST4_PART_D => GM%WORK7D
            IL_UP_D      => GM%WORK8D
            SCAEFF_G_D   => GM%WORK9D
            OUTRAD_W_D   => GM%RADI%WALL_WORK1_D
            INRAD_W_D    => GM%RADI%WALL_WORK2_D
            ! Ratio of solid angle used in scattering.

            NRA     = NUMBER_RADIATION_ANGLES
            RSA_RAT = 1._EB/(1._EB-1._EB/NRA)

            ! Determine whether this call updates the intensity field.
            IF ( MOD(RAD_CALL_COUNTER,TIME_STEP_INCREMENT)==0 .OR. INITIALIZATION_PHASE .OR. ICYC==1 .OR. UPDATE_ALL_ANGLES) THEN
                UPDATE_INTENSITY   = .TRUE.
                EXCHANGE_RADIATION = .TRUE.
            ELSE
                UPDATE_INTENSITY   = .FALSE.
                EXCHANGE_RADIATION = .FALSE.
            ENDIF

            ! Advance the radiation-call counter after the final iteration of a time step.
            IF (GM%RADI%RAD_ITER==RADIATION_ITERATIONS) RAD_CALL_COUNTER  = RAD_CALL_COUNTER + 1

            IF (GPU_DEBUG_H2D_D2H) THEN
                GM%ZZD = ZZ
                GM%TMPD = TMP
                GM%QD = Q
                GM%CHI_RD = CHI_R
                GM%RADI%UII_D = UII
                GM%RADI%UIID_D = UIID
                CALL RADI_SYNC_WALL_RADIATION_H2D
                CALL RADI_SYNC_BOUNDARY_RADIA_ILW_H2D
                IF (N_LP_ARRAY_INDICES>0 .AND. ASSOCIATED(QR_W)) GM%RADI%QR_W_D = QR_W
                IF (STORE_RADIATION_TERMS) THEN
                    GM%RADI%EMISSION_D = RADIATION_EMISSION
                    GM%RADI%ABSORPTION_D = RADIATION_ABSORPTION
                ENDIF
            ENDIF

            IF (INIT_HRRPUV) THEN
                IF (GM%RADI%RAD_ITER==1 .AND. N_REACTIONS==0) THEN
                    GM%QD = 0._EB
                    GM%CHI_RD = 0._EB
                ELSEIF (GM%RADI%RAD_ITER>1) THEN
                    CALL RADIATION_APPLY_INIT_HRRPUV_GPU(GM%RADI%T,2)
                ENDIF
            ENDIF

            IF (UPDATE_INTENSITY) THEN
                GM%WALL_Q_RAD_IN = 0._EB
            ENDIF

            ! IF (WIDE_BAND_MODEL .OR. WSGG_MODEL) THEN
            !    QRD = 0._EB
            !    IF (STORE_RADIATION_TERMS) THEN
            !       RADIATION_EMISSION_D = 0._EB
            !       RADIATION_ABSORPTION_D = 0._EB
            !    ENDIF
            !    IF (N_LP_ARRAY_INDICES>0 .AND. PARTICLES_EXISTED .AND. ASSOCIATED(QR_W)) RADI_QR_W_D = 0._EB
            ! ENDIF

            ! Loop over spectral bands
            BAND_LOOP: DO IBND = 1,NUMBER_SPECTRAL_BANDS

                KAPPA_PART = 0._EB
                KFST4_GAS  = 0._EB
                KFST4_PART = 0._EB
                SCAEFF = 0._EB
                SCAEFF_G = 0._EB

                BBFA = 1._EB

                CALL RADIATION_ZERO_BAND_FIELDS_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1>>>( &
                    GM%RADI%KAPPA_GAS_D,KAPPA_PART_D,KFST4_GAS_D,KFST4_PART_D,SCAEFF_D,SCAEFF_G_D,IBP1,JBP1,KBP1)

                ! Set or interpolate gas absorption coefficient KAPPA_GAS.

                IF (KAPPA0>=0._EB) THEN

                    KAPPA_GAS = KAPPA0
                    GM%RADI%KAPPA_GAS_D = KAPPA0

                ELSEIF (N_RADCAL_ARRAY_SIZE>0) THEN
                    CALL RADIATION_KAPPA_GAS_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1>>>( &
                        GM%RADI%KAPPA_GAS_D,GM%ZZD,GM%TMPD,GM%CELL_INDEXD,GM%CELL_SOLID,MWR_ZD, &
                        RADI_Z2RADCAL_SPECIES_D,RADI_RADCAL_SPECIES2KAPPA_D,RADI_RADCAL_SPECIES_INDEX_D, &
                        IBND,N_TRACKED_SPECIES,N_RADCAL_ARRAY_SIZE,N_KAPPA_T,N_KAPPA_Y,IBAR,JBAR,KBAR,RTMPMIN,RTMPMAX,TWO_EPSILON_EB)
                ELSE
                    GM%RADI%KAPPA_GAS_D = KAPPA_GAS
                ENDIF

                ! Assemble KFST4_GAS and its reduction sums for the gray-gas source.
                CALL RADIATION_KFST4_SUM_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1>>>( &
                    KFST4_GAS_D,GM%RADI%KAPPA_GAS_D,GM%TMPD,GM%QD,GM%CHI_RD,GM%RADI%UII_D,GM%CELL_INDEXD, &
                    GM%CELL_SOLID,GM%RD,GM%DXD,GM%DYD,GM%DZD,GM%RADI%RAD_Q_SUM_D,GM%RADI%KFST4_SUM_D, &
                    IBAR,JBAR,KBAR,FOUR_SIGMA,QR_CLIP)

                ! Apply the REAC radiative-fraction correction to the RTE source.
                CALL RADIATION_RTE_SOURCE_CORRECTION_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1>>>( &
                    KFST4_GAS_D,GM%QD,GM%CHI_RD,GM%CELL_INDEXD,GM%CELL_SOLID,IBAR,JBAR,KBAR,QR_CLIP, &
                    RTE_SOURCE_CORRECTION_FACTOR)

                IF (INIT_HRRPUV) CALL RADIATION_APPLY_INIT_HRRPUV_GPU(GM%RADI%T,1,KFST4_GAS_D)

                ! Assemble EXTCOE from gas, particle, and scattering coefficients.
                CALL RADIATION_EXTCOE_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1>>>( &
                    EXTCOE_D,GM%RADI%KAPPA_GAS_D,KAPPA_PART_D,SCAEFF_D,SCAEFF_G_D,GM%CELL_INDEXD,GM%CELL_SOLID, &
                    IBAR,JBAR,KBAR,RSA_RAT)

                INTENSITY_UPDATE: IF (UPDATE_INTENSITY) THEN

                    ! Preserve UII in UIIOLD and reset the full-mesh UII accumulator.
                    CALL RADIATION_UIIOLD_ZERO_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1>>>( &
                        UIIOLD_D,GM%RADI%UII_D,IBP1,JBP1,KBP1)

                    ! Assemble wall and open-boundary outgoing intensity.
                    BBF = 1.0_EB
                    CALL RADIATION_OUTRAD_W_GPU<<<GM%GRDD_TOTAL_WALL_1D,GPU_BLKD_1D>>>( &
                        OUTRAD_W_D,GM%WALL_TOPO_AOS,GM%WALL_Q_RAD_OUT, &
                        N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS, BBFA, RPI, SIGMA, TMPA, NULL_BOUNDARY, OPEN_BOUNDARY)
                    ! Initialize the solid-wall incoming-intensity term.
                    INRAD_W_D = 0._EB
                    CALL RADIATION_INRAD_W_GPU<<<GM%GRDD_TOTAL_WALL_1D,GPU_BLKD_1D>>>( &
                        INRAD_W_D,GM%WALL_TOPO_AOS,GM%RADI%WALL_BR_INDEX_D,DLN_D, &
                        GM%RADI%BOUNDARY_RADIA_ILW_D,IBND,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS,NUMBER_RADIATION_ANGLES, &
                        SOLID_BOUNDARY,SIZE(GM%RADI%BOUNDARY_RADIA_ILW_D,1),SIZE(GM%RADI%BOUNDARY_RADIA_ILW_D,2))

                    ! If updating intensities first time, sweep ALL angles

                    N_UPDATES = 1
                    IF (INITIALIZATION_PHASE .OR. ICYC==1 .OR. UPDATE_ALL_ANGLES) N_UPDATES = ANGLE_INCREMENT

                    UPDATE_LOOP: DO I_UIID = 1,N_UPDATES

                        ! Advance the host-side angle-group counter.
                        ANGLE_INC_COUNTER = MOD(ANGLE_INC_COUNTER,ANGLE_INCREMENT) + 1

                        IF (ANGLE_INC_COUNTER==ANGLE_INCREMENT) RADIATION_COMPLETED(GM%NM) = .TRUE.

                        ! Clear current-angle UIID and open-wall ILW values.
                        CALL RADIATION_ZERO_UIID_ANGLE_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1>>>( &
                            GM%RADI%UIID_D,ANGLE_INC_COUNTER,IBP1,JBP1,KBP1,UIIDIM)
                        CALL RADIATION_ZERO_OPEN_ILW_ANGLE_GPU<<<GM%GRDD_TOTAL_WALL_1D,GPU_BLKD_1D>>>( &
                            GM%RADI%BOUNDARY_RADIA_ILW_D,GM%WALL_TOPO_AOS,GM%RADI%WALL_BR_INDEX_D, &
                            IBND,ANGLE_INC_COUNTER,N_EXTERNAL_WALL_CELLS,OPEN_BOUNDARY, &
                            SIZE(GM%RADI%BOUNDARY_RADIA_ILW_D,1),SIZE(GM%RADI%BOUNDARY_RADIA_ILW_D,2))

                        ! Sweep angles in descending order for cylindrical boundary coupling.
                        NSTART    = NRA - ANGLE_INC_COUNTER + 1
                        NEND      = 1
                        NSTEP     = -ANGLE_INCREMENT

                        ! Initialize IL to ambient incoming intensity.
                        CALL RADIATION_INIT_IL_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1>>>( &
                            IL_D,BBFA*RPI_SIGMA*TMPA4,IBP1,JBP1,KBP1)

                        ! Sweep the selected discrete directions.
                        ANGLE_LOOP: DO N = NSTART,NEND,NSTEP

                            CELL_ILW = -HUGE(EB)

                            ! Apply wall conditions first, then dispatch interpolated walls by donor mesh.

                            ! Boundary conditions: Intensities leaving the boundaries.

                            IF (.NOT.CYLINDRICAL) THEN
                                GM%RADI%CELL_ILW_D = -HUGE(EB)
                                ! Apply wall leaving-boundary intensities for this direction.
                                CALL RADIATION_WALL_LOOP1_GPU<<<GM%GRDD_TOTAL_WALL_1D,GPU_BLKD_1D>>>( &
                                    IL_D,GM%RADI%CELL_ILW_D,GM%RADI%BOUNDARY_RADIA_ILW_D,OUTRAD_W_D,INRAD_W_D, &
                                    GM%WALL_TOPO_AOS,GM%WALL_EMISSIVITY,GM%RADI%WALL_BR_INDEX_D,GM%CELL_INDEXD, &
                                    DLN_D, DLM_D, N, IBND, N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS, &
                                    BBFA*RPI_SIGMA*TMPA4, RPI, TWO_D, NULL_BOUNDARY, OPEN_BOUNDARY, MIRROR_BOUNDARY, &
                                    INTERPOLATED_BOUNDARY,IBP1,JBP1,SIZE(GM%RADI%CELL_ILW_D,1), &
                                    SIZE(GM%RADI%BOUNDARY_RADIA_ILW_D,1),SIZE(GM%RADI%BOUNDARY_RADIA_ILW_D,2), &
                                    NUMBER_RADIATION_ANGLES)

                                DO NOM=1,NMESHES
                                    IF (.NOT.ALLOCATED(GM%RADI%OMESH(NOM)%IL_R_D)) CYCLE
                                    CALL RADIATION_INTERPOLATED_WALL_GPU<<<GM%GRDD_EXTERNAL_WALL_1D,GPU_BLKD_1D>>>( &
                                        IL_D,GM%WALL_TOPO_AOS,GM%EXTERNAL_WALL_NOMD,DLN_D,GM%RADI%WALL_NIC_MIN_D, &
                                        GM%RADI%WALL_NIC_MAX_D,GM%RADI%OMESH(NOM)%IL_R_D,N,IBND, &
                                        N_EXTERNAL_WALL_CELLS,NOM,SIZE(GM%RADI%OMESH(NOM)%IL_R_D,1),INTERPOLATED_BOUNDARY, &
                                        SIZE(IL_D,1),SIZE(IL_D,2),LBOUND(IL_D,1),LBOUND(IL_D,2),LBOUND(IL_D,3), &
                                        SIZE(GM%RADI%OMESH(NOM)%IL_R_D,2),LBOUND(GM%RADI%OMESH(NOM)%IL_R_D,1), &
                                        LBOUND(GM%RADI%OMESH(NOM)%IL_R_D,2),LBOUND(GM%RADI%OMESH(NOM)%IL_R_D,3), &
                                        SIZE(DLN_D,1),LBOUND(DLN_D,1),LBOUND(DLN_D,2),GM%WALL_LB)
                                ENDDO
                            ENDIF


                            ! Determine sweep direction in physical space
                            ! Select I/J/K sweep direction and bounds for one of eight octants.
                            ISTART = 1
                            JSTART = 1
                            KSTART = 1
                            IEND   = IBAR
                            JEND   = JBAR
                            KEND   = KBAR
                            ISTEP  = 1
                            JSTEP  = 1
                            KSTEP  = 1
                            IMIN = ISTART
                            JMIN = JSTART
                            KMIN = KSTART
                            IMAX = IEND
                            JMAX = JEND
                            KMAX = KEND
                            IF (DLX(N) < 0._EB) THEN
                                ISTART = IBAR
                                IEND   = 1
                                ISTEP  = -1
                                IMIN = IEND
                                IMAX = ISTART
                            ENDIF
                            IF (DLY(N) < 0._EB) THEN
                                JSTART = JBAR
                                JEND   = 1
                                JSTEP  = -1
                                JMIN = JEND
                                JMAX = JSTART
                            ENDIF
                            IF (DLZ(N) < 0._EB) THEN
                                KSTART = KBAR
                                KEND   = 1
                                KSTEP  = -1
                                KMIN = KEND
                                KMAX = KSTART
                            ENDIF


                            ! Forward sweep maps are reused by coordinate reflection for all eight directions.
                            CALL RADI_RUN_SWEEP(IL_D,IL_UP_D,GM%RADI%CELL_ILW_D,EXTCOE_D,KFST4_GAS_D,KFST4_PART_D, &
                                SCAEFF_D,SCAEFF_G_D,UIIOLD_D,N,ISTEP,JSTEP,KSTEP,RSA_RAT)


                            ! Copy the Y-downwind intensities to Y-upwind in cylindrical case

                            ! Boundary values: Incoming radiation
                            ! Update solid/interpolated walls and accumulate open-boundary intensity.
                            CALL RADIATION_WALL_INCOMING_GPU<<<GM%GRDD_TOTAL_WALL_1D,GPU_BLKD_1D>>>( &
                                INRAD_W_D,GM%RADI%BOUNDARY_RADIA_ILW_D,IL_D,GM%WALL_TOPO_AOS,GM%RADI%WALL_BR_INDEX_D,DLN_D, &
                                N,IBND,ANGLE_INC_COUNTER,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS,TWO_D,CYLINDRICAL, &
                                NULL_BOUNDARY,OPEN_BOUNDARY,IBP1,JBP1, &
                                SIZE(GM%RADI%BOUNDARY_RADIA_ILW_D,1),SIZE(GM%RADI%BOUNDARY_RADIA_ILW_D,2))

                            ! ENDIF

                            ! Calculate integrated intensity UIID

                            ! Accumulate angle-weighted intensity.
                            IF (WIDE_BAND_MODEL .OR. WSGG_MODEL) THEN
                                CALL RADIATION_ACCUMULATE_UIID_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1>>>( &
                                    GM%RADI%UIID_D,IL_D,GM%CELL_INDEXD,GM%CELL_SOLID,IBND,N,IBP1,JBP1,KBP1,UIIDIM,WEIGH_CYL,RSA_D)
                            ELSE
                                CALL RADIATION_ACCUMULATE_UIID_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1>>>( &
                                    GM%RADI%UIID_D,IL_D,GM%CELL_INDEXD,GM%CELL_SOLID,ANGLE_INC_COUNTER,N, &
                                    IBP1,JBP1,KBP1,UIIDIM,WEIGH_CYL,RSA_D)
                            ENDIF

                            DO NOM=1,NMESHES
                                IF (.NOT.ALLOCATED(GM%RADI%OMESH(NOM)%IL_S_D)) CYCLE
                                GRDD_1D = (GM%OMESH(NOM)%NIC_S + GPU_BLKD_1D - 1)/GPU_BLKD_1D
                                CALL RADIATION_PACK_IL_S_GPU<<<GRDD_1D,GPU_BLKD_1D>>>( &
                                    IL_D,GM%RADI%OMESH(NOM)%IL_S_D,GM%OMESH(NOM)%IIO_S_D,GM%OMESH(NOM)%JJO_S_D, &
                                    GM%OMESH(NOM)%KKO_S_D,N,IBND,GM%OMESH(NOM)%NIC_S,IBP1,JBP1,KBP1, &
                                    SIZE(GM%RADI%OMESH(NOM)%IL_S_D,1),SIZE(GM%RADI%OMESH(NOM)%IL_S_D,2))
                            ENDDO


                        ENDDO ANGLE_LOOP

                    ENDDO UPDATE_LOOP

                    ! Compute incoming flux on walls and particles
                    ! Convert the incoming intensity integral to wall heat flux.
                    CALL RADIATION_SOLID_WALL_Q_RAD_IN_GPU<<<GM%GRDD_TOTAL_WALL_1D,GPU_BLKD_1D>>>( &
                        GM%WALL_Q_RAD_IN,INRAD_W_D,GM%WALL_TOPO_AOS,GM%WALL_SF_INDEX,GM%WALL_EMISSIVITY, &
                        SF_EXTERNAL_FLUX,SF_RAMP_EFLUX_INDEX,SF_RAMP_EFLUX_TAU,GM%WALL_T_IGN, &
                        RAMP_T_MIN_D,RAMP_SPAN_D,RAMP_RDT_D,RAMP_NPOINTS_D,RAMP_OFFSET_D,RAMP_INTERP_D, &
                        GM%RADI%T,T_BEGIN,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS,BBFA,SOLID_BOUNDARY)

                ENDIF INTENSITY_UPDATE


            ENDDO BAND_LOOP

            ! Sum direction groups and update open-boundary incoming heat flux.

            IF (UPDATE_INTENSITY) THEN
                CALL RADIATION_SUM_UIID_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1>>>( &
                    GM%RADI%UII_D,GM%RADI%UIID_D,IBP1,JBP1,KBP1,UIIDIM)
                CALL RADIATION_OPEN_Q_RAD_IN_GPU<<<GM%GRDD_TOTAL_WALL_1D,GPU_BLKD_1D>>>( &
                    GM%WALL_Q_RAD_IN,GM%RADI%BOUNDARY_RADIA_ILW_D,GM%WALL_TOPO_AOS,GM%RADI%WALL_BR_INDEX_D, &
                    N_EXTERNAL_WALL_CELLS,NUMBER_SPECTRAL_BANDS, &
                    NUMBER_RADIATION_ANGLES,OPEN_BOUNDARY,SIZE(GM%RADI%BOUNDARY_RADIA_ILW_D,1))

            ENDIF

            ! Assemble the gray-gas energy source QR; wide-band QR is handled per band.
            IF (.NOT. (WIDE_BAND_MODEL .OR. WSGG_MODEL)) THEN
                CALL RADIATION_QR_GRAY_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1>>>( &
                    GM%QRD,GM%RADI%QR_W_D,GM%RADI%EMISSION_D,GM%RADI%ABSORPTION_D,GM%RADI%KAPPA_GAS_D,KFST4_GAS_D, &
                    KAPPA_PART_D,KFST4_PART_D,GM%RADI%UII_D,GM%CELL_INDEXD,GM%CELL_SOLID,IBP1,JBP1,KBP1, &
                    STORE_RADIATION_TERMS,N_LP_ARRAY_INDICES>0)
            ENDIF

            IF (GPU_DEBUG_H2D_D2H) THEN
                Q     = GM%QD
                QR    = GM%QRD
                UII   = GM%RADI%UII_D
                UIID  = GM%RADI%UIID_D
                CALL RADI_SYNC_WALL_RADIATION_D2H
                CALL RADI_SYNC_BOUNDARY_RADIA_ILW_D2H
                IF (N_LP_ARRAY_INDICES>0 .AND. ASSOCIATED(QR_W)) QR_W = GM%RADI%QR_W_D
                IF (STORE_RADIATION_TERMS) THEN
                    RADIATION_EMISSION   = GM%RADI%EMISSION_D
                    RADIATION_ABSORPTION = GM%RADI%ABSORPTION_D
                ENDIF
            ENDIF

            IF (.NOT.GPU_DEBUG_H2D_D2H .AND. UPDATE_INTENSITY .AND. RADI_NEEDS_UII_D2H()) UII = GM%RADI%UII_D

            ! UPDATE_ALL_ANGLES is a one-time only logical. If needed again, it will be reset again.

            UPDATE_ALL_ANGLES = .FALSE.

        END SUBROUTINE RADIATION_FVM_GPU

        SUBROUTINE RADI_RUN_SWEEP(IL_D,IL_UP_D,CELL_ILW_D,EXTCOE_D,KFST4_GAS_D,KFST4_PART_D, &
            SCAEFF_D,SCAEFF_G_D,UIIOLD_D,N_IN,ISTEP_IN,JSTEP_IN,KSTEP_IN,RSA_RAT_IN)
            IMPLICIT NONE
            REAL(EB), DEVICE, DIMENSION(0:,0:,0:) :: IL_D,IL_UP_D,EXTCOE_D,KFST4_GAS_D,KFST4_PART_D,SCAEFF_D,SCAEFF_G_D,UIIOLD_D
            REAL(EB), DEVICE, DIMENSION(:,:) :: CELL_ILW_D
            INTEGER, INTENT(IN) :: N_IN,ISTEP_IN,JSTEP_IN,KSTEP_IN
            REAL(EB), INTENT(IN) :: RSA_RAT_IN
            INTEGER :: SID,BASE,M_IJK,GRDD_1D

            DO SID=1,GM%RADI%NSLICES
                BASE = GM%RADI%SLICE_OFFSET_H(SID)
                M_IJK = GM%RADI%SLICE_OFFSET_H(SID+1) - BASE
                IF (M_IJK <= 0) CYCLE
                GRDD_1D = (M_IJK + GPU_BLKD_1D - 1) / GPU_BLKD_1D
                CALL RADIATION_SWEEP_3D_SLICE_GPU<<<GRDD_1D,GPU_BLKD_1D>>>( &
                    IL_D,IL_UP_D,CELL_ILW_D,EXTCOE_D,KFST4_GAS_D,KFST4_PART_D,SCAEFF_D,SCAEFF_G_D,UIIOLD_D, &
                    GM%CELL_INDEXD,GM%CELL_SOLID,GM%RADI%IJK_SLICE_MAP_D,GM%DXD,GM%DYD,GM%DZD,DLX_D,DLY_D,DLZ_D,RSA_D, &
                    N_IN,BASE,M_IJK,ISTEP_IN,JSTEP_IN,KSTEP_IN,RSA_RAT_IN,RFPI,SOLID_PARTICLES, &
                    IBP1,JBP1,KBP1,GM%CELL_LB,GM%CELL_UB,SIZE(CELL_ILW_D,1),SIZE(GM%RADI%IJK_SLICE_MAP_D,1))
            ENDDO

        END SUBROUTINE RADI_RUN_SWEEP

        !--------------------------------------------------------------------------
        ! Prescribed HRRPUV source from INIT records
        !--------------------------------------------------------------------------
        SUBROUTINE RADIATION_APPLY_INIT_HRRPUV_GPU(T_IN,MODE,KFST4_GAS_D)
            USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
            IMPLICIT NONE

            REAL(EB), INTENT(IN) :: T_IN
            INTEGER, INTENT(IN) :: MODE
            REAL(EB), DEVICE, OPTIONAL, DIMENSION(0:,0:,0:) :: KFST4_GAS_D
            REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: HRRPUV_SOURCE,HRRPUV_RADI_SOURCE
            TYPE(INITIALIZATION_TYPE), POINTER :: IN
            REAL(EB) :: TIME_RAMP_FACTOR,RR,SOURCE_VALUE
            INTEGER :: N,I,J,K

            ALLOCATE(HRRPUV_SOURCE(0:IBP1,0:JBP1,0:KBP1))
            HRRPUV_SOURCE = 0._EB
            IF (MODE==1) THEN
                ALLOCATE(HRRPUV_RADI_SOURCE(0:IBP1,0:JBP1,0:KBP1))
                HRRPUV_RADI_SOURCE = 0._EB
            ENDIF

            DO N=1,N_INIT
                IN => INITIALIZATION(N)
                IF (IN%HRRPUV<=0._EB) CYCLE
                TIME_RAMP_FACTOR = 1._EB
                IF (IN%RAMP_Q_INDEX>0) TIME_RAMP_FACTOR = EVALUATE_RAMP(T_IN,IN%RAMP_Q_INDEX)
                SOURCE_VALUE = TIME_RAMP_FACTOR*IN%HRRPUV*IN%VOLUME_ADJUST(GM%NM)
                DO K=0,KBP1
                    DO J=0,JBP1
                        DO I=0,IBP1
                            SELECT CASE(IN%SHAPE)
                                CASE DEFAULT
                                    IF (XC(I)<IN%X1 .OR. XC(I)>IN%X2 .OR. &
                                        YC(J)<IN%Y1 .OR. YC(J)>IN%Y2 .OR. &
                                        ZC(K)<IN%Z1 .OR. ZC(K)>IN%Z2) CYCLE
                                CASE('CONE')
                                    RR = MAX(0._EB,IN%RADIUS*(1._EB-(ZC(K)-IN%Z0)/IN%HEIGHT))
                                    IF ((XC(I)-IN%X0)**2+(YC(J)-IN%Y0)**2>RR**2) CYCLE
                            END SELECT
                            HRRPUV_SOURCE(I,J,K) = HRRPUV_SOURCE(I,J,K) + SOURCE_VALUE
                            IF (MODE==1) HRRPUV_RADI_SOURCE(I,J,K) = HRRPUV_RADI_SOURCE(I,J,K) + IN%CHI_R*SOURCE_VALUE
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO

            GM%WORK11D = HRRPUV_SOURCE
            IF (MODE==1) THEN
                GM%WORK12D = HRRPUV_RADI_SOURCE
                CALL RADIATION_APPLY_INIT_HRRPUV_KERNEL<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1>>>( &
                    GM%QD,GM%WORK11D,KFST4_GAS_D,GM%WORK12D,MODE,IBP1,JBP1,KBP1)
                DEALLOCATE(HRRPUV_RADI_SOURCE)
            ELSE
                CALL RADIATION_APPLY_INIT_HRRPUV_KERNEL<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1>>>( &
                    GM%QD,GM%WORK11D,GM%WORK12D,GM%WORK12D,MODE,IBP1,JBP1,KBP1)
            ENDIF
            DEALLOCATE(HRRPUV_SOURCE)

        END SUBROUTINE RADIATION_APPLY_INIT_HRRPUV_GPU

        !--------------------------------------------------------------------------
        ! Radiation device helper and CUDA kernels
        !--------------------------------------------------------------------------
        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_APPLY_INIT_HRRPUV_KERNEL( &
            Q_D,HRRPUV_SOURCE_D,KFST4_GAS_D,HRRPUV_RADI_SOURCE_D,MODE,IBP1_IN,JBP1_IN,KBP1_IN)

            INTEGER, VALUE :: MODE,IBP1_IN,JBP1_IN,KBP1_IN
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: Q_D,HRRPUV_SOURCE_D,KFST4_GAS_D,HRRPUV_RADI_SOURCE_D
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1)*BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1)*BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1)*BLOCKDIM%Z + THREADIDX%Z - 1
            IF (I<0 .OR. I>IBP1_IN .OR. J<0 .OR. J>JBP1_IN .OR. K<0 .OR. K>KBP1_IN) RETURN

            SELECT CASE(MODE)
                CASE(1)
                    Q_D(I,J,K) = Q_D(I,J,K) + HRRPUV_SOURCE_D(I,J,K)
                    KFST4_GAS_D(I,J,K) = KFST4_GAS_D(I,J,K) + HRRPUV_RADI_SOURCE_D(I,J,K)
                CASE(2)
                    Q_D(I,J,K) = Q_D(I,J,K) - HRRPUV_SOURCE_D(I,J,K)
                CASE DEFAULT
                    Q_D(I,J,K) = Q_D(I,J,K) + HRRPUV_SOURCE_D(I,J,K)
            END SELECT

        END SUBROUTINE RADIATION_APPLY_INIT_HRRPUV_KERNEL

        ATTRIBUTES(DEVICE) REAL(EB) FUNCTION GET_KAPPA_GPU_DEVICE(ZZ_D,I,J,K,IBAR_IN,JBAR_IN,KBAR_IN,TMP_IN,IBND,N_TRACKED_SPECIES_IN, &
            N_RADCAL_ARRAY_SIZE_IN,N_KAPPA_T_IN,N_KAPPA_Y_IN, &
            RTMPMIN_IN,RTMPMAX_IN,TWO_EPS_IN, &
            MWR_Z_D,Z2RADCAL_SPECIES_D,RADCAL_SPECIES2KAPPA_D, &
            RADCAL_SPECIES_INDEX_D)

            IMPLICIT NONE
            INTEGER, VALUE :: I,J,K,IBAR_IN,JBAR_IN,KBAR_IN,IBND,N_TRACKED_SPECIES_IN,N_RADCAL_ARRAY_SIZE_IN,N_KAPPA_T_IN,N_KAPPA_Y_IN
            REAL(EB), VALUE :: TMP_IN,RTMPMIN_IN,RTMPMAX_IN,TWO_EPS_IN
            REAL(EB), DEVICE, DIMENSION(0:IBAR_IN+1,0:JBAR_IN+1,0:KBAR_IN+1,1:*) :: ZZ_D
            REAL(EB), DEVICE, DIMENSION(1:*) :: MWR_Z_D
            REAL(EB), DEVICE, DIMENSION(1:N_RADCAL_ARRAY_SIZE_IN,1:*) :: Z2RADCAL_SPECIES_D
            REAL(EB), DEVICE, DIMENSION(1:N_RADCAL_ARRAY_SIZE_IN,0:N_KAPPA_Y_IN,0:N_KAPPA_T_IN,1:*) :: RADCAL_SPECIES2KAPPA_D
            INTEGER, DEVICE, DIMENSION(1:*) :: RADCAL_SPECIES_INDEX_D
            REAL(EB) :: KAPPA_TEMP,INT_FAC,KAPPA_SUM,SCALED_Y_RADCAL_SPECIES,MWA,TYY_FAC_LOCAL
            INTEGER :: LBND,UBND,N,NS,TYY

            KAPPA_SUM = 0._EB
            TYY_FAC_LOCAL = REAL(N_KAPPA_T_IN,EB)/(RTMPMAX_IN-RTMPMIN_IN)
            TYY = MAX(0,MIN(N_KAPPA_T_IN,INT((TMP_IN-RTMPMIN_IN)*TYY_FAC_LOCAL)))

            MWA = 0._EB
            DO NS = 1,N_TRACKED_SPECIES_IN
                MWA = MWA + MWR_Z_D(NS)*ZZ_D(I,J,K,NS)
            ENDDO
            MWA = 1._EB/MWA

            DO N = 1,N_RADCAL_ARRAY_SIZE_IN
                SCALED_Y_RADCAL_SPECIES = 0._EB
                DO NS = 1,N_TRACKED_SPECIES_IN
                    SCALED_Y_RADCAL_SPECIES = SCALED_Y_RADCAL_SPECIES + Z2RADCAL_SPECIES_D(N,NS)*ZZ_D(I,J,K,NS)
                ENDDO
                IF (SCALED_Y_RADCAL_SPECIES<TWO_EPS_IN) CYCLE
                IF (RADCAL_SPECIES_INDEX_D(N)==16) THEN
                    INT_FAC = MAX(0._EB,SCALED_Y_RADCAL_SPECIES)**0.25_EB
                ELSE
                    INT_FAC = MAX(0._EB,SCALED_Y_RADCAL_SPECIES*MWA)**0.25_EB
                ENDIF
                LBND = INT(INT_FAC)
                INT_FAC = INT_FAC - LBND
                LBND = MIN(LBND,N_KAPPA_Y_IN)
                UBND = MIN(LBND+1,N_KAPPA_Y_IN)
                KAPPA_TEMP = RADCAL_SPECIES2KAPPA_D(N,LBND,TYY,IBND)
                KAPPA_SUM = KAPPA_SUM + KAPPA_TEMP + INT_FAC*(RADCAL_SPECIES2KAPPA_D(N,UBND,TYY,IBND)-KAPPA_TEMP)
            ENDDO

            GET_KAPPA_GPU_DEVICE = KAPPA_SUM

        END FUNCTION GET_KAPPA_GPU_DEVICE

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_ZERO_BAND_FIELDS_GPU( &
            KAPPA_GAS_D,KAPPA_PART_D,KFST4_GAS_D,KFST4_PART_D,SCAEFF_D,SCAEFF_G_D,IBP1_IN,JBP1_IN,KBP1_IN)
            INTEGER, VALUE :: IBP1_IN,JBP1_IN,KBP1_IN
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: KAPPA_GAS_D,KAPPA_PART_D
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: KFST4_GAS_D,KFST4_PART_D,SCAEFF_D,SCAEFF_G_D
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1)*BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1)*BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1)*BLOCKDIM%Z + THREADIDX%Z - 1
            IF (I<0 .OR. I>IBP1_IN .OR. J<0 .OR. J>JBP1_IN .OR. K<0 .OR. K>KBP1_IN) RETURN

            KAPPA_GAS_D(I,J,K) = 0._EB
            KAPPA_PART_D(I,J,K) = 0._EB
            KFST4_GAS_D(I,J,K) = 0._EB
            KFST4_PART_D(I,J,K) = 0._EB
            SCAEFF_D(I,J,K) = 0._EB
            SCAEFF_G_D(I,J,K) = 0._EB
        END SUBROUTINE RADIATION_ZERO_BAND_FIELDS_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_KAPPA_GAS_GPU(KAPPA_GAS_D,ZZ_D,TMP_D,CELL_INDEX_D,CELL_SOLID_D,MWR_Z_D, &
            Z2RADCAL_SPECIES_D,RADCAL_SPECIES2KAPPA_D,RADCAL_SPECIES_INDEX_D, &
            IBND,N_TRACKED_SPECIES_IN,N_RADCAL_ARRAY_SIZE_IN,N_KAPPA_T_IN,N_KAPPA_Y_IN, &
            IBAR_IN,JBAR_IN,KBAR_IN,RTMPMIN_IN,RTMPMAX_IN,TWO_EPS_IN)

            IMPLICIT NONE
            INTEGER, VALUE :: IBND,N_TRACKED_SPECIES_IN,N_RADCAL_ARRAY_SIZE_IN,N_KAPPA_T_IN,N_KAPPA_Y_IN
            INTEGER, VALUE :: IBAR_IN,JBAR_IN,KBAR_IN
            REAL(EB), VALUE :: RTMPMIN_IN,RTMPMAX_IN,TWO_EPS_IN
            REAL(EB), DEVICE, DIMENSION(0:IBAR_IN+1,0:JBAR_IN+1,0:*) :: KAPPA_GAS_D,TMP_D
            REAL(EB), DEVICE, DIMENSION(0:IBAR_IN+1,0:JBAR_IN+1,0:KBAR_IN+1,1:*) :: ZZ_D
            INTEGER, DEVICE, DIMENSION(0:IBAR_IN+1,0:JBAR_IN+1,0:*) :: CELL_INDEX_D
            LOGICAL, DEVICE, DIMENSION(0:*) :: CELL_SOLID_D
            REAL(EB), DEVICE, DIMENSION(1:*) :: MWR_Z_D
            REAL(EB), DEVICE, DIMENSION(1:N_RADCAL_ARRAY_SIZE_IN,1:*) :: Z2RADCAL_SPECIES_D
            REAL(EB), DEVICE, DIMENSION(1:N_RADCAL_ARRAY_SIZE_IN,0:N_KAPPA_Y_IN,0:N_KAPPA_T_IN,1:*) :: RADCAL_SPECIES2KAPPA_D
            INTEGER, DEVICE, DIMENSION(1:*) :: RADCAL_SPECIES_INDEX_D
            INTEGER :: I,J,K,IC

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR_IN) RETURN
            IF (J < 1 .OR. J > JBAR_IN) RETURN
            IF (K < 1 .OR. K > KBAR_IN) RETURN

            IC = CELL_INDEX_D(I,J,K)
            IF (CELL_SOLID_D(IC)) RETURN

            KAPPA_GAS_D(I,J,K) = GET_KAPPA_GPU_DEVICE(ZZ_D,I,J,K,IBAR_IN,JBAR_IN,KBAR_IN,TMP_D(I,J,K),IBND,N_TRACKED_SPECIES_IN, &
                N_RADCAL_ARRAY_SIZE_IN,N_KAPPA_T_IN,N_KAPPA_Y_IN, &
                RTMPMIN_IN,RTMPMAX_IN,TWO_EPS_IN, &
                MWR_Z_D,Z2RADCAL_SPECIES_D,RADCAL_SPECIES2KAPPA_D, &
                RADCAL_SPECIES_INDEX_D)

        END SUBROUTINE RADIATION_KAPPA_GAS_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_KFST4_SUM_GPU(KFST4_GAS_D,KAPPA_GAS_D,TMP_D,Q_D,CHI_R_D,UII_D, &
            CELL_INDEX_D,CELL_SOLID_D,R_D,DX_D,DY_D,DZ_D, &
            RAD_Q_SUM_D,KFST4_SUM_D,IBAR_IN,JBAR_IN,KBAR_IN, &
            FOUR_SIGMA_IN,QR_CLIP_IN)

            IMPLICIT NONE
            INTEGER, VALUE :: IBAR_IN,JBAR_IN,KBAR_IN
            REAL(EB), VALUE :: FOUR_SIGMA_IN,QR_CLIP_IN
            REAL(EB), DEVICE, DIMENSION(0:IBAR_IN+1,0:JBAR_IN+1,0:*) :: KFST4_GAS_D,KAPPA_GAS_D,TMP_D,Q_D,CHI_R_D,UII_D
            INTEGER, DEVICE, DIMENSION(0:IBAR_IN+1,0:JBAR_IN+1,0:*) :: CELL_INDEX_D
            LOGICAL, DEVICE, DIMENSION(0:*) :: CELL_SOLID_D
            REAL(EB), DEVICE, DIMENSION(0:*) :: R_D,DX_D,DY_D,DZ_D
            REAL(EB), DEVICE, DIMENSION(1:*) :: RAD_Q_SUM_D,KFST4_SUM_D
            INTEGER :: I,J,K,IC,IERR
            REAL(EB) :: VOL,Q_RAD_CELL,KFST4_CELL

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR_IN) RETURN
            IF (J < 1 .OR. J > JBAR_IN) RETURN
            IF (K < 1 .OR. K > KBAR_IN) RETURN

            IC = CELL_INDEX_D(I,J,K)
            IF (CELL_SOLID_D(IC)) RETURN

            KFST4_CELL = KAPPA_GAS_D(I,J,K)*FOUR_SIGMA_IN*TMP_D(I,J,K)**4
            KFST4_GAS_D(I,J,K) = KFST4_CELL

            Q_RAD_CELL = CHI_R_D(I,J,K)*Q_D(I,J,K)
            IF (Q_RAD_CELL > QR_CLIP_IN) THEN
                VOL = R_D(I)*DX_D(I)*DY_D(J)*DZ_D(K)
                IERR = atomicAdd(RAD_Q_SUM_D(1),(Q_RAD_CELL+KAPPA_GAS_D(I,J,K)*UII_D(I,J,K))*VOL)
                IERR = atomicAdd(KFST4_SUM_D(1),KFST4_CELL*VOL)
            ENDIF

        END SUBROUTINE RADIATION_KFST4_SUM_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_RTE_SOURCE_CORRECTION_GPU(KFST4_GAS_D,Q_D,CHI_R_D,CELL_INDEX_D,CELL_SOLID_D, &
            IBAR_IN,JBAR_IN,KBAR_IN,QR_CLIP_IN, &
            RTE_SOURCE_CORRECTION_FACTOR_IN)

            IMPLICIT NONE
            INTEGER, VALUE :: IBAR_IN,JBAR_IN,KBAR_IN
            REAL(EB), VALUE :: QR_CLIP_IN,RTE_SOURCE_CORRECTION_FACTOR_IN
            REAL(EB), DEVICE, DIMENSION(0:IBAR_IN+1,0:JBAR_IN+1,0:*) :: KFST4_GAS_D,Q_D,CHI_R_D
            INTEGER, DEVICE, DIMENSION(0:IBAR_IN+1,0:JBAR_IN+1,0:*) :: CELL_INDEX_D
            LOGICAL, DEVICE, DIMENSION(0:*) :: CELL_SOLID_D
            INTEGER :: I,J,K,IC

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR_IN) RETURN
            IF (J < 1 .OR. J > JBAR_IN) RETURN
            IF (K < 1 .OR. K > KBAR_IN) RETURN

            IC = CELL_INDEX_D(I,J,K)
            IF (CELL_SOLID_D(IC)) RETURN

            IF (CHI_R_D(I,J,K)*Q_D(I,J,K)>QR_CLIP_IN) &
                KFST4_GAS_D(I,J,K) = KFST4_GAS_D(I,J,K)*RTE_SOURCE_CORRECTION_FACTOR_IN

        END SUBROUTINE RADIATION_RTE_SOURCE_CORRECTION_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_EXTCOE_GPU(EXTCOE_D,KAPPA_GAS_D,KAPPA_PART_D,SCAEFF_D,SCAEFF_G_D, &
            CELL_INDEX_D,CELL_SOLID_D,IBAR_IN,JBAR_IN,KBAR_IN,RSA_RAT_IN)

            IMPLICIT NONE
            INTEGER, VALUE :: IBAR_IN,JBAR_IN,KBAR_IN
            REAL(EB), VALUE :: RSA_RAT_IN
            REAL(EB), DEVICE, DIMENSION(0:IBAR_IN+1,0:JBAR_IN+1,0:*) :: EXTCOE_D,KAPPA_GAS_D,KAPPA_PART_D,SCAEFF_D,SCAEFF_G_D
            INTEGER, DEVICE, DIMENSION(0:IBAR_IN+1,0:JBAR_IN+1,0:*) :: CELL_INDEX_D
            LOGICAL, DEVICE, DIMENSION(0:*) :: CELL_SOLID_D
            INTEGER :: I,J,K,IC

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR_IN) RETURN
            IF (J < 1 .OR. J > JBAR_IN) RETURN
            IF (K < 1 .OR. K > KBAR_IN) RETURN

            IC = CELL_INDEX_D(I,J,K)
            IF (CELL_SOLID_D(IC)) RETURN

            EXTCOE_D(I,J,K) = KAPPA_GAS_D(I,J,K) + KAPPA_PART_D(I,J,K) + &
                (SCAEFF_D(I,J,K)+SCAEFF_G_D(I,J,K))*RSA_RAT_IN

        END SUBROUTINE RADIATION_EXTCOE_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_UIIOLD_ZERO_GPU(UIIOLD_D,UII_D,IBP1_IN,JBP1_IN,KBP1_IN)

            IMPLICIT NONE
            INTEGER, VALUE :: IBP1_IN,JBP1_IN,KBP1_IN
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: UIIOLD_D,UII_D
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBP1_IN) RETURN
            IF (J < 0 .OR. J > JBP1_IN) RETURN
            IF (K < 0 .OR. K > KBP1_IN) RETURN

            UIIOLD_D(I,J,K) = UII_D(I,J,K)
            UII_D(I,J,K) = 0._EB

        END SUBROUTINE RADIATION_UIIOLD_ZERO_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_OUTRAD_W_GPU(OUTRAD_W_D, WALL_TOPO, WALL_Q_RAD_OUT_D, N_WALL_TOTAL, BBFA_IN, &
            RPI_IN, SIGMA_IN, TMPA_IN, NULL_BOUNDARY_IN, OPEN_BOUNDARY_IN)
            IMPLICIT NONE
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_WALL_TOTAL,NULL_BOUNDARY_IN,OPEN_BOUNDARY_IN
            REAL(EB), VALUE :: BBFA_IN,RPI_IN,SIGMA_IN,TMPA_IN
            REAL(EB), DEVICE, DIMENSION(*) :: OUTRAD_W_D
            REAL(EB), DEVICE, DIMENSION(0:*) :: WALL_Q_RAD_OUT_D
            INTEGER :: IW

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X

            IF (IW < 1 .OR. IW > N_WALL_TOTAL) RETURN
            IF (WALL_TOPO(IW)%BOUNDARY_TYPE==NULL_BOUNDARY_IN) RETURN

            IF (WALL_TOPO(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY_IN) THEN
                OUTRAD_W_D(IW) = BBFA_IN*RPI_IN*SIGMA_IN*TMPA_IN**4
            ELSE
                OUTRAD_W_D(IW) = RPI_IN*WALL_Q_RAD_OUT_D(IW)
            ENDIF

        END SUBROUTINE RADIATION_OUTRAD_W_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_INRAD_W_GPU(INRAD_W_D, WALL_TOPO, WALL_BR_INDEX_D, DLN_D, &
            BOUNDARY_RADIA_ILW_D, IBND_IN, N_WALL_TOTAL, NRA_IN, SOLID_BOUNDARY_IN, NBR_IN, NBAND_IN)
            IMPLICIT NONE
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: IBND_IN,N_WALL_TOTAL,NRA_IN,SOLID_BOUNDARY_IN,NBR_IN,NBAND_IN
            REAL(EB), DEVICE, DIMENSION(*) :: INRAD_W_D
            INTEGER, DEVICE, DIMENSION(0:*) :: WALL_BR_INDEX_D
            REAL(EB), DEVICE, DIMENSION(-3:3,1:*) :: DLN_D
            REAL(EB), DEVICE, DIMENSION(1:NBR_IN,1:NBAND_IN,1:*) :: BOUNDARY_RADIA_ILW_D
            INTEGER :: IW,N,IOR,IBR
            REAL(EB) :: SUM_IN,DLN_VAL

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X

            IF (IW < 1 .OR. IW > N_WALL_TOTAL) RETURN
            IF (WALL_TOPO(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY_IN) RETURN

            IOR = WALL_TOPO(IW)%IOR
            IBR = WALL_BR_INDEX_D(IW)
            IF (IBR <= 0) RETURN

            SUM_IN = 0._EB
            DO N=1,NRA_IN
                DLN_VAL = DLN_D(IOR,N)
                IF (DLN_VAL < 0._EB) SUM_IN = SUM_IN - DLN_VAL*BOUNDARY_RADIA_ILW_D(IBR,IBND_IN,N)
            ENDDO
            INRAD_W_D(IW) = SUM_IN

        END SUBROUTINE RADIATION_INRAD_W_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_INTERPOLATED_WALL_GPU(IL_D,WALL_TOPO,EXTERNAL_WALL_NOM_D,DLN_D, &
            WALL_NIC_MIN_D,WALL_NIC_MAX_D,OM_IL_R_D,N_IN,IBND_IN,N_EXTERNAL_WALL_CELLS_IN,NOM_FILTER,NIC_R_IN, &
            INTERPOLATED_BOUNDARY_IN,IL_NI,IL_NJ,IL_LB_I,IL_LB_J,IL_LB_K,OM_N2,OM_LB_1,OM_LB_2,OM_LB_3, &
            DLN_N1,DLN_LB_1,DLN_LB_2,WALL_LB)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(*)
            INTEGER, VALUE :: N_IN,IBND_IN,N_EXTERNAL_WALL_CELLS_IN,NOM_FILTER,NIC_R_IN,INTERPOLATED_BOUNDARY_IN
            INTEGER, VALUE :: IL_NI,IL_NJ,IL_LB_I,IL_LB_J,IL_LB_K,OM_N2,OM_LB_1,OM_LB_2,OM_LB_3
            INTEGER, VALUE :: DLN_N1,DLN_LB_1,DLN_LB_2,WALL_LB
            REAL(EB), DEVICE, DIMENSION(*) :: IL_D,OM_IL_R_D,DLN_D
            INTEGER, DEVICE, DIMENSION(*) :: EXTERNAL_WALL_NOM_D,WALL_NIC_MIN_D,WALL_NIC_MAX_D
            INTEGER :: IW,IW_WALL,LL,NIC_MIN,NIC_MAX,II,JJ,KK,IOR,IL_IDX,OM_IDX,DLN_IDX
            REAL(EB) :: IL_SUM

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS_IN) RETURN
            IW_WALL = IW - WALL_LB + 1
            IF (WALL_TOPO(IW_WALL)%BOUNDARY_TYPE /= INTERPOLATED_BOUNDARY_IN) RETURN
            IF (EXTERNAL_WALL_NOM_D(IW) /= NOM_FILTER) RETURN
            IOR = WALL_TOPO(IW_WALL)%IOR
            DLN_IDX = 1 + (IOR-DLN_LB_1) + DLN_N1*(N_IN-DLN_LB_2)
            IF (DLN_D(DLN_IDX) < 0._EB) RETURN

            NIC_MIN = WALL_NIC_MIN_D(IW)
            NIC_MAX = WALL_NIC_MAX_D(IW)
            IF (NIC_MIN < 1 .OR. NIC_MAX > NIC_R_IN .OR. NIC_MIN > NIC_MAX) RETURN

            IL_SUM = 0._EB
            DO LL=NIC_MIN,NIC_MAX
                OM_IDX = 1 + (LL-OM_LB_1) + NIC_R_IN*((N_IN-OM_LB_2) + OM_N2*(IBND_IN-OM_LB_3))
                IL_SUM = IL_SUM + OM_IL_R_D(OM_IDX)
            ENDDO
            II = WALL_TOPO(IW_WALL)%II
            JJ = WALL_TOPO(IW_WALL)%JJ
            KK = WALL_TOPO(IW_WALL)%KK
            IL_IDX = 1 + (II-IL_LB_I) + IL_NI*((JJ-IL_LB_J) + IL_NJ*(KK-IL_LB_K))
            IL_D(IL_IDX) = IL_SUM/REAL(NIC_MAX-NIC_MIN+1,EB)

        END SUBROUTINE RADIATION_INTERPOLATED_WALL_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_PACK_IL_S_GPU(IL_D,OM_IL_S_D,IIO_S_D,JJO_S_D,KKO_S_D, &
            N_IN,IBND_IN,NIC_S_IN,IBP1_IN,JBP1_IN,KBP1_IN,OM_N1,OM_N2)
            INTEGER, VALUE :: N_IN,IBND_IN,NIC_S_IN,IBP1_IN,JBP1_IN,KBP1_IN,OM_N1,OM_N2
            REAL(EB), DEVICE :: IL_D(0:IBP1_IN,0:JBP1_IN,0:KBP1_IN)
            REAL(EB), DEVICE, DIMENSION(*) :: OM_IL_S_D
            INTEGER, DEVICE, DIMENSION(*) :: IIO_S_D,JJO_S_D,KKO_S_D
            INTEGER :: LL,OM_IDX

            LL = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (LL < 1 .OR. LL > NIC_S_IN) RETURN
            OM_IDX = 1 + (LL-1) + OM_N1*((N_IN-1) + OM_N2*(IBND_IN-1))
            OM_IL_S_D(OM_IDX) = IL_D(IIO_S_D(LL),JJO_S_D(LL),KKO_S_D(LL))

        END SUBROUTINE RADIATION_PACK_IL_S_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_ZERO_UIID_ANGLE_GPU(UIID_D,ANGLE_INC_COUNTER_IN,IBP1_IN,JBP1_IN,KBP1_IN,UIIDIM_IN)

            IMPLICIT NONE
            INTEGER, VALUE :: ANGLE_INC_COUNTER_IN,IBP1_IN,JBP1_IN,KBP1_IN,UIIDIM_IN
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:KBP1_IN,1:*) :: UIID_D
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBP1_IN) RETURN
            IF (J < 0 .OR. J > JBP1_IN) RETURN
            IF (K < 0 .OR. K > KBP1_IN) RETURN
            IF (ANGLE_INC_COUNTER_IN < 1 .OR. ANGLE_INC_COUNTER_IN > UIIDIM_IN) RETURN

            UIID_D(I,J,K,ANGLE_INC_COUNTER_IN) = 0._EB

        END SUBROUTINE RADIATION_ZERO_UIID_ANGLE_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_ZERO_OPEN_ILW_ANGLE_GPU(BOUNDARY_RADIA_ILW_D, WALL_TOPO, WALL_BR_INDEX_D, &
            IBND_IN, ANGLE_INC_COUNTER_IN, N_EXTERNAL_WALL_CELLS_IN, OPEN_BOUNDARY_IN, NBR_IN, NBAND_IN)
            IMPLICIT NONE
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: IBND_IN,ANGLE_INC_COUNTER_IN,N_EXTERNAL_WALL_CELLS_IN,OPEN_BOUNDARY_IN,NBR_IN,NBAND_IN
            REAL(EB), DEVICE, DIMENSION(1:NBR_IN,1:NBAND_IN,1:*) :: BOUNDARY_RADIA_ILW_D
            INTEGER, DEVICE, DIMENSION(0:*) :: WALL_BR_INDEX_D
            INTEGER :: IW,IBR

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X

            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS_IN) RETURN
            IF (WALL_TOPO(IW)%BOUNDARY_TYPE/=OPEN_BOUNDARY_IN) RETURN

            IBR = WALL_BR_INDEX_D(IW)
            IF (IBR <= 0) RETURN

            BOUNDARY_RADIA_ILW_D(IBR,IBND_IN,ANGLE_INC_COUNTER_IN) = 0._EB

        END SUBROUTINE RADIATION_ZERO_OPEN_ILW_ANGLE_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_INIT_IL_GPU(IL_D,IL_INIT,IBP1_IN,JBP1_IN,KBP1_IN)

            IMPLICIT NONE
            INTEGER, VALUE :: IBP1_IN,JBP1_IN,KBP1_IN
            REAL(EB), VALUE :: IL_INIT
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: IL_D
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBP1_IN) RETURN
            IF (J < 0 .OR. J > JBP1_IN) RETURN
            IF (K < 0 .OR. K > KBP1_IN) RETURN

            IL_D(I,J,K) = IL_INIT

        END SUBROUTINE RADIATION_INIT_IL_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_WALL_LOOP1_GPU(IL_D, CELL_ILW_D, BOUNDARY_RADIA_ILW_D, OUTRAD_W_D, INRAD_W_D, &
            WALL_TOPO, WALL_EMISSIVITY_D, WALL_BR_INDEX_D, CELL_INDEX_D, DLN_D, DLM_D, N_IN, IBND_IN, N_WALL_TOTAL, &
            IL_OPEN_IN, RPI_IN, TWO_D_IN, NULL_BOUNDARY_IN, OPEN_BOUNDARY_IN, MIRROR_BOUNDARY_IN, &
            INTERPOLATED_BOUNDARY_IN, IBP1_IN, JBP1_IN, CELL_ILW_SIZE_IN, NBR_IN, NBAND_IN, NRA_IN)
            IMPLICIT NONE
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_IN,IBND_IN,N_WALL_TOTAL,IBP1_IN,JBP1_IN,CELL_ILW_SIZE_IN,NBR_IN,NBAND_IN,NRA_IN
            INTEGER, VALUE :: NULL_BOUNDARY_IN,OPEN_BOUNDARY_IN,MIRROR_BOUNDARY_IN,INTERPOLATED_BOUNDARY_IN
            REAL(EB), VALUE :: IL_OPEN_IN,RPI_IN
            LOGICAL, VALUE :: TWO_D_IN
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: IL_D
            REAL(EB), DEVICE, DIMENSION(1:CELL_ILW_SIZE_IN,1:*) :: CELL_ILW_D
            REAL(EB), DEVICE, DIMENSION(1:NBR_IN,1:NBAND_IN,1:*) :: BOUNDARY_RADIA_ILW_D
            REAL(EB), DEVICE, DIMENSION(*) :: OUTRAD_W_D,INRAD_W_D
            REAL(EB), DEVICE, DIMENSION(0:*) :: WALL_EMISSIVITY_D
            INTEGER, DEVICE, DIMENSION(0:*) :: WALL_BR_INDEX_D
            INTEGER, DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: CELL_INDEX_D
            REAL(EB), DEVICE, DIMENSION(-3:3,1:*) :: DLN_D
            INTEGER, DEVICE, DIMENSION(1:NRA_IN,1:*) :: DLM_D
            INTEGER :: IW,BOUNDARY_TYPE,IOR,II,JJ,KK,IIG,JJG,KKG,IBR,IC,ABS_IOR
            REAL(EB) :: ILW_VAL

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X

            IF (IW < 1 .OR. IW > N_WALL_TOTAL) RETURN

            BOUNDARY_TYPE = WALL_TOPO(IW)%BOUNDARY_TYPE
            IF (BOUNDARY_TYPE==NULL_BOUNDARY_IN) RETURN

            IOR = WALL_TOPO(IW)%IOR
            IF (DLN_D(IOR,N_IN) < 0._EB) RETURN
            IF (TWO_D_IN .AND. ABS(IOR)==2) RETURN

            II  = WALL_TOPO(IW)%II
            JJ  = WALL_TOPO(IW)%JJ
            KK  = WALL_TOPO(IW)%KK
            IBR = WALL_BR_INDEX_D(IW)
            IF (IBR <= 0) RETURN

            IF (BOUNDARY_TYPE==OPEN_BOUNDARY_IN) THEN
                IL_D(II,JJ,KK) = IL_OPEN_IN
            ELSEIF (BOUNDARY_TYPE==MIRROR_BOUNDARY_IN) THEN
                ABS_IOR = ABS(IOR)
                ILW_VAL = BOUNDARY_RADIA_ILW_D(IBR,IBND_IN,DLM_D(N_IN,ABS_IOR))
                BOUNDARY_RADIA_ILW_D(IBR,IBND_IN,N_IN) = ILW_VAL
                IL_D(II,JJ,KK) = ILW_VAL
            ELSEIF (BOUNDARY_TYPE==INTERPOLATED_BOUNDARY_IN) THEN
                CONTINUE
            ELSE
                IIG = WALL_TOPO(IW)%IIG
                JJG = WALL_TOPO(IW)%JJG
                KKG = WALL_TOPO(IW)%KKG
                ABS_IOR = ABS(IOR)
                ILW_VAL = OUTRAD_W_D(IW) + RPI_IN*(1._EB-WALL_EMISSIVITY_D(IW))*INRAD_W_D(IW)
                BOUNDARY_RADIA_ILW_D(IBR,IBND_IN,N_IN) = ILW_VAL
                IC = CELL_INDEX_D(IIG,JJG,KKG)
                CELL_ILW_D(IC,ABS_IOR) = ILW_VAL
            ENDIF

        END SUBROUTINE RADIATION_WALL_LOOP1_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_SWEEP_3D_SLICE_GPU( &
            IL_D,IL_UP_D,CELL_ILW_D,EXTCOE_D,KFST4_GAS_D,KFST4_PART_D,SCAEFF_D,SCAEFF_G_D,UIIOLD_D, &
            CELL_INDEX_D,CELL_SOLID_D,IJK_SLICE_MAP_D,DX_D,DY_D,DZ_D,DLX_D,DLY_D,DLZ_D,RSA_D, &
            N_IN,BASE_IN,M_IJK_IN,ISTEP_IN,JSTEP_IN,KSTEP_IN,RSA_RAT_IN,RFPI_IN,SOLID_PARTICLES_IN, &
            IBP1_IN,JBP1_IN,KBP1_IN,CELL_LB_IN,CELL_UB_IN,CELL_ILW_SIZE_IN,NCELL_IN)

            IMPLICIT NONE
            INTEGER, VALUE :: N_IN,BASE_IN,M_IJK_IN,ISTEP_IN,JSTEP_IN,KSTEP_IN
            INTEGER, VALUE :: IBP1_IN,JBP1_IN,KBP1_IN,CELL_LB_IN,CELL_UB_IN,CELL_ILW_SIZE_IN,NCELL_IN
            REAL(EB), VALUE :: RSA_RAT_IN,RFPI_IN
            LOGICAL, VALUE :: SOLID_PARTICLES_IN
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: IL_D,IL_UP_D,EXTCOE_D,KFST4_GAS_D,KFST4_PART_D, &
                SCAEFF_D,SCAEFF_G_D,UIIOLD_D
            REAL(EB), DEVICE, DIMENSION(1:CELL_ILW_SIZE_IN,1:*) :: CELL_ILW_D
            INTEGER, DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: CELL_INDEX_D
            LOGICAL, DEVICE, DIMENSION(CELL_LB_IN:*) :: CELL_SOLID_D
            INTEGER, DEVICE, DIMENSION(1:NCELL_IN,1:3) :: IJK_SLICE_MAP_D
            REAL(EB), DEVICE, DIMENSION(0:*) :: DX_D,DY_D,DZ_D
            REAL(EB), DEVICE, DIMENSION(*) :: DLX_D,DLY_D,DLZ_D,RSA_D
            INTEGER :: IJK,IDX,I,J,K,IC
            REAL(EB) :: AY1,AX,VC1,AZ1,ILXU,ILYU,ILZU,VC,AY,AZ,AXU,AYU,AZU,AXD,AYD,AZD
            REAL(EB) :: AILFU,AFD,A_SUM,AIU_SUM,RAP

            IJK = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X

            IF (IJK < 1 .OR. IJK > M_IJK_IN) RETURN

            IDX = BASE_IN + IJK - 1
            I = IJK_SLICE_MAP_D(IDX,1)
            J = IJK_SLICE_MAP_D(IDX,2)
            K = IJK_SLICE_MAP_D(IDX,3)
            IF (ISTEP_IN < 0) I = IBP1_IN - I
            IF (JSTEP_IN < 0) J = JBP1_IN - J
            IF (KSTEP_IN < 0) K = KBP1_IN - K

            AY1 = DZ_D(K) * ABS(DLY_D(N_IN))
            AX  = DY_D(J) * DZ_D(K) * ABS(DLX_D(N_IN))
            VC1 = DY_D(J) * DZ_D(K)
            AZ1 = DY_D(J) * ABS(DLZ_D(N_IN))

            ILXU = IL_D(I-ISTEP_IN,J,K)
            ILYU = IL_D(I,J-JSTEP_IN,K)
            ILZU = IL_D(I,J,K-KSTEP_IN)

            IC = CELL_INDEX_D(I,J,K)
            IF (IC /= 0) THEN
                IF (CELL_SOLID_D(IC)) RETURN
                IF (CELL_ILW_D(IC,1)>-1.E6_EB) ILXU = CELL_ILW_D(IC,1)
                IF (CELL_ILW_D(IC,2)>-1.E6_EB) ILYU = CELL_ILW_D(IC,2)
                IF (CELL_ILW_D(IC,3)>-1.E6_EB) ILZU = CELL_ILW_D(IC,3)
            ENDIF

            VC = DX_D(I) * VC1
            AY = DX_D(I) * AY1
            AZ = DX_D(I) * AZ1

            AXU = AX; AXD = AX
            AYU = AY; AYD = AY
            AZU = AZ; AZD = AZ

            AILFU = 0._EB
            AFD = 0._EB

            A_SUM = AXD + AYD + AZD + AFD
            AIU_SUM = AXU*ILXU + AYU*ILYU + AZU*ILZU + AILFU
            IF (SOLID_PARTICLES_IN) IL_UP_D(I,J,K) = MAX(0._EB,AIU_SUM/A_SUM)

            RAP = 1._EB/(A_SUM + EXTCOE_D(I,J,K)*VC*RSA_D(N_IN))
            IL_D(I,J,K) = MAX(0._EB, RAP * (AIU_SUM + VC*RSA_D(N_IN)*RFPI_IN * &
                ( KFST4_GAS_D(I,J,K) + KFST4_PART_D(I,J,K) + RSA_RAT_IN * &
                (SCAEFF_D(I,J,K)+SCAEFF_G_D(I,J,K))*UIIOLD_D(I,J,K) ) ) )

        END SUBROUTINE RADIATION_SWEEP_3D_SLICE_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_WALL_INCOMING_GPU(INRAD_W_D,BOUNDARY_RADIA_ILW_D,IL_D,WALL_TOPO, &
            WALL_BR_INDEX_D,DLN_D,N_IN,IBND_IN,ANGLE_INC_COUNTER_IN,N_WALL_TOTAL,TWO_D_IN,CYLINDRICAL_IN, &
            NULL_BOUNDARY_IN,OPEN_BOUNDARY_IN,IBP1_IN,JBP1_IN,NBR_IN,NBAND_IN)
            IMPLICIT NONE
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_IN,IBND_IN,ANGLE_INC_COUNTER_IN,N_WALL_TOTAL,NULL_BOUNDARY_IN,OPEN_BOUNDARY_IN
            INTEGER, VALUE :: IBP1_IN,JBP1_IN,NBR_IN,NBAND_IN
            LOGICAL, VALUE :: TWO_D_IN,CYLINDRICAL_IN
            REAL(EB), DEVICE, DIMENSION(*) :: INRAD_W_D
            REAL(EB), DEVICE, DIMENSION(1:NBR_IN,1:NBAND_IN,1:*) :: BOUNDARY_RADIA_ILW_D
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: IL_D
            INTEGER, DEVICE, DIMENSION(0:*) :: WALL_BR_INDEX_D
            REAL(EB), DEVICE, DIMENSION(-3:3,1:*) :: DLN_D
            INTEGER :: IW,IOR,IBR,IIG,JJG,KKG
            REAL(EB) :: OLD_ILW,NEW_ILW,DLN_VAL

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X

            IF (IW < 1 .OR. IW > N_WALL_TOTAL) RETURN
            IF (WALL_TOPO(IW)%BOUNDARY_TYPE==NULL_BOUNDARY_IN) RETURN

            IOR = WALL_TOPO(IW)%IOR
            DLN_VAL = DLN_D(IOR,N_IN)
            IF (DLN_VAL >= 0._EB) RETURN

            IBR = WALL_BR_INDEX_D(IW)
            IF (IBR <= 0) RETURN

            IIG = WALL_TOPO(IW)%IIG
            JJG = WALL_TOPO(IW)%JJG
            KKG = WALL_TOPO(IW)%KKG
            IF (WALL_TOPO(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY_IN) THEN
                BOUNDARY_RADIA_ILW_D(IBR,IBND_IN,ANGLE_INC_COUNTER_IN) = &
                    BOUNDARY_RADIA_ILW_D(IBR,IBND_IN,ANGLE_INC_COUNTER_IN) - DLN_VAL*IL_D(IIG,JJG,KKG)
            ELSE
                IF (TWO_D_IN .AND. (.NOT.CYLINDRICAL_IN) .AND. ABS(IOR)==2) RETURN
                OLD_ILW = BOUNDARY_RADIA_ILW_D(IBR,IBND_IN,N_IN)
                NEW_ILW = IL_D(IIG,JJG,KKG)
                INRAD_W_D(IW) = INRAD_W_D(IW) + DLN_VAL*OLD_ILW - DLN_VAL*NEW_ILW
                BOUNDARY_RADIA_ILW_D(IBR,IBND_IN,N_IN) = NEW_ILW
            ENDIF

        END SUBROUTINE RADIATION_WALL_INCOMING_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_ACCUMULATE_UIID_GPU(UIID_D,IL_D,CELL_INDEX_D,CELL_SOLID_D, &
            UIID_INDEX_IN,N_IN,IBP1_IN,JBP1_IN,KBP1_IN, &
            UIIDIM_IN,WEIGH_CYL_IN,RSA_D)

            IMPLICIT NONE
            INTEGER, VALUE :: UIID_INDEX_IN,N_IN,IBP1_IN,JBP1_IN,KBP1_IN,UIIDIM_IN
            REAL(EB), VALUE :: WEIGH_CYL_IN
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:KBP1_IN,1:*) :: UIID_D
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: IL_D
            INTEGER, DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: CELL_INDEX_D
            LOGICAL, DEVICE, DIMENSION(0:*) :: CELL_SOLID_D
            REAL(EB), DEVICE, DIMENSION(*) :: RSA_D
            INTEGER :: I,J,K,IC

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBP1_IN) RETURN
            IF (J < 0 .OR. J > JBP1_IN) RETURN
            IF (K < 0 .OR. K > KBP1_IN) RETURN
            IF (UIID_INDEX_IN < 1 .OR. UIID_INDEX_IN > UIIDIM_IN) RETURN

            IC = CELL_INDEX_D(I,J,K)
            IF (IC /= 0) THEN
                IF (CELL_SOLID_D(IC)) RETURN
            ENDIF

            UIID_D(I,J,K,UIID_INDEX_IN) = UIID_D(I,J,K,UIID_INDEX_IN) + WEIGH_CYL_IN*RSA_D(N_IN)*IL_D(I,J,K)

        END SUBROUTINE RADIATION_ACCUMULATE_UIID_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_SOLID_WALL_Q_RAD_IN_GPU(WALL_Q_RAD_IN_D, INRAD_W_D, WALL_TOPO, WALL_SF_INDEX_D, &
            WALL_EMISSIVITY_D,SF_EXTERNAL_FLUX_D,SF_RAMP_EFLUX_INDEX_D,SF_RAMP_EFLUX_TAU_D,WALL_T_IGN_D, &
            RAMP_T_MIN_D,RAMP_SPAN_D,RAMP_RDT_D,RAMP_NPOINTS_D,RAMP_OFFSET_D,RAMP_INTERP_D,T_IN,T_BEGIN_IN, &
            N_WALL_TOTAL, BBFA_IN, SOLID_BOUNDARY_IN)
            IMPLICIT NONE
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_WALL_TOTAL,SOLID_BOUNDARY_IN
            REAL(EB), VALUE :: BBFA_IN,T_IN,T_BEGIN_IN
            INTEGER, DEVICE, DIMENSION(0:*) :: WALL_SF_INDEX_D,SF_RAMP_EFLUX_INDEX_D
            REAL(EB), DEVICE, DIMENSION(0:*) :: WALL_Q_RAD_IN_D,WALL_EMISSIVITY_D,SF_EXTERNAL_FLUX_D, &
                SF_RAMP_EFLUX_TAU_D,WALL_T_IGN_D
            REAL(EB), DEVICE, DIMENSION(0:*) :: RAMP_T_MIN_D,RAMP_SPAN_D,RAMP_RDT_D
            INTEGER, DEVICE, DIMENSION(0:*) :: RAMP_NPOINTS_D,RAMP_OFFSET_D
            REAL(EB), DEVICE, DIMENSION(1:*) :: RAMP_INTERP_D
            REAL(EB), DEVICE, DIMENSION(*) :: INRAD_W_D
            INTEGER :: IW,SF_IDX,RAMP_INDEX
            REAL(EB) :: TSI,EFLUX

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X

            IF (IW < 1 .OR. IW > N_WALL_TOTAL) RETURN
            IF (WALL_TOPO(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY_IN) RETURN

            SF_IDX = WALL_SF_INDEX_D(IW)
            IF (ABS(T_BEGIN_IN) <= SPACING(WALL_T_IGN_D(IW))) THEN
                TSI = T_IN
            ELSE
                TSI = T_IN-WALL_T_IGN_D(IW)
            ENDIF
            RAMP_INDEX = SF_RAMP_EFLUX_INDEX_D(SF_IDX)
            EFLUX = EVALUATE_RAMP_DEVICE(TSI,RAMP_INDEX,SF_RAMP_EFLUX_TAU_D(SF_IDX),RAMP_T_MIN_D,RAMP_SPAN_D, &
                                          RAMP_RDT_D,RAMP_NPOINTS_D,RAMP_OFFSET_D,RAMP_INTERP_D) * SF_EXTERNAL_FLUX_D(SF_IDX)
            WALL_Q_RAD_IN_D(IW) = WALL_Q_RAD_IN_D(IW) + WALL_EMISSIVITY_D(IW)*(INRAD_W_D(IW)+BBFA_IN*EFLUX)

        END SUBROUTINE RADIATION_SOLID_WALL_Q_RAD_IN_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_SUM_UIID_GPU(UII_D,UIID_D,IBP1_IN,JBP1_IN,KBP1_IN,UIIDIM_IN)

            IMPLICIT NONE
            INTEGER, VALUE :: IBP1_IN,JBP1_IN,KBP1_IN,UIIDIM_IN
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: UII_D
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:KBP1_IN,1:*) :: UIID_D
            INTEGER :: I,J,K,N
            REAL(EB) :: SUM_UIID

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBP1_IN) RETURN
            IF (J < 0 .OR. J > JBP1_IN) RETURN
            IF (K < 0 .OR. K > KBP1_IN) RETURN

            SUM_UIID = 0._EB
            DO N=1,UIIDIM_IN
                SUM_UIID = SUM_UIID + UIID_D(I,J,K,N)
            ENDDO
            UII_D(I,J,K) = SUM_UIID

        END SUBROUTINE RADIATION_SUM_UIID_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_OPEN_Q_RAD_IN_GPU(WALL_Q_RAD_IN_D, BOUNDARY_RADIA_ILW_D, WALL_TOPO, &
            WALL_BR_INDEX_D, N_EXTERNAL_WALL_CELLS_IN, N_BANDS_IN, NRA_IN, OPEN_BOUNDARY_IN, NBR_IN)
            IMPLICIT NONE
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS_IN,N_BANDS_IN,NRA_IN,OPEN_BOUNDARY_IN,NBR_IN
            REAL(EB), DEVICE, DIMENSION(0:*) :: WALL_Q_RAD_IN_D
            REAL(EB), DEVICE, DIMENSION(1:NBR_IN,1:N_BANDS_IN,1:*) :: BOUNDARY_RADIA_ILW_D
            INTEGER, DEVICE, DIMENSION(0:*) :: WALL_BR_INDEX_D
            INTEGER :: IW,IBR,IB,N
            REAL(EB) :: SUM_ILW

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X

            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS_IN) RETURN
            IF (WALL_TOPO(IW)%BOUNDARY_TYPE/=OPEN_BOUNDARY_IN) RETURN

            IBR = WALL_BR_INDEX_D(IW)
            IF (IBR <= 0) RETURN

            SUM_ILW = 0._EB
            DO IB=1,N_BANDS_IN
                DO N=1,NRA_IN
                    SUM_ILW = SUM_ILW + BOUNDARY_RADIA_ILW_D(IBR,IB,N)
                ENDDO
            ENDDO
            WALL_Q_RAD_IN_D(IW) = SUM_ILW

        END SUBROUTINE RADIATION_OPEN_Q_RAD_IN_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_QR_WIDEBAND_GPU(QR_D,QR_W_D,RADIATION_EMISSION_D,RADIATION_ABSORPTION_D, &
            KAPPA_GAS_D,KFST4_GAS_D,KAPPA_PART_D,KFST4_PART_D,UIID_D, &
            CELL_INDEX_D,CELL_SOLID_D,IBND_IN,IBP1_IN,JBP1_IN,KBP1_IN, &
            UIIDIM_IN,STORE_RADIATION_TERMS_IN,DO_QR_W_IN)

            IMPLICIT NONE
            INTEGER, VALUE :: IBND_IN,IBP1_IN,JBP1_IN,KBP1_IN,UIIDIM_IN
            LOGICAL, VALUE :: STORE_RADIATION_TERMS_IN,DO_QR_W_IN
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: QR_D,QR_W_D,RADIATION_EMISSION_D,RADIATION_ABSORPTION_D
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: KAPPA_GAS_D,KFST4_GAS_D,KAPPA_PART_D,KFST4_PART_D
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:KBP1_IN,1:*) :: UIID_D
            INTEGER, DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: CELL_INDEX_D
            LOGICAL, DEVICE, DIMENSION(0:*) :: CELL_SOLID_D
            INTEGER :: I,J,K,IC
            REAL(EB) :: ABS_TERM

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBP1_IN) RETURN
            IF (J < 0 .OR. J > JBP1_IN) RETURN
            IF (K < 0 .OR. K > KBP1_IN) RETURN
            IF (IBND_IN < 1 .OR. IBND_IN > UIIDIM_IN) RETURN

            IC = CELL_INDEX_D(I,J,K)
            IF (IC /= 0) THEN
                IF (CELL_SOLID_D(IC)) RETURN
            ENDIF

            ABS_TERM = KAPPA_GAS_D(I,J,K)*UIID_D(I,J,K,IBND_IN)
            QR_D(I,J,K) = QR_D(I,J,K) + ABS_TERM - KFST4_GAS_D(I,J,K)
            IF (STORE_RADIATION_TERMS_IN) THEN
                RADIATION_EMISSION_D(I,J,K) = RADIATION_EMISSION_D(I,J,K) + KFST4_GAS_D(I,J,K)
                RADIATION_ABSORPTION_D(I,J,K) = RADIATION_ABSORPTION_D(I,J,K) + ABS_TERM
            ENDIF
            IF (DO_QR_W_IN) QR_W_D(I,J,K) = QR_W_D(I,J,K) + KAPPA_PART_D(I,J,K)*UIID_D(I,J,K,IBND_IN) - KFST4_PART_D(I,J,K)

        END SUBROUTINE RADIATION_QR_WIDEBAND_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_QR_GRAY_GPU(QR_D,QR_W_D,RADIATION_EMISSION_D,RADIATION_ABSORPTION_D, &
            KAPPA_GAS_D,KFST4_GAS_D,KAPPA_PART_D,KFST4_PART_D,UII_D, &
            CELL_INDEX_D,CELL_SOLID_D,IBP1_IN,JBP1_IN,KBP1_IN, &
            STORE_RADIATION_TERMS_IN,DO_QR_W_IN)

            IMPLICIT NONE
            INTEGER, VALUE :: IBP1_IN,JBP1_IN,KBP1_IN
            LOGICAL, VALUE :: STORE_RADIATION_TERMS_IN,DO_QR_W_IN
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: QR_D,QR_W_D,RADIATION_EMISSION_D,RADIATION_ABSORPTION_D
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: KAPPA_GAS_D,KFST4_GAS_D,KAPPA_PART_D,KFST4_PART_D,UII_D
            INTEGER, DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: CELL_INDEX_D
            LOGICAL, DEVICE, DIMENSION(0:*) :: CELL_SOLID_D
            INTEGER :: I,J,K,IC
            REAL(EB) :: ABS_TERM

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBP1_IN) RETURN
            IF (J < 0 .OR. J > JBP1_IN) RETURN
            IF (K < 0 .OR. K > KBP1_IN) RETURN

            IC = CELL_INDEX_D(I,J,K)
            IF (IC /= 0) THEN
                IF (CELL_SOLID_D(IC)) RETURN
            ENDIF

            ABS_TERM = KAPPA_GAS_D(I,J,K)*UII_D(I,J,K)
            QR_D(I,J,K) = ABS_TERM - KFST4_GAS_D(I,J,K)
            IF (DO_QR_W_IN) QR_W_D(I,J,K) = QR_W_D(I,J,K) + KAPPA_PART_D(I,J,K)*UII_D(I,J,K) - KFST4_PART_D(I,J,K)
            IF (STORE_RADIATION_TERMS_IN) THEN
                RADIATION_EMISSION_D(I,J,K) = KFST4_GAS_D(I,J,K)
                RADIATION_ABSORPTION_D(I,J,K) = ABS_TERM
            ENDIF

        END SUBROUTINE RADIATION_QR_GRAY_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE RADIATION_OFF_GPU(QR_D,CHI_R_D,Q_D,IBP1_IN,JBP1_IN,KBP1_IN)

            INTEGER, VALUE :: IBP1_IN,JBP1_IN,KBP1_IN
            REAL(EB), DEVICE, DIMENSION(0:IBP1_IN,0:JBP1_IN,0:*) :: QR_D,CHI_R_D,Q_D
            INTEGER :: I,J,K

            I = (blockIdx%x - 1) * blockDim%x + threadIdx%x - 1
            J = (blockIdx%y - 1) * blockDim%y + threadIdx%y - 1
            K = (blockIdx%z - 1) * blockDim%z + threadIdx%z - 1

            IF (I < 0 .OR. I > IBP1_IN) RETURN
            IF (J < 0 .OR. J > JBP1_IN) RETURN
            IF (K < 0 .OR. K > KBP1_IN) RETURN

            QR_D(I,J,K) = -CHI_R_D(I,J,K)*Q_D(I,J,K)

        END SUBROUTINE RADIATION_OFF_GPU

END MODULE G_RADI
