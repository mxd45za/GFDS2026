!------------------------------------------------------------------------------
! Module: G_EXCH
! Purpose: Exchange GPU-resident fields between local and MPI-connected meshes.
!------------------------------------------------------------------------------
MODULE G_EXCH
    USE G_MESH_POINTERS
    USE G_CONS
    USE G_MESH, ONLY: GPU_MESH_TYPE,GPU_OMESH_TYPE,GPU_MESHES,HT1D_BACKSIDE_FIELD_COUNT
    USE MPI_F08, ONLY: MPI_REQUEST
    USE G_ERRO, ONLY: GPU_WAIT,GPU_WAIT_ALL,GPU_MPI_IRECV,GPU_MPI_ISEND,GPU_MPI_WAITALL
    USE CUDAFOR
    INTEGER, PARAMETER :: HT1D_EXCHANGE_THREADS = 128
    CONTAINS

        !--------------------------------------------------------------------------
        ! GPU mesh-exchange orchestration
        !--------------------------------------------------------------------------
        SUBROUTINE GPU_MESH_EXCHANGE(CODE)
            IMPLICIT NONE

            INTEGER, INTENT(IN) :: CODE
            INTEGER :: NSRC
            REAL(EB) :: TNOW,EXCHANGE_TIME

            TNOW = CURRENT_TIME()

            IF (GPU_DEBUG_H2D_D2H) CALL SYNC_GPU_EXCHANGE_INPUTS_FROM_CPU(CODE)

            IF (N_MPI_PROCESSES>1) THEN
                CALL GPU_MESH_EXCHANGE_MPI(CODE)
            ELSE
                DO NSRC=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                    CALL COPY_GPU_EXCHANGE_SOURCE(CODE,NSRC)
                ENDDO
            ENDIF

            IF (N_MPI_PROCESSES==1) CALL GPU_WAIT_ALL
            IF (GPU_DEBUG_H2D_D2H) CALL SYNC_GPU_EXCHANGE_OUTPUTS_TO_CPU(CODE)
            IF (CODE==6) CALL GPU_MESH_EXCHANGE_HT1D_BACKSIDE

            EXCHANGE_TIME = CURRENT_TIME() - TNOW
            IF (ALLOCATED(T_USED)) T_USED(11) = T_USED(11) + EXCHANGE_TIME
            MESH_EXCHANGE_WALL_TIME = MESH_EXCHANGE_WALL_TIME + EXCHANGE_TIME

        END SUBROUTINE GPU_MESH_EXCHANGE

        ! Seed exchange-owned wall data during GPU initialization.
        SUBROUTINE GPU_MESH_EXCHANGE_INIT
            IMPLICIT NONE

            CALL GPU_MESH_EXCHANGE_HT1D_BACKSIDE
        END SUBROUTINE GPU_MESH_EXCHANGE_INIT

        !--------------------------------------------------------------------------
        ! MPI exchange of GPU-resident mesh fields
        !--------------------------------------------------------------------------
        SUBROUTINE GPU_MESH_EXCHANGE_MPI(CODE)
            IMPLICIT NONE

            INTEGER, INTENT(IN) :: CODE
            INTEGER, PARAMETER :: THREADS = 128
            INTEGER :: NM,NOM,NN,N_REQ,N_REQ_MAX,N_BLOCKS,N_ITEMS
            LOGICAL :: SAME_RANK_COPY_PENDING
            TYPE(MPI_REQUEST), ALLOCATABLE :: REQUEST(:)
            TYPE(GPU_OMESH_TYPE), POINTER :: GOM_LOCAL

            N_REQ_MAX = 0
            DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                IF (.NOT.GPU_MESHES(NM)%INITIALIZED) CYCLE
                IF (.NOT.ALLOCATED(MESHES(NM)%NEIGHBORING_MESH)) CYCLE
                DO NN=1,MESHES(NM)%N_NEIGHBORING_MESHES
                    NOM = MESHES(NM)%NEIGHBORING_MESH(NN)
                    IF (PROCESS(NOM)==MY_RANK) CYCLE
                    GOM_LOCAL => GPU_MESHES(NM)%OMESH(NOM)
                    IF (GOM_LOCAL%NIC_R>0) N_REQ_MAX = N_REQ_MAX + 1
                    IF (GOM_LOCAL%NIC_S>0) N_REQ_MAX = N_REQ_MAX + 1
                ENDDO
            ENDDO
            IF (N_REQ_MAX>0) ALLOCATE(REQUEST(N_REQ_MAX))

            ! Post every receive before any rank starts sending.
            N_REQ = 0
            DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                IF (.NOT.GPU_MESHES(NM)%INITIALIZED) CYCLE
                IF (.NOT.ALLOCATED(MESHES(NM)%NEIGHBORING_MESH)) CYCLE
                DO NN=1,MESHES(NM)%N_NEIGHBORING_MESHES
                    NOM = MESHES(NM)%NEIGHBORING_MESH(NN)
                    IF (PROCESS(NOM)==MY_RANK) CYCLE
                    GOM_LOCAL => GPU_MESHES(NM)%OMESH(NOM)
                    IF (GOM_LOCAL%NIC_R<=0) CYCLE
                    N_REQ = N_REQ + 1
                    SELECT CASE(CODE)
                        CASE(1,4)
                            CALL GPU_MPI_IRECV(GOM_LOCAL%RECV_PKG1_D,SIZE(GOM_LOCAL%RECV_PKG1_D), &
                                           PROCESS(NOM),NOM,REQUEST(N_REQ))
                        CASE(3,6)
                            CALL GPU_MPI_IRECV(GOM_LOCAL%RECV_PKG3_D,SIZE(GOM_LOCAL%RECV_PKG3_D), &
                                           PROCESS(NOM),NOM,REQUEST(N_REQ))
                        CASE(5)
                            CALL GPU_MPI_IRECV(GOM_LOCAL%RECV_PKG7_D,SIZE(GOM_LOCAL%RECV_PKG7_D), &
                                           PROCESS(NOM),NOM,REQUEST(N_REQ))
                    END SELECT
                ENDDO
            ENDDO

            ! Keep same-rank neighbors on the existing D2D path and pack remote neighbors on the GPU.
            SAME_RANK_COPY_PENDING = .FALSE.
            DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                IF (.NOT.GPU_MESHES(NM)%INITIALIZED) CYCLE
                CALL COPY_GPU_EXCHANGE_SOURCE(CODE,NM)
                IF (.NOT.ALLOCATED(MESHES(NM)%NEIGHBORING_MESH)) CYCLE
                DO NN=1,MESHES(NM)%N_NEIGHBORING_MESHES
                    NOM = MESHES(NM)%NEIGHBORING_MESH(NN)
                    GOM_LOCAL => GPU_MESHES(NM)%OMESH(NOM)
                    IF (PROCESS(NOM)==MY_RANK) THEN
                        IF (GOM_LOCAL%NIC_S>0) SAME_RANK_COPY_PENDING = .TRUE.
                        CYCLE
                    ENDIF
                    IF (GOM_LOCAL%NIC_S<=0) CYCLE
                    SELECT CASE(CODE)
                        CASE(1)
                            N_BLOCKS = (GOM_LOCAL%NIC_S+THREADS-1)/THREADS
                            CALL PACK_GPU_EXCHANGE_CODE_1_KERNEL<<<N_BLOCKS,THREADS,0,GPU_MESHES(NM)%STRM>>>( &
                                GPU_MESHES(NM)%RHOSD,GPU_MESHES(NM)%D_D,GPU_MESHES(NM)%MUD, &
                                GPU_MESHES(NM)%KRESD,GPU_MESHES(NM)%QD,GPU_MESHES(NM)%ZZSD, &
                                GOM_LOCAL%IIO_S_D,GOM_LOCAL%JJO_S_D,GOM_LOCAL%KKO_S_D,GOM_LOCAL%IOR_S_D, &
                                GOM_LOCAL%SEND_PKG1_D,GOM_LOCAL%NIC_S,MESHES(NM)%IBP1+1, &
                                MESHES(NM)%JBP1+1,MESHES(NM)%KBP1+1,N_TOTAL_SCALARS)
                        CASE(4)
                            N_BLOCKS = (GOM_LOCAL%NIC_S+THREADS-1)/THREADS
                            CALL PACK_GPU_EXCHANGE_CODE_1_KERNEL<<<N_BLOCKS,THREADS,0,GPU_MESHES(NM)%STRM>>>( &
                                GPU_MESHES(NM)%RHOD,GPU_MESHES(NM)%DS_D,GPU_MESHES(NM)%MUD, &
                                GPU_MESHES(NM)%KRESD,GPU_MESHES(NM)%QD,GPU_MESHES(NM)%ZZD, &
                                GOM_LOCAL%IIO_S_D,GOM_LOCAL%JJO_S_D,GOM_LOCAL%KKO_S_D,GOM_LOCAL%IOR_S_D, &
                                GOM_LOCAL%SEND_PKG1_D,GOM_LOCAL%NIC_S,MESHES(NM)%IBP1+1, &
                                MESHES(NM)%JBP1+1,MESHES(NM)%KBP1+1,N_TOTAL_SCALARS)
                        CASE(3)
                            N_ITEMS = (GOM_LOCAL%I_MAX_S-GOM_LOCAL%I_MIN_S+1)* &
                                      (GOM_LOCAL%J_MAX_S-GOM_LOCAL%J_MIN_S+1)* &
                                      (GOM_LOCAL%K_MAX_S-GOM_LOCAL%K_MIN_S+1)
                            N_BLOCKS = (N_ITEMS+THREADS-1)/THREADS
                            CALL PACK_GPU_EXCHANGE_CODE_3_KERNEL<<<N_BLOCKS,THREADS,0,GPU_MESHES(NM)%STRM>>>( &
                                GPU_MESHES(NM)%HSD,GPU_MESHES(NM)%USD,GPU_MESHES(NM)%VSD,GPU_MESHES(NM)%WSD, &
                                GOM_LOCAL%SEND_PKG3_D,GOM_LOCAL%I_MIN_S,GOM_LOCAL%I_MAX_S, &
                                GOM_LOCAL%J_MIN_S,GOM_LOCAL%J_MAX_S,GOM_LOCAL%K_MIN_S,GOM_LOCAL%K_MAX_S, &
                                MESHES(NM)%IBP1+1,MESHES(NM)%JBP1+1,N_ITEMS)
                        CASE(6)
                            N_ITEMS = (GOM_LOCAL%I_MAX_S-GOM_LOCAL%I_MIN_S+1)* &
                                      (GOM_LOCAL%J_MAX_S-GOM_LOCAL%J_MIN_S+1)* &
                                      (GOM_LOCAL%K_MAX_S-GOM_LOCAL%K_MIN_S+1)
                            N_BLOCKS = (N_ITEMS+THREADS-1)/THREADS
                            CALL PACK_GPU_EXCHANGE_CODE_3_KERNEL<<<N_BLOCKS,THREADS,0,GPU_MESHES(NM)%STRM>>>( &
                                GPU_MESHES(NM)%HD,GPU_MESHES(NM)%UD,GPU_MESHES(NM)%VD,GPU_MESHES(NM)%WD, &
                                GOM_LOCAL%SEND_PKG3_D,GOM_LOCAL%I_MIN_S,GOM_LOCAL%I_MAX_S, &
                                GOM_LOCAL%J_MIN_S,GOM_LOCAL%J_MAX_S,GOM_LOCAL%K_MIN_S,GOM_LOCAL%K_MAX_S, &
                                MESHES(NM)%IBP1+1,MESHES(NM)%JBP1+1,N_ITEMS)
                        CASE(5)
                            N_BLOCKS = (GOM_LOCAL%NIC_S+THREADS-1)/THREADS
                            IF (PREDICTOR) THEN
                                CALL PACK_GPU_EXCHANGE_CODE_5_KERNEL<<<N_BLOCKS,THREADS,0,GPU_MESHES(NM)%STRM>>>( &
                                    GPU_MESHES(NM)%FVXD,GPU_MESHES(NM)%FVYD,GPU_MESHES(NM)%FVZD,GPU_MESHES(NM)%HD, &
                                    GOM_LOCAL%IIO_S_D,GOM_LOCAL%JJO_S_D,GOM_LOCAL%KKO_S_D,GOM_LOCAL%IOR_S_D, &
                                    GOM_LOCAL%SEND_PKG7_D,GOM_LOCAL%NIC_S,MESHES(NM)%IBP1+1,MESHES(NM)%JBP1+1)
                            ELSE
                                CALL PACK_GPU_EXCHANGE_CODE_5_KERNEL<<<N_BLOCKS,THREADS,0,GPU_MESHES(NM)%STRM>>>( &
                                    GPU_MESHES(NM)%FVXD,GPU_MESHES(NM)%FVYD,GPU_MESHES(NM)%FVZD,GPU_MESHES(NM)%HSD, &
                                    GOM_LOCAL%IIO_S_D,GOM_LOCAL%JJO_S_D,GOM_LOCAL%KKO_S_D,GOM_LOCAL%IOR_S_D, &
                                    GOM_LOCAL%SEND_PKG7_D,GOM_LOCAL%NIC_S,MESHES(NM)%IBP1+1,MESHES(NM)%JBP1+1)
                            ENDIF
                    END SELECT
                ENDDO
            ENDDO

            IF (SAME_RANK_COPY_PENDING) THEN
                ! Same-rank D2D copies are launched on the source mesh stream,
                ! while the next stage reads the target mesh mirror on another
                ! stream. Complete those copies before MPI reads send packages
                ! and before the following velocity-matching stage.
                DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                    IF (.NOT.GPU_MESHES(NM)%INITIALIZED) CYCLE
                    CALL GPU_WAIT(GPU_MESHES(NM)%STRM)
                ENDDO
            ENDIF

            CALL GPU_WAIT_ALL

            ! Device send packages are complete and can now be read by CUDA-aware MPI.
            DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                IF (.NOT.GPU_MESHES(NM)%INITIALIZED) CYCLE
                IF (.NOT.ALLOCATED(MESHES(NM)%NEIGHBORING_MESH)) CYCLE
                DO NN=1,MESHES(NM)%N_NEIGHBORING_MESHES
                    NOM = MESHES(NM)%NEIGHBORING_MESH(NN)
                    IF (PROCESS(NOM)==MY_RANK) CYCLE
                    GOM_LOCAL => GPU_MESHES(NM)%OMESH(NOM)
                    IF (GOM_LOCAL%NIC_S<=0) CYCLE
                    N_REQ = N_REQ + 1
                    SELECT CASE(CODE)
                        CASE(1,4)
                            CALL GPU_MPI_ISEND(GOM_LOCAL%SEND_PKG1_D,SIZE(GOM_LOCAL%SEND_PKG1_D), &
                                           PROCESS(NOM),NM,REQUEST(N_REQ))
                        CASE(3,6)
                            CALL GPU_MPI_ISEND(GOM_LOCAL%SEND_PKG3_D,SIZE(GOM_LOCAL%SEND_PKG3_D), &
                                           PROCESS(NOM),NM,REQUEST(N_REQ))
                        CASE(5)
                            CALL GPU_MPI_ISEND(GOM_LOCAL%SEND_PKG7_D,SIZE(GOM_LOCAL%SEND_PKG7_D), &
                                           PROCESS(NOM),NM,REQUEST(N_REQ))
                    END SELECT
                ENDDO
            ENDDO

            IF (N_REQ>0) THEN
                CALL GPU_MPI_WAITALL(REQUEST(1:N_REQ))
            ENDIF

            ! Unpack each remote source package into the local target mesh's GPU OMESH mirror.
            DO NM=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                IF (.NOT.GPU_MESHES(NM)%INITIALIZED) CYCLE
                IF (.NOT.ALLOCATED(MESHES(NM)%NEIGHBORING_MESH)) CYCLE
                DO NN=1,MESHES(NM)%N_NEIGHBORING_MESHES
                    NOM = MESHES(NM)%NEIGHBORING_MESH(NN)
                    IF (PROCESS(NOM)==MY_RANK) CYCLE
                    GOM_LOCAL => GPU_MESHES(NM)%OMESH(NOM)
                    IF (GOM_LOCAL%NIC_R<=0) CYCLE
                    SELECT CASE(CODE)
                        CASE(1)
                            N_BLOCKS = (GOM_LOCAL%NIC_R+THREADS-1)/THREADS
                            CALL UNPACK_GPU_EXCHANGE_CODE_1_KERNEL<<<N_BLOCKS,THREADS,0,GPU_MESHES(NM)%STRM>>>( &
                                GOM_LOCAL%RECV_PKG1_D,GOM_LOCAL%IIO_R_D,GOM_LOCAL%JJO_R_D, &
                                GOM_LOCAL%KKO_R_D,GOM_LOCAL%IOR_R_D,GOM_LOCAL%RHOS_D,GOM_LOCAL%D_D, &
                                GOM_LOCAL%MU_D,GOM_LOCAL%KRES_D,GOM_LOCAL%Q_D,GOM_LOCAL%ZZS_D, &
                                GOM_LOCAL%NIC_R,LBOUND(GOM_LOCAL%RHOS_D,1),LBOUND(GOM_LOCAL%RHOS_D,2), &
                                LBOUND(GOM_LOCAL%RHOS_D,3),SIZE(GOM_LOCAL%RHOS_D,1),SIZE(GOM_LOCAL%RHOS_D,2), &
                                SIZE(GOM_LOCAL%RHOS_D,3),N_TOTAL_SCALARS)
                        CASE(4)
                            N_BLOCKS = (GOM_LOCAL%NIC_R+THREADS-1)/THREADS
                            CALL UNPACK_GPU_EXCHANGE_CODE_1_KERNEL<<<N_BLOCKS,THREADS,0,GPU_MESHES(NM)%STRM>>>( &
                                GOM_LOCAL%RECV_PKG1_D,GOM_LOCAL%IIO_R_D,GOM_LOCAL%JJO_R_D, &
                                GOM_LOCAL%KKO_R_D,GOM_LOCAL%IOR_R_D,GOM_LOCAL%RHO_D,GOM_LOCAL%DS_D, &
                                GOM_LOCAL%MU_D,GOM_LOCAL%KRES_D,GOM_LOCAL%Q_D,GOM_LOCAL%ZZ_D, &
                                GOM_LOCAL%NIC_R,LBOUND(GOM_LOCAL%RHO_D,1),LBOUND(GOM_LOCAL%RHO_D,2), &
                                LBOUND(GOM_LOCAL%RHO_D,3),SIZE(GOM_LOCAL%RHO_D,1),SIZE(GOM_LOCAL%RHO_D,2), &
                                SIZE(GOM_LOCAL%RHO_D,3),N_TOTAL_SCALARS)
                        CASE(3)
                            N_ITEMS = (GOM_LOCAL%I_MAX_R-GOM_LOCAL%I_MIN_R+1)* &
                                      (GOM_LOCAL%J_MAX_R-GOM_LOCAL%J_MIN_R+1)* &
                                      (GOM_LOCAL%K_MAX_R-GOM_LOCAL%K_MIN_R+1)
                            N_BLOCKS = (N_ITEMS+THREADS-1)/THREADS
                            CALL UNPACK_GPU_EXCHANGE_CODE_3_KERNEL<<<N_BLOCKS,THREADS,0,GPU_MESHES(NM)%STRM>>>( &
                                GOM_LOCAL%RECV_PKG3_D,GOM_LOCAL%HS_D,GOM_LOCAL%US_D,GOM_LOCAL%VS_D,GOM_LOCAL%WS_D, &
                                GOM_LOCAL%I_MIN_R,GOM_LOCAL%I_MAX_R,GOM_LOCAL%J_MIN_R,GOM_LOCAL%J_MAX_R, &
                                GOM_LOCAL%K_MIN_R,GOM_LOCAL%K_MAX_R,LBOUND(GOM_LOCAL%HS_D,1), &
                                LBOUND(GOM_LOCAL%HS_D,2),LBOUND(GOM_LOCAL%HS_D,3),SIZE(GOM_LOCAL%HS_D,1), &
                                SIZE(GOM_LOCAL%HS_D,2),N_ITEMS)
                        CASE(6)
                            N_ITEMS = (GOM_LOCAL%I_MAX_R-GOM_LOCAL%I_MIN_R+1)* &
                                      (GOM_LOCAL%J_MAX_R-GOM_LOCAL%J_MIN_R+1)* &
                                      (GOM_LOCAL%K_MAX_R-GOM_LOCAL%K_MIN_R+1)
                            N_BLOCKS = (N_ITEMS+THREADS-1)/THREADS
                            CALL UNPACK_GPU_EXCHANGE_CODE_3_KERNEL<<<N_BLOCKS,THREADS,0,GPU_MESHES(NM)%STRM>>>( &
                                GOM_LOCAL%RECV_PKG3_D,GOM_LOCAL%H_D,GOM_LOCAL%U_D,GOM_LOCAL%V_D,GOM_LOCAL%W_D, &
                                GOM_LOCAL%I_MIN_R,GOM_LOCAL%I_MAX_R,GOM_LOCAL%J_MIN_R,GOM_LOCAL%J_MAX_R, &
                                GOM_LOCAL%K_MIN_R,GOM_LOCAL%K_MAX_R,LBOUND(GOM_LOCAL%H_D,1), &
                                LBOUND(GOM_LOCAL%H_D,2),LBOUND(GOM_LOCAL%H_D,3),SIZE(GOM_LOCAL%H_D,1), &
                                SIZE(GOM_LOCAL%H_D,2),N_ITEMS)
                        CASE(5)
                            N_BLOCKS = (GOM_LOCAL%NIC_R+THREADS-1)/THREADS
                            IF (PREDICTOR) THEN
                                CALL UNPACK_GPU_EXCHANGE_CODE_5_KERNEL<<<N_BLOCKS,THREADS,0,GPU_MESHES(NM)%STRM>>>( &
                                    GOM_LOCAL%RECV_PKG7_D,GOM_LOCAL%IIO_R_D,GOM_LOCAL%JJO_R_D,GOM_LOCAL%KKO_R_D, &
                                    GOM_LOCAL%IOR_R_D,GOM_LOCAL%FVX_D,GOM_LOCAL%FVY_D,GOM_LOCAL%FVZ_D,GOM_LOCAL%H_D, &
                                    GOM_LOCAL%NIC_R,LBOUND(GOM_LOCAL%H_D,1),LBOUND(GOM_LOCAL%H_D,2),LBOUND(GOM_LOCAL%H_D,3), &
                                    SIZE(GOM_LOCAL%H_D,1),SIZE(GOM_LOCAL%H_D,2))
                            ELSE
                                CALL UNPACK_GPU_EXCHANGE_CODE_5_KERNEL<<<N_BLOCKS,THREADS,0,GPU_MESHES(NM)%STRM>>>( &
                                    GOM_LOCAL%RECV_PKG7_D,GOM_LOCAL%IIO_R_D,GOM_LOCAL%JJO_R_D,GOM_LOCAL%KKO_R_D, &
                                    GOM_LOCAL%IOR_R_D,GOM_LOCAL%FVX_D,GOM_LOCAL%FVY_D,GOM_LOCAL%FVZ_D,GOM_LOCAL%HS_D, &
                                    GOM_LOCAL%NIC_R,LBOUND(GOM_LOCAL%HS_D,1),LBOUND(GOM_LOCAL%HS_D,2),LBOUND(GOM_LOCAL%HS_D,3), &
                                    SIZE(GOM_LOCAL%HS_D,1),SIZE(GOM_LOCAL%HS_D,2))
                            ENDIF
                    END SELECT
                ENDDO
            ENDDO

            CALL GPU_WAIT_ALL

            IF (ALLOCATED(REQUEST)) DEALLOCATE(REQUEST)

        END SUBROUTINE GPU_MESH_EXCHANGE_MPI

        ATTRIBUTES(GLOBAL) SUBROUTINE PACK_GPU_EXCHANGE_CODE_1_KERNEL( &
            RHOS,D_FIELD,MU,KRES,Q,ZZS,IIO,JJO,KKO,IOR,SEND_PKG,NIC,NI,NJ,NK,N_SCALARS)
            IMPLICIT NONE

            INTEGER, VALUE :: NIC,NI,NJ,NK,N_SCALARS
            INTEGER, DEVICE, DIMENSION(*) :: IIO,JJO,KKO,IOR
            REAL(EB), DEVICE, DIMENSION(*) :: RHOS,D_FIELD,MU,KRES,Q,ZZS,SEND_PKG
            INTEGER :: LL,II1,II2,JJ1,JJ2,KK1,KK2,IDX1,IDX2,NN,NQT,PKG_BASE,FIELD_SIZE

            LL = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
            IF (LL>NIC) RETURN

            II1 = IIO(LL); II2 = II1
            JJ1 = JJO(LL); JJ2 = JJ1
            KK1 = KKO(LL); KK2 = KK1
            SELECT CASE(IOR(LL))
                CASE(-1); II2 = II1 + 1
                CASE( 1); II2 = II1 - 1
                CASE(-2); JJ2 = JJ1 + 1
                CASE( 2); JJ2 = JJ1 - 1
                CASE(-3); KK2 = KK1 + 1
                CASE( 3); KK2 = KK1 - 1
            END SELECT

            IDX1 = 1 + II1 + NI*(JJ1 + NJ*KK1)
            IDX2 = 1 + II2 + NI*(JJ2 + NJ*KK2)
            NQT = 6 + 2*N_SCALARS
            PKG_BASE = NQT*(LL-1)
            FIELD_SIZE = NI*NJ*NK

            SEND_PKG(PKG_BASE+1) = RHOS(IDX1)
            SEND_PKG(PKG_BASE+2) = RHOS(IDX2)
            SEND_PKG(PKG_BASE+3) = MU(IDX1)
            SEND_PKG(PKG_BASE+4) = KRES(IDX1)
            SEND_PKG(PKG_BASE+5) = D_FIELD(IDX1)
            SEND_PKG(PKG_BASE+6) = Q(IDX1)
            DO NN=1,N_SCALARS
                SEND_PKG(PKG_BASE+6+2*NN-1) = ZZS(IDX1+FIELD_SIZE*(NN-1))
                SEND_PKG(PKG_BASE+6+2*NN  ) = ZZS(IDX2+FIELD_SIZE*(NN-1))
            ENDDO

        END SUBROUTINE PACK_GPU_EXCHANGE_CODE_1_KERNEL

        ATTRIBUTES(GLOBAL) SUBROUTINE UNPACK_GPU_EXCHANGE_CODE_1_KERNEL( &
            RECV_PKG,IIO,JJO,KKO,IOR,RHOS,D_FIELD,MU,KRES,Q,ZZS,NIC,LB_I,LB_J,LB_K,NI,NJ,NK,N_SCALARS)
            IMPLICIT NONE

            INTEGER, VALUE :: NIC,LB_I,LB_J,LB_K,NI,NJ,NK,N_SCALARS
            INTEGER, DEVICE, DIMENSION(*) :: IIO,JJO,KKO,IOR
            REAL(EB), DEVICE, DIMENSION(*) :: RECV_PKG,RHOS,D_FIELD,MU,KRES,Q,ZZS
            INTEGER :: LL,II1,II2,JJ1,JJ2,KK1,KK2,IDX1,IDX2,NN,NQT,PKG_BASE,FIELD_SIZE

            LL = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
            IF (LL>NIC) RETURN

            II1 = IIO(LL); II2 = II1
            JJ1 = JJO(LL); JJ2 = JJ1
            KK1 = KKO(LL); KK2 = KK1
            SELECT CASE(IOR(LL))
                CASE(-1); II2 = II1 + 1
                CASE( 1); II2 = II1 - 1
                CASE(-2); JJ2 = JJ1 + 1
                CASE( 2); JJ2 = JJ1 - 1
                CASE(-3); KK2 = KK1 + 1
                CASE( 3); KK2 = KK1 - 1
            END SELECT

            IDX1 = 1 + (II1-LB_I) + NI*((JJ1-LB_J) + NJ*(KK1-LB_K))
            IDX2 = 1 + (II2-LB_I) + NI*((JJ2-LB_J) + NJ*(KK2-LB_K))
            NQT = 6 + 2*N_SCALARS
            PKG_BASE = NQT*(LL-1)
            FIELD_SIZE = NI*NJ*NK

            RHOS(IDX1) = RECV_PKG(PKG_BASE+1)
            RHOS(IDX2) = RECV_PKG(PKG_BASE+2)
            MU(IDX1) = RECV_PKG(PKG_BASE+3)
            KRES(IDX1) = RECV_PKG(PKG_BASE+4)
            D_FIELD(IDX1) = RECV_PKG(PKG_BASE+5)
            Q(IDX1) = RECV_PKG(PKG_BASE+6)
            DO NN=1,N_SCALARS
                ZZS(IDX1+FIELD_SIZE*(NN-1)) = RECV_PKG(PKG_BASE+6+2*NN-1)
                ZZS(IDX2+FIELD_SIZE*(NN-1)) = RECV_PKG(PKG_BASE+6+2*NN  )
            ENDDO

        END SUBROUTINE UNPACK_GPU_EXCHANGE_CODE_1_KERNEL

        ATTRIBUTES(GLOBAL) SUBROUTINE PACK_GPU_EXCHANGE_CODE_3_KERNEL( &
            H,U,V,W,SEND_PKG,I_MIN,I_MAX,J_MIN,J_MAX,K_MIN,K_MAX,NI,NJ,N_ITEMS)
            IMPLICIT NONE

            INTEGER, VALUE :: I_MIN,I_MAX,J_MIN,J_MAX,K_MIN,K_MAX,NI,NJ,N_ITEMS
            REAL(EB), DEVICE, DIMENSION(*) :: H,U,V,W,SEND_PKG
            INTEGER :: LL,I,J,K,NX,NY,SRC_IDX,PKG_BASE

            LL = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X - 1
            IF (LL>=N_ITEMS) RETURN

            NX = I_MAX-I_MIN+1
            NY = J_MAX-J_MIN+1
            I = I_MIN + MOD(LL,NX)
            J = J_MIN + MOD(LL/NX,NY)
            K = K_MIN + LL/(NX*NY)
            SRC_IDX = 1 + I + NI*(J + NJ*K)
            PKG_BASE = 4*LL

            SEND_PKG(PKG_BASE+1) = H(SRC_IDX)
            SEND_PKG(PKG_BASE+2) = U(SRC_IDX)
            SEND_PKG(PKG_BASE+3) = V(SRC_IDX)
            SEND_PKG(PKG_BASE+4) = W(SRC_IDX)

        END SUBROUTINE PACK_GPU_EXCHANGE_CODE_3_KERNEL

        ATTRIBUTES(GLOBAL) SUBROUTINE UNPACK_GPU_EXCHANGE_CODE_3_KERNEL( &
            RECV_PKG,H,U,V,W,I_MIN,I_MAX,J_MIN,J_MAX,K_MIN,K_MAX,LB_I,LB_J,LB_K,NI,NJ,N_ITEMS)
            IMPLICIT NONE

            INTEGER, VALUE :: I_MIN,I_MAX,J_MIN,J_MAX,K_MIN,K_MAX,LB_I,LB_J,LB_K,NI,NJ,N_ITEMS
            REAL(EB), DEVICE, DIMENSION(*) :: RECV_PKG,H,U,V,W
            INTEGER :: LL,I,J,K,NX,NY,DST_IDX,PKG_BASE

            LL = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X - 1
            IF (LL>=N_ITEMS) RETURN

            NX = I_MAX-I_MIN+1
            NY = J_MAX-J_MIN+1
            I = I_MIN + MOD(LL,NX)
            J = J_MIN + MOD(LL/NX,NY)
            K = K_MIN + LL/(NX*NY)
            DST_IDX = 1 + (I-LB_I) + NI*((J-LB_J) + NJ*(K-LB_K))
            PKG_BASE = 4*LL

            H(DST_IDX) = RECV_PKG(PKG_BASE+1)
            U(DST_IDX) = RECV_PKG(PKG_BASE+2)
            V(DST_IDX) = RECV_PKG(PKG_BASE+3)
            W(DST_IDX) = RECV_PKG(PKG_BASE+4)

        END SUBROUTINE UNPACK_GPU_EXCHANGE_CODE_3_KERNEL

        ATTRIBUTES(GLOBAL) SUBROUTINE PACK_GPU_EXCHANGE_CODE_5_KERNEL( &
            FVX,FVY,FVZ,H,IIO,JJO,KKO,IOR,SEND_PKG,NIC,NI,NJ)
            IMPLICIT NONE

            INTEGER, VALUE :: NIC,NI,NJ
            INTEGER, DEVICE, DIMENSION(*) :: IIO,JJO,KKO,IOR
            REAL(EB), DEVICE, DIMENSION(*) :: FVX,FVY,FVZ,H,SEND_PKG
            INTEGER :: LL,I,J,K,FACE_IDX,H1_IDX,H2_IDX,PKG_BASE
            REAL(EB) :: FACE_VALUE

            LL = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
            IF (LL>NIC) RETURN

            I = IIO(LL)
            J = JJO(LL)
            K = KKO(LL)
            SELECT CASE(IOR(LL))
                CASE(-1)
                    FACE_IDX = 1 + (I-1) + NI*(J + NJ*K)
                    H1_IDX = FACE_IDX
                    H2_IDX = 1 + I + NI*(J + NJ*K)
                    FACE_VALUE = FVX(FACE_IDX)
                CASE(1)
                    FACE_IDX = 1 + I + NI*(J + NJ*K)
                    H1_IDX = FACE_IDX
                    H2_IDX = 1 + (I+1) + NI*(J + NJ*K)
                    FACE_VALUE = FVX(FACE_IDX)
                CASE(-2)
                    FACE_IDX = 1 + I + NI*((J-1) + NJ*K)
                    H1_IDX = FACE_IDX
                    H2_IDX = 1 + I + NI*(J + NJ*K)
                    FACE_VALUE = FVY(FACE_IDX)
                CASE(2)
                    FACE_IDX = 1 + I + NI*(J + NJ*K)
                    H1_IDX = FACE_IDX
                    H2_IDX = 1 + I + NI*((J+1) + NJ*K)
                    FACE_VALUE = FVY(FACE_IDX)
                CASE(-3)
                    FACE_IDX = 1 + I + NI*(J + NJ*(K-1))
                    H1_IDX = FACE_IDX
                    H2_IDX = 1 + I + NI*(J + NJ*K)
                    FACE_VALUE = FVZ(FACE_IDX)
                CASE(3)
                    FACE_IDX = 1 + I + NI*(J + NJ*K)
                    H1_IDX = FACE_IDX
                    H2_IDX = 1 + I + NI*(J + NJ*(K+1))
                    FACE_VALUE = FVZ(FACE_IDX)
                CASE DEFAULT
                    RETURN
            END SELECT

            PKG_BASE = 3*(LL-1)
            SEND_PKG(PKG_BASE+1) = FACE_VALUE
            SEND_PKG(PKG_BASE+2) = H(H1_IDX)
            SEND_PKG(PKG_BASE+3) = H(H2_IDX)

        END SUBROUTINE PACK_GPU_EXCHANGE_CODE_5_KERNEL

        ATTRIBUTES(GLOBAL) SUBROUTINE UNPACK_GPU_EXCHANGE_CODE_5_KERNEL( &
            RECV_PKG,IIO,JJO,KKO,IOR,FVX,FVY,FVZ,H,NIC,LB_I,LB_J,LB_K,NI,NJ)
            IMPLICIT NONE

            INTEGER, VALUE :: NIC,LB_I,LB_J,LB_K,NI,NJ
            INTEGER, DEVICE, DIMENSION(*) :: IIO,JJO,KKO,IOR
            REAL(EB), DEVICE, DIMENSION(*) :: RECV_PKG,FVX,FVY,FVZ,H
            INTEGER :: LL,I,J,K,FACE_IDX,H1_IDX,H2_IDX,PKG_BASE

            LL = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
            IF (LL>NIC) RETURN

            I = IIO(LL)
            J = JJO(LL)
            K = KKO(LL)
            SELECT CASE(IOR(LL))
                CASE(-1)
                    FACE_IDX = 1 + (I-1-LB_I) + NI*((J-LB_J) + NJ*(K-LB_K))
                    H1_IDX = FACE_IDX
                    H2_IDX = 1 + (I-LB_I) + NI*((J-LB_J) + NJ*(K-LB_K))
                CASE(1)
                    FACE_IDX = 1 + (I-LB_I) + NI*((J-LB_J) + NJ*(K-LB_K))
                    H1_IDX = FACE_IDX
                    H2_IDX = 1 + (I+1-LB_I) + NI*((J-LB_J) + NJ*(K-LB_K))
                CASE(-2)
                    FACE_IDX = 1 + (I-LB_I) + NI*((J-1-LB_J) + NJ*(K-LB_K))
                    H1_IDX = FACE_IDX
                    H2_IDX = 1 + (I-LB_I) + NI*((J-LB_J) + NJ*(K-LB_K))
                CASE(2)
                    FACE_IDX = 1 + (I-LB_I) + NI*((J-LB_J) + NJ*(K-LB_K))
                    H1_IDX = FACE_IDX
                    H2_IDX = 1 + (I-LB_I) + NI*((J+1-LB_J) + NJ*(K-LB_K))
                CASE(-3)
                    FACE_IDX = 1 + (I-LB_I) + NI*((J-LB_J) + NJ*(K-1-LB_K))
                    H1_IDX = FACE_IDX
                    H2_IDX = 1 + (I-LB_I) + NI*((J-LB_J) + NJ*(K-LB_K))
                CASE(3)
                    FACE_IDX = 1 + (I-LB_I) + NI*((J-LB_J) + NJ*(K-LB_K))
                    H1_IDX = FACE_IDX
                    H2_IDX = 1 + (I-LB_I) + NI*((J-LB_J) + NJ*(K+1-LB_K))
                CASE DEFAULT
                    RETURN
            END SELECT

            PKG_BASE = 3*(LL-1)
            SELECT CASE(ABS(IOR(LL)))
                CASE(1); FVX(FACE_IDX) = RECV_PKG(PKG_BASE+1)
                CASE(2); FVY(FACE_IDX) = RECV_PKG(PKG_BASE+1)
                CASE(3); FVZ(FACE_IDX) = RECV_PKG(PKG_BASE+1)
            END SELECT
            H(H1_IDX) = RECV_PKG(PKG_BASE+2)
            H(H2_IDX) = RECV_PKG(PKG_BASE+3)

        END SUBROUTINE UNPACK_GPU_EXCHANGE_CODE_5_KERNEL

        !--------------------------------------------------------------------------
        ! Same-rank device-to-device source copy
        !--------------------------------------------------------------------------
        SUBROUTINE COPY_GPU_EXCHANGE_SOURCE(CODE,NSRC)
            IMPLICIT NONE

            INTEGER, INTENT(IN) :: CODE,NSRC
            INTEGER, PARAMETER :: BS = 8
            INTEGER :: NN,NDST,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX
            INTEGER :: SRC_NI,SRC_NJ,DST_NI,DST_NJ,DST_LB_I,DST_LB_J,DST_LB_K
            TYPE(OMESH_TYPE), POINTER :: SEND_OM
            TYPE(GPU_OMESH_TYPE), POINTER :: RECV_GOM
            TYPE(DIM3) :: BLKD,GRDD
            REAL(EB), DEVICE, POINTER, CONTIGUOUS :: SRC3_1(:,:,:),SRC3_2(:,:,:),SRC3_3(:,:,:),SRC3_4(:,:,:),SRC3_5(:,:,:)
            REAL(EB), DEVICE, POINTER, CONTIGUOUS :: DST3_1(:,:,:),DST3_2(:,:,:),DST3_3(:,:,:),DST3_4(:,:,:),DST3_5(:,:,:)
            REAL(EB), DEVICE, POINTER, CONTIGUOUS :: SRC4(:,:,:,:),DST4(:,:,:,:)

            IF (.NOT.ALLOCATED(GPU_MESHES)) RETURN
            IF (NSRC < LBOUND(GPU_MESHES,1) .OR. NSRC > UBOUND(GPU_MESHES,1)) RETURN
            IF (.NOT.GPU_MESHES(NSRC)%INITIALIZED) RETURN
            IF (.NOT.ALLOCATED(MESHES(NSRC)%NEIGHBORING_MESH)) RETURN

            DO NN=1,MESHES(NSRC)%N_NEIGHBORING_MESHES
                NDST = MESHES(NSRC)%NEIGHBORING_MESH(NN)
                IF (NDST < LBOUND(GPU_MESHES,1) .OR. NDST > UBOUND(GPU_MESHES,1)) CYCLE
                IF (.NOT.GPU_MESHES(NDST)%INITIALIZED) CYCLE
                IF (.NOT.ALLOCATED(GPU_MESHES(NDST)%OMESH)) CYCLE

                SEND_OM => MESHES(NSRC)%OMESH(NDST)
                IF (SEND_OM%NIC_S <= 0) CYCLE

                RECV_GOM => GPU_MESHES(NDST)%OMESH(NSRC)

                IMIN = SEND_OM%I_MIN_S
                IMAX = SEND_OM%I_MAX_S
                JMIN = SEND_OM%J_MIN_S
                JMAX = SEND_OM%J_MAX_S
                KMIN = SEND_OM%K_MIN_S
                KMAX = SEND_OM%K_MAX_S
                IF (IMIN>IMAX .OR. JMIN>JMAX .OR. KMIN>KMAX) CYCLE

                SRC_NI = MESHES(NSRC)%IBP1 + 1
                SRC_NJ = MESHES(NSRC)%JBP1 + 1
                DST_NI = SIZE(RECV_GOM%RHO_D,1)
                DST_NJ = SIZE(RECV_GOM%RHO_D,2)
                DST_LB_I = LBOUND(RECV_GOM%RHO_D,1)
                DST_LB_J = LBOUND(RECV_GOM%RHO_D,2)
                DST_LB_K = LBOUND(RECV_GOM%RHO_D,3)

                BLKD = DIM3(BS,BS,BS)
                GRDD = DIM3((IMAX-IMIN+BS)/BS,(JMAX-JMIN+BS)/BS,(KMAX-KMIN+BS)/BS)

                SELECT CASE(CODE)
                    CASE(1)
                        SRC3_1 => GPU_MESHES(NSRC)%RHOSD; DST3_1 => RECV_GOM%RHOS_D
                        SRC3_2 => GPU_MESHES(NSRC)%D_D;   DST3_2 => RECV_GOM%D_D
                        SRC3_3 => GPU_MESHES(NSRC)%MUD;   DST3_3 => RECV_GOM%MU_D
                        SRC3_4 => GPU_MESHES(NSRC)%KRESD; DST3_4 => RECV_GOM%KRES_D
                        SRC3_5 => GPU_MESHES(NSRC)%QD;    DST3_5 => RECV_GOM%Q_D
                        CALL COPY_GPU_FIELDS_5X3D_KERNEL<<<GRDD,BLKD,0,GPU_MESHES(NSRC)%STRM>>>(SRC3_1,SRC3_2,SRC3_3,SRC3_4,SRC3_5, &
                            DST3_1,DST3_2,DST3_3,DST3_4,DST3_5,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX, &
                            DST_LB_I,DST_LB_J,DST_LB_K,SRC_NI,SRC_NJ,DST_NI,DST_NJ)
                        SRC4 => GPU_MESHES(NSRC)%ZZSD
                        DST4 => RECV_GOM%ZZS_D
                        CALL COPY_GPU_FIELD_4D_KERNEL<<<GRDD,BLKD,0,GPU_MESHES(NSRC)%STRM>>>(SRC4,DST4,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX, &
                                                                     LBOUND(DST4,1),LBOUND(DST4,2),LBOUND(DST4,3), &
                                                                     SIZE(SRC4,1),SIZE(SRC4,2),SIZE(SRC4,3), &
                                                                     SIZE(DST4,1),SIZE(DST4,2),SIZE(DST4,3),N_TOTAL_SCALARS)
                    CASE(3)
                        SRC3_1 => GPU_MESHES(NSRC)%HSD; DST3_1 => RECV_GOM%HS_D
                        SRC3_2 => GPU_MESHES(NSRC)%USD; DST3_2 => RECV_GOM%US_D
                        SRC3_3 => GPU_MESHES(NSRC)%VSD; DST3_3 => RECV_GOM%VS_D
                        SRC3_4 => GPU_MESHES(NSRC)%WSD; DST3_4 => RECV_GOM%WS_D
                        CALL COPY_GPU_FIELDS_4X3D_KERNEL<<<GRDD,BLKD,0,GPU_MESHES(NSRC)%STRM>>>(SRC3_1,SRC3_2,SRC3_3,SRC3_4, &
                            DST3_1,DST3_2,DST3_3,DST3_4,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX, &
                            DST_LB_I,DST_LB_J,DST_LB_K,SRC_NI,SRC_NJ,DST_NI,DST_NJ)
                    CASE(4)
                        SRC3_1 => GPU_MESHES(NSRC)%RHOD;  DST3_1 => RECV_GOM%RHO_D
                        SRC3_2 => GPU_MESHES(NSRC)%DS_D;  DST3_2 => RECV_GOM%DS_D
                        SRC3_3 => GPU_MESHES(NSRC)%MUD;   DST3_3 => RECV_GOM%MU_D
                        SRC3_4 => GPU_MESHES(NSRC)%KRESD; DST3_4 => RECV_GOM%KRES_D
                        SRC3_5 => GPU_MESHES(NSRC)%QD;    DST3_5 => RECV_GOM%Q_D
                        CALL COPY_GPU_FIELDS_5X3D_KERNEL<<<GRDD,BLKD,0,GPU_MESHES(NSRC)%STRM>>>(SRC3_1,SRC3_2,SRC3_3,SRC3_4,SRC3_5, &
                            DST3_1,DST3_2,DST3_3,DST3_4,DST3_5,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX, &
                            DST_LB_I,DST_LB_J,DST_LB_K,SRC_NI,SRC_NJ,DST_NI,DST_NJ)
                        SRC4 => GPU_MESHES(NSRC)%ZZD
                        DST4 => RECV_GOM%ZZ_D
                        CALL COPY_GPU_FIELD_4D_KERNEL<<<GRDD,BLKD,0,GPU_MESHES(NSRC)%STRM>>>(SRC4,DST4,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX, &
                                                                     LBOUND(DST4,1),LBOUND(DST4,2),LBOUND(DST4,3), &
                                                                     SIZE(SRC4,1),SIZE(SRC4,2),SIZE(SRC4,3), &
                                                                     SIZE(DST4,1),SIZE(DST4,2),SIZE(DST4,3),N_TOTAL_SCALARS)
                    CASE(5)
                        SRC3_1 => GPU_MESHES(NSRC)%FVXD; DST3_1 => RECV_GOM%FVX_D
                        SRC3_2 => GPU_MESHES(NSRC)%FVYD; DST3_2 => RECV_GOM%FVY_D
                        SRC3_3 => GPU_MESHES(NSRC)%FVZD; DST3_3 => RECV_GOM%FVZ_D
                        IF (PREDICTOR) THEN
                            SRC3_4 => GPU_MESHES(NSRC)%HD
                            DST3_4 => RECV_GOM%H_D
                        ELSE
                            SRC3_4 => GPU_MESHES(NSRC)%HSD
                            DST3_4 => RECV_GOM%HS_D
                        ENDIF
                        CALL COPY_GPU_FIELDS_4X3D_KERNEL<<<GRDD,BLKD,0,GPU_MESHES(NSRC)%STRM>>>(SRC3_1,SRC3_2,SRC3_3,SRC3_4, &
                            DST3_1,DST3_2,DST3_3,DST3_4,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX, &
                            DST_LB_I,DST_LB_J,DST_LB_K,SRC_NI,SRC_NJ,DST_NI,DST_NJ)
                    CASE(6)
                        SRC3_1 => GPU_MESHES(NSRC)%HD; DST3_1 => RECV_GOM%H_D
                        SRC3_2 => GPU_MESHES(NSRC)%UD; DST3_2 => RECV_GOM%U_D
                        SRC3_3 => GPU_MESHES(NSRC)%VD; DST3_3 => RECV_GOM%V_D
                        SRC3_4 => GPU_MESHES(NSRC)%WD; DST3_4 => RECV_GOM%W_D
                        CALL COPY_GPU_FIELDS_4X3D_KERNEL<<<GRDD,BLKD,0,GPU_MESHES(NSRC)%STRM>>>(SRC3_1,SRC3_2,SRC3_3,SRC3_4, &
                            DST3_1,DST3_2,DST3_3,DST3_4,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX, &
                            DST_LB_I,DST_LB_J,DST_LB_K,SRC_NI,SRC_NJ,DST_NI,DST_NJ)
                END SELECT
            ENDDO

        END SUBROUTINE COPY_GPU_EXCHANGE_SOURCE

        ATTRIBUTES(GLOBAL) SUBROUTINE COPY_GPU_FIELDS_4X3D_KERNEL(SRC1,SRC2,SRC3,SRC4,DST1,DST2,DST3,DST4, &
            IMIN,IMAX,JMIN,JMAX,KMIN,KMAX,DST_LB_I,DST_LB_J,DST_LB_K,SRC_NI,SRC_NJ,DST_NI,DST_NJ)
            IMPLICIT NONE

            INTEGER, VALUE :: IMIN,IMAX,JMIN,JMAX,KMIN,KMAX,DST_LB_I,DST_LB_J,DST_LB_K
            INTEGER, VALUE :: SRC_NI,SRC_NJ,DST_NI,DST_NJ
            REAL(EB), DEVICE, DIMENSION(*) :: SRC1,SRC2,SRC3,SRC4,DST1,DST2,DST3,DST4
            INTEGER :: I,J,K,SRC_IDX,DST_IDX

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X + IMIN - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y + JMIN - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z + KMIN - 1

            IF (I>IMAX .OR. J>JMAX .OR. K>KMAX) RETURN

            SRC_IDX = 1 + I + SRC_NI*(J + SRC_NJ*K)
            DST_IDX = 1 + (I-DST_LB_I) + DST_NI*((J-DST_LB_J) + DST_NJ*(K-DST_LB_K))
            DST1(DST_IDX) = SRC1(SRC_IDX)
            DST2(DST_IDX) = SRC2(SRC_IDX)
            DST3(DST_IDX) = SRC3(SRC_IDX)
            DST4(DST_IDX) = SRC4(SRC_IDX)

        END SUBROUTINE COPY_GPU_FIELDS_4X3D_KERNEL

        ATTRIBUTES(GLOBAL) SUBROUTINE COPY_GPU_FIELDS_5X3D_KERNEL(SRC1,SRC2,SRC3,SRC4,SRC5,DST1,DST2,DST3,DST4,DST5, &
            IMIN,IMAX,JMIN,JMAX,KMIN,KMAX,DST_LB_I,DST_LB_J,DST_LB_K,SRC_NI,SRC_NJ,DST_NI,DST_NJ)
            IMPLICIT NONE

            INTEGER, VALUE :: IMIN,IMAX,JMIN,JMAX,KMIN,KMAX,DST_LB_I,DST_LB_J,DST_LB_K
            INTEGER, VALUE :: SRC_NI,SRC_NJ,DST_NI,DST_NJ
            REAL(EB), DEVICE, DIMENSION(*) :: SRC1,SRC2,SRC3,SRC4,SRC5,DST1,DST2,DST3,DST4,DST5
            INTEGER :: I,J,K,SRC_IDX,DST_IDX

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X + IMIN - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y + JMIN - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z + KMIN - 1

            IF (I>IMAX .OR. J>JMAX .OR. K>KMAX) RETURN

            SRC_IDX = 1 + I + SRC_NI*(J + SRC_NJ*K)
            DST_IDX = 1 + (I-DST_LB_I) + DST_NI*((J-DST_LB_J) + DST_NJ*(K-DST_LB_K))
            DST1(DST_IDX) = SRC1(SRC_IDX)
            DST2(DST_IDX) = SRC2(SRC_IDX)
            DST3(DST_IDX) = SRC3(SRC_IDX)
            DST4(DST_IDX) = SRC4(SRC_IDX)
            DST5(DST_IDX) = SRC5(SRC_IDX)

        END SUBROUTINE COPY_GPU_FIELDS_5X3D_KERNEL

        ATTRIBUTES(GLOBAL) SUBROUTINE COPY_GPU_FIELD_4D_KERNEL( &
            SRC,DST,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX,DST_LB_I,DST_LB_J,DST_LB_K, &
            SRC_NI,SRC_NJ,SRC_NK,DST_NI,DST_NJ,DST_NK,N_SCALARS)
            IMPLICIT NONE

            INTEGER, VALUE :: IMIN,IMAX,JMIN,JMAX,KMIN,KMAX,DST_LB_I,DST_LB_J,DST_LB_K
            INTEGER, VALUE :: SRC_NI,SRC_NJ,SRC_NK,DST_NI,DST_NJ,DST_NK,N_SCALARS
            REAL(EB), DEVICE, DIMENSION(*) :: SRC,DST
            INTEGER :: I,J,K,N,SRC_IDX,DST_IDX

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X + IMIN - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y + JMIN - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z + KMIN - 1

            IF (I>IMAX .OR. J>JMAX .OR. K>KMAX) RETURN

            DO N=1,N_SCALARS
                SRC_IDX = 1 + I + SRC_NI*(J + SRC_NJ*(K + SRC_NK*(N-1)))
                DST_IDX = 1 + (I-DST_LB_I) + DST_NI*((J-DST_LB_J) + DST_NJ*((K-DST_LB_K) + DST_NK*(N-1)))
                DST(DST_IDX) = SRC(SRC_IDX)
            ENDDO

        END SUBROUTINE COPY_GPU_FIELD_4D_KERNEL

        !--------------------------------------------------------------------------
        ! CPU/GPU debug synchronization for exchange inputs
        !--------------------------------------------------------------------------
        SUBROUTINE SYNC_GPU_EXCHANGE_INPUTS_FROM_CPU(CODE)
            IMPLICIT NONE

            INTEGER, INTENT(IN) :: CODE
            INTEGER :: NSRC

            DO NSRC=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                IF (.NOT.GPU_MESHES(NSRC)%INITIALIZED) CYCLE
                SELECT CASE(CODE)
                    CASE(1)
                        GPU_MESHES(NSRC)%RHOSD = MESHES(NSRC)%RHOS
                        GPU_MESHES(NSRC)%D_D = MESHES(NSRC)%D
                        GPU_MESHES(NSRC)%MUD = MESHES(NSRC)%MU
                        GPU_MESHES(NSRC)%KRESD = MESHES(NSRC)%KRES
                        GPU_MESHES(NSRC)%QD = MESHES(NSRC)%Q
                        GPU_MESHES(NSRC)%ZZSD = MESHES(NSRC)%ZZS
                    CASE(3)
                        GPU_MESHES(NSRC)%HSD = MESHES(NSRC)%HS
                        GPU_MESHES(NSRC)%USD = MESHES(NSRC)%US
                        GPU_MESHES(NSRC)%VSD = MESHES(NSRC)%VS
                        GPU_MESHES(NSRC)%WSD = MESHES(NSRC)%WS
                    CASE(4)
                        GPU_MESHES(NSRC)%RHOD = MESHES(NSRC)%RHO
                        GPU_MESHES(NSRC)%DS_D = MESHES(NSRC)%DS
                        GPU_MESHES(NSRC)%MUD = MESHES(NSRC)%MU
                        GPU_MESHES(NSRC)%KRESD = MESHES(NSRC)%KRES
                        GPU_MESHES(NSRC)%QD = MESHES(NSRC)%Q
                        GPU_MESHES(NSRC)%ZZD = MESHES(NSRC)%ZZ
                    CASE(5)
                        GPU_MESHES(NSRC)%FVXD = MESHES(NSRC)%FVX
                        GPU_MESHES(NSRC)%FVYD = MESHES(NSRC)%FVY
                        GPU_MESHES(NSRC)%FVZD = MESHES(NSRC)%FVZ
                        IF (PREDICTOR) THEN
                            GPU_MESHES(NSRC)%HD = MESHES(NSRC)%H
                        ELSE
                            GPU_MESHES(NSRC)%HSD = MESHES(NSRC)%HS
                        ENDIF
                    CASE(6)
                        GPU_MESHES(NSRC)%HD = MESHES(NSRC)%H
                        GPU_MESHES(NSRC)%UD = MESHES(NSRC)%U
                        GPU_MESHES(NSRC)%VD = MESHES(NSRC)%V
                        GPU_MESHES(NSRC)%WD = MESHES(NSRC)%W
                END SELECT
            ENDDO

        END SUBROUTINE SYNC_GPU_EXCHANGE_INPUTS_FROM_CPU

        !--------------------------------------------------------------------------
        ! CPU/GPU debug synchronization for exchange outputs
        !--------------------------------------------------------------------------
        SUBROUTINE SYNC_GPU_EXCHANGE_OUTPUTS_TO_CPU(CODE)
            IMPLICIT NONE

            INTEGER, INTENT(IN) :: CODE
            INTEGER :: NSRC,NN,NDST,IMIN,IMAX,JMIN,JMAX,KMIN,KMAX
            TYPE(OMESH_TYPE), POINTER :: SEND_OM,RECV_OM
            TYPE(GPU_OMESH_TYPE), POINTER :: RECV_GOM

            DO NSRC=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                IF (.NOT.ALLOCATED(MESHES(NSRC)%NEIGHBORING_MESH)) CYCLE
                DO NN=1,MESHES(NSRC)%N_NEIGHBORING_MESHES
                    NDST = MESHES(NSRC)%NEIGHBORING_MESH(NN)
                    IF (NDST < LBOUND(GPU_MESHES,1) .OR. NDST > UBOUND(GPU_MESHES,1)) CYCLE
                    IF (.NOT.GPU_MESHES(NDST)%INITIALIZED) CYCLE
                    SEND_OM => MESHES(NSRC)%OMESH(NDST)
                    IF (SEND_OM%NIC_S <= 0) CYCLE

                    RECV_OM => MESHES(NDST)%OMESH(NSRC)
                    RECV_GOM => GPU_MESHES(NDST)%OMESH(NSRC)
                    IMIN = SEND_OM%I_MIN_S
                    IMAX = SEND_OM%I_MAX_S
                    JMIN = SEND_OM%J_MIN_S
                    JMAX = SEND_OM%J_MAX_S
                    KMIN = SEND_OM%K_MIN_S
                    KMAX = SEND_OM%K_MAX_S

                    SELECT CASE(CODE)
                        CASE(1)
                            RECV_OM%RHOS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%RHOS_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%D_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%MU(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%MU_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%KRES(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%KRES_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%Q(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%Q_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%ZZS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,1:N_TOTAL_SCALARS) = &
                                RECV_GOM%ZZS_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,1:N_TOTAL_SCALARS)
                        CASE(3)
                            RECV_OM%HS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%HS_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%US(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%US_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%VS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%VS_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%WS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%WS_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                        CASE(4)
                            RECV_OM%RHO(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%RHO_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%DS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%DS_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%MU(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%MU_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%KRES(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%KRES_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%Q(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%Q_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%ZZ(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,1:N_TOTAL_SCALARS) = &
                                RECV_GOM%ZZ_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX,1:N_TOTAL_SCALARS)
                        CASE(5)
                            RECV_OM%FVX(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%FVX_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%FVY(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%FVY_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%FVZ(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%FVZ_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            IF (PREDICTOR) THEN
                                RECV_OM%H(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                    RECV_GOM%H_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            ELSE
                                RECV_OM%HS(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                    RECV_GOM%HS_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            ENDIF
                        CASE(6)
                            RECV_OM%H(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%H_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%U(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%U_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%V(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%V_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                            RECV_OM%W(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) = &
                                RECV_GOM%W_D(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX)
                    END SELECT
                ENDDO
            ENDDO

            ! In --gt mode, expose remotely received mirrors to the following CPU stages.
            IF (N_MPI_PROCESSES>1) THEN
                DO NDST=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                    IF (.NOT.GPU_MESHES(NDST)%INITIALIZED) CYCLE
                    IF (.NOT.ALLOCATED(MESHES(NDST)%NEIGHBORING_MESH)) CYCLE
                    DO NN=1,MESHES(NDST)%N_NEIGHBORING_MESHES
                        NSRC = MESHES(NDST)%NEIGHBORING_MESH(NN)
                        IF (PROCESS(NSRC)==MY_RANK) CYCLE
                        RECV_OM => MESHES(NDST)%OMESH(NSRC)
                        RECV_GOM => GPU_MESHES(NDST)%OMESH(NSRC)
                        IF (RECV_GOM%NIC_R<=0 .OR. .NOT.RECV_GOM%ACTIVE) CYCLE
                        SELECT CASE(CODE)
                            CASE(1)
                                RECV_OM%RHOS = RECV_GOM%RHOS_D
                                RECV_OM%D = RECV_GOM%D_D
                                RECV_OM%MU = RECV_GOM%MU_D
                                RECV_OM%KRES = RECV_GOM%KRES_D
                                RECV_OM%Q = RECV_GOM%Q_D
                                RECV_OM%ZZS = RECV_GOM%ZZS_D
                            CASE(3)
                                RECV_OM%HS = RECV_GOM%HS_D
                                RECV_OM%US = RECV_GOM%US_D
                                RECV_OM%VS = RECV_GOM%VS_D
                                RECV_OM%WS = RECV_GOM%WS_D
                            CASE(4)
                                RECV_OM%RHO = RECV_GOM%RHO_D
                                RECV_OM%DS = RECV_GOM%DS_D
                                RECV_OM%MU = RECV_GOM%MU_D
                                RECV_OM%KRES = RECV_GOM%KRES_D
                                RECV_OM%Q = RECV_GOM%Q_D
                                RECV_OM%ZZ = RECV_GOM%ZZ_D
                            CASE(5)
                                RECV_OM%FVX = RECV_GOM%FVX_D
                                RECV_OM%FVY = RECV_GOM%FVY_D
                                RECV_OM%FVZ = RECV_GOM%FVZ_D
                                IF (PREDICTOR) THEN
                                    RECV_OM%H = RECV_GOM%H_D
                                ELSE
                                    RECV_OM%HS = RECV_GOM%HS_D
                                ENDIF
                            CASE(6)
                                RECV_OM%H = RECV_GOM%H_D
                                RECV_OM%U = RECV_GOM%U_D
                                RECV_OM%V = RECV_GOM%V_D
                                RECV_OM%W = RECV_GOM%W_D
                        END SELECT
                    ENDDO
                ENDDO
            ENDIF

        END SUBROUTINE SYNC_GPU_EXCHANGE_OUTPUTS_TO_CPU


        SUBROUTINE GPU_MESH_EXCHANGE_HT1D_BACKSIDE
            INTEGER :: NSRC,NDST,GRDD,N_REQ,N_REQ_MAX,N_ITEMS,MPI_TAG
            TYPE(GPU_MESH_TYPE), POINTER :: GM_TARGET
            TYPE(GPU_OMESH_TYPE), POINTER :: GOM_LOCAL
            TYPE(MPI_REQUEST), ALLOCATABLE :: REQUEST(:)
            REAL(EB), DEVICE, POINTER :: TMP_G_SOURCE(:),HTC_SOURCE(:),Q_RAD_SOURCE(:),EMISS_SOURCE(:)

            IF (.NOT.ALLOCATED(GPU_MESHES) .OR. NMESHES<2) RETURN

            ! Finish wall preparation before reading another mesh's wall data.
            CALL GPU_WAIT_ALL

            N_REQ_MAX = 0
            IF (N_MPI_PROCESSES>1) THEN
                DO NDST=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                    IF (.NOT.ALLOCATED(GPU_MESHES(NDST)%OMESH)) CYCLE
                    DO NSRC=1,NMESHES
                        IF (PROCESS(NSRC)==MY_RANK) CYCLE
                        IF (GPU_MESHES(NDST)%OMESH(NSRC)%HT1D_N_RECV>0) N_REQ_MAX = N_REQ_MAX + 1
                    ENDDO
                ENDDO
                DO NSRC=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                    IF (.NOT.ALLOCATED(GPU_MESHES(NSRC)%OMESH)) CYCLE
                    DO NDST=1,NMESHES
                        IF (PROCESS(NDST)==MY_RANK) CYCLE
                        IF (GPU_MESHES(NSRC)%OMESH(NDST)%HT1D_N_SEND>0) N_REQ_MAX = N_REQ_MAX + 1
                    ENDDO
                ENDDO
            ENDIF
            IF (N_REQ_MAX>0) ALLOCATE(REQUEST(N_REQ_MAX))

            ! Receive packages use the same source-WALL ordering established
            ! by INITIALIZE_BACK_WALL_EXCHANGE.
            N_REQ = 0
            DO NDST=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                IF (.NOT.ALLOCATED(GPU_MESHES(NDST)%OMESH)) CYCLE
                DO NSRC=1,NMESHES
                    IF (PROCESS(NSRC)==MY_RANK) CYCLE
                    GOM_LOCAL => GPU_MESHES(NDST)%OMESH(NSRC)
                    N_ITEMS = GOM_LOCAL%HT1D_N_RECV
                    IF (N_ITEMS<=0) CYCLE
                    N_REQ = N_REQ + 1
                    MPI_TAG = NSRC + NMESHES*(NDST-1)
                    CALL GPU_MPI_IRECV(GOM_LOCAL%HT1D_RECV_PKG_D,HT1D_BACKSIDE_FIELD_COUNT*N_ITEMS, &
                        PROCESS(NSRC),MPI_TAG,REQUEST(N_REQ))
                ENDDO
            ENDDO

            ! Same-rank backsides remain direct device-to-device copies.
            DO NDST=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                GM_TARGET => GPU_MESHES(NDST)
                IF (.NOT.GM_TARGET%INITIALIZED .OR. GM_TARGET%HT1D_N_LINES<=0) CYCLE
                IF (.NOT.ALLOCATED(GM_TARGET%HT1D_BACKSIDE_VALID_D)) CYCLE
                GRDD = (GM_TARGET%HT1D_N_LINES+HT1D_EXCHANGE_THREADS-1)/HT1D_EXCHANGE_THREADS

                DO NSRC=1,NMESHES
                    IF (NSRC==NDST) CYCLE
                    IF (PROCESS(NSRC)/=MY_RANK) CYCLE
                    IF (.NOT.GPU_MESHES(NSRC)%INITIALIZED) CYCLE
                    TMP_G_SOURCE => GPU_MESHES(NSRC)%WALL_TMP_G
                    HTC_SOURCE   => GPU_MESHES(NSRC)%WALL_HEAT_TRANS_COEF
                    Q_RAD_SOURCE => GPU_MESHES(NSRC)%WALL_Q_RAD_IN
                    EMISS_SOURCE => GPU_MESHES(NSRC)%WALL_EMISSIVITY
                    CALL GPU_MESH_EXCHANGE_COPY_HT1D_BACKSIDE_GPU<<<GRDD,HT1D_EXCHANGE_THREADS,0,GM_TARGET%STRM>>>(GM_TARGET%HT1D_N_LINES,NSRC, &
                        GM_TARGET%HT1D_BACK_MESH_D,GM_TARGET%HT1D_BACK_INDEX_D,GM_TARGET%HT1D_BACKSIDE_VALID_D, &
                        TMP_G_SOURCE,HTC_SOURCE,Q_RAD_SOURCE,EMISS_SOURCE,GPU_MESHES(NSRC)%WALL_UB, &
                        GM_TARGET%HT1D_TMP_G_B_D,GM_TARGET%HT1D_HTCB_D, &
                        GM_TARGET%HT1D_Q_RAD_IN_B_D,GM_TARGET%HT1D_EMISSIVITY_B_D)
                ENDDO
            ENDDO

            ! Cross-rank sources are packed on the GPU from the same WALL
            ! arrays used by the local D2D path. TMP_F is deliberately not
            ! exchanged: CPU SOLID_HEAT_TRANSFER uses the target line's own
            ! TMP_B as its old solid backside temperature.
            DO NSRC=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                IF (.NOT.GPU_MESHES(NSRC)%INITIALIZED) CYCLE
                IF (.NOT.ALLOCATED(GPU_MESHES(NSRC)%OMESH)) CYCLE
                DO NDST=1,NMESHES
                    IF (PROCESS(NDST)==MY_RANK) CYCLE
                    GOM_LOCAL => GPU_MESHES(NSRC)%OMESH(NDST)
                    N_ITEMS = GOM_LOCAL%HT1D_N_SEND
                    IF (N_ITEMS<=0) CYCLE
                    GRDD = (N_ITEMS+HT1D_EXCHANGE_THREADS-1)/HT1D_EXCHANGE_THREADS
                    CALL GPU_MESH_EXCHANGE_PACK_HT1D_BACKSIDE_GPU<<<GRDD,HT1D_EXCHANGE_THREADS,0,GPU_MESHES(NSRC)%STRM>>>(N_ITEMS, &
                        GOM_LOCAL%HT1D_WALL_INDEX_S_D,GPU_MESHES(NSRC)%WALL_TMP_G,GPU_MESHES(NSRC)%WALL_HEAT_TRANS_COEF, &
                        GPU_MESHES(NSRC)%WALL_Q_RAD_IN, &
                        GPU_MESHES(NSRC)%WALL_EMISSIVITY,GPU_MESHES(NSRC)%WALL_UB,GOM_LOCAL%HT1D_SEND_PKG_D)
                ENDDO
            ENDDO

            CALL GPU_WAIT_ALL

            ! CUDA-aware MPI sends the completed device packages directly.
            DO NSRC=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                IF (.NOT.ALLOCATED(GPU_MESHES(NSRC)%OMESH)) CYCLE
                DO NDST=1,NMESHES
                    IF (PROCESS(NDST)==MY_RANK) CYCLE
                    GOM_LOCAL => GPU_MESHES(NSRC)%OMESH(NDST)
                    N_ITEMS = GOM_LOCAL%HT1D_N_SEND
                    IF (N_ITEMS<=0) CYCLE
                    N_REQ = N_REQ + 1
                    MPI_TAG = NSRC + NMESHES*(NDST-1)
                    CALL GPU_MPI_ISEND(GOM_LOCAL%HT1D_SEND_PKG_D,HT1D_BACKSIDE_FIELD_COUNT*N_ITEMS, &
                        PROCESS(NDST),MPI_TAG,REQUEST(N_REQ))
                ENDDO
            ENDDO

            IF (N_REQ>0) THEN
                CALL GPU_MPI_WAITALL(REQUEST(1:N_REQ))
            ENDIF

            DO NDST=LOWER_MESH_INDEX,UPPER_MESH_INDEX
                GM_TARGET => GPU_MESHES(NDST)
                IF (.NOT.GM_TARGET%INITIALIZED .OR. GM_TARGET%HT1D_N_LINES<=0) CYCLE
                IF (.NOT.ALLOCATED(GM_TARGET%OMESH)) CYCLE
                GRDD = (GM_TARGET%HT1D_N_LINES+HT1D_EXCHANGE_THREADS-1)/HT1D_EXCHANGE_THREADS
                DO NSRC=1,NMESHES
                    IF (PROCESS(NSRC)==MY_RANK) CYCLE
                    GOM_LOCAL => GM_TARGET%OMESH(NSRC)
                    N_ITEMS = GOM_LOCAL%HT1D_N_RECV
                    IF (N_ITEMS<=0) CYCLE
                    CALL GPU_MESH_EXCHANGE_UNPACK_HT1D_BACKSIDE_GPU<<<GRDD,HT1D_EXCHANGE_THREADS,0,GM_TARGET%STRM>>>(GM_TARGET%HT1D_N_LINES, &
                        NSRC,N_ITEMS,GM_TARGET%HT1D_BACK_MESH_D,GM_TARGET%HT1D_BACK_SLOT_D, &
                        GM_TARGET%HT1D_BACKSIDE_VALID_D,GOM_LOCAL%HT1D_RECV_PKG_D,GM_TARGET%HT1D_TMP_G_B_D, &
                        GM_TARGET%HT1D_HTCB_D,GM_TARGET%HT1D_Q_RAD_IN_B_D, &
                        GM_TARGET%HT1D_EMISSIVITY_B_D)
                ENDDO
            ENDDO

            CALL GPU_WAIT_ALL
            IF (ALLOCATED(REQUEST)) DEALLOCATE(REQUEST)
        END SUBROUTINE GPU_MESH_EXCHANGE_HT1D_BACKSIDE

        ATTRIBUTES(GLOBAL) SUBROUTINE GPU_MESH_EXCHANGE_PACK_HT1D_BACKSIDE_GPU(N_ITEMS,WALL_INDEX_D,WALL_TMP_G_D, &
            WALL_HTC_D,WALL_Q_RAD_IN_D,WALL_EMISSIVITY_D,WALL_UB,SEND_PKG_D)
            INTEGER, VALUE :: N_ITEMS,WALL_UB
            INTEGER, DEVICE :: WALL_INDEX_D(*)
            REAL(EB), DEVICE :: WALL_TMP_G_D(0:*),WALL_HTC_D(0:*),WALL_Q_RAD_IN_D(0:*), &
                                WALL_EMISSIVITY_D(0:*),SEND_PKG_D(*)

            INTEGER :: ITEM,IW

            ITEM = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
            IF (ITEM<1 .OR. ITEM>N_ITEMS) RETURN
            IW = WALL_INDEX_D(ITEM)
            IF (IW<1 .OR. IW>WALL_UB) RETURN
            SEND_PKG_D(ITEM) = WALL_TMP_G_D(IW)
            SEND_PKG_D(N_ITEMS+ITEM) = MAX(0._EB,WALL_HTC_D(IW))
            SEND_PKG_D(2*N_ITEMS+ITEM) = WALL_Q_RAD_IN_D(IW)
            SEND_PKG_D(3*N_ITEMS+ITEM) = MAX(0._EB,WALL_EMISSIVITY_D(IW))
        END SUBROUTINE GPU_MESH_EXCHANGE_PACK_HT1D_BACKSIDE_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE GPU_MESH_EXCHANGE_UNPACK_HT1D_BACKSIDE_GPU(N_LINES,SOURCE_MESH,N_ITEMS,BACK_MESH_D, &
            BACK_SLOT_D,VALID_D,RECV_PKG_D,TMP_G_B_D,HTCB_D,Q_RAD_IN_B_D,EMISSIVITY_B_D)
            INTEGER, VALUE :: N_LINES,SOURCE_MESH,N_ITEMS
            INTEGER, DEVICE :: BACK_MESH_D(*),BACK_SLOT_D(*),VALID_D(*)
            REAL(EB), DEVICE :: RECV_PKG_D(*),TMP_G_B_D(*),HTCB_D(*),Q_RAD_IN_B_D(*),EMISSIVITY_B_D(*)

            INTEGER :: LINE,SLOT

            LINE = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
            IF (LINE<1 .OR. LINE>N_LINES) RETURN
            IF (BACK_MESH_D(LINE)/=SOURCE_MESH .OR. VALID_D(LINE)==0) RETURN
            SLOT = BACK_SLOT_D(LINE)
            IF (SLOT<1 .OR. SLOT>N_ITEMS) THEN
                VALID_D(LINE) = 0
                RETURN
            ENDIF
            TMP_G_B_D(LINE) = RECV_PKG_D(SLOT)
            HTCB_D(LINE) = MAX(0._EB,RECV_PKG_D(N_ITEMS+SLOT))
            Q_RAD_IN_B_D(LINE) = RECV_PKG_D(2*N_ITEMS+SLOT)
            EMISSIVITY_B_D(LINE) = MAX(0._EB,RECV_PKG_D(3*N_ITEMS+SLOT))
        END SUBROUTINE GPU_MESH_EXCHANGE_UNPACK_HT1D_BACKSIDE_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE GPU_MESH_EXCHANGE_COPY_HT1D_BACKSIDE_GPU(N_LINES,SOURCE_MESH,BACK_MESH_D,BACK_INDEX_D,VALID_D, &
            TMP_G_SOURCE,HTC_SOURCE,Q_RAD_SOURCE,EMISS_SOURCE,SOURCE_WALL_UB,TMP_G_B_D,HTCB_D, &
            Q_RAD_IN_B_D,EMISSIVITY_B_D)
            INTEGER, VALUE :: N_LINES,SOURCE_MESH,SOURCE_WALL_UB
            INTEGER, DEVICE :: BACK_MESH_D(*),BACK_INDEX_D(*),VALID_D(*)
            REAL(EB), DEVICE :: TMP_G_SOURCE(0:*),HTC_SOURCE(0:*),Q_RAD_SOURCE(0:*),EMISS_SOURCE(0:*)
            REAL(EB), DEVICE :: TMP_G_B_D(*),HTCB_D(*),Q_RAD_IN_B_D(*),EMISSIVITY_B_D(*)

            INTEGER :: LINE,IW

            LINE = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
            IF (LINE<1 .OR. LINE>N_LINES) RETURN
            IF (BACK_MESH_D(LINE)/=SOURCE_MESH .OR. VALID_D(LINE)==0) RETURN

            IW = BACK_INDEX_D(LINE)
            IF (IW<1 .OR. IW>SOURCE_WALL_UB) THEN
                VALID_D(LINE) = 0
                RETURN
            ENDIF

            TMP_G_B_D(LINE) = TMP_G_SOURCE(IW)
            HTCB_D(LINE) = MAX(0._EB,HTC_SOURCE(IW))
            Q_RAD_IN_B_D(LINE) = Q_RAD_SOURCE(IW)
            EMISSIVITY_B_D(LINE) = MAX(0._EB,EMISS_SOURCE(IW))
        END SUBROUTINE GPU_MESH_EXCHANGE_COPY_HT1D_BACKSIDE_GPU

END MODULE G_EXCH
