!------------------------------------------------------------------------------
! Module: G_POIS
! Purpose: GPU Fishpak pressure-solver setup and FCTRD/transform solution.
! Scope: Supports Dirichlet, Neumann, mixed, and periodic pressure boundaries.
!------------------------------------------------------------------------------
MODULE G_POIS

    USE PRECISION_PARAMETERS, ONLY: EB
    USE CUDAFOR
    USE MESH_POINTERS, ONLY: IBAR, JBAR, KBAR, IBP1, JBP1, KBP1
    USE G_MESH_POINTERS
    USE G_ERRO, ONLY: GPU_WAIT
    IMPLICIT NONE
    PRIVATE

    PUBLIC :: FISHPAK_SOLVER_GPU_SAVE1_FUSED_YZ_GLOBAL_FT
    PUBLIC :: CHECK_GPU_POISSON_SETUP
    PUBLIC :: INITIALIZE_GPU_POISSON_SOLVER
    PUBLIC :: RELEASE_FISHPAK_SOLVER_GPU
    PUBLIC :: STRETCHED_POISSON_SOLVER_GPU,COPY_PRESSURE_RHS_GPU
    PUBLIC :: INITIALIZE_STRETCHED_POISSON,SOLVE_STRETCHED_POISSON

    INTEGER, PARAMETER :: BC_PERIODIC = 0
    INTEGER, PARAMETER :: BC_DD       = 1
    INTEGER, PARAMETER :: BC_DN       = 2
    INTEGER, PARAMETER :: BC_NN       = 3
    INTEGER, PARAMETER :: BC_ND       = 4
    INTEGER, PARAMETER :: PRESSURE_BC_CLASS_DIRICHLET_ONLY = 1
    INTEGER, PARAMETER :: PRESSURE_BC_CLASS_MIXED_DN       = 2
    INTEGER, PARAMETER :: PRESSURE_BC_CLASS_PURE_NEUMANN   = 3
    INTEGER, PARAMETER :: PRESSURE_BC_CLASS_PERIODIC       = 4

    INTEGER, PARAMETER :: MAX_FISHPAK_LINE = 32768
    INTEGER, PARAMETER :: REDUCE_BLOCK = 256
    INTEGER, PARAMETER :: FISHPAK_X_LINES_PER_BLOCK_DEFAULT = 16

    CONTAINS

      !----------------------------------------------------------------------------
      ! Pressure-boundary classification helpers
      !----------------------------------------------------------------------------
      INTEGER FUNCTION PRESSURE_BC_CLASS(LBC, MBC, NBC)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: LBC, MBC, NBC
        LOGICAL :: HAS_PERIODIC, HAS_DIRICHLET, HAS_NEUMANN

        HAS_PERIODIC = (LBC == BC_PERIODIC .OR. MBC == BC_PERIODIC .OR. NBC == BC_PERIODIC)
        HAS_DIRICHLET = (LBC == BC_DD .OR. LBC == BC_DN .OR. LBC == BC_ND .OR. &
                         MBC == BC_DD .OR. MBC == BC_DN .OR. MBC == BC_ND .OR. &
                         NBC == BC_DD .OR. NBC == BC_DN .OR. NBC == BC_ND)
        HAS_NEUMANN = (LBC == BC_NN .OR. LBC == BC_DN .OR. LBC == BC_ND .OR. &
                       MBC == BC_NN .OR. MBC == BC_DN .OR. MBC == BC_ND .OR. &
                       NBC == BC_NN .OR. NBC == BC_DN .OR. NBC == BC_ND)

        IF (HAS_PERIODIC) THEN
           PRESSURE_BC_CLASS = PRESSURE_BC_CLASS_PERIODIC
        ELSEIF (LBC == BC_NN .AND. MBC == BC_NN .AND. NBC == BC_NN) THEN
           PRESSURE_BC_CLASS = PRESSURE_BC_CLASS_PURE_NEUMANN
        ELSEIF (HAS_DIRICHLET .AND. HAS_NEUMANN) THEN
           PRESSURE_BC_CLASS = PRESSURE_BC_CLASS_MIXED_DN
        ELSE
           PRESSURE_BC_CLASS = PRESSURE_BC_CLASS_DIRICHLET_ONLY
        ENDIF
      END FUNCTION PRESSURE_BC_CLASS

      LOGICAL FUNCTION PRESSURE_BC_IS_SINGULAR(LBC, MBC, NBC)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: LBC, MBC, NBC
        PRESSURE_BC_IS_SINGULAR = ((LBC == BC_NN .OR. LBC == BC_PERIODIC) .AND. &
                                   (MBC == BC_NN .OR. MBC == BC_PERIODIC) .AND. &
                                   (NBC == BC_NN .OR. NBC == BC_PERIODIC))
      END FUNCTION PRESSURE_BC_IS_SINGULAR


      !----------------------------------------------------------------------------
      ! Per-mesh Fishpak setup and pressure-solver host orchestration
      !----------------------------------------------------------------------------
      SUBROUTINE CHECK_GPU_POISSON_SETUP(NM)
        USE MESH_VARIABLES, ONLY: MESHES
        INTEGER, INTENT(IN) :: NM

        IF (MIN(MESHES(NM)%LBC,MESHES(NM)%MBC,MESHES(NM)%NBC)<0 .OR. &
            MAX(MESHES(NM)%LBC,MESHES(NM)%MBC,MESHES(NM)%NBC)>4) &
            ERROR STOP 'GPU Fishpak supports pressure boundary codes 0:4'
        IF (ALLOCATED(GPU_MESHES(NM)%POIS%STRETCH)) RETURN
        IF (MESHES(NM)%LBC==0 .AND. MESHES(NM)%IBAR>MAX_FISHPAK_LINE) &
            ERROR STOP 'GPU Fishpak periodic x line exceeds MAX_FISHPAK_LINE'
      END SUBROUTINE CHECK_GPU_POISSON_SETUP

      SUBROUTINE INITIALIZE_GPU_POISSON_SOLVER(NM,DX,DY,DZ,LBC,MBC,NBC,SAVE_F)
        USE MESH_VARIABLES, ONLY: MESHES
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NM,LBC,MBC,NBC
        REAL(EB), INTENT(IN) :: DX,DY,DZ
        REAL(EB), INTENT(IN), DIMENSION(-3:) :: SAVE_F

        GM => GPU_MESHES(NM)
        IF (.NOT.ALLOCATED(GM%POIS)) ALLOCATE(GM%POIS)
        IF (MESHES(NM)%IPS>=1 .AND. MESHES(NM)%IPS<=6) THEN
           ALLOCATE(GM%POIS%STRETCH)
           CALL INITIALIZE_STRETCHED_POISSON(GM%POIS%STRETCH,[IBAR,JBAR,KBAR],MESHES(NM)%IPS, &
                [MESHES(NM)%DXI,MESHES(NM)%DETA,MESHES(NM)%DZETA],[LBC,MBC,NBC], &
                MESHES(NM)%HX,MESHES(NM)%HY,MESHES(NM)%HZ)
           RETURN
        ENDIF
        CALL ENSURE_FISHPAK_SAVE1_SETUP(NM,DX,DY,DZ,LBC,MBC,NBC,SAVE_F)
        CALL ENSURE_FISHPAK_FT_WORKSPACE(NM)
      END SUBROUTINE INITIALIZE_GPU_POISSON_SOLVER

      SUBROUTINE STRETCHED_POISSON_SOLVER_GPU(NM,H,PRHS,BXS,BXF,BYS,BYF,BZS,BZF)
        INTEGER, INTENT(IN) :: NM
        REAL(EB), DEVICE, INTENT(INOUT) :: H(0:,0:,0:)
        REAL(EB), DEVICE, INTENT(IN) :: PRHS(:,:,:),BXS(:,:),BXF(:,:),BYS(:,:),BYF(:,:),BZS(:,:),BZF(:,:)
        REAL(EB) :: DS(3)
        INTEGER :: BC(3)
        TYPE(DIM3) :: BLK,GRD

        GM => GPU_MESHES(NM)
        CALL SOLVE_STRETCHED_POISSON(GM%POIS%STRETCH,H,PRHS,BXS,BXF,BYS,BYF,BZS,BZF,GM%STRM)
        DS=GM%POIS%STRETCH%SPACING; BC=GM%POIS%STRETCH%BC
        BLK=DIM3(8,4,4); GRD=DIM3((IBP1+8)/8,(JBP1+4)/4,(KBP1+4)/4)
        CALL FILL_PRESSURE_GHOST_KERNEL<<<GRD,BLK,0,GM%STRM>>>(H,BXS,BXF,BYS,BYF,BZS,BZF, &
             IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,DS(1),DS(2),DS(3),BC(1),BC(2),BC(3))
      END SUBROUTINE STRETCHED_POISSON_SOLVER_GPU

      ! Map CPU PRHS storage order to the GPU's Cartesian array order.
      SUBROUTINE COPY_PRESSURE_RHS_GPU(NM,TO_DEVICE)
        USE MESH_VARIABLES, ONLY: MESHES
        INTEGER, INTENT(IN) :: NM
        LOGICAL, INTENT(IN) :: TO_DEVICE
        REAL(EB), ALLOCATABLE :: RHS_H(:,:,:)
        INTEGER :: I,J,K,XYZ(3),ORDER(3)

        GM => GPU_MESHES(NM)
        ORDER=[1,2,3]
        SELECT CASE(MESHES(NM)%IPS)
           CASE(2); ORDER=[2,1,3]
           CASE(3,6); ORDER=[3,2,1]
           CASE(5); ORDER=[1,3,2]
        END SELECT
        IF (ALL(ORDER==[1,2,3])) THEN
           IF (TO_DEVICE) THEN
              GM%PRHSD(1:IBAR,1:JBAR,1:KBAR)=MESHES(NM)%PRHS(1:IBAR,1:JBAR,1:KBAR)
           ELSE
              MESHES(NM)%PRHS(1:IBAR,1:JBAR,1:KBAR)=GM%PRHSD(1:IBAR,1:JBAR,1:KBAR)
           ENDIF
           RETURN
        ENDIF
        ALLOCATE(RHS_H(IBAR,JBAR,KBAR))
        IF (.NOT.TO_DEVICE) RHS_H=GM%PRHSD(1:IBAR,1:JBAR,1:KBAR)
        DO K=1,KBAR
           DO J=1,JBAR
              DO I=1,IBAR
                 XYZ=[I,J,K]
                 IF (TO_DEVICE) THEN
                    RHS_H(I,J,K)=MESHES(NM)%PRHS(XYZ(ORDER(1)),XYZ(ORDER(2)),XYZ(ORDER(3)))
                 ELSE
                    MESHES(NM)%PRHS(XYZ(ORDER(1)),XYZ(ORDER(2)),XYZ(ORDER(3)))=RHS_H(I,J,K)
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        IF (TO_DEVICE) GM%PRHSD(1:IBAR,1:JBAR,1:KBAR)=RHS_H
      END SUBROUTINE COPY_PRESSURE_RHS_GPU


      SUBROUTINE FISHPAK_SOLVER_GPU_SAVE1_FUSED_YZ_GLOBAL_FT(NM,H,PRHS,BXS,BXF,BYS,BYF,BZS,BZF, &
                                                              DX,DY,DZ,LBC,MBC,NBC,SAVE_F)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NM
        REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:KBP1), INTENT(INOUT) :: H
        REAL(EB), DEVICE, DIMENSION(1:IBP1,1:JBP1,1:KBP1), INTENT(IN) :: PRHS
        REAL(EB), DEVICE, DIMENSION(:,:), INTENT(IN) :: BXS, BXF, BYS, BYF, BZS, BZF
        REAL(EB), INTENT(IN) :: DX, DY, DZ
        INTEGER, INTENT(IN) :: LBC, MBC, NBC
        REAL(EB), INTENT(IN), DIMENSION(-3:) :: SAVE_F

        TYPE(DIM3) :: BLK, GRD, RBLK, RGRD, XBLK, XGRD
        INTEGER :: N_LINES,X_LINES_PER_BLOCK
        LOGICAL :: SINGULAR_POISSON
        INTEGER :: BC_CLASS

        GM => GPU_MESHES(NM)
        BC_CLASS = PRESSURE_BC_CLASS(LBC,MBC,NBC)
        CALL ENSURE_FISHPAK_SAVE1_SETUP(NM,DX,DY,DZ,LBC,MBC,NBC,SAVE_F)
        CALL ENSURE_FISHPAK_FT_WORKSPACE(NM)

        BLK = DIM3(8,4,4)
        GRD = DIM3((IBAR+7)/8, (JBAR+3)/4, (KBAR+3)/4)
        CALL FISHPAK_APPLY_BC_KERNEL<<<GRD,BLK,0,GM%STRM>>>(PRHS, GM%POIS%RHS_D, BXS, BXF, BYS, BYF, BZS, BZF, &
                                                   GM%POIS%AX_D, GM%POIS%CX_D, IBAR, JBAR, KBAR, IBP1, JBP1, KBP1, &
                                                   DX, DY, DZ, LBC, MBC, NBC)

        SINGULAR_POISSON = (BC_CLASS == PRESSURE_BC_CLASS_PURE_NEUMANN .OR. PRESSURE_BC_IS_SINGULAR(LBC,MBC,NBC))
        IF (SINGULAR_POISSON) THEN
           RBLK = DIM3(REDUCE_BLOCK,1,1)
           RGRD = DIM3(GM%POIS%SUM_BLOCKS,1,1)
           CALL REDUCE_INTERIOR_SUM_KERNEL<<<RGRD,RBLK,0,GM%STRM>>>(GM%POIS%RHS_D, GM%POIS%SUM_PARTIAL_D, IBAR, JBAR, KBAR)
           CALL FINALIZE_MEAN_KERNEL<<<DIM3(1,1,1),RBLK,0,GM%STRM>>>(GM%POIS%SUM_PARTIAL_D, GM%POIS%MEAN_D, GM%POIS%SUM_BLOCKS, GM%POIS%N)
           CALL SUBTRACT_MEAN_KERNEL<<<GRD,BLK,0,GM%STRM>>>(GM%POIS%RHS_D, GM%POIS%MEAN_D, IBAR, JBAR, KBAR)
        ENDIF

        CALL FUSED_FORWARD_Z_KERNEL<<<GRD,BLK,0,GM%STRM>>>(GM%POIS%RHS_D, GM%POIS%W1_D, GM%POIS%UZ_D, IBAR, JBAR, KBAR)
        CALL FUSED_FORWARD_Y_KERNEL<<<GRD,BLK,0,GM%STRM>>>(GM%POIS%W1_D, GM%POIS%RHS_D, GM%POIS%UY_D, IBAR, JBAR, KBAR)

        N_LINES = JBAR*KBAR
        X_LINES_PER_BLOCK = FISHPAK_X_LINES_PER_BLOCK_DEFAULT
        XBLK = DIM3(X_LINES_PER_BLOCK,1,1)
        XGRD = DIM3((N_LINES+X_LINES_PER_BLOCK-1)/X_LINES_PER_BLOCK,1,1)
        IF (LBC == BC_PERIODIC) THEN
           CALL SOLVE_X_LINES_FCTRD_PERIODIC_GLOBAL_FT_KERNEL<<<XGRD,XBLK,0,GM%STRM>>>(GM%POIS%RHS_D, GM%POIS%AX_D, GM%POIS%CX_D, GM%POIS%FCTRD_D, GM%POIS%FT_LINE_D, &
                                                                             IBAR, JBAR, KBAR, N_LINES)
        ELSE
           CALL SOLVE_X_LINES_FCTRD_GLOBAL_FT_KERNEL<<<XGRD,XBLK,0,GM%STRM>>>(GM%POIS%RHS_D, GM%POIS%AX_D, GM%POIS%CX_D, GM%POIS%FCTRD_D, GM%POIS%FT_LINE_D, &
                                                                    IBAR, JBAR, KBAR, N_LINES)
        ENDIF

        CALL FUSED_INVERSE_Y_KERNEL<<<GRD,BLK,0,GM%STRM>>>(GM%POIS%RHS_D, GM%POIS%W2_D, GM%POIS%UY_D, IBAR, JBAR, KBAR)
        CALL FUSED_INVERSE_Z_KERNEL<<<GRD,BLK,0,GM%STRM>>>(GM%POIS%W2_D, GM%POIS%RHS_D, GM%POIS%UZ_D, IBAR, JBAR, KBAR)

        IF (SINGULAR_POISSON) THEN
           RBLK = DIM3(REDUCE_BLOCK,1,1)
           RGRD = DIM3(GM%POIS%SUM_BLOCKS,1,1)
           CALL REDUCE_INTERIOR_SUM_KERNEL<<<RGRD,RBLK,0,GM%STRM>>>(GM%POIS%RHS_D, GM%POIS%SUM_PARTIAL_D, IBAR, JBAR, KBAR)
           CALL FINALIZE_MEAN_KERNEL<<<DIM3(1,1,1),RBLK,0,GM%STRM>>>(GM%POIS%SUM_PARTIAL_D, GM%POIS%MEAN_D, GM%POIS%SUM_BLOCKS, GM%POIS%N)
           CALL SUBTRACT_MEAN_KERNEL<<<GRD,BLK,0,GM%STRM>>>(GM%POIS%RHS_D, GM%POIS%MEAN_D, IBAR, JBAR, KBAR)
        ENDIF
        CALL SCATTER_INTERIOR_KERNEL<<<GRD,BLK,0,GM%STRM>>>(GM%POIS%RHS_D, H, IBAR, JBAR, KBAR, IBP1, JBP1, KBP1)

        GRD = DIM3((IBP1+8)/8, (JBP1+4)/4, (KBP1+4)/4)
        CALL FILL_PRESSURE_GHOST_KERNEL<<<GRD,BLK,0,GM%STRM>>>(H, BXS, BXF, BYS, BYF, BZS, BZF, &
                                                      IBAR, JBAR, KBAR, IBP1, JBP1, KBP1, DX, DY, DZ, LBC, MBC, NBC)
      END SUBROUTINE FISHPAK_SOLVER_GPU_SAVE1_FUSED_YZ_GLOBAL_FT

      SUBROUTINE ENSURE_FISHPAK_SETUP(NM,DX,DY,DZ,LBC,MBC,NBC)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NM
        REAL(EB), INTENT(IN) :: DX, DY, DZ
        INTEGER, INTENT(IN) :: LBC, MBC, NBC
        LOGICAL :: NEED_REBUILD
        REAL(EB), ALLOCATABLE :: UY_H(:,:), UZ_H(:,:), LAMY_H(:), LAMZ_H(:)
        REAL(EB), ALLOCATABLE :: AX_H(:), BX_H(:), CX_H(:)

        GM => GPU_MESHES(NM)
        NEED_REBUILD = .FALSE.
        IF (.NOT. GM%POIS%FISHPAK_READY) NEED_REBUILD = .TRUE.
        IF (GM%POIS%NX /= IBAR .OR. GM%POIS%NY /= JBAR .OR. GM%POIS%NZ /= KBAR) NEED_REBUILD = .TRUE.
        IF (DX /= GM%POIS%DX .OR. DY /= GM%POIS%DY .OR. DZ /= GM%POIS%DZ) NEED_REBUILD = .TRUE.
        IF (LBC /= GM%POIS%LBC .OR. MBC /= GM%POIS%MBC .OR. NBC /= GM%POIS%NBC) NEED_REBUILD = .TRUE.
        IF (.NOT. NEED_REBUILD) RETURN

        CALL RELEASE_FISHPAK_SOLVER_GPU(NM,RELEASE_CONTAINER=.FALSE.)

        GM%POIS%NX = IBAR
        GM%POIS%NY = JBAR
        GM%POIS%NZ = KBAR
        GM%POIS%N = IBAR*JBAR*KBAR
        GM%POIS%DX = DX
        GM%POIS%DY = DY
        GM%POIS%DZ = DZ
        GM%POIS%LBC = LBC
        GM%POIS%MBC = MBC
        GM%POIS%NBC = NBC
        GM%POIS%SUM_BLOCKS = (GM%POIS%N + REDUCE_BLOCK - 1) / REDUCE_BLOCK
        IF (GM%POIS%SUM_BLOCKS < 1) GM%POIS%SUM_BLOCKS = 1

        ALLOCATE(UY_H(JBAR,JBAR), UZ_H(KBAR,KBAR))
        ALLOCATE(LAMY_H(JBAR), LAMZ_H(KBAR))
        ALLOCATE(AX_H(IBAR), BX_H(IBAR), CX_H(IBAR))

        CALL BUILD_1D_BASIS(JBAR, DY, MBC, UY_H, LAMY_H)
        CALL BUILD_1D_BASIS(KBAR, DZ, NBC, UZ_H, LAMZ_H)
        CALL BUILD_X_COEFFICIENTS(IBAR, DX, LBC, AX_H, BX_H, CX_H)

        ALLOCATE(GM%POIS%UY_D(JBAR,JBAR),GM%POIS%UZ_D(KBAR,KBAR))
        ALLOCATE(GM%POIS%AX_D(IBAR),GM%POIS%BX_D(IBAR),GM%POIS%CX_D(IBAR))

        ALLOCATE(GM%POIS%RHS_D(IBAR,JBAR,KBAR))
        ALLOCATE(GM%POIS%W1_D(IBAR,JBAR,KBAR),GM%POIS%W2_D(IBAR,JBAR,KBAR))
        ALLOCATE(GM%POIS%SUM_PARTIAL_D(GM%POIS%SUM_BLOCKS),GM%POIS%MEAN_D(1))

        GM%POIS%UY_D = UY_H
        GM%POIS%UZ_D = UZ_H
        GM%POIS%AX_D = AX_H
        GM%POIS%BX_D = BX_H
        GM%POIS%CX_D = CX_H

        GM%POIS%FISHPAK_READY = .TRUE.
        GM%POIS%FISHPAK_SAVE1_READY = .FALSE.
        DEALLOCATE(UY_H, UZ_H, LAMY_H, LAMZ_H, AX_H, BX_H, CX_H)
      END SUBROUTINE ENSURE_FISHPAK_SETUP


      SUBROUTINE ENSURE_FISHPAK_SAVE1_SETUP(NM,DX,DY,DZ,LBC,MBC,NBC,SAVE_F)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NM
        REAL(EB), INTENT(IN) :: DX, DY, DZ
        INTEGER, INTENT(IN) :: LBC, MBC, NBC
        REAL(EB), INTENT(IN), DIMENSION(-3:) :: SAVE_F
        REAL(EB), ALLOCATABLE :: AX_H(:), CX_H(:), FCTRD_H(:,:,:)
        REAL(EB), ALLOCATABLE :: UY_H(:,:), UZ_H(:,:), LAMY_H(:), LAMZ_H(:)
        INTEGER :: I, J, K, IS_TOP, IA, IB, IC, ICFY, ICFZ, IFCTRD, LDIMFC, OFF

        CALL ENSURE_FISHPAK_SETUP(NM,DX,DY,DZ,LBC,MBC,NBC)
        IF (GM%POIS%FISHPAK_SAVE1_READY .AND. ALLOCATED(GM%POIS%FCTRD_D)) RETURN

        IF (ALLOCATED(GM%POIS%FCTRD_D)) DEALLOCATE(GM%POIS%FCTRD_D)
        ALLOCATE(AX_H(IBAR), CX_H(IBAR), FCTRD_H(IBAR,JBAR,KBAR))
        ALLOCATE(UY_H(JBAR,JBAR), UZ_H(KBAR,KBAR))
        ALLOCATE(LAMY_H(JBAR), LAMZ_H(KBAR))

        IS_TOP = 12 + 4*IBAR
        IA = 12
        IB = IA + IBAR
        IC = IB + IBAR
        ICFY = IC + IBAR
        ICFZ = ICFY + 4*JBAR
        IFCTRD = ICFZ + 4*KBAR
        LDIMFC = IBAR
        IF (IBP1 > IBAR .AND. MOD(IBAR,2) == 0) LDIMFC = IBAR + 1

        DO I = 1, IBAR
           AX_H(I) = SAVE_F(IS_TOP + IA + (I-1) + 3)
           CX_H(I) = SAVE_F(IS_TOP + IC + (I-1) + 3)
        ENDDO

        DO K = 1, KBAR
           DO J = 1, JBAR
              DO I = 1, IBAR
                 OFF = (I-1) + (K-1)*LDIMFC + (J-1)*LDIMFC*KBAR
                 FCTRD_H(I,J,K) = SAVE_F(IS_TOP + IFCTRD + OFF + 3)
              ENDDO
           ENDDO
        ENDDO

        CALL BUILD_1D_BASIS(JBAR, DY, MBC, UY_H, LAMY_H)
        CALL BUILD_1D_BASIS(KBAR, DZ, NBC, UZ_H, LAMZ_H)

        GM%POIS%UY_D = UY_H
        GM%POIS%UZ_D = UZ_H
        GM%POIS%AX_D = AX_H
        GM%POIS%CX_D = CX_H
        ALLOCATE(GM%POIS%FCTRD_D(IBAR,JBAR,KBAR))
        GM%POIS%FCTRD_D = FCTRD_H
        GM%POIS%FISHPAK_SAVE1_READY = .TRUE.

        DEALLOCATE(AX_H, CX_H, FCTRD_H, UY_H, UZ_H, LAMY_H, LAMZ_H)
      END SUBROUTINE ENSURE_FISHPAK_SAVE1_SETUP


      SUBROUTINE BUILD_X_COEFFICIENTS(N, D, BC, A, B, C)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N, BC
        REAL(EB), INTENT(IN) :: D
        REAL(EB), INTENT(OUT) :: A(N), B(N), C(N)
        REAL(EB) :: D2I

        D2I = 1._EB/(D*D)
        A = D2I
        B = -2._EB*D2I
        C = D2I

        SELECT CASE(BC)
           CASE(BC_DD)
              B(1) = B(1) - A(1)
              B(N) = B(N) - C(N)
           CASE(BC_DN)
              B(1) = B(1) - A(1)
              B(N) = B(N) + C(N)
           CASE(BC_NN)
              B(1) = B(1) + A(1)
              B(N) = B(N) + C(N)
           CASE(BC_ND)
              B(1) = B(1) + A(1)
              B(N) = B(N) - C(N)
        END SELECT
      END SUBROUTINE BUILD_X_COEFFICIENTS


      SUBROUTINE BUILD_1D_BASIS(N, D, BC, U, LAM)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N, BC
        REAL(EB), INTENT(IN) :: D
        REAL(EB), INTENT(OUT) :: U(N,N), LAM(N)
        INTEGER :: I, M
        REAL(EB) :: PI, THETA, ANG, NORM0, NORM1, NORMP

        PI = ACOS(-1._EB)
        NORM0 = 1._EB/SQRT(REAL(N,EB))
        NORM1 = SQRT(2._EB/REAL(N,EB))
        NORMP = NORM1
        U = 0._EB
        LAM = 0._EB

        DO M = 1, N
           SELECT CASE(BC)
              CASE(BC_PERIODIC)
                 IF (M == 1) THEN
                    DO I = 1, N
                       U(I,M) = NORM0
                    ENDDO
                    THETA = 0._EB
                 ELSEIF (MOD(N,2) == 0 .AND. M == N) THEN
                    THETA = PI
                    DO I = 1, N
                       ANG = THETA*REAL(I-1,EB)
                       U(I,M) = NORM0*COS(ANG)
                    ENDDO
                 ELSEIF (MOD(M,2) == 0) THEN
                    THETA = 2._EB*PI*REAL(M/2,EB)/REAL(N,EB)
                    DO I = 1, N
                       ANG = THETA*REAL(I-1,EB)
                       U(I,M) = NORMP*SIN(ANG)
                    ENDDO
                 ELSE
                    THETA = 2._EB*PI*REAL((M-1)/2,EB)/REAL(N,EB)
                    DO I = 1, N
                       ANG = THETA*REAL(I-1,EB)
                       U(I,M) = NORMP*COS(ANG)
                    ENDDO
                 ENDIF
              CASE(BC_NN)
                 THETA = PI*REAL(M-1,EB)/REAL(N,EB)
                 IF (M == 1) THEN
                    DO I = 1, N
                       U(I,M) = NORM0
                    ENDDO
                 ELSE
                    DO I = 1, N
                       ANG = THETA*(REAL(I,EB)-0.5_EB)
                       U(I,M) = NORM1*COS(ANG)
                    ENDDO
                 ENDIF
              CASE(BC_DD)
                 THETA = PI*REAL(M,EB)/REAL(N,EB)
                 DO I = 1, N
                    ANG = THETA*(REAL(I,EB)-0.5_EB)
                    IF (M == N) THEN
                       U(I,M) = NORM0*SIN(ANG)
                    ELSE
                       U(I,M) = NORM1*SIN(ANG)
                    ENDIF
                 ENDDO
              CASE(BC_DN)
                 THETA = PI*(REAL(M,EB)-0.5_EB)/REAL(N,EB)
                 DO I = 1, N
                    ANG = THETA*(REAL(I,EB)-0.5_EB)
                    U(I,M) = NORM1*SIN(ANG)
                 ENDDO
              CASE(BC_ND)
                 THETA = PI*(REAL(M,EB)-0.5_EB)/REAL(N,EB)
                 DO I = 1, N
                    ANG = THETA*(REAL(I,EB)-0.5_EB)
                    U(I,M) = NORM1*COS(ANG)
                 ENDDO
           END SELECT
           LAM(M) = -(4._EB/(D*D))*SIN(0.5_EB*THETA)**2
        ENDDO
      END SUBROUTINE BUILD_1D_BASIS


      !----------------------------------------------------------------------------
      ! Fishpak CUDA kernels: boundary application, transforms, line solve, and copy
      !----------------------------------------------------------------------------
      ATTRIBUTES(GLOBAL) SUBROUTINE FISHPAK_APPLY_BC_KERNEL(PRHS, RHS, BXS, BXF, BYS, BYF, BZS, BZF, &
                                                             AX, CX, NX, NY, NZ, PX, PY, PZ, DX, DY, DZ, LBC, MBC, NBC)
        IMPLICIT NONE
        INTEGER, VALUE :: NX, NY, NZ, PX, PY, PZ, LBC, MBC, NBC
        REAL(EB), VALUE :: DX, DY, DZ
        REAL(EB), DEVICE :: PRHS(1:PX,1:PY,1:PZ), RHS(NX,NY,NZ)
        REAL(EB), DEVICE :: BXS(1:MERGE(PY,1,NY>1),1:*), BXF(1:MERGE(PY,1,NY>1),1:*)
        REAL(EB), DEVICE :: BYS(1:PX,1:*), BYF(1:PX,1:*), BZS(1:PX,1:*), BZF(1:PX,1:*)
        REAL(EB), DEVICE :: AX(NX), CX(NX)
        INTEGER :: I, J, K
        REAL(EB) :: VAL, DYRCP, DZRCP, TWDYSQ, TWDZSQ

        I = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
        J = (BLOCKIDX%Y-1)*BLOCKDIM%Y + THREADIDX%Y
        K = (BLOCKIDX%Z-1)*BLOCKDIM%Z + THREADIDX%Z
        IF (I < 1 .OR. I > NX .OR. J < 1 .OR. J > NY .OR. K < 1 .OR. K > NZ) RETURN

        DYRCP = 1._EB/DY
        DZRCP = 1._EB/DZ
        TWDYSQ = 2._EB/(DY*DY)
        TWDZSQ = 2._EB/(DZ*DZ)
        VAL = PRHS(I,J,K)

        IF (I == 1) THEN
           IF (LBC == BC_DD .OR. LBC == BC_DN) VAL = VAL - 2._EB*BXS(J,K)*AX(1)
           IF (LBC == BC_NN .OR. LBC == BC_ND) VAL = VAL + AX(1)*DX*BXS(J,K)
        ENDIF
        IF (I == NX) THEN
           IF (LBC == BC_DD .OR. LBC == BC_ND) VAL = VAL - 2._EB*BXF(J,K)*CX(NX)
           IF (LBC == BC_DN .OR. LBC == BC_NN) VAL = VAL - CX(NX)*DX*BXF(J,K)
        ENDIF

        IF (J == 1) THEN
           IF (MBC == BC_DD .OR. MBC == BC_DN) VAL = VAL - BYS(I,K)*TWDYSQ
           IF (MBC == BC_NN .OR. MBC == BC_ND) VAL = VAL + BYS(I,K)*DYRCP
        ENDIF
        IF (J == NY) THEN
           IF (MBC == BC_DD .OR. MBC == BC_ND) VAL = VAL - BYF(I,K)*TWDYSQ
           IF (MBC == BC_DN .OR. MBC == BC_NN) VAL = VAL - BYF(I,K)*DYRCP
        ENDIF

        IF (K == 1) THEN
           IF (NBC == BC_DD .OR. NBC == BC_DN) VAL = VAL - BZS(I,J)*TWDZSQ
           IF (NBC == BC_NN .OR. NBC == BC_ND) VAL = VAL + BZS(I,J)*DZRCP
        ENDIF
        IF (K == NZ) THEN
           IF (NBC == BC_DD .OR. NBC == BC_ND) VAL = VAL - BZF(I,J)*TWDZSQ
           IF (NBC == BC_DN .OR. NBC == BC_NN) VAL = VAL - BZF(I,J)*DZRCP
        ENDIF

        RHS(I,J,K) = VAL
      END SUBROUTINE FISHPAK_APPLY_BC_KERNEL



      ATTRIBUTES(GLOBAL) SUBROUTINE SOLVE_X_LINES_FCTRD_GLOBAL_FT_KERNEL(F, AX, CX, FCTRD, FT, NX, NY, NZ, N_LINES)
        IMPLICIT NONE
        INTEGER, VALUE :: NX, NY, NZ, N_LINES
        REAL(EB), DEVICE :: F(NX,NY,NZ), AX(NX), CX(NX), FCTRD(NX,NY,NZ), FT(NX,NY,NZ)
        INTEGER :: I, J, K, LINE

        LINE = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
        IF (LINE < 1 .OR. LINE > N_LINES) RETURN

        J = MOD(LINE-1,NY) + 1
        K = (LINE-1)/NY + 1

        FT(1,J,K) = F(1,J,K)*FCTRD(1,J,K)
        DO I = 2, NX
           FT(I,J,K) = (F(I,J,K) - AX(I)*FT(I-1,J,K))*FCTRD(I,J,K)
        ENDDO

        DO I = NX-1, 1, -1
           FT(I,J,K) = FT(I,J,K) - CX(I)*FCTRD(I,J,K)*FT(I+1,J,K)
        ENDDO

        DO I = 1, NX
           F(I,J,K) = FT(I,J,K)
        ENDDO
      END SUBROUTINE SOLVE_X_LINES_FCTRD_GLOBAL_FT_KERNEL


      ATTRIBUTES(GLOBAL) SUBROUTINE SOLVE_X_LINES_FCTRD_PERIODIC_GLOBAL_FT_KERNEL(F, AX, CX, FCTRD, FT, NX, NY, NZ, N_LINES)
        IMPLICIT NONE
        INTEGER, VALUE :: NX, NY, NZ, N_LINES
        REAL(EB), DEVICE :: F(NX,NY,NZ), AX(NX), CX(NX), FCTRD(NX,NY,NZ), FT(NX,NY,NZ)
        REAL(EB) :: WORK(MAX_FISHPAK_LINE)
        INTEGER :: I, J, K, LINE, LH, LQ
        REAL(EB) :: T1, T2

        LINE = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
        IF (LINE < 1 .OR. LINE > N_LINES) RETURN
        IF (NX > MAX_FISHPAK_LINE) RETURN

        J = MOD(LINE-1,NY) + 1
        K = (LINE-1)/NY + 1
        LH = (NX+1)/2
        LQ = (LH-1)/2

        DO I = 1, NX
           WORK(I) = F(I,J,K)
        ENDDO

        DO I = 1, LQ
           WORK(I)       = WORK(I)       - WORK(2*LH-I)
           WORK(LH-I)    = WORK(LH-I)    - WORK(LH+I)
           WORK(LH+I)    = WORK(LH-I)    + 2._EB*WORK(LH+I)
           WORK(2*LH-I)  = WORK(I)       + 2._EB*WORK(2*LH-I)
           T1 = WORK(I)
           WORK(I) = WORK(LH-I)
           WORK(LH-I) = T1
        ENDDO
        WORK(LH) = 2._EB*WORK(LH)
        IF (MOD(NX,2) == 0) WORK(NX) = 2._EB*WORK(NX)
        IF (MOD(LH,2) == 0) THEN
           T1 = WORK(LH/2)
           T2 = WORK(3*LH/2)
           WORK(LH/2)   = T1 - T2
           WORK(3*LH/2) = T1 + 2._EB*T2
        ENDIF

        FT(1,J,K) = 0.5_EB*WORK(1)*FCTRD(1,J,K)
        DO I = 2, NX
           FT(I,J,K) = (0.5_EB*WORK(I) - AX(I)*FT(I-1,J,K))*FCTRD(I,J,K)
        ENDDO

        DO I = NX-1, 1, -1
           FT(I,J,K) = FT(I,J,K) - CX(I)*FCTRD(I,J,K)*FT(I+1,J,K)
        ENDDO

        DO I = 1, LQ
           FT(I,J,K)       = FT(LH+I,J,K)   + FT(I,J,K)
           FT(LH-I,J,K)    = FT(2*LH-I,J,K) + FT(LH-I,J,K)
           FT(2*LH-I,J,K)  = 2._EB*FT(2*LH-I,J,K) - FT(LH-I,J,K)
           FT(LH+I,J,K)    = 2._EB*FT(LH+I,J,K)   - FT(I,J,K)
           T1 = FT(I,J,K)
           FT(I,J,K) = FT(LH-I,J,K)
           FT(LH-I,J,K) = T1
        ENDDO
        IF (MOD(LH,2) == 0) THEN
           T1 = FT(LH/2,J,K)
           T2 = FT(3*LH/2,J,K)
           FT(LH/2,J,K)   = T2 + T1
           FT(3*LH/2,J,K) = 2._EB*T2 - FT(LH/2,J,K)
        ENDIF

        DO I = 1, NX
           F(I,J,K) = FT(I,J,K)
        ENDDO
      END SUBROUTINE SOLVE_X_LINES_FCTRD_PERIODIC_GLOBAL_FT_KERNEL



      ATTRIBUTES(GLOBAL) SUBROUTINE REDUCE_INTERIOR_SUM_KERNEL(F, PARTIAL, NX, NY, NZ)
        IMPLICIT NONE
        INTEGER, VALUE :: NX, NY, NZ
        REAL(EB), DEVICE :: F(NX,NY,NZ), PARTIAL(*)
        REAL(EB), SHARED :: SDATA(REDUCE_BLOCK)
        INTEGER :: TID, BID, GTID, STRIDE, IDX, OFS, TMP, I, J, K, NALL
        REAL(EB) :: LOCAL

        TID = THREADIDX%X
        BID = BLOCKIDX%X
        GTID = (BID-1)*BLOCKDIM%X + TID
        STRIDE = GRIDDIM%X*BLOCKDIM%X
        NALL = NX*NY*NZ
        LOCAL = 0._EB
        IDX = GTID
        DO WHILE (IDX <= NALL)
           OFS = IDX - 1
           I = MOD(OFS,NX) + 1
           TMP = OFS/NX
           J = MOD(TMP,NY) + 1
           K = TMP/NY + 1
           LOCAL = LOCAL + F(I,J,K)
           IDX = IDX + STRIDE
        ENDDO
        SDATA(TID) = LOCAL
        CALL SYNCTHREADS()
        OFS = BLOCKDIM%X/2
        DO WHILE (OFS >= 1)
           IF (TID <= OFS) SDATA(TID) = SDATA(TID) + SDATA(TID+OFS)
           CALL SYNCTHREADS()
           OFS = OFS/2
        ENDDO
        IF (TID == 1) PARTIAL(BID) = SDATA(1)
      END SUBROUTINE REDUCE_INTERIOR_SUM_KERNEL


      ATTRIBUTES(GLOBAL) SUBROUTINE FINALIZE_MEAN_KERNEL(PARTIAL, MEAN, NBLOCKS, NTOTAL)
        IMPLICIT NONE
        INTEGER, VALUE :: NBLOCKS, NTOTAL
        REAL(EB), DEVICE :: PARTIAL(NBLOCKS), MEAN(1)
        REAL(EB), SHARED :: SDATA(REDUCE_BLOCK)
        INTEGER :: TID, IDX, OFS
        REAL(EB) :: LOCAL

        TID = THREADIDX%X
        LOCAL = 0._EB
        IDX = TID
        DO WHILE (IDX <= NBLOCKS)
           LOCAL = LOCAL + PARTIAL(IDX)
           IDX = IDX + BLOCKDIM%X
        ENDDO
        SDATA(TID) = LOCAL
        CALL SYNCTHREADS()
        OFS = BLOCKDIM%X/2
        DO WHILE (OFS >= 1)
           IF (TID <= OFS) SDATA(TID) = SDATA(TID) + SDATA(TID+OFS)
           CALL SYNCTHREADS()
           OFS = OFS/2
        ENDDO
        IF (TID == 1) MEAN(1) = SDATA(1)/REAL(NTOTAL,EB)
      END SUBROUTINE FINALIZE_MEAN_KERNEL


      ATTRIBUTES(GLOBAL) SUBROUTINE SUBTRACT_MEAN_KERNEL(F, MEAN, NX, NY, NZ)
        IMPLICIT NONE
        INTEGER, VALUE :: NX, NY, NZ
        REAL(EB), DEVICE :: F(NX,NY,NZ), MEAN(1)
        INTEGER :: I, J, K
        I = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
        J = (BLOCKIDX%Y-1)*BLOCKDIM%Y + THREADIDX%Y
        K = (BLOCKIDX%Z-1)*BLOCKDIM%Z + THREADIDX%Z
        IF (I < 1 .OR. I > NX .OR. J < 1 .OR. J > NY .OR. K < 1 .OR. K > NZ) RETURN
        F(I,J,K) = F(I,J,K) - MEAN(1)
      END SUBROUTINE SUBTRACT_MEAN_KERNEL


      ATTRIBUTES(GLOBAL) SUBROUTINE FUSED_FORWARD_Z_KERNEL(IN_D, OUT_D, UZ, NX, NY, NZ)
        IMPLICIT NONE
        INTEGER, VALUE :: NX, NY, NZ
        REAL(EB), DEVICE :: IN_D(NX,NY,NZ), OUT_D(NX,NY,NZ), UZ(NZ,NZ)
        INTEGER :: I, J, PZ, K
        REAL(EB) :: S

        I  = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
        J  = (BLOCKIDX%Y-1)*BLOCKDIM%Y + THREADIDX%Y
        PZ = (BLOCKIDX%Z-1)*BLOCKDIM%Z + THREADIDX%Z
        IF (I < 1 .OR. I > NX .OR. J < 1 .OR. J > NY .OR. PZ < 1 .OR. PZ > NZ) RETURN

        S = 0._EB
        DO K = 1, NZ
           S = S + UZ(K,PZ)*IN_D(I,J,K)
        ENDDO
        OUT_D(I,J,PZ) = S
      END SUBROUTINE FUSED_FORWARD_Z_KERNEL


      ATTRIBUTES(GLOBAL) SUBROUTINE FUSED_FORWARD_Y_KERNEL(IN_D, OUT_D, UY, NX, NY, NZ)
        IMPLICIT NONE
        INTEGER, VALUE :: NX, NY, NZ
        REAL(EB), DEVICE :: IN_D(NX,NY,NZ), OUT_D(NX,NY,NZ), UY(NY,NY)
        INTEGER :: I, PY, K, J
        REAL(EB) :: S

        I  = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
        PY = (BLOCKIDX%Y-1)*BLOCKDIM%Y + THREADIDX%Y
        K  = (BLOCKIDX%Z-1)*BLOCKDIM%Z + THREADIDX%Z
        IF (I < 1 .OR. I > NX .OR. PY < 1 .OR. PY > NY .OR. K < 1 .OR. K > NZ) RETURN

        S = 0._EB
        DO J = 1, NY
           S = S + UY(J,PY)*IN_D(I,J,K)
        ENDDO
        OUT_D(I,PY,K) = S
      END SUBROUTINE FUSED_FORWARD_Y_KERNEL


      ATTRIBUTES(GLOBAL) SUBROUTINE FUSED_INVERSE_Z_KERNEL(IN_D, OUT_D, UZ, NX, NY, NZ)
        IMPLICIT NONE
        INTEGER, VALUE :: NX, NY, NZ
        REAL(EB), DEVICE :: IN_D(NX,NY,NZ), OUT_D(NX,NY,NZ), UZ(NZ,NZ)
        INTEGER :: I, J, K, PZ
        REAL(EB) :: S

        I = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
        J = (BLOCKIDX%Y-1)*BLOCKDIM%Y + THREADIDX%Y
        K = (BLOCKIDX%Z-1)*BLOCKDIM%Z + THREADIDX%Z
        IF (I < 1 .OR. I > NX .OR. J < 1 .OR. J > NY .OR. K < 1 .OR. K > NZ) RETURN

        S = 0._EB
        DO PZ = 1, NZ
           S = S + UZ(K,PZ)*IN_D(I,J,PZ)
        ENDDO
        OUT_D(I,J,K) = S
      END SUBROUTINE FUSED_INVERSE_Z_KERNEL


      ATTRIBUTES(GLOBAL) SUBROUTINE FUSED_INVERSE_Y_KERNEL(IN_D, OUT_D, UY, NX, NY, NZ)
        IMPLICIT NONE
        INTEGER, VALUE :: NX, NY, NZ
        REAL(EB), DEVICE :: IN_D(NX,NY,NZ), OUT_D(NX,NY,NZ), UY(NY,NY)
        INTEGER :: I, J, K, PY
        REAL(EB) :: S

        I = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
        J = (BLOCKIDX%Y-1)*BLOCKDIM%Y + THREADIDX%Y
        K = (BLOCKIDX%Z-1)*BLOCKDIM%Z + THREADIDX%Z
        IF (I < 1 .OR. I > NX .OR. J < 1 .OR. J > NY .OR. K < 1 .OR. K > NZ) RETURN

        S = 0._EB
        DO PY = 1, NY
           S = S + UY(J,PY)*IN_D(I,PY,K)
        ENDDO
        OUT_D(I,J,K) = S
      END SUBROUTINE FUSED_INVERSE_Y_KERNEL


      SUBROUTINE ENSURE_FISHPAK_FT_WORKSPACE(NM)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NM

        GM => GPU_MESHES(NM)
        IF (ALLOCATED(GM%POIS%FT_LINE_D)) THEN
           IF (SIZE(GM%POIS%FT_LINE_D,1) == IBAR .AND. SIZE(GM%POIS%FT_LINE_D,2) == JBAR .AND. &
               SIZE(GM%POIS%FT_LINE_D,3) == KBAR) RETURN
           DEALLOCATE(GM%POIS%FT_LINE_D)
        ENDIF
        ALLOCATE(GM%POIS%FT_LINE_D(IBAR,JBAR,KBAR))
      END SUBROUTINE ENSURE_FISHPAK_FT_WORKSPACE


      ATTRIBUTES(GLOBAL) SUBROUTINE SCATTER_INTERIOR_KERNEL(SOL, H, NX, NY, NZ, GX, GY, GZ)
        IMPLICIT NONE
        INTEGER, VALUE :: NX, NY, NZ, GX, GY, GZ
        REAL(EB), DEVICE :: SOL(NX,NY,NZ), H(0:GX,0:GY,0:GZ)
        INTEGER :: I, J, K
        I = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
        J = (BLOCKIDX%Y-1)*BLOCKDIM%Y + THREADIDX%Y
        K = (BLOCKIDX%Z-1)*BLOCKDIM%Z + THREADIDX%Z
        IF (I < 1 .OR. I > NX .OR. J < 1 .OR. J > NY .OR. K < 1 .OR. K > NZ) RETURN
        H(I,J,K) = SOL(I,J,K)
      END SUBROUTINE SCATTER_INTERIOR_KERNEL


      ATTRIBUTES(GLOBAL) SUBROUTINE FILL_PRESSURE_GHOST_KERNEL(H, BXS, BXF, BYS, BYF, BZS, BZF, &
                                                                NX, NY, NZ, GX, GY, GZ, DX, DY, DZ, LBC, MBC, NBC)
        IMPLICIT NONE
        INTEGER, VALUE :: NX, NY, NZ, GX, GY, GZ, LBC, MBC, NBC
        REAL(EB), VALUE :: DX, DY, DZ
        REAL(EB), DEVICE :: H(0:GX,0:GY,0:GZ)
        REAL(EB), DEVICE :: BXS(1:MERGE(GY,1,NY>1),1:*), BXF(1:MERGE(GY,1,NY>1),1:*)
        REAL(EB), DEVICE :: BYS(1:GX,1:*), BYF(1:GX,1:*), BZS(1:GX,1:*), BZF(1:GX,1:*)
        INTEGER :: I, J, K

        I = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X - 1
        J = (BLOCKIDX%Y-1)*BLOCKDIM%Y + THREADIDX%Y - 1
        K = (BLOCKIDX%Z-1)*BLOCKDIM%Z + THREADIDX%Z - 1
        IF (I < 0 .OR. I > GX .OR. J < 0 .OR. J > GY .OR. K < 0 .OR. K > GZ) RETURN

        IF (I == 0 .AND. J >= 1 .AND. J <= NY .AND. K >= 1 .AND. K <= NZ) THEN
           IF (LBC == BC_PERIODIC) H(0,J,K) = H(NX,J,K)
           IF (LBC == BC_NN .OR. LBC == BC_ND) H(0,J,K) = H(1,J,K) - DX*BXS(J,K)
           IF (LBC == BC_DD .OR. LBC == BC_DN) H(0,J,K) = -H(1,J,K) + 2._EB*BXS(J,K)
        ENDIF
        IF (I == GX .AND. J >= 1 .AND. J <= NY .AND. K >= 1 .AND. K <= NZ) THEN
           IF (LBC == BC_PERIODIC) H(GX,J,K) = H(1,J,K)
           IF (LBC == BC_NN .OR. LBC == BC_DN) H(GX,J,K) = H(NX,J,K) + DX*BXF(J,K)
           IF (LBC == BC_DD .OR. LBC == BC_ND) H(GX,J,K) = -H(NX,J,K) + 2._EB*BXF(J,K)
        ENDIF

        IF (J == 0 .AND. I >= 1 .AND. I <= NX .AND. K >= 1 .AND. K <= NZ) THEN
           IF (MBC == BC_PERIODIC) H(I,0,K) = H(I,NY,K)
           IF (MBC == BC_NN .OR. MBC == BC_ND) H(I,0,K) = H(I,1,K) - DY*BYS(I,K)
           IF (MBC == BC_DD .OR. MBC == BC_DN) H(I,0,K) = -H(I,1,K) + 2._EB*BYS(I,K)
        ENDIF
        IF (J == GY .AND. I >= 1 .AND. I <= NX .AND. K >= 1 .AND. K <= NZ) THEN
           IF (MBC == BC_PERIODIC) H(I,GY,K) = H(I,1,K)
           IF (MBC == BC_NN .OR. MBC == BC_DN) H(I,GY,K) = H(I,NY,K) + DY*BYF(I,K)
           IF (MBC == BC_DD .OR. MBC == BC_ND) H(I,GY,K) = -H(I,NY,K) + 2._EB*BYF(I,K)
        ENDIF

        IF (K == 0 .AND. I >= 1 .AND. I <= NX .AND. J >= 1 .AND. J <= NY) THEN
           IF (NBC == BC_PERIODIC) H(I,J,0) = H(I,J,NZ)
           IF (NBC == BC_NN .OR. NBC == BC_ND) H(I,J,0) = H(I,J,1) - DZ*BZS(I,J)
           IF (NBC == BC_DD .OR. NBC == BC_DN) H(I,J,0) = -H(I,J,1) + 2._EB*BZS(I,J)
        ENDIF
        IF (K == GZ .AND. I >= 1 .AND. I <= NX .AND. J >= 1 .AND. J <= NY) THEN
           IF (NBC == BC_PERIODIC) H(I,J,GZ) = H(I,J,1)
           IF (NBC == BC_NN .OR. NBC == BC_DN) H(I,J,GZ) = H(I,J,NZ) + DZ*BZF(I,J)
           IF (NBC == BC_DD .OR. NBC == BC_ND) H(I,J,GZ) = -H(I,J,NZ) + 2._EB*BZF(I,J)
        ENDIF
      END SUBROUTINE FILL_PRESSURE_GHOST_KERNEL


      !----------------------------------------------------------------------------
      ! Fishpak-state release
      !----------------------------------------------------------------------------
      SUBROUTINE RELEASE_FISHPAK_SOLVER_GPU(NM,RELEASE_CONTAINER)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NM
        LOGICAL, INTENT(IN), OPTIONAL :: RELEASE_CONTAINER
        LOGICAL :: RELEASE_PARENT

        RELEASE_PARENT = .TRUE.
        IF (PRESENT(RELEASE_CONTAINER)) RELEASE_PARENT = RELEASE_CONTAINER

        GM => GPU_MESHES(NM)
        IF (.NOT.ALLOCATED(GM%POIS)) RETURN

        IF (ALLOCATED(GM%POIS%STRETCH)) DEALLOCATE(GM%POIS%STRETCH)
        IF (ALLOCATED(GM%POIS%UY_D)) DEALLOCATE(GM%POIS%UY_D)
        IF (ALLOCATED(GM%POIS%UZ_D)) DEALLOCATE(GM%POIS%UZ_D)
        IF (ALLOCATED(GM%POIS%AX_D)) DEALLOCATE(GM%POIS%AX_D)
        IF (ALLOCATED(GM%POIS%BX_D)) DEALLOCATE(GM%POIS%BX_D)
        IF (ALLOCATED(GM%POIS%CX_D)) DEALLOCATE(GM%POIS%CX_D)
        IF (ALLOCATED(GM%POIS%FCTRD_D)) DEALLOCATE(GM%POIS%FCTRD_D)
        IF (ALLOCATED(GM%POIS%FT_LINE_D)) DEALLOCATE(GM%POIS%FT_LINE_D)
        IF (ALLOCATED(GM%POIS%RHS_D)) DEALLOCATE(GM%POIS%RHS_D)
        IF (ALLOCATED(GM%POIS%W1_D)) DEALLOCATE(GM%POIS%W1_D)
        IF (ALLOCATED(GM%POIS%W2_D)) DEALLOCATE(GM%POIS%W2_D)
        IF (ALLOCATED(GM%POIS%SUM_PARTIAL_D)) DEALLOCATE(GM%POIS%SUM_PARTIAL_D)
        IF (ALLOCATED(GM%POIS%MEAN_D)) DEALLOCATE(GM%POIS%MEAN_D)

        GM%POIS%FISHPAK_READY = .FALSE.
        GM%POIS%FISHPAK_SAVE1_READY = .FALSE.
        GM%POIS%NX = -1; GM%POIS%NY = -1; GM%POIS%NZ = -1
        GM%POIS%LBC = -1; GM%POIS%MBC = -1; GM%POIS%NBC = -1
        GM%POIS%N = 0; GM%POIS%SUM_BLOCKS = 0
        GM%POIS%DX = -1._EB; GM%POIS%DY = -1._EB; GM%POIS%DZ = -1._EB

        IF (RELEASE_PARENT) DEALLOCATE(GM%POIS)
      END SUBROUTINE RELEASE_FISHPAK_SOLVER_GPU

      !----------------------------------------------------------------------------
      ! Stretched-grid transforms and tridiagonal line solves
      !----------------------------------------------------------------------------
      SUBROUTINE INITIALIZE_STRETCHED_POISSON(S,XYZ,IPS,DS,BC,HX,HY,HZ)
         TYPE(GPU_STRETCH_POIS_TYPE), INTENT(OUT) :: S
         INTEGER, INTENT(IN) :: XYZ(3),IPS,BC(3)
         REAL(EB), INTENT(IN) :: DS(3),HX(0:),HY(0:),HZ(0:)
         REAL(EB), ALLOCATABLE :: METRIC(:,:),A(:,:),B(:,:),C(:,:),F(:,:),V(:,:),LAM2(:),LAM3(:)
         REAL(EB), DEVICE, ALLOCATABLE :: B_D(:),LAM2_D(:),LAM3_D(:)
         REAL(EB) :: FAC(2,3)
         INTEGER :: ORDER(3),N1,N2,N3,ND,D

         ! The first axis supplies independent tridiagonal lines. The other two are transformed.
         SELECT CASE(IPS)
            CASE(1,4); ORDER=[1,2,3]
            CASE(2);   ORDER=[2,1,3]
            CASE(3,6); ORDER=[3,2,1]
            CASE(5);   ORDER=[1,3,2]
            CASE DEFAULT; ERROR STOP 'Stretched GPU Poisson: invalid IPS'
         END SELECT
         IF (BC(ORDER(1))==0) ERROR STOP 'Stretched GPU Poisson: periodic line axis'
         S%AXIS=ORDER; S%N=XYZ(ORDER); S%XYZ=XYZ; S%BC=BC; S%SPACING=DS
         S%SINGULAR=ALL(BC==0 .OR. BC==3)
         N1=S%N(1); N2=S%N(2); N3=S%N(3)
         ND=MAXVAL(XYZ)
         ALLOCATE(METRIC(0:ND+1,3),A(ND,3),B(ND,3),C(ND,3))
         METRIC=1._EB
         METRIC(0:XYZ(1)+1,1)=HX(0:XYZ(1)+1)
         METRIC(0:XYZ(2)+1,2)=HY(0:XYZ(2)+1)
         METRIC(0:XYZ(3)+1,3)=HZ(0:XYZ(3)+1)
         DO D=1,3
            CALL AXIS_COEFFICIENTS(XYZ(D),DS(D),BC(D),METRIC(:,D),A(:,D),B(:,D),C(:,D))
            FAC(:,D)=0._EB
            IF (BC(D)==1 .OR. BC(D)==2) FAC(1,D)=-2._EB*A(1,D)
            IF (BC(D)==3 .OR. BC(D)==4) FAC(1,D)=DS(D)*A(1,D)
            IF (BC(D)==1 .OR. BC(D)==4) FAC(2,D)=-2._EB*C(XYZ(D),D)
            IF (BC(D)==2 .OR. BC(D)==3) FAC(2,D)=-DS(D)*C(XYZ(D),D)
         ENDDO

         ALLOCATE(S%FORWARD2(N2,N2),S%INVERSE2(N2,N2),S%FORWARD3(N3,N3),S%INVERSE3(N3,N3))
         ALLOCATE(F(N2,N2),V(N2,N2),LAM2(N2),LAM3(N3))
         D=ORDER(2)
         CALL WEIGHTED_BASIS(N2,DS(D),BC(D),METRIC(:,D),A(:,D),B(:,D),C(:,D),F,V,LAM2,S%ZERO2)
         S%FORWARD2=F; S%INVERSE2=V
         DEALLOCATE(F,V)
         ALLOCATE(F(N3,N3),V(N3,N3))
         D=ORDER(3)
         CALL WEIGHTED_BASIS(N3,DS(D),BC(D),METRIC(:,D),A(:,D),B(:,D),C(:,D),F,V,LAM3,S%ZERO3)
         S%FORWARD3=F; S%INVERSE3=V

         ALLOCATE(S%A(N1),S%C(N1),S%PIVOT(N1,N2,N3),S%RHS(N1,N2,N3),S%WORK(N1,N2,N3))
         D=ORDER(1)
         S%A=A(1:N1,D); S%C=C(1:N1,D)
         ! Cache inverse Thomas pivots once. Pin the last cell of the single null mode.
         ALLOCATE(B_D(N1),LAM2_D(N2),LAM3_D(N3))
         B_D=B(1:N1,D); LAM2_D=LAM2; LAM3_D=LAM3
         CALL FACTOR_LINES<<<(N2*N3+31)/32,32>>>(S%A,B_D,S%C,LAM2_D,LAM3_D,S%PIVOT,N1,N2,N3, &
              S%SINGULAR,S%ZERO2,S%ZERO3)
         CALL GPU_WAIT(0_CUDA_STREAM_KIND)
         ALLOCATE(S%WEIGHT1(N1),S%WEIGHT2(N2),S%WEIGHT3(N3),S%BC_FACTOR(2,3))
         S%WEIGHT1=METRIC(1:N1,ORDER(1)); S%WEIGHT2=METRIC(1:N2,ORDER(2)); S%WEIGHT3=METRIC(1:N3,ORDER(3))
         S%TOTAL_WEIGHT=SUM(METRIC(1:N1,ORDER(1)))*SUM(METRIC(1:N2,ORDER(2)))*SUM(METRIC(1:N3,ORDER(3)))
         S%BC_FACTOR=FAC
         S%SUM_BLOCKS=MIN(1024,(N1*N2*N3+REDUCE_BLOCK-1)/REDUCE_BLOCK)
         ALLOCATE(S%PARTIAL(S%SUM_BLOCKS),S%MEAN(1))
      END SUBROUTINE INITIALIZE_STRETCHED_POISSON

      SUBROUTINE AXIS_COEFFICIENTS(N,DS,BC,H,A,B,C)
         INTEGER, INTENT(IN) :: N,BC
         REAL(EB), INTENT(IN) :: DS,H(0:)
         REAL(EB), INTENT(OUT) :: A(:),B(:),C(:)
         INTEGER :: I
         DO I=1,N
            A(I)=2._EB/(DS**2*H(I)*(H(I-1)+H(I)))
            C(I)=2._EB/(DS**2*H(I)*(H(I)+H(I+1)))
            B(I)=-A(I)-C(I)
         ENDDO
         IF (BC==1 .OR. BC==2) B(1)=B(1)-A(1)
         IF (BC==3 .OR. BC==4) B(1)=B(1)+A(1)
         IF (BC==1 .OR. BC==4) B(N)=B(N)-C(N)
         IF (BC==2 .OR. BC==3) B(N)=B(N)+C(N)
      END SUBROUTINE AXIS_COEFFICIENTS

      SUBROUTINE WEIGHTED_BASIS(N,DS,BC,H,A,B,C,F,V,LAM,ZERO)
         INTEGER, INTENT(IN) :: N,BC
         REAL(EB), INTENT(IN) :: DS,H(0:),A(:),B(:),C(:)
         REAL(EB), INTENT(OUT) :: F(N,N),V(N,N),LAM(N)
         INTEGER, INTENT(OUT) :: ZERO
         ! Setup in double precision also serves GFDS_SINGLE_PRECISION builds.
         DOUBLE PRECISION :: Q(N,N),E(MAX(1,N-1)),D(N),WORK(MAX(1,2*N-2)),THETA,PI,NORM
         INTEGER :: I,P,INFO
         ZERO=0
         IF (BC==0) THEN
            ! CPU Fishpak admits periodicity on uniform transform axes.
            PI=ACOS(-1.D0)
            DO P=1,N
               IF (P==1) THEN
                  THETA=0.D0; Q(:,P)=1.D0/SQRT(DBLE(N))
               ELSEIF (MOD(N,2)==0 .AND. P==N) THEN
                  THETA=PI
                  DO I=1,N
                     Q(I,P)=COS(THETA*(I-1))/SQRT(DBLE(N))
                  ENDDO
               ELSE
                  THETA=2.D0*PI*(P/2)/N
                  NORM=SQRT(2.D0/N)
                  DO I=1,N
                     IF (MOD(P,2)==0) Q(I,P)=NORM*SIN(THETA*(I-1))
                     IF (MOD(P,2)==1) Q(I,P)=NORM*COS(THETA*(I-1))
                  ENDDO
               ENDIF
               D(P)=-4.D0*SIN(THETA/2.D0)**2/DBLE(DS)**2
            ENDDO
            ZERO=1
         ELSE
            ! sqrt(H) A / sqrt(H) is symmetric tridiagonal; DSTEV returns orthonormal Q.
            D=REAL(B(1:N),KIND(1.D0))
            DO I=1,N-1
               E(I)=SQRT(DBLE(C(I))*DBLE(A(I+1)))
            ENDDO
            CALL DSTEV('V',N,D,E,Q,N,WORK,INFO)
            IF (INFO/=0) ERROR STOP 'Stretched GPU Poisson: eigensystem setup failed'
            IF (BC==3) THEN
               ZERO=N; D(N)=0.D0
               Q(:,N)=SQRT(REAL(H(1:N),KIND(1.D0))/SUM(REAL(H(1:N),KIND(1.D0))))
            ENDIF
         ENDIF
         LAM=REAL(D,EB)
         DO P=1,N
            DO I=1,N
               V(I,P)=REAL(Q(I,P),EB)/SQRT(H(I))
               F(I,P)=REAL(Q(I,P),EB)*SQRT(H(I))
            ENDDO
         ENDDO
      END SUBROUTINE WEIGHTED_BASIS

      SUBROUTINE SOLVE_STRETCHED_POISSON(S,H,PRHS,BXS,BXF,BYS,BYF,BZS,BZF,STREAM)
         TYPE(GPU_STRETCH_POIS_TYPE), INTENT(INOUT) :: S
         REAL(EB), DEVICE, INTENT(INOUT) :: H(0:,0:,0:)
         REAL(EB), DEVICE, INTENT(IN) :: PRHS(:,:,:),BXS(:,:),BXF(:,:),BYS(:,:),BYF(:,:),BZS(:,:),BZF(:,:)
         INTEGER(KIND=CUDA_STREAM_KIND), INTENT(IN) :: STREAM
         TYPE(DIM3) :: BLOCK,GRID
         INTEGER :: N1,N2,N3,NX,NY,NZ,P1,P2,P3
         N1=S%N(1); N2=S%N(2); N3=S%N(3)
         NX=S%XYZ(1); NY=S%XYZ(2); NZ=S%XYZ(3)
         P1=S%AXIS(1); P2=S%AXIS(2); P3=S%AXIS(3)
         BLOCK=DIM3(8,4,4); GRID=DIM3((N1+7)/8,(N2+3)/4,(N3+3)/4)
         CALL PACK_RHS<<<GRID,BLOCK,0,STREAM>>>(PRHS,S%RHS,BXS,BXF,BYS,BYF,BZS,BZF,S%BC_FACTOR, &
              NX,NY,NZ,N1,N2,N3,P1,P2,P3)
         IF (S%SINGULAR) CALL REMOVE_WEIGHTED_MEAN(S,STREAM)
         CALL TRANSFORM_AXIS3<<<GRID,BLOCK,0,STREAM>>>(S%RHS,S%WORK,S%FORWARD3,N1,N2,N3,.TRUE.)
         CALL TRANSFORM_AXIS2<<<GRID,BLOCK,0,STREAM>>>(S%WORK,S%RHS,S%FORWARD2,N1,N2,N3,.TRUE.)
         CALL SOLVE_LINES<<<(N2*N3+31)/32,32,0,STREAM>>>(S%RHS,S%A,S%C,S%PIVOT,N1,N2,N3)
         CALL TRANSFORM_AXIS2<<<GRID,BLOCK,0,STREAM>>>(S%RHS,S%WORK,S%INVERSE2,N1,N2,N3,.FALSE.)
         CALL TRANSFORM_AXIS3<<<GRID,BLOCK,0,STREAM>>>(S%WORK,S%RHS,S%INVERSE3,N1,N2,N3,.FALSE.)
         IF (S%SINGULAR) CALL REMOVE_WEIGHTED_MEAN(S,STREAM)
         CALL UNPACK_H<<<GRID,BLOCK,0,STREAM>>>(S%RHS,H,NX,NY,NZ,N1,N2,N3,P1,P2,P3)
      END SUBROUTINE SOLVE_STRETCHED_POISSON

      SUBROUTINE REMOVE_WEIGHTED_MEAN(S,STREAM)
         TYPE(GPU_STRETCH_POIS_TYPE), INTENT(INOUT) :: S
         INTEGER(KIND=CUDA_STREAM_KIND), INTENT(IN) :: STREAM
         CALL WEIGHTED_SUM<<<S%SUM_BLOCKS,REDUCE_BLOCK,0,STREAM>>>(S%RHS,S%WEIGHT1,S%WEIGHT2,S%WEIGHT3, &
              S%PARTIAL,S%N(1),S%N(2),S%N(3))
         CALL SUM_MEAN<<<1,REDUCE_BLOCK,0,STREAM>>>(S%PARTIAL,S%MEAN,S%SUM_BLOCKS,S%TOTAL_WEIGHT)
         CALL SUBTRACT_MEAN<<<S%SUM_BLOCKS,REDUCE_BLOCK,0,STREAM>>>(S%RHS,S%MEAN,PRODUCT(S%N))
      END SUBROUTINE REMOVE_WEIGHTED_MEAN

      ATTRIBUTES(GLOBAL) SUBROUTINE PACK_RHS(PRHS,RHS,BXS,BXF,BYS,BYF,BZS,BZF,FAC,NX,NY,NZ,N1,N2,N3,P1,P2,P3)
         INTEGER, VALUE :: NX,NY,NZ,N1,N2,N3,P1,P2,P3
         REAL(EB), DEVICE :: PRHS(NX+1,NY+1,NZ+1),RHS(N1,N2,N3),FAC(2,3)
         REAL(EB), DEVICE :: BXS(MERGE(NY+1,1,NY>1),*),BXF(MERGE(NY+1,1,NY>1),*)
         REAL(EB), DEVICE :: BYS(NX+1,*),BYF(NX+1,*),BZS(NX+1,*),BZF(NX+1,*)
         INTEGER :: L,M,N,XYZ(3),I,J,K
         REAL(EB) :: V
         L=(BLOCKIDX%X-1)*BLOCKDIM%X+THREADIDX%X
         M=(BLOCKIDX%Y-1)*BLOCKDIM%Y+THREADIDX%Y
         N=(BLOCKIDX%Z-1)*BLOCKDIM%Z+THREADIDX%Z
         IF (L>N1 .OR. M>N2 .OR. N>N3) RETURN
         XYZ(P1)=L; XYZ(P2)=M; XYZ(P3)=N
         I=XYZ(1); J=XYZ(2); K=XYZ(3)
         V=PRHS(I,J,K)
         IF (I==1)  V=V+FAC(1,1)*BXS(J,K)
         IF (I==NX) V=V+FAC(2,1)*BXF(J,K)
         IF (J==1)  V=V+FAC(1,2)*BYS(I,K)
         IF (J==NY) V=V+FAC(2,2)*BYF(I,K)
         IF (K==1)  V=V+FAC(1,3)*BZS(I,J)
         IF (K==NZ) V=V+FAC(2,3)*BZF(I,J)
         RHS(L,M,N)=V
      END SUBROUTINE PACK_RHS

      ATTRIBUTES(GLOBAL) SUBROUTINE TRANSFORM_AXIS3(F,G,U,N1,N2,N3,FORWARD)
         INTEGER, VALUE :: N1,N2,N3
         LOGICAL, VALUE :: FORWARD
         REAL(EB), DEVICE :: F(N1,N2,N3),G(N1,N2,N3),U(N3,N3)
         INTEGER :: I,J,K,P
         REAL(EB) :: V
         I=(BLOCKIDX%X-1)*BLOCKDIM%X+THREADIDX%X
         J=(BLOCKIDX%Y-1)*BLOCKDIM%Y+THREADIDX%Y
         K=(BLOCKIDX%Z-1)*BLOCKDIM%Z+THREADIDX%Z
         IF (I>N1 .OR. J>N2 .OR. K>N3) RETURN
         V=0._EB
         DO P=1,N3
            IF (FORWARD) THEN
               V=V+U(P,K)*F(I,J,P)
            ELSE
               V=V+U(K,P)*F(I,J,P)
            ENDIF
         ENDDO
         G(I,J,K)=V
      END SUBROUTINE TRANSFORM_AXIS3

      ATTRIBUTES(GLOBAL) SUBROUTINE TRANSFORM_AXIS2(F,G,U,N1,N2,N3,FORWARD)
         INTEGER, VALUE :: N1,N2,N3
         LOGICAL, VALUE :: FORWARD
         REAL(EB), DEVICE :: F(N1,N2,N3),G(N1,N2,N3),U(N2,N2)
         INTEGER :: I,J,K,P
         REAL(EB) :: V
         I=(BLOCKIDX%X-1)*BLOCKDIM%X+THREADIDX%X
         J=(BLOCKIDX%Y-1)*BLOCKDIM%Y+THREADIDX%Y
         K=(BLOCKIDX%Z-1)*BLOCKDIM%Z+THREADIDX%Z
         IF (I>N1 .OR. J>N2 .OR. K>N3) RETURN
         V=0._EB
         DO P=1,N2
            IF (FORWARD) THEN
               V=V+U(P,J)*F(I,P,K)
            ELSE
               V=V+U(J,P)*F(I,P,K)
            ENDIF
         ENDDO
         G(I,J,K)=V
      END SUBROUTINE TRANSFORM_AXIS2

      ATTRIBUTES(GLOBAL) SUBROUTINE FACTOR_LINES(A,B,C,LAM2,LAM3,PIV,N1,N2,N3,SINGULAR,ZERO2,ZERO3)
         INTEGER, VALUE :: N1,N2,N3,ZERO2,ZERO3
         LOGICAL, VALUE :: SINGULAR
         REAL(EB), DEVICE :: A(N1),B(N1),C(N1),LAM2(N2),LAM3(N3),PIV(N1,N2,N3)
         INTEGER :: LINE,I,J,K,LAST
         REAL(EB) :: SHIFT
         LINE=(BLOCKIDX%X-1)*BLOCKDIM%X+THREADIDX%X
         IF (LINE>N2*N3) RETURN
         J=MOD(LINE-1,N2)+1; K=(LINE-1)/N2+1
         LAST=N1
         IF (SINGULAR .AND. J==ZERO2 .AND. K==ZERO3) LAST=N1-1
         SHIFT=LAM2(J)+LAM3(K)
         PIV(1,J,K)=1._EB/(B(1)+SHIFT)
         DO I=2,LAST
            PIV(I,J,K)=1._EB/(B(I)+SHIFT-A(I)*C(I-1)*PIV(I-1,J,K))
         ENDDO
         IF (LAST<N1) PIV(N1,J,K)=0._EB
      END SUBROUTINE FACTOR_LINES

      ATTRIBUTES(GLOBAL) SUBROUTINE SOLVE_LINES(F,A,C,PIV,N1,N2,N3)
         INTEGER, VALUE :: N1,N2,N3
         REAL(EB), DEVICE :: F(N1,N2,N3),A(N1),C(N1),PIV(N1,N2,N3)
         INTEGER :: LINE,I,J,K
         LINE=(BLOCKIDX%X-1)*BLOCKDIM%X+THREADIDX%X
         IF (LINE>N2*N3) RETURN
         J=MOD(LINE-1,N2)+1; K=(LINE-1)/N2+1
         F(1,J,K)=F(1,J,K)*PIV(1,J,K)
         DO I=2,N1
            F(I,J,K)=(F(I,J,K)-A(I)*F(I-1,J,K))*PIV(I,J,K)
         ENDDO
         DO I=N1-1,1,-1
            F(I,J,K)=F(I,J,K)-C(I)*PIV(I,J,K)*F(I+1,J,K)
         ENDDO
      END SUBROUTINE SOLVE_LINES

      ATTRIBUTES(GLOBAL) SUBROUTINE WEIGHTED_SUM(F,W1,W2,W3,PARTIAL,N1,N2,N3)
         INTEGER, VALUE :: N1,N2,N3
         REAL(EB), DEVICE :: F(N1,N2,N3),W1(N1),W2(N2),W3(N3),PARTIAL(*)
         REAL(EB), SHARED :: BUF(REDUCE_BLOCK)
         INTEGER :: T,L,I,J,K,S
         REAL(EB) :: V
         T=THREADIDX%X; V=0._EB
         DO L=(BLOCKIDX%X-1)*BLOCKDIM%X+T,N1*N2*N3,BLOCKDIM%X*GRIDDIM%X
            I=MOD(L-1,N1)+1; J=MOD((L-1)/N1,N2)+1; K=(L-1)/(N1*N2)+1
            V=V+F(I,J,K)*W1(I)*W2(J)*W3(K)
         ENDDO
         BUF(T)=V
         CALL SYNCTHREADS()
         S=REDUCE_BLOCK/2
         DO WHILE(S>0)
            IF (T<=S) BUF(T)=BUF(T)+BUF(T+S)
            CALL SYNCTHREADS()
            S=S/2
         ENDDO
         IF (T==1) PARTIAL(BLOCKIDX%X)=BUF(1)
      END SUBROUTINE WEIGHTED_SUM

      ATTRIBUTES(GLOBAL) SUBROUTINE SUM_MEAN(PARTIAL,MEAN,N,TOTAL)
         INTEGER, VALUE :: N
         REAL(EB), VALUE :: TOTAL
         REAL(EB), DEVICE :: PARTIAL(N),MEAN(1)
         REAL(EB), SHARED :: BUF(REDUCE_BLOCK)
         REAL(EB) :: V
         INTEGER :: T,I,S
         T=THREADIDX%X; V=0._EB
         DO I=T,N,REDUCE_BLOCK
            V=V+PARTIAL(I)
         ENDDO
         BUF(T)=V
         CALL SYNCTHREADS()
         S=REDUCE_BLOCK/2
         DO WHILE(S>0)
            IF (T<=S) BUF(T)=BUF(T)+BUF(T+S)
            CALL SYNCTHREADS()
            S=S/2
         ENDDO
         IF (T==1) MEAN(1)=BUF(1)/TOTAL
      END SUBROUTINE SUM_MEAN

      ATTRIBUTES(GLOBAL) SUBROUTINE SUBTRACT_MEAN(F,MEAN,N)
         INTEGER, VALUE :: N
         REAL(EB), DEVICE :: F(N),MEAN(1)
         INTEGER :: I
         DO I=(BLOCKIDX%X-1)*BLOCKDIM%X+THREADIDX%X,N,BLOCKDIM%X*GRIDDIM%X
            F(I)=F(I)-MEAN(1)
         ENDDO
      END SUBROUTINE SUBTRACT_MEAN

      ATTRIBUTES(GLOBAL) SUBROUTINE UNPACK_H(F,H,NX,NY,NZ,N1,N2,N3,P1,P2,P3)
         INTEGER, VALUE :: NX,NY,NZ,N1,N2,N3,P1,P2,P3
         REAL(EB), DEVICE :: F(N1,N2,N3),H(0:NX+1,0:NY+1,0:NZ+1)
         INTEGER :: I,J,K,XYZ(3)
         I=(BLOCKIDX%X-1)*BLOCKDIM%X+THREADIDX%X
         J=(BLOCKIDX%Y-1)*BLOCKDIM%Y+THREADIDX%Y
         K=(BLOCKIDX%Z-1)*BLOCKDIM%Z+THREADIDX%Z
         IF (I>N1 .OR. J>N2 .OR. K>N3) RETURN
         XYZ(P1)=I; XYZ(P2)=J; XYZ(P3)=K
         H(XYZ(1),XYZ(2),XYZ(3))=F(I,J,K)
      END SUBROUTINE UNPACK_H

END MODULE G_POIS
