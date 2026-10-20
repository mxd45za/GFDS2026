!------------------------------------------------------------------------------
! Module: G_DIVG
! Purpose: GPU divergence-source, enthalpy-transport, and density-derivative assembly.
!------------------------------------------------------------------------------
MODULE G_DIVG
    USE G_ERRO, ONLY: GPU_WAIT
    USE G_MESH_POINTERS
    USE G_CONS
    USE G_FUNC
    USE CUDAFOR
    CONTAINS

        !--------------------------------------------------------------------------
        ! Scalar-face interpolation kernels
        !--------------------------------------------------------------------------

        ATTRIBUTES(GLOBAL) SUBROUTINE DIVG_GET_SCALAR_FACE_VALUE_KERNEL_XYZ(AX,AY,AZ,U,FX,FY,FZ, &
            IBAR,JBAR,KBAR,N,TWO_EPSILON_EB)
            INTEGER, VALUE :: IBAR, JBAR, KBAR, N
            REAL(EB), VALUE :: TWO_EPSILON_EB
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: AX,AY,AZ,U
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:KBAR+1,1:*) :: FX,FY,FZ
            INTEGER :: I,J,K,S,NEG,IUP,IDN,IUP2,JUP,JDN,JUP2,KUP,KDN,KUP2
            REAL(EB) :: R,B,DU_UP,DU_LOC,U_UP,U_DN,U_UP2

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I>=1 .AND. I<=IBAR-1 .AND. J>=1 .AND. J<=JBAR .AND. K>=1 .AND. K<=KBAR) THEN
                NEG = MERGE(0,1,AX(I,J,K)>0._EB)
                S = 1-2*NEG
                IUP = I+NEG
                IDN = I+1-NEG
                IUP2 = I-1+3*NEG
                U_UP = U(IUP,J,K)
                U_DN = U(IDN,J,K)
                U_UP2 = U(IUP2,J,K)
                DU_LOC = S*(U_DN-U_UP)
                DU_UP = S*(U_UP-U_UP2)
                R = DU_UP/(DU_LOC+SIGN(TWO_EPSILON_EB,DU_LOC))
                B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
                FX(I,J,K,N) = U_UP+0.5_EB*S*B*DU_LOC
            ENDIF

            IF (I>=1 .AND. I<=IBAR .AND. J>=1 .AND. J<=JBAR-1 .AND. K>=1 .AND. K<=KBAR) THEN
                NEG = MERGE(0,1,AY(I,J,K)>0._EB)
                S = 1-2*NEG
                JUP = J+NEG
                JDN = J+1-NEG
                JUP2 = J-1+3*NEG
                U_UP = U(I,JUP,K)
                U_DN = U(I,JDN,K)
                U_UP2 = U(I,JUP2,K)
                DU_LOC = S*(U_DN-U_UP)
                DU_UP = S*(U_UP-U_UP2)
                R = DU_UP/(DU_LOC+SIGN(TWO_EPSILON_EB,DU_LOC))
                B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
                FY(I,J,K,N) = U_UP+0.5_EB*S*B*DU_LOC
            ENDIF

            IF (I>=1 .AND. I<=IBAR .AND. J>=1 .AND. J<=JBAR .AND. K>=1 .AND. K<=KBAR-1) THEN
                NEG = MERGE(0,1,AZ(I,J,K)>0._EB)
                S = 1-2*NEG
                KUP = K+NEG
                KDN = K+1-NEG
                KUP2 = K-1+3*NEG
                U_UP = U(I,J,KUP)
                U_DN = U(I,J,KDN)
                U_UP2 = U(I,J,KUP2)
                DU_LOC = S*(U_DN-U_UP)
                DU_UP = S*(U_UP-U_UP2)
                R = DU_UP/(DU_LOC+SIGN(TWO_EPSILON_EB,DU_LOC))
                B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
                FZ(I,J,K,N) = U_UP+0.5_EB*S*B*DU_LOC
            ENDIF
        END SUBROUTINE DIVG_GET_SCALAR_FACE_VALUE_KERNEL_XYZ

        ATTRIBUTES(GLOBAL) SUBROUTINE DIVG_GET_SCALAR_FACE_VALUE_KERNEL_BOUNDARY(WALL_TOPO,CELL_INDEXD,CELL_WALL_INDEX, &
            CELL_SOLID,CELL_EXTERIOR,WALL_RHO_F,RHOD,UD,VD,WD,FXD,FYD,FZD, &
            N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,N,TWO_EPSILON_EB,IBP1,JBP1,KBP1,CELL_LB,CELL_UB)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,N,IBP1,JBP1,KBP1,CELL_LB,CELL_UB
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: RHOD,UD,VD,WD
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:KBP1,1:*) :: FXD,FYD,FZD
            REAL(EB), DEVICE, DIMENSION(0:*) :: WALL_RHO_F
            INTEGER, DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: CELL_INDEXD
            INTEGER, DEVICE, DIMENSION(CELL_LB:CELL_UB,-3:*) :: CELL_WALL_INDEX
            LOGICAL, DEVICE, DIMENSION(0:*) :: CELL_SOLID,CELL_EXTERIOR
            REAL(EB), DEVICE :: TWO_EPSILON_EB
            REAL(EB) :: R,B,DU_UP,DU_LOC
            INTEGER :: IW,II,JJ,KK,IOR,S,VEL_INDEX,ADJ_INDEX
            REAL(EB) :: Q0,Q1,Q2

            IW = (BLOCKIDX%X - 1)*BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS) RETURN
            IF (WALL_TOPO(IW)%BOUNDARY_TYPE == 0) RETURN

            II = WALL_TOPO(IW)%II
            JJ = WALL_TOPO(IW)%JJ
            KK = WALL_TOPO(IW)%KK
            IOR = WALL_TOPO(IW)%IOR

            S = MERGE(1,-1,IOR>0)

            SELECT CASE(ABS(IOR))
                CASE(1)
                    VEL_INDEX = II + MERGE(1,-2,S>0)
                    ADJ_INDEX = II + MERGE(1,-1,S>0)
                    IF (WALL_TOPO(IW)%BOUNDARY_TYPE /= INTERPOLATED_BOUNDARY .AND. &
                        WALL_TOPO(IW)%BOUNDARY_TYPE /= OPEN_BOUNDARY .AND. S*UD(VEL_INDEX,JJ,KK)>0._EB .AND. &
                        .NOT.(CELL_WALL_INDEX(CELL_INDEXD(ADJ_INDEX,JJ,KK),IOR)>0)) THEN
                        IF (S > 0) THEN
                            Q0=RHOD(II+1,JJ,KK); Q1=RHOD(II+1,JJ,KK); Q2=RHOD(II+2,JJ,KK)
                            DU_LOC=Q2-Q1; DU_UP=Q1-Q0
                            R=DU_UP/(DU_LOC+SIGN(TWO_EPSILON_EB,DU_LOC)); B=MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
                            FXD(VEL_INDEX,JJ,KK,N)=Q1+0.5_EB*B*DU_LOC
                        ELSE
                            Q0=RHOD(II-2,JJ,KK); Q1=RHOD(II-1,JJ,KK); Q2=RHOD(II-1,JJ,KK)
                            DU_LOC=Q1-Q0; DU_UP=Q2-Q1
                            R=DU_UP/(DU_LOC+SIGN(TWO_EPSILON_EB,DU_LOC)); B=MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
                            FXD(VEL_INDEX,JJ,KK,N)=Q1-0.5_EB*B*DU_LOC
                        ENDIF
                    ENDIF
                CASE(2)
                    VEL_INDEX = JJ + MERGE(1,-2,S>0)
                    ADJ_INDEX = JJ + MERGE(1,-1,S>0)
                    IF (WALL_TOPO(IW)%BOUNDARY_TYPE /= INTERPOLATED_BOUNDARY .AND. &
                        WALL_TOPO(IW)%BOUNDARY_TYPE /= OPEN_BOUNDARY .AND. S*VD(II,VEL_INDEX,KK)>0._EB .AND. &
                        .NOT.(CELL_WALL_INDEX(CELL_INDEXD(II,ADJ_INDEX,KK),IOR)>0)) THEN
                        IF (S > 0) THEN
                            Q0=RHOD(II,JJ+1,KK); Q1=RHOD(II,JJ+1,KK); Q2=RHOD(II,JJ+2,KK)
                            DU_LOC=Q2-Q1; DU_UP=Q1-Q0
                            R=DU_UP/(DU_LOC+SIGN(TWO_EPSILON_EB,DU_LOC)); B=MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
                            FYD(II,VEL_INDEX,KK,N)=Q1+0.5_EB*B*DU_LOC
                        ELSE
                            Q0=RHOD(II,JJ-2,KK); Q1=RHOD(II,JJ-1,KK); Q2=RHOD(II,JJ-1,KK)
                            DU_LOC=Q1-Q0; DU_UP=Q2-Q1
                            R=DU_UP/(DU_LOC+SIGN(TWO_EPSILON_EB,DU_LOC)); B=MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
                            FYD(II,VEL_INDEX,KK,N)=Q1-0.5_EB*B*DU_LOC
                        ENDIF
                    ENDIF
                CASE(3)
                    VEL_INDEX = KK + MERGE(1,-2,S>0)
                    ADJ_INDEX = KK + MERGE(1,-1,S>0)
                    IF (WALL_TOPO(IW)%BOUNDARY_TYPE /= INTERPOLATED_BOUNDARY .AND. &
                        WALL_TOPO(IW)%BOUNDARY_TYPE /= OPEN_BOUNDARY .AND. S*WD(II,JJ,VEL_INDEX)>0._EB .AND. &
                        .NOT.(CELL_WALL_INDEX(CELL_INDEXD(II,JJ,ADJ_INDEX),IOR)>0)) THEN
                        IF (S > 0) THEN
                            Q0=RHOD(II,JJ,KK+1); Q1=RHOD(II,JJ,KK+1); Q2=RHOD(II,JJ,KK+2)
                            DU_LOC=Q2-Q1; DU_UP=Q1-Q0
                            R=DU_UP/(DU_LOC+SIGN(TWO_EPSILON_EB,DU_LOC)); B=MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
                            FZD(II,JJ,VEL_INDEX,N)=Q1+0.5_EB*B*DU_LOC
                        ELSE
                            Q0=RHOD(II,JJ,KK-2); Q1=RHOD(II,JJ,KK-1); Q2=RHOD(II,JJ,KK-1)
                            DU_LOC=Q1-Q0; DU_UP=Q2-Q1
                            R=DU_UP/(DU_LOC+SIGN(TWO_EPSILON_EB,DU_LOC)); B=MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
                            FZD(II,JJ,VEL_INDEX,N)=Q1-0.5_EB*B*DU_LOC
                        ENDIF
                    ENDIF
            END SELECT
        END SUBROUTINE DIVG_GET_SCALAR_FACE_VALUE_KERNEL_BOUNDARY


        !--------------------------------------------------------------------------
        ! Divergence part 1 host orchestration
        !--------------------------------------------------------------------------
        SUBROUTINE DIVERGENCE_PART_1_GPU(T,DT,NM)
            USE PHYSICAL_FUNCTIONS, ONLY: GET_CONDUCTIVITY,GET_SPECIFIC_HEAT,GET_SENSIBLE_ENTHALPY_Z,GET_SENSIBLE_ENTHALPY,&
                                            GET_VISCOSITY,GET_MOLECULAR_WEIGHT
            USE MANUFACTURED_SOLUTIONS, ONLY: DIFF_MMS,UF_MMS,WF_MMS,VD2D_MMS_Z_SRC
            USE COMPLEX_GEOMETRY, ONLY : CC_CGSC, CC_UNKZ, CC_SOLID, CC_CUTCFE
            USE CC_SCALARS, ONLY : ADD_CUTCELL_PSUM,ADD_LINKEDCELL_PSUM,SET_EXIMDIFFLX_3D,SET_EXIMRHOHSLIM_3D,&
                                    SET_EXIMRHOZZLIM_3D,CC_DIVERGENCE_PART_1,CC_VELOCITY_FLUX,CFACE_PREDICT_NORMAL_VELOCITY
            USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
            USE MATH_FUNCTIONS, ONLY: INTERPOLATE1D_UNIFORM

            USE TURBULENCE, ONLY: TENSOR_DIFFUSIVITY_MODEL
            IMPLICIT NONE

            REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW,RHOP
            REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZP

            REAL(EB) :: SM_RCON_1
            INTEGER, INTENT(IN) :: NM
            REAL(EB), INTENT(IN) :: T,DT
            REAL(EB), POINTER, DIMENSION(:,:) :: PBAR_P
            REAL(EB), POINTER, DIMENSION(:,:,:) :: KDTDX,KDTDY,KDTDZ,DP,KP,CP, &
                        RHO_D,H_RHO_D_DZDX,H_RHO_D_DZDY,H_RHO_D_DZDZ,RTRM, &
                        U_DOT_DEL_RHO_H_S,U_DOT_DEL_RHO_Z,RHO_D_TURB,R_H_G
            REAL(EB), POINTER, DIMENSION(:,:,:,:) :: RHO_D_DZDX,RHO_D_DZDY,RHO_D_DZDZ
            REAL(EB) :: DELKDELT,VC,VC1,DTDX,DTDY,DTDZ,TNOW,DZDX,DZDY,DZDZ,TMP_G,DIV_DIFF_HEAT_FLUX,H_S,XHAT,ZHAT,TT,Q_Z,&
                        D_Z_TEMP,D_Z_N(0:I_MAX_TEMP),RHO_D_DZDN_GET(1:N_TRACKED_SPECIES),UN_P,TMP_F_GAS,R_PFCT,RHO_D_DZDN
            INTEGER :: IW,N,I,J,K,IPZ,N_ZZ_MAX,ICF
            REAL(EB), ALLOCATABLE, DIMENSION(:) :: WALL_U_NORMAL_H_TMP, WALL_U_NORMAL_S_H_TMP
            REAL(EB), ALLOCATABLE, DIMENSION(:) :: WALL_U_NORMAL_0_H_TMP, WALL_AREA_H_TMP, WALL_AREA_ADJUST_H_TMP
            REAL(EB), ALLOCATABLE, DIMENSION(:) :: WALL_K_G_H_TMP, WALL_DUNDT_H_TMP
            REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: WALL_ZZ_F_H_TMP, WALL_RHO_D_F_H_TMP, WALL_RHO_D_DZDN_F_H_TMP
            REAL(EB), ALLOCATABLE, DIMENSION(:) :: WALL_Q_CON_F_H_TMP, WALL_RHO_F_H_TMP, WALL_TMP_F_H_TMP
            REAL(EB), ALLOCATABLE, DIMENSION(:) :: DSUM_ZONE_H, PSUM_ZONE_H, USUM_ZONE_H
            INTEGER, ALLOCATABLE, DIMENSION(:) :: WALL_PRESSURE_ZONE_H_TMP
            TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1


            REAL(EB), DEVICE, POINTER, DIMENSION(:,:,:)   :: UUD,VVD,WWD,RHOPD
            REAL(EB), DEVICE, POINTER, DIMENSION(:,:,:,:) :: ZZPD
            REAL(EB), DEVICE, POINTER, DIMENSION(:,:,:)   :: KDTDXD,KDTDYD,KDTDZD,DPD,KPD,CPD,RTRMD,R_H_GD
            REAL(EB), DEVICE, POINTER, DIMENSION(:,:)     :: PBAR_PD
            REAL(EB), DEVICE, POINTER, DIMENSION(:,:,:) :: FX_H_SD,FY_H_SD,FZ_H_SD,RHO_H_S_PD,U_DOT_DEL_RHO_H_SD
            REAL(EB), DEVICE, POINTER, DIMENSION(:,:,:) :: RHO_Z_PD,U_DOT_DEL_RHO_ZD

            REAL(EB), DEVICE, POINTER, CONTIGUOUS :: WALL_K_G_D_TMP(:)
            INTEGER :: BLKD_1D,GRDD_1D
            INTEGER :: BLKD_SUM, GRDD_SUM, NCELL_SUM
            TYPE(DIM3) :: BLKD, GRDD
            INTEGER, PARAMETER :: BS = 8
            LOGICAL :: HAS_D_SOURCE, HAS_M_DOT_PPP

            IF (SOLID_PHASE_ONLY) RETURN

            TNOW=CURRENT_TIME()
            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)
            HAS_D_SOURCE  = ALLOCATED(MESHES(NM)%D_SOURCE)
            HAS_M_DOT_PPP = ALLOCATED(MESHES(NM)%M_DOT_PPP)

            BLKD_1D = GPU_BLKD_1D
            GRDD_1D = GM%GRDD_TOTAL_WALL_1D

            BLKD = GM%BLKD_3D_BP1
            GRDD = GM%GRDD_3D_BP1
            SELECT CASE(PREDICTOR)
            CASE(.TRUE.)
                DP => DS
                PBAR_P => PBAR_S
                RHOP => RHOS
                UU=>U
                VV=>V
                WW=>W
                ZZP=>ZZS

                DPD => GM%DS_D
                PBAR_PD => GM%PBAR_SD
                RHOPD => GM%RHOSD
                UUD => GM%UD
                VVD => GM%VD
                WWD => GM%WD
                ZZPD => GM%ZZSD

            CASE(.FALSE.)
                DP => D
                PBAR_P => PBAR
                RHOP => RHO
                UU=>US
                VV=>VS
                WW=>WS
                ZZP=>ZZ

                DPD => GM%D_D
                PBAR_PD => GM%PBARD
                RHOPD => GM%RHOD
                UUD => GM%USD
                VVD => GM%VSD
                WWD => GM%WSD
                ZZPD => GM%ZZD
            END SELECT

            DPD = 0._EB

            R_PBAR = 1._EB/PBAR_P
            SM_RCON_1 = SPECIES_MIXTURE(1)%RCON

            KDTDXD => GM%WORK1D
            KDTDYD => GM%WORK2D
            KDTDZD => GM%WORK3D
            KPD    => GM%WORK4D
            CPD    => GM%WORK5D
            R_H_GD => GM%WORK9D

            IF (.NOT.ALLOCATED(GM%WALL_K_G_WORKD)) ALLOCATE(GM%WALL_K_G_WORKD(LBOUND(WALL,1):UBOUND(WALL,1)))
            WALL_K_G_D_TMP => GM%WALL_K_G_WORKD
            WALL_K_G_D_TMP = 0._EB
            IF (CHECK_VN) GM%D_Z_MAXD = 0._EB

            IF (GPU_DEBUG_H2D_D2H) THEN
                GM%TMPD      = TMP
                GM%MUD       = MU
                PBAR_PD   = PBAR_P
                GM%UVW_SAVED = UVW_SAVE
                UUD       = UU
                VVD       = VV
                WWD       = WW
                RHOPD     = RHOP
                ZZPD      = ZZP
                GM%RSUMD     = RSUM
                GM%QD        = Q
                GM%QRD       = QR
                GM%D_SOURCED = 0._EB
                GM%M_DOT_PPPD = 0._EB
                IF (HAS_D_SOURCE) GM%D_SOURCED = D_SOURCE
                IF (HAS_M_DOT_PPP) GM%M_DOT_PPPD = M_DOT_PPP
                ALLOCATE(WALL_Q_CON_F_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)))
                ALLOCATE(WALL_RHO_F_H_TMP  (LBOUND(WALL,1):UBOUND(WALL,1)))
                ALLOCATE(WALL_TMP_F_H_TMP  (LBOUND(WALL,1):UBOUND(WALL,1)))
                ALLOCATE(WALL_U_NORMAL_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)))
                ALLOCATE(WALL_U_NORMAL_S_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)))
                ALLOCATE(WALL_U_NORMAL_0_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)))
                ALLOCATE(WALL_AREA_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)))
                ALLOCATE(WALL_AREA_ADJUST_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)))
                ALLOCATE(WALL_K_G_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)))
                ALLOCATE(WALL_DUNDT_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)))
                ALLOCATE(WALL_ZZ_F_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1),1:N_TOTAL_SCALARS))
                ALLOCATE(WALL_RHO_D_F_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1),1:N_TRACKED_SPECIES))
                ALLOCATE(WALL_RHO_D_DZDN_F_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1),1:N_TRACKED_SPECIES))
                ALLOCATE(WALL_PRESSURE_ZONE_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)))
                WALL_Q_CON_F_H_TMP      = 0._EB
                WALL_RHO_F_H_TMP        = 0._EB
                WALL_TMP_F_H_TMP        = 0._EB
                WALL_U_NORMAL_H_TMP     = 0._EB
                WALL_U_NORMAL_S_H_TMP   = 0._EB
                WALL_U_NORMAL_0_H_TMP   = 0._EB
                WALL_AREA_H_TMP         = 0._EB
                WALL_AREA_ADJUST_H_TMP  = 1._EB
                WALL_K_G_H_TMP          = 0._EB
                WALL_DUNDT_H_TMP        = 0._EB
                WALL_ZZ_F_H_TMP         = 0._EB
                WALL_RHO_D_F_H_TMP      = 0._EB
                WALL_RHO_D_DZDN_F_H_TMP = 0._EB
                IF (N_TOTAL_SCALARS >= 1) WALL_ZZ_F_H_TMP(:,1) = 1._EB
                WALL_PRESSURE_ZONE_H_TMP = 0
                DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
                    IF (WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE
                    B1 => BOUNDARY_PROP1(WALL(IW)%B1_INDEX)
                    WALL_Q_CON_F_H_TMP(IW)      = B1%Q_CON_F
                    WALL_RHO_F_H_TMP(IW)        = B1%RHO_F
                    WALL_TMP_F_H_TMP(IW)        = B1%TMP_F
                    WALL_U_NORMAL_H_TMP(IW)     = B1%U_NORMAL
                    WALL_U_NORMAL_S_H_TMP(IW)   = B1%U_NORMAL_S
                    WALL_U_NORMAL_0_H_TMP(IW)   = B1%U_NORMAL_0
                    WALL_AREA_H_TMP(IW)         = B1%AREA
                    WALL_AREA_ADJUST_H_TMP(IW)  = B1%AREA_ADJUST
                    WALL_K_G_H_TMP(IW)          = B1%K_G
                    IF (IW <= N_EXTERNAL_WALL_CELLS) WALL_DUNDT_H_TMP(IW) = EXTERNAL_WALL(IW)%DUNDT
                    IF (ALLOCATED(B1%ZZ_F)) WALL_ZZ_F_H_TMP(IW,1:N_TOTAL_SCALARS) = B1%ZZ_F(1:N_TOTAL_SCALARS)
                    IF (ALLOCATED(B1%RHO_D_F)) WALL_RHO_D_F_H_TMP(IW,1:N_TRACKED_SPECIES) = B1%RHO_D_F(1:N_TRACKED_SPECIES)
                    IF (ALLOCATED(B1%RHO_D_DZDN_F)) WALL_RHO_D_DZDN_F_H_TMP(IW,1:N_TRACKED_SPECIES) = B1%RHO_D_DZDN_F(1:N_TRACKED_SPECIES)
                    WALL_PRESSURE_ZONE_H_TMP(IW)= B1%PRESSURE_ZONE
                ENDDO
                GM%WALL_Q_CON_F      = WALL_Q_CON_F_H_TMP
                GM%WALL_RHO_F        = WALL_RHO_F_H_TMP
                GM%WALL_TMP_F        = WALL_TMP_F_H_TMP
                GM%WALL_U_NORMAL     = WALL_U_NORMAL_H_TMP
                GM%WALL_U_NORMAL_S   = WALL_U_NORMAL_S_H_TMP
                GM%WALL_U_NORMAL_0   = WALL_U_NORMAL_0_H_TMP
                GM%WALL_AREA         = WALL_AREA_H_TMP
                GM%WALL_AREA_ADJUST  = WALL_AREA_ADJUST_H_TMP
                GM%WALL_K_G          = WALL_K_G_H_TMP
                GM%WALL_DUNDTD       = WALL_DUNDT_H_TMP
                DO N=1,N_TOTAL_SCALARS
                    GM%WALL_ZZ_F(:,N) = WALL_ZZ_F_H_TMP(:,N)
                ENDDO
                DO N=1,N_TRACKED_SPECIES
                    GM%WALL_RHO_D_F(:,N)      = WALL_RHO_D_F_H_TMP(:,N)
                    GM%WALL_RHO_D_DZDN_F(:,N) = WALL_RHO_D_DZDN_F_H_TMP(:,N)
                ENDDO
                GM%WALL_PRESSURE_ZONE = WALL_PRESSURE_ZONE_H_TMP
            ELSE
                IF (.NOT. WALL_BC_ON) CALL SYNC_DIV1_WALL_INPUTS_TO_GPU(NM)
                IF (.NOT. COMPUTE_RADIATION_ON) GM%QRD = QR
                IF (.NOT. COMBUSTION_LOAD_BALANCED_ON) THEN
                    GM%QD = Q
                    GM%D_SOURCED = 0._EB
                    GM%M_DOT_PPPD = 0._EB
                    IF (HAS_D_SOURCE) GM%D_SOURCED = D_SOURCE
                    IF (HAS_M_DOT_PPP) GM%M_DOT_PPPD = M_DOT_PPP
                ENDIF
            ENDIF

            CALL PREDICT_NORMAL_VELOCITY_GPU<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>( &
                GM%WALL_U_NORMAL,GM%WALL_U_NORMAL_S,GM%WALL_U_NORMAL_0,GM%WALL_RHO_F,GM%WALL_T_IGN,GM%WALL_TOPO_AOS, &
                GM%WALL_SF_INDEX,SF_SPECIES_BC_INDEX,SF_MASS_FLUX_TOTAL,SF_RAMP_VELO_INDEX,SF_RAMP_VELO_TAU, &
                RAMP_T_MIN_D,RAMP_SPAN_D,RAMP_RDT_D,RAMP_NPOINTS_D,RAMP_OFFSET_D,RAMP_INTERP_D, &
                UUD,VVD,WWD,N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,MERGE(1,0,PREDICTOR),T,DT,T_BEGIN,RHOA,IBP1,JBP1)
            CALL UPDATE_WALL_DUNDT_GPU<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>( &
                GM%WALL_DUNDTD,GM%WALL_U_NORMAL,GM%WALL_U_NORMAL_S,N_EXTERNAL_WALL_CELLS,DT,MERGE(1,0,PREDICTOR))

            IF (N_TOTAL_SCALARS > 1) THEN
                CALL ZERO_DIFFUSIVE_WORK_GPU<<<GRDD,BLKD,0,GM%STRM>>>( &
                    GM%DEL_RHO_D_DEL_ZD,GM%SWORK1D,GM%SWORK2D,GM%SWORK3D,IBP1,JBP1,KBP1,N_TOTAL_SCALARS)

                CALL SPECIES_RHO_D_VLES_GPU<<<GRDD,BLKD,0,GM%STRM>>>( &
                    GM%MUD,GM%WORK4D,RSC,GM%D_Z_MAXD,RHOPD,IBP1,JBP1,KBP1,MERGE(1,0,CHECK_VN))

                DO N=1,N_TRACKED_SPECIES
                    CALL SPECIES_DIFFUSIVE_FLUX_GPU<<<GRDD,BLKD,0,GM%STRM>>>( &
                        ZZPD,GM%WORK4D,GM%SWORK1D,GM%SWORK2D,GM%SWORK3D,GM%RDXND,GM%RDYND,GM%RDZND,IBAR,JBAR,KBAR,N)

                    CALL SPECIES_DIFFUSIVE_FLUX_BOUNDARY_GPU<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>(GM%SWORK1D, GM%SWORK2D, GM%SWORK3D, GM%WALL_RHO_D_DZDN_F, &
                        GM%WALL_TOPO_AOS,GM%EXTERNAL_WALL_NICD,N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,N, &
                        IBP1,JBP1,KBP1,GM%WALL_LB,GM%WALL_UB)

                    CALL SPECIES_DIFF_HEAT_FACE_GPU<<<GRDD,BLKD,0,GM%STRM>>>( &
                        GM%TMPD,H_SENS_ZD,GM%SWORK1D,GM%SWORK2D,GM%SWORK3D,GM%WORK5D,GM%WORK6D,GM%WORK7D,IBAR,JBAR,KBAR,N)

                    CALL SPECIES_DIFF_HEAT_BOUNDARY_GPU<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>(GM%TMPD, H_SENS_ZD, ZZPD, GM%SWORK1D, GM%SWORK2D, GM%SWORK3D, GM%WORK5D, &
                        GM%WORK6D, GM%WORK7D, GM%WALL_ZZ_F, GM%WALL_RHO_D_F, GM%WALL_RHO_D_DZDN_F, GM%WALL_TMP_F, GM%WALL_U_NORMAL, &
                        GM%WALL_U_NORMAL_S, GM%WALL_RDN, GM%WALL_TOPO_AOS, N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS, &
                        N_TRACKED_SPECIES,N,MERGE(1,0,PREDICTOR),IBP1,JBP1,KBP1,GM%WALL_LB,GM%WALL_UB)

                    CALL SPECIES_DIFF_HEAT_DIV_GPU<<<GRDD,BLKD,0,GM%STRM>>>( &
                        DPD,GM%DEL_RHO_D_DEL_ZD,GM%SWORK1D,GM%SWORK2D,GM%SWORK3D,GM%WORK5D,GM%WORK6D,GM%WORK7D,GM%RDXD,GM%RDYD,GM%RDZD,IBAR,JBAR,KBAR,N)
                ENDDO
            ELSE
                GM%DEL_RHO_D_DEL_ZD = 0._EB
            ENDIF

            ! ---- CPU block: cp / r_h_g ----
            CALL GET_SPECIFIC_HEAT_GPU<<<GRDD,BLKD,0,GM%STRM>>>(ZZPD,CPD,CP_ZD,R_H_GD,GM%TMPD,N_TRACKED_SPECIES,IBAR,JBAR,KBAR)

            ! ---- CPU block: KP and D_Z_MAX thermal contribution ----
            CALL KDTDN_GPU_1<<<GRDD,BLKD,0,GM%STRM>>>(GM%MUD,CPOPR,KPD,IBAR,JBAR,KBAR)

            IF (CHECK_VN .AND. .NOT.CONSTANT_SPECIFIC_HEAT_RATIO) THEN
                CALL D_Z_MAX_THERMAL_GPU<<<GRDD,BLKD,0,GM%STRM>>>(GM%D_Z_MAXD,KPD,CPD,RHOPD,GM%CELL_INDEXD,GM%CELL_SOLID,IBAR,JBAR,KBAR)
            ENDIF

            ! ---- CPU block: k*dT/dx, k*dT/dy, k*dT/dz ----
            CALL KDTDN_GPU_2<<<GRDD,BLKD,0,GM%STRM>>>(GM%TMPD,KPD,KDTDXD,KDTDYD,KDTDZD,IBAR,JBAR,KBAR,GM%RDXND,GM%RDYND,GM%RDZND)

            ! ---- CPU block: boundary correction on thermal gradients and wall conductive source ----
            CALL KDTDN_GPU_BOUNDARY<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>(DPD, GM%WALL_Q_CON_F, GM%WALL_AREA_ADJUST, GM%WALL_TOPO_AOS, GM%WALL_RDN, KPD, &
                WALL_K_G_D_TMP,KDTDXD,KDTDYD,KDTDZD,N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,IBP1,JBP1)

            ! ---- CPU block: add del dot k del T ----
            CALL DP_ADD_DELKDELT_GPU<<<GRDD,BLKD,0,GM%STRM>>>(DPD,KDTDXD,KDTDYD,KDTDZD,GM%QD,GM%QRD,GM%RDXD,GM%RDYD,GM%RDZD,IBAR,JBAR,KBAR)

            ! ---- CPU block: ENTHALPY_ADVECTION ----
            U_DOT_DEL_RHO_H_SD => GM%WORK6D
            FX_H_SD            => GM%WORK1D
            FY_H_SD            => GM%WORK2D
            FZ_H_SD            => GM%WORK3D
            RHO_H_S_PD         => GM%WORK4D

            CALL ENTHALPY_ADVECTION_GPU_INIT<<<GRDD,BLKD,0,GM%STRM>>>( &
                H_SENS_ZD, U_DOT_DEL_RHO_H_SD, FX_H_SD, FY_H_SD, FZ_H_SD, RHO_H_S_PD, &
                RHOPD, ZZPD, GM%TMPD, N_TRACKED_SPECIES, IBP1, JBP1, KBP1)

            CALL GET_SCALAR_FACE_VALUE_KERNEL_XYZ<<<GRDD,BLKD,0,GM%STRM>>>( &
                UUD,VVD,WWD,RHO_H_S_PD,FX_H_SD,FY_H_SD,FZ_H_SD,IBAR,JBAR,KBAR,TWO_EPSILON_EB)

            CALL GET_SCALAR_FACE_VALUE_KERNEL_BOUNDARY_H<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>( &
                GM%WALL_TOPO_AOS,GM%CELL_INDEXD,GM%CELL_WALL_INDEX,GM%CELL_SOLID,GM%CELL_EXTERIOR,GM%WALL_RHO_F, &
                RHO_H_S_PD,UUD,VVD,WWD,FX_H_SD,FY_H_SD,FZ_H_SD,N_EXTERNAL_WALL_CELLS, &
                N_INTERNAL_WALL_CELLS,TWO_EPSILON_EBD,IBP1,JBP1,KBP1,GM%CELL_LB,GM%CELL_UB)

            CALL ENTHALPY_ADVECTION_BOUNDARY_GPU<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>(U_DOT_DEL_RHO_H_SD, RHO_H_S_PD, UUD, VVD, WWD, GM%TMPD, &
                H_SENS_ZD, GM%WALL_TMP_F, GM%WALL_RHO_F, GM%WALL_U_NORMAL, GM%WALL_U_NORMAL_S, GM%WALL_ZZ_F, GM%WALL_RDN, GM%UVW_SAVED, &
                GM%WALL_TOPO_AOS,N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,N_TRACKED_SPECIES,MERGE(1,0,PREDICTOR), &
                IBP1,JBP1,GM%WALL_LB,GM%WALL_UB)

            CALL ENTHALPY_ADVECTION_GPU_END<<<GRDD,BLKD,0,GM%STRM>>>(U_DOT_DEL_RHO_H_SD,FX_H_SD,FY_H_SD,FZ_H_SD,RHO_H_S_PD,UUD,VVD,WWD, &
                    GM%CELL_INDEXD,GM%CELL_SOLID,GM%CELL_WALL_INDEX,GM%RDXD,GM%RDYD,GM%RDZD,IBAR,JBAR,KBAR,GM%CELL_LB,GM%CELL_UB)

            CALL DP_SUB_ENTHALPY_ADV_GPU<<<GRDD,BLKD,0,GM%STRM>>>(DPD,U_DOT_DEL_RHO_H_SD,IBAR,JBAR,KBAR)

            ! ---- CPU block: RTRM = R_H_G / RHOP; DP = RTRM * DP ----
            RTRMD => GM%WORK1D
            RTRM => WORK1
            CALL DP_MUL_RTRM_GPU<<<GRDD,BLKD,0,GM%STRM>>>(DPD,RTRMD,R_H_GD,RHOPD,IBAR,JBAR,KBAR)
            ! ---- CPU block: SPECIES_ADVECTION_PART_1 / PART_2 ----
            IF (N_TRACKED_SPECIES >= 1) THEN
                RHO_Z_PD         => GM%WORK4D
                U_DOT_DEL_RHO_ZD => GM%WORK5D

                DO N=1,N_TRACKED_SPECIES
                    CALL SINGLE_SPECIES_ADVECTION_INIT_GPU<<<GRDD,BLKD,0,GM%STRM>>>( &
                        RHOPD,ZZPD,RHO_Z_PD,U_DOT_DEL_RHO_ZD,GM%SWORK1D,GM%SWORK2D,GM%SWORK3D,IBP1,JBP1,KBP1,N)

                    CALL DIVG_GET_SCALAR_FACE_VALUE_KERNEL_XYZ<<<GRDD,BLKD,0,GM%STRM>>>( &
                        UUD,VVD,WWD,RHO_Z_PD,GM%SWORK1D,GM%SWORK2D,GM%SWORK3D, &
                        IBAR,JBAR,KBAR,N,TWO_EPSILON_EB)

                    CALL DIVG_GET_SCALAR_FACE_VALUE_KERNEL_BOUNDARY<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>( &
                        GM%WALL_TOPO_AOS,GM%CELL_INDEXD,GM%CELL_WALL_INDEX,GM%CELL_SOLID,GM%CELL_EXTERIOR,GM%WALL_RHO_F, &
                        RHO_Z_PD,UUD,VVD,WWD,GM%SWORK1D,GM%SWORK2D,GM%SWORK3D,N_EXTERNAL_WALL_CELLS, &
                        N_INTERNAL_WALL_CELLS,N,TWO_EPSILON_EBD,IBP1,JBP1,KBP1,GM%CELL_LB,GM%CELL_UB)

                    CALL SINGLE_SPECIES_ADVECTION_BOUNDARY_GPU<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>(U_DOT_DEL_RHO_ZD, RHOPD, ZZPD, UUD, VVD, WWD, &
                        GM%WALL_RHO_F, GM%WALL_ZZ_F, GM%WALL_U_NORMAL, GM%WALL_U_NORMAL_S, GM%WALL_RDN, GM%UVW_SAVED, GM%WALL_TOPO_AOS, &
                        N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,MERGE(1,0,PREDICTOR),N,IBP1,JBP1,KBP1,GM%WALL_LB,GM%WALL_UB)

                    CALL SINGLE_SPECIES_ADVECTION_CORRECTION_GPU<<<GRDD,BLKD,0,GM%STRM>>>( &
                        DPD,U_DOT_DEL_RHO_ZD,GM%DEL_RHO_D_DEL_ZD,GM%SWORK1D,GM%SWORK2D,GM%SWORK3D,RHO_Z_PD, &
                        UUD,VVD,WWD,RHOPD,GM%TMPD,GM%RSUMD,R_H_GD,H_SENS_ZD,SPECIES_RCOND,GM%CELL_INDEXD,GM%CELL_SOLID,GM%CELL_WALL_INDEX, &
                        GM%RDXD,GM%RDYD,GM%RDZD,R0,IBAR,JBAR,KBAR,N,GM%CELL_LB,GM%CELL_UB)
                ENDDO
            ENDIF

            IF (N_REACTIONS > 0 .OR. N_LP_ARRAY_INDICES > 0 .OR. &
                ANY(SPECIES_MIXTURE%DEPOSITING) .OR. &
                ANY(SPECIES_MIXTURE%CONDENSATION_SMIX_INDEX > 0)) THEN
                CALL DP_ADD_D_SOURCE_GPU<<<GRDD,BLKD,0,GM%STRM>>>(DPD,GM%D_SOURCED,IBAR,JBAR,KBAR)
            ENDIF

            IF (STRATIFICATION) THEN
                CALL STRATIFICATION_GPU<<<GRDD,BLKD,0,GM%STRM>>>(DPD,RTRMD,WWD,GM%RHO_0D,GVEC(3),IBAR,JBAR,KBAR)
            ENDIF

            IF (N_ZONE > 0) THEN
                ! Pressure-zone reduction bridge: accumulate zone integrals on the GPU,
                ! then copy the compact zone arrays back for the existing CPU zone algebra.
                BLKD_SUM  = 256
                NCELL_SUM = IBAR*JBAR*KBAR
                GRDD_SUM  = (NCELL_SUM + BLKD_SUM - 1) / BLKD_SUM

                IF (.NOT.ALLOCATED(GM%DSUM_BLOCK_D)) THEN
                    ALLOCATE(GM%DSUM_BLOCK_D(1:N_ZONE),GM%PSUM_BLOCK_D(1:N_ZONE),GM%USUM_BLOCK_D(1:N_ZONE))
                ELSEIF (SIZE(GM%DSUM_BLOCK_D) < N_ZONE) THEN
                    DEALLOCATE(GM%DSUM_BLOCK_D,GM%PSUM_BLOCK_D,GM%USUM_BLOCK_D)
                    ALLOCATE(GM%DSUM_BLOCK_D(1:N_ZONE),GM%PSUM_BLOCK_D(1:N_ZONE),GM%USUM_BLOCK_D(1:N_ZONE))
                ENDIF
                IF (ALLOCATED(GM%DSUM_ZONE_PARTIAL_D)) THEN
                    IF (SIZE(GM%DSUM_ZONE_PARTIAL_D,1) < GRDD_SUM .OR. SIZE(GM%DSUM_ZONE_PARTIAL_D,2) < N_ZONE) THEN
                        DEALLOCATE(GM%DSUM_ZONE_PARTIAL_D,GM%PSUM_ZONE_PARTIAL_D)
                    ENDIF
                ENDIF
                IF (.NOT.ALLOCATED(GM%DSUM_ZONE_PARTIAL_D)) THEN
                    ALLOCATE(GM%DSUM_ZONE_PARTIAL_D(1:GRDD_SUM,1:N_ZONE),GM%PSUM_ZONE_PARTIAL_D(1:GRDD_SUM,1:N_ZONE))
                ENDIF
                GM%DSUM_BLOCK_D(1:N_ZONE) = 0._EB
                GM%PSUM_BLOCK_D(1:N_ZONE) = 0._EB
                GM%USUM_BLOCK_D(1:N_ZONE) = 0._EB

                CALL D_SUM_P_SUM_ZONES_GPU<<<GRDD_SUM,BLKD_SUM,0,GM%STRM>>>( &
                    GM%DSUM_ZONE_PARTIAL_D,GM%PSUM_ZONE_PARTIAL_D,DPD,PBAR_PD,RTRMD,GM%DXD,GM%DYD,GM%DZD,IBAR,JBAR,KBAR, &
                    GM%CELL_INDEXD,GM%CELL_SOLID,GM%PRESSURE_ZONED,N_ZONE,GRDD_SUM)
                CALL FINALIZE_D_SUM_P_SUM_ZONES_GPU<<<N_ZONE,BLKD_SUM,0,GM%STRM>>>( &
                    GM%DSUM_BLOCK_D,GM%PSUM_BLOCK_D,GM%DSUM_ZONE_PARTIAL_D,GM%PSUM_ZONE_PARTIAL_D,GRDD_SUM,N_ZONE)

                BLKD_SUM = 256
                GRDD_SUM = (N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS + BLKD_SUM - 1) / BLKD_SUM

                CALL USUM_WALL_ZONES_GPU<<<GRDD_SUM,BLKD_SUM,0,GM%STRM>>>(GM%USUM_BLOCK_D, GM%WALL_U_NORMAL, GM%WALL_U_NORMAL_S, GM%WALL_AREA, &
                    GM%WALL_PRESSURE_ZONE, GM%WALL_TOPO_AOS, N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS, N_ZONE, &
                    MERGE(1,0,PREDICTOR))

                ALLOCATE(DSUM_ZONE_H(1:N_ZONE),PSUM_ZONE_H(1:N_ZONE),USUM_ZONE_H(1:N_ZONE))
                CALL GPU_WAIT(GM%STRM)
                DSUM_ZONE_H = GM%DSUM_BLOCK_D(1:N_ZONE)
                PSUM_ZONE_H = GM%PSUM_BLOCK_D(1:N_ZONE)
                USUM_ZONE_H = GM%USUM_BLOCK_D(1:N_ZONE)
                DSUM(1:N_ZONE) = DSUM(1:N_ZONE) + DSUM_ZONE_H(1:N_ZONE)
                PSUM(1:N_ZONE) = PSUM(1:N_ZONE) + PSUM_ZONE_H(1:N_ZONE)
                USUM(1:N_ZONE) = USUM(1:N_ZONE) + USUM_ZONE_H(1:N_ZONE)
                DEALLOCATE(DSUM_ZONE_H,PSUM_ZONE_H,USUM_ZONE_H)

            ENDIF

            IF (GPU_DEBUG_H2D_D2H) THEN
                CALL GPU_WAIT(GM%STRM)
                IF (CHECK_VN) D_Z_MAX = GM%D_Z_MAXD
                WORK1 = GM%WORK1D
                DP = DPD
                DEL_RHO_D_DEL_Z = GM%DEL_RHO_D_DEL_ZD
                WALL_U_NORMAL_H_TMP   = GM%WALL_U_NORMAL
                WALL_U_NORMAL_S_H_TMP = GM%WALL_U_NORMAL_S
                WALL_K_G_H_TMP        = WALL_K_G_D_TMP
                DO N=1,N_TRACKED_SPECIES
                    WALL_RHO_D_DZDN_F_H_TMP(:,N) = GM%WALL_RHO_D_DZDN_F(:,N)
                ENDDO
                DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
                    IF (WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE
                    B1 => BOUNDARY_PROP1(WALL(IW)%B1_INDEX)
                    B1%K_G = WALL_K_G_H_TMP(IW)
                    IF (ALLOCATED(B1%RHO_D_DZDN_F)) B1%RHO_D_DZDN_F(1:N_TRACKED_SPECIES) = &
                        WALL_RHO_D_DZDN_F_H_TMP(IW,1:N_TRACKED_SPECIES)
                    IF (PREDICTOR) THEN
                        B1%U_NORMAL_S = WALL_U_NORMAL_S_H_TMP(IW)
                        IF (IW<=N_EXTERNAL_WALL_CELLS) EXTERNAL_WALL(IW)%DUNDT = (B1%U_NORMAL_S-B1%U_NORMAL)/DT
                    ELSE
                        B1%U_NORMAL = WALL_U_NORMAL_H_TMP(IW)
                        IF (IW<=N_EXTERNAL_WALL_CELLS) EXTERNAL_WALL(IW)%DUNDT = EXTERNAL_WALL(IW)%DUNDT + &
                            2._EB*(B1%U_NORMAL-B1%U_NORMAL_S)/DT
                    ENDIF
                ENDDO
            ENDIF
            IF (GPU_DEBUG_H2D_D2H) THEN
                DEALLOCATE(WALL_U_NORMAL_H_TMP,WALL_U_NORMAL_S_H_TMP,WALL_U_NORMAL_0_H_TMP,WALL_AREA_H_TMP, &
                    WALL_AREA_ADJUST_H_TMP, &
                    WALL_K_G_H_TMP,WALL_DUNDT_H_TMP,WALL_ZZ_F_H_TMP,WALL_RHO_D_F_H_TMP,WALL_RHO_D_DZDN_F_H_TMP, &
                    WALL_PRESSURE_ZONE_H_TMP,WALL_Q_CON_F_H_TMP,WALL_RHO_F_H_TMP,WALL_TMP_F_H_TMP)
            ENDIF

        END SUBROUTINE DIVERGENCE_PART_1_GPU

        !--------------------------------------------------------------------------
        ! Divergence part 1 wall-input synchronization
        !--------------------------------------------------------------------------
        SUBROUTINE SYNC_DIV1_WALL_INPUTS_TO_GPU(NM)
            INTEGER, INTENT(IN) :: NM
            INTEGER :: IW, N, ISTAT
            REAL(EB), ALLOCATABLE, DIMENSION(:) :: WALL_Q_CON_F_H, WALL_RHO_F_H, WALL_TMP_F_H
            REAL(EB), ALLOCATABLE, DIMENSION(:) :: WALL_U_NORMAL_H, WALL_U_NORMAL_S_H, WALL_U_NORMAL_0_H
            REAL(EB), ALLOCATABLE, DIMENSION(:) :: WALL_AREA_H, WALL_AREA_ADJUST_H, WALL_K_G_H, WALL_DUNDT_H
            REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: WALL_ZZ_F_H, WALL_RHO_D_F_H, WALL_RHO_D_DZDN_F_H
            INTEGER, ALLOCATABLE, DIMENSION(:) :: WALL_PRESSURE_ZONE_H
            TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)

            ALLOCATE(WALL_Q_CON_F_H(LBOUND(WALL,1):UBOUND(WALL,1)))
            ALLOCATE(WALL_RHO_F_H  (LBOUND(WALL,1):UBOUND(WALL,1)))
            ALLOCATE(WALL_TMP_F_H  (LBOUND(WALL,1):UBOUND(WALL,1)))
            ALLOCATE(WALL_U_NORMAL_H  (LBOUND(WALL,1):UBOUND(WALL,1)))
            ALLOCATE(WALL_U_NORMAL_S_H(LBOUND(WALL,1):UBOUND(WALL,1)))
            ALLOCATE(WALL_U_NORMAL_0_H(LBOUND(WALL,1):UBOUND(WALL,1)))
            ALLOCATE(WALL_AREA_H       (LBOUND(WALL,1):UBOUND(WALL,1)))
            ALLOCATE(WALL_AREA_ADJUST_H(LBOUND(WALL,1):UBOUND(WALL,1)))
            ALLOCATE(WALL_K_G_H        (LBOUND(WALL,1):UBOUND(WALL,1)))
            ALLOCATE(WALL_DUNDT_H      (LBOUND(WALL,1):UBOUND(WALL,1)))
            ALLOCATE(WALL_ZZ_F_H        (LBOUND(WALL,1):UBOUND(WALL,1),1:N_TOTAL_SCALARS))
            ALLOCATE(WALL_RHO_D_F_H     (LBOUND(WALL,1):UBOUND(WALL,1),1:N_TRACKED_SPECIES))
            ALLOCATE(WALL_RHO_D_DZDN_F_H(LBOUND(WALL,1):UBOUND(WALL,1),1:N_TRACKED_SPECIES))
            ALLOCATE(WALL_PRESSURE_ZONE_H(LBOUND(WALL,1):UBOUND(WALL,1)))

            WALL_Q_CON_F_H = 0._EB
            WALL_RHO_F_H   = 0._EB
            WALL_TMP_F_H   = 0._EB
            WALL_U_NORMAL_H   = 0._EB
            WALL_U_NORMAL_S_H = 0._EB
            WALL_U_NORMAL_0_H = 0._EB
            WALL_AREA_H        = 0._EB
            WALL_AREA_ADJUST_H = 1._EB
            WALL_K_G_H         = 0._EB
            WALL_DUNDT_H       = 0._EB
            WALL_ZZ_F_H         = 0._EB
            WALL_RHO_D_F_H      = 0._EB
            WALL_RHO_D_DZDN_F_H = 0._EB
            IF (N_TOTAL_SCALARS >= 1) WALL_ZZ_F_H(:,1) = 1._EB
            WALL_PRESSURE_ZONE_H = 0

            DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
                IF (WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE
                B1 => BOUNDARY_PROP1(WALL(IW)%B1_INDEX)
                WALL_Q_CON_F_H(IW)      = B1%Q_CON_F
                WALL_RHO_F_H(IW)        = B1%RHO_F
                WALL_TMP_F_H(IW)        = B1%TMP_F
                WALL_U_NORMAL_H(IW)     = B1%U_NORMAL
                WALL_U_NORMAL_S_H(IW)   = B1%U_NORMAL_S
                WALL_U_NORMAL_0_H(IW)   = B1%U_NORMAL_0
                WALL_AREA_H(IW)         = B1%AREA
                WALL_AREA_ADJUST_H(IW)  = B1%AREA_ADJUST
                WALL_K_G_H(IW)          = B1%K_G
                IF (IW <= N_EXTERNAL_WALL_CELLS) WALL_DUNDT_H(IW) = EXTERNAL_WALL(IW)%DUNDT
                IF (ALLOCATED(B1%ZZ_F)) WALL_ZZ_F_H(IW,1:N_TOTAL_SCALARS) = B1%ZZ_F(1:N_TOTAL_SCALARS)
                IF (ALLOCATED(B1%RHO_D_F)) WALL_RHO_D_F_H(IW,1:N_TRACKED_SPECIES) = B1%RHO_D_F(1:N_TRACKED_SPECIES)
                IF (ALLOCATED(B1%RHO_D_DZDN_F)) WALL_RHO_D_DZDN_F_H(IW,1:N_TRACKED_SPECIES) = &
                    B1%RHO_D_DZDN_F(1:N_TRACKED_SPECIES)
                WALL_PRESSURE_ZONE_H(IW)= B1%PRESSURE_ZONE
            ENDDO

            GM%WALL_Q_CON_F      = WALL_Q_CON_F_H
            GM%WALL_RHO_F        = WALL_RHO_F_H
            GM%WALL_TMP_F        = WALL_TMP_F_H
            GM%WALL_U_NORMAL     = WALL_U_NORMAL_H
            GM%WALL_U_NORMAL_S   = WALL_U_NORMAL_S_H
            GM%WALL_U_NORMAL_0   = WALL_U_NORMAL_0_H
            GM%WALL_AREA         = WALL_AREA_H
            GM%WALL_AREA_ADJUST  = WALL_AREA_ADJUST_H
            GM%WALL_K_G          = WALL_K_G_H
            GM%WALL_DUNDTD       = WALL_DUNDT_H
            DO N=1,N_TOTAL_SCALARS
                GM%WALL_ZZ_F(:,N) = WALL_ZZ_F_H(:,N)
            ENDDO
            DO N=1,N_TRACKED_SPECIES
                GM%WALL_RHO_D_F(:,N)      = WALL_RHO_D_F_H(:,N)
                GM%WALL_RHO_D_DZDN_F(:,N) = WALL_RHO_D_DZDN_F_H(:,N)
            ENDDO
            GM%WALL_PRESSURE_ZONE = WALL_PRESSURE_ZONE_H

            DEALLOCATE(WALL_Q_CON_F_H,WALL_RHO_F_H,WALL_TMP_F_H,STAT=ISTAT)
            DEALLOCATE(WALL_U_NORMAL_H,WALL_U_NORMAL_S_H,WALL_U_NORMAL_0_H,STAT=ISTAT)
            DEALLOCATE(WALL_AREA_H,WALL_AREA_ADJUST_H,WALL_K_G_H,WALL_DUNDT_H,STAT=ISTAT)
            DEALLOCATE(WALL_ZZ_F_H,WALL_RHO_D_F_H,WALL_RHO_D_DZDN_F_H,WALL_PRESSURE_ZONE_H,STAT=ISTAT)
        END SUBROUTINE SYNC_DIV1_WALL_INPUTS_TO_GPU
        !--------------------------------------------------------------------------
        ! Wall, diffusion, enthalpy, and divergence-source kernels
        !--------------------------------------------------------------------------
        ATTRIBUTES(GLOBAL) SUBROUTINE PREDICT_NORMAL_VELOCITY_GPU( &
            WALL_U_NORMAL,WALL_U_NORMAL_S,WALL_U_NORMAL_0,WALL_RHO_F,WALL_T_IGN,WALL_TOPO,WALL_SF_INDEX, &
            SF_SPECIES_BC_INDEX,SF_MASS_FLUX_TOTAL,SF_RAMP_VELO_INDEX,SF_RAMP_VELO_TAU, &
            RAMP_T_MIN,RAMP_SPAN,RAMP_RDT,RAMP_NPOINTS,RAMP_OFFSET,RAMP_INTERP, &
            UD,VD,WD,N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,PREDICTOR_FLAG,T_IN,DT_IN,T_BEGIN_IN,RHOA_IN,IBP1,JBP1)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS, PREDICTOR_FLAG,IBP1,JBP1
            REAL(EB), VALUE :: T_IN,DT_IN,T_BEGIN_IN,RHOA_IN
            REAL(EB), DEVICE, DIMENSION(0:*) :: WALL_U_NORMAL,WALL_U_NORMAL_S,WALL_U_NORMAL_0,WALL_RHO_F,WALL_T_IGN
            INTEGER, DEVICE, DIMENSION(0:*) :: WALL_SF_INDEX,SF_SPECIES_BC_INDEX,SF_RAMP_VELO_INDEX
            INTEGER, DEVICE, DIMENSION(0:*) :: RAMP_NPOINTS,RAMP_OFFSET
            REAL(EB), DEVICE, DIMENSION(0:*) :: SF_MASS_FLUX_TOTAL,SF_RAMP_VELO_TAU,RAMP_T_MIN,RAMP_SPAN,RAMP_RDT
            REAL(EB), DEVICE, DIMENSION(1:*) :: RAMP_INTERP
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: UD, VD, WD
            INTEGER :: IW, II, JJ, KK, IOR, BTYPE, SF_IDX
            REAL(EB) :: TSI,TIME_RAMP_FACTOR

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X

            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS) RETURN

            BTYPE = WALL_TOPO(IW)%BOUNDARY_TYPE
            IF (PREDICTOR_FLAG == 1) THEN
                IF (BTYPE == NULL_BOUNDARY) THEN
                    WALL_U_NORMAL_S(IW) = 0._EB
                    RETURN
                ENDIF

                IF (BTYPE == SOLID_BOUNDARY) THEN
                    SF_IDX = WALL_SF_INDEX(IW)
                    IF (SF_SPECIES_BC_INDEX(SF_IDX) == SPECIFIED_MASS_FLUX) RETURN
                    IF (ABS(WALL_T_IGN(IW)-T_BEGIN_IN) < SPACING(WALL_T_IGN(IW)) .AND. &
                        SF_RAMP_VELO_INDEX(SF_IDX) >= 1) THEN
                        TSI = T_IN + DT_IN
                    ELSE
                        TSI = T_IN + DT_IN - WALL_T_IGN(IW)
                    ENDIF
                    IF (TSI < 0._EB) THEN
                        WALL_U_NORMAL_S(IW) = 0._EB
                        RETURN
                    ENDIF
                    TIME_RAMP_FACTOR = EVALUATE_RAMP_DEVICE(TSI,SF_RAMP_VELO_INDEX(SF_IDX),SF_RAMP_VELO_TAU(SF_IDX), &
                        RAMP_T_MIN,RAMP_SPAN,RAMP_RDT,RAMP_NPOINTS,RAMP_OFFSET,RAMP_INTERP)
                    WALL_U_NORMAL_S(IW) = TIME_RAMP_FACTOR*WALL_U_NORMAL_0(IW)
                    IF (ABS(SF_MASS_FLUX_TOTAL(SF_IDX)) >= TWO_EPSILON_EB) THEN
                        WALL_U_NORMAL_S(IW) = WALL_U_NORMAL_S(IW)*RHOA_IN/WALL_RHO_F(IW)
                    ENDIF
                    RETURN
                ENDIF

                IF (BTYPE == OPEN_BOUNDARY .OR. BTYPE == INTERPOLATED_BOUNDARY) THEN
                    II  = WALL_TOPO(IW)%II
                    JJ  = WALL_TOPO(IW)%JJ
                    KK  = WALL_TOPO(IW)%KK
                    IOR = WALL_TOPO(IW)%IOR
                    SELECT CASE(IOR)
                        CASE( 1); WALL_U_NORMAL_S(IW) = -UD(II,JJ,KK)
                        CASE(-1); WALL_U_NORMAL_S(IW) =  UD(II-1,JJ,KK)
                        CASE( 2); WALL_U_NORMAL_S(IW) = -VD(II,JJ,KK)
                        CASE(-2); WALL_U_NORMAL_S(IW) =  VD(II,JJ-1,KK)
                        CASE( 3); WALL_U_NORMAL_S(IW) = -WD(II,JJ,KK)
                        CASE(-3); WALL_U_NORMAL_S(IW) =  WD(II,JJ,KK-1)
                    END SELECT
                ENDIF
            ELSE
                IF (BTYPE == SOLID_BOUNDARY) THEN
                    SF_IDX = WALL_SF_INDEX(IW)
                    IF (SF_SPECIES_BC_INDEX(SF_IDX) == SPECIFIED_MASS_FLUX) RETURN
                ENDIF
                WALL_U_NORMAL(IW) = WALL_U_NORMAL_S(IW)
            ENDIF
        END SUBROUTINE PREDICT_NORMAL_VELOCITY_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE UPDATE_WALL_DUNDT_GPU(WALL_DUNDT,WALL_U_NORMAL,WALL_U_NORMAL_S, &
            N_EXTERNAL_WALL_CELLS,DT,PREDICTOR_FLAG)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS, PREDICTOR_FLAG
            REAL(EB), VALUE :: DT
            REAL(EB), DEVICE, DIMENSION(0:*) :: WALL_DUNDT, WALL_U_NORMAL, WALL_U_NORMAL_S
            INTEGER :: IW

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS) RETURN

            IF (PREDICTOR_FLAG == 1) THEN
                WALL_DUNDT(IW) = (WALL_U_NORMAL_S(IW)-WALL_U_NORMAL(IW))/DT
            ELSE
                WALL_DUNDT(IW) = WALL_DUNDT(IW) + 2._EB*(WALL_U_NORMAL(IW)-WALL_U_NORMAL_S(IW))/DT
            ENDIF

        END SUBROUTINE UPDATE_WALL_DUNDT_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE SPECIES_RHO_D_VLES_GPU( &
            MUD,RHO_D,RSC_IN,D_Z_MAX,RHOP,IBP1,JBP1,KBP1,CHECK_VN_FLAG)
            INTEGER, VALUE :: IBP1, JBP1, KBP1, CHECK_VN_FLAG
            REAL(EB), VALUE :: RSC_IN
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: MUD, RHO_D, D_Z_MAX, RHOP
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBP1) RETURN
            IF (J < 0 .OR. J > JBP1) RETURN
            IF (K < 0 .OR. K > KBP1) RETURN

            RHO_D(I,J,K) = MAX(0._EB,MUD(I,J,K))*RSC_IN
            IF (CHECK_VN_FLAG == 1) D_Z_MAX(I,J,K) = MAX(D_Z_MAX(I,J,K),RHO_D(I,J,K)/(RHOP(I,J,K)+TWO_EPSILON_EB))
        END SUBROUTINE SPECIES_RHO_D_VLES_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE SPECIES_DIFFUSIVE_FLUX_GPU( &
            ZZPD,RHO_D,RHO_D_DZDX,RHO_D_DZDY,RHO_D_DZDZ,RDXN,RDYN,RDZN,IBAR,JBAR,KBAR,N)
            INTEGER, VALUE :: IBAR, JBAR, KBAR, N
            REAL(EB), DEVICE, DIMENSION(0:*) :: RDXN, RDYN, RDZN
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: RHO_D
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:KBAR+1,1:*) :: ZZPD, RHO_D_DZDX, RHO_D_DZDY, RHO_D_DZDZ
            INTEGER :: I,J,K
            REAL(EB) :: DZDX, DZDY, DZDZ

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBAR) RETURN
            IF (J < 0 .OR. J > JBAR) RETURN
            IF (K < 0 .OR. K > KBAR) RETURN

            DZDX = (ZZPD(I+1,J,K,N)-ZZPD(I,J,K,N))*RDXN(I)
            DZDY = (ZZPD(I,J+1,K,N)-ZZPD(I,J,K,N))*RDYN(J)
            DZDZ = (ZZPD(I,J,K+1,N)-ZZPD(I,J,K,N))*RDZN(K)

            RHO_D_DZDX(I,J,K,N) = 0.5_EB*(RHO_D(I+1,J,K)+RHO_D(I,J,K))*DZDX
            RHO_D_DZDY(I,J,K,N) = 0.5_EB*(RHO_D(I,J+1,K)+RHO_D(I,J,K))*DZDY
            RHO_D_DZDZ(I,J,K,N) = 0.5_EB*(RHO_D(I,J,K+1)+RHO_D(I,J,K))*DZDZ
        END SUBROUTINE SPECIES_DIFFUSIVE_FLUX_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE SPECIES_DIFFUSIVE_FLUX_BOUNDARY_GPU(RHO_D_DZDX, RHO_D_DZDY, RHO_D_DZDZ, &
            WALL_RHO_D_DZDN_F,WALL_TOPO,EXTERNAL_WALL_NICD,N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,N, &
            IBP1,JBP1,KBP1,WALL_LB,WALL_UB)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,N,IBP1,JBP1,KBP1,WALL_LB,WALL_UB
            INTEGER, DEVICE, DIMENSION(1:*) :: EXTERNAL_WALL_NICD
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:KBP1,1:*) :: RHO_D_DZDX,RHO_D_DZDY,RHO_D_DZDZ
            REAL(EB), DEVICE, DIMENSION(WALL_LB:WALL_UB,1:*) :: WALL_RHO_D_DZDN_F
            INTEGER :: IW, IOR, IIG, JJG, KKG

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS) RETURN
            IF (WALL_TOPO(IW)%BOUNDARY_TYPE == NULL_BOUNDARY) RETURN

            IOR = WALL_TOPO(IW)%IOR
            IIG = WALL_TOPO(IW)%IIG
            JJG = WALL_TOPO(IW)%JJG
            KKG = WALL_TOPO(IW)%KKG

            IF (WALL_TOPO(IW)%BOUNDARY_TYPE == OPEN_BOUNDARY .OR. WALL_TOPO(IW)%BOUNDARY_TYPE == INTERPOLATED_BOUNDARY) THEN
                IF (IW <= N_EXTERNAL_WALL_CELLS .AND. EXTERNAL_WALL_NICD(IW) > 1) THEN
                    SELECT CASE(IOR)
                        CASE( 1); RHO_D_DZDX(IIG-1,JJG  ,KKG  ,N) =  WALL_RHO_D_DZDN_F(IW,N)
                        CASE(-1); RHO_D_DZDX(IIG  ,JJG  ,KKG  ,N) = -WALL_RHO_D_DZDN_F(IW,N)
                        CASE( 2); RHO_D_DZDY(IIG  ,JJG-1,KKG  ,N) =  WALL_RHO_D_DZDN_F(IW,N)
                        CASE(-2); RHO_D_DZDY(IIG  ,JJG  ,KKG  ,N) = -WALL_RHO_D_DZDN_F(IW,N)
                        CASE( 3); RHO_D_DZDZ(IIG  ,JJG  ,KKG-1,N) =  WALL_RHO_D_DZDN_F(IW,N)
                        CASE(-3); RHO_D_DZDZ(IIG  ,JJG  ,KKG  ,N) = -WALL_RHO_D_DZDN_F(IW,N)
                    END SELECT
                ELSE
                    SELECT CASE(IOR)
                        CASE( 1); WALL_RHO_D_DZDN_F(IW,N) =  RHO_D_DZDX(IIG-1,JJG  ,KKG  ,N)
                        CASE(-1); WALL_RHO_D_DZDN_F(IW,N) = -RHO_D_DZDX(IIG  ,JJG  ,KKG  ,N)
                        CASE( 2); WALL_RHO_D_DZDN_F(IW,N) =  RHO_D_DZDY(IIG  ,JJG-1,KKG  ,N)
                        CASE(-2); WALL_RHO_D_DZDN_F(IW,N) = -RHO_D_DZDY(IIG  ,JJG  ,KKG  ,N)
                        CASE( 3); WALL_RHO_D_DZDN_F(IW,N) =  RHO_D_DZDZ(IIG  ,JJG  ,KKG-1,N)
                        CASE(-3); WALL_RHO_D_DZDN_F(IW,N) = -RHO_D_DZDZ(IIG  ,JJG  ,KKG  ,N)
                    END SELECT
                ENDIF
                RETURN
            ENDIF

            SELECT CASE(IOR)
                CASE( 1); RHO_D_DZDX(IIG-1,JJG  ,KKG  ,N) = 0._EB
                CASE(-1); RHO_D_DZDX(IIG  ,JJG  ,KKG  ,N) = 0._EB
                CASE( 2); RHO_D_DZDY(IIG  ,JJG-1,KKG  ,N) = 0._EB
                CASE(-2); RHO_D_DZDY(IIG  ,JJG  ,KKG  ,N) = 0._EB
                CASE( 3); RHO_D_DZDZ(IIG  ,JJG  ,KKG-1,N) = 0._EB
                CASE(-3); RHO_D_DZDZ(IIG  ,JJG  ,KKG  ,N) = 0._EB
            END SELECT
        END SUBROUTINE SPECIES_DIFFUSIVE_FLUX_BOUNDARY_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE SPECIES_DIFF_HEAT_FACE_GPU( &
            TMPD,H_SENS_ZD,RHO_D_DZDX,RHO_D_DZDY,RHO_D_DZDZ,H_RHO_D_DZDX,H_RHO_D_DZDY,H_RHO_D_DZDZ,IBAR,JBAR,KBAR,N)
            INTEGER, VALUE :: IBAR, JBAR, KBAR, N
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: TMPD, H_RHO_D_DZDX, H_RHO_D_DZDY, H_RHO_D_DZDZ
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:KBAR+1,1:*) :: RHO_D_DZDX, RHO_D_DZDY, RHO_D_DZDZ
            REAL(EB), DEVICE, DIMENSION(0:5000,1:*) :: H_SENS_ZD
            INTEGER :: I,J,K,ITMP
            REAL(EB) :: TMP_G, FRAC, H_S

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBAR) RETURN
            IF (J < 0 .OR. J > JBAR) RETURN
            IF (K < 0 .OR. K > KBAR) RETURN

            TMP_G = 0.5_EB*(TMPD(I+1,J,K)+TMPD(I,J,K))
            ITMP = MAX(0,MIN(4999,INT(TMP_G)))
            FRAC = TMP_G - REAL(ITMP,EB)
            H_S = H_SENS_ZD(ITMP,N) + FRAC*(H_SENS_ZD(ITMP+1,N)-H_SENS_ZD(ITMP,N))
            H_RHO_D_DZDX(I,J,K) = H_S*RHO_D_DZDX(I,J,K,N)

            TMP_G = 0.5_EB*(TMPD(I,J+1,K)+TMPD(I,J,K))
            ITMP = MAX(0,MIN(4999,INT(TMP_G)))
            FRAC = TMP_G - REAL(ITMP,EB)
            H_S = H_SENS_ZD(ITMP,N) + FRAC*(H_SENS_ZD(ITMP+1,N)-H_SENS_ZD(ITMP,N))
            H_RHO_D_DZDY(I,J,K) = H_S*RHO_D_DZDY(I,J,K,N)

            TMP_G = 0.5_EB*(TMPD(I,J,K+1)+TMPD(I,J,K))
            ITMP = MAX(0,MIN(4999,INT(TMP_G)))
            FRAC = TMP_G - REAL(ITMP,EB)
            H_S = H_SENS_ZD(ITMP,N) + FRAC*(H_SENS_ZD(ITMP+1,N)-H_SENS_ZD(ITMP,N))
            H_RHO_D_DZDZ(I,J,K) = H_S*RHO_D_DZDZ(I,J,K,N)
        END SUBROUTINE SPECIES_DIFF_HEAT_FACE_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE SPECIES_DIFF_HEAT_BOUNDARY_GPU(TMPD, H_SENS_ZD, ZZPD, RHO_D_DZDX, RHO_D_DZDY, &
            RHO_D_DZDZ, H_RHO_D_DZDX, H_RHO_D_DZDY, H_RHO_D_DZDZ, WALL_ZZ_F, WALL_RHO_D_F, WALL_RHO_D_DZDN_F, WALL_TMP_F, &
            WALL_U_NORMAL, WALL_U_NORMAL_S, WALL_RDN, WALL_TOPO, N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS, &
            N_TRACKED_SPECIES,N,PREDICTOR_FLAG,IBP1,JBP1,KBP1,WALL_LB,WALL_UB)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,N_TRACKED_SPECIES,N,PREDICTOR_FLAG
            INTEGER, VALUE :: IBP1,JBP1,KBP1,WALL_LB,WALL_UB
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: TMPD,H_RHO_D_DZDX,H_RHO_D_DZDY,H_RHO_D_DZDZ
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:KBP1,1:*) :: ZZPD,RHO_D_DZDX,RHO_D_DZDY,RHO_D_DZDZ
            REAL(EB), DEVICE, DIMENSION(0:5000,1:*) :: H_SENS_ZD
            REAL(EB), DEVICE, DIMENSION(WALL_LB:WALL_UB,1:*) :: WALL_ZZ_F,WALL_RHO_D_F,WALL_RHO_D_DZDN_F
            REAL(EB), DEVICE, DIMENSION(0:*) :: WALL_TMP_F, WALL_U_NORMAL, WALL_U_NORMAL_S, WALL_RDN
            INTEGER :: IW, IOR, IIG, JJG, KKG, NN, N_ZZ_MAX, ITMP
            REAL(EB) :: RHO_D_DZDN, SUM_GET, TMP_F_GAS, H_S, FRAC, UN_P

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS) RETURN
            IF (WALL_TOPO(IW)%BOUNDARY_TYPE == NULL_BOUNDARY) RETURN
            IF (WALL_TOPO(IW)%BOUNDARY_TYPE == OPEN_BOUNDARY .OR. WALL_TOPO(IW)%BOUNDARY_TYPE == INTERPOLATED_BOUNDARY) RETURN

            IOR = WALL_TOPO(IW)%IOR
            IIG = WALL_TOPO(IW)%IIG
            JJG = WALL_TOPO(IW)%JJG
            KKG = WALL_TOPO(IW)%KKG

            N_ZZ_MAX = MAXLOC(WALL_ZZ_F(IW,1:N_TRACKED_SPECIES),1)

            RHO_D_DZDN = 2._EB*WALL_RHO_D_F(IW,N)*(ZZPD(IIG,JJG,KKG,N)-WALL_ZZ_F(IW,N))*WALL_RDN(IW)
            IF (N == N_ZZ_MAX) THEN
                SUM_GET = 0._EB
                DO NN=1,N_TRACKED_SPECIES
                    SUM_GET = SUM_GET + 2._EB*WALL_RHO_D_F(IW,NN)*(ZZPD(IIG,JJG,KKG,NN)-WALL_ZZ_F(IW,NN))*WALL_RDN(IW)
                ENDDO
                RHO_D_DZDN = -(SUM_GET-RHO_D_DZDN)
            ENDIF
            WALL_RHO_D_DZDN_F(IW,N) = RHO_D_DZDN

            IF (PREDICTOR_FLAG == 1) THEN
                UN_P = WALL_U_NORMAL_S(IW)
            ELSE
                UN_P = WALL_U_NORMAL(IW)
            ENDIF

            IF (WALL_TOPO(IW)%BOUNDARY_TYPE == SOLID_BOUNDARY .AND. UN_P > 0._EB) THEN
                TMP_F_GAS = TMPD(IIG,JJG,KKG)
            ELSE
                TMP_F_GAS = WALL_TMP_F(IW)
            ENDIF

            ITMP = MAX(0,MIN(4999,INT(TMP_F_GAS)))
            FRAC = TMP_F_GAS - REAL(ITMP,EB)
            H_S = H_SENS_ZD(ITMP,N) + FRAC*(H_SENS_ZD(ITMP+1,N)-H_SENS_ZD(ITMP,N))

            SELECT CASE(IOR)
                CASE( 1)
                    RHO_D_DZDX(IIG-1,JJG  ,KKG  ,N) =  RHO_D_DZDN
                    H_RHO_D_DZDX(IIG-1,JJG  ,KKG    ) =  RHO_D_DZDN*H_S
                CASE(-1)
                    RHO_D_DZDX(IIG  ,JJG  ,KKG  ,N) = -RHO_D_DZDN
                    H_RHO_D_DZDX(IIG  ,JJG  ,KKG    ) = -RHO_D_DZDN*H_S
                CASE( 2)
                    RHO_D_DZDY(IIG  ,JJG-1,KKG  ,N) =  RHO_D_DZDN
                    H_RHO_D_DZDY(IIG  ,JJG-1,KKG    ) =  RHO_D_DZDN*H_S
                CASE(-2)
                    RHO_D_DZDY(IIG  ,JJG  ,KKG  ,N) = -RHO_D_DZDN
                    H_RHO_D_DZDY(IIG  ,JJG  ,KKG    ) = -RHO_D_DZDN*H_S
                CASE( 3)
                    RHO_D_DZDZ(IIG  ,JJG  ,KKG-1,N) =  RHO_D_DZDN
                    H_RHO_D_DZDZ(IIG  ,JJG  ,KKG-1  ) =  RHO_D_DZDN*H_S
                CASE(-3)
                    RHO_D_DZDZ(IIG  ,JJG  ,KKG  ,N) = -RHO_D_DZDN
                    H_RHO_D_DZDZ(IIG  ,JJG  ,KKG    ) = -RHO_D_DZDN*H_S
            END SELECT
        END SUBROUTINE SPECIES_DIFF_HEAT_BOUNDARY_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE SPECIES_DIFF_HEAT_DIV_GPU( &
            DP,DEL_RHO_D_DEL_Z,RHO_D_DZDX,RHO_D_DZDY,RHO_D_DZDZ,H_RHO_D_DZDX,H_RHO_D_DZDY,H_RHO_D_DZDZ,RDX,RDY,RDZ,IBAR,JBAR,KBAR,N)
            INTEGER, VALUE :: IBAR, JBAR, KBAR, N
            REAL(EB), DEVICE, DIMENSION(0:*) :: RDX, RDY, RDZ
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: DP, H_RHO_D_DZDX, H_RHO_D_DZDY, H_RHO_D_DZDZ
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:KBAR+1,1:*) :: RHO_D_DZDX, RHO_D_DZDY, RHO_D_DZDZ
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:KBAR+1,1:*) :: DEL_RHO_D_DEL_Z
            INTEGER :: I,J,K
            REAL(EB) :: DIV_DIFF_HEAT_FLUX, DEL_RHO_D_DEL_Z_LOCAL

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            DIV_DIFF_HEAT_FLUX = (H_RHO_D_DZDX(I,J,K)-H_RHO_D_DZDX(I-1,J,K))*RDX(I) + &
                                 (H_RHO_D_DZDY(I,J,K)-H_RHO_D_DZDY(I,J-1,K))*RDY(J) + &
                                 (H_RHO_D_DZDZ(I,J,K)-H_RHO_D_DZDZ(I,J,K-1))*RDZ(K)
            DP(I,J,K) = DP(I,J,K) + DIV_DIFF_HEAT_FLUX

            DEL_RHO_D_DEL_Z_LOCAL = (RHO_D_DZDX(I,J,K,N)-RHO_D_DZDX(I-1,J,K,N))*RDX(I) + &
                                    (RHO_D_DZDY(I,J,K,N)-RHO_D_DZDY(I,J-1,K,N))*RDY(J) + &
                                    (RHO_D_DZDZ(I,J,K,N)-RHO_D_DZDZ(I,J,K-1,N))*RDZ(K)
            DEL_RHO_D_DEL_Z(I,J,K,N) = DEL_RHO_D_DEL_Z_LOCAL
        END SUBROUTINE SPECIES_DIFF_HEAT_DIV_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE USUM_WALL_ZONES_GPU(USUM_BLOCK, WALL_U_NORMAL, WALL_U_NORMAL_S, WALL_AREA, &
            WALL_PRESSURE_ZONE, WALL_TOPO, N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS, N_ZONE_LOCAL, PREDICTOR_FLAG)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS, N_ZONE_LOCAL, PREDICTOR_FLAG
            REAL(EB), DEVICE, DIMENSION(1:*) :: USUM_BLOCK
            REAL(EB), DEVICE, DIMENSION(0:*) :: WALL_U_NORMAL, WALL_U_NORMAL_S, WALL_AREA
            INTEGER, DEVICE, DIMENSION(0:*) :: WALL_PRESSURE_ZONE
            INTEGER :: IW, IPZ, N_WALL, IERR
            REAL(EB) :: SUM_LOCAL

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            N_WALL = N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS

            IF (IW >= 1 .AND. IW <= N_WALL) THEN
                IPZ = WALL_PRESSURE_ZONE(IW)
                IF (WALL_TOPO(IW)%BOUNDARY_TYPE == SOLID_BOUNDARY .AND. IPZ >= 1 .AND. IPZ <= N_ZONE_LOCAL) THEN
                    IF (PREDICTOR_FLAG == 1) THEN
                        SUM_LOCAL = WALL_U_NORMAL_S(IW) * WALL_AREA(IW)
                    ELSE
                        SUM_LOCAL = WALL_U_NORMAL(IW) * WALL_AREA(IW)
                    ENDIF
                    IERR = atomicAdd(USUM_BLOCK(IPZ), SUM_LOCAL)
                ENDIF
            ENDIF

        END SUBROUTINE USUM_WALL_ZONES_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE ZERO_DIFFUSIVE_WORK_GPU(DEL_RHO_D_DEL_Z,SWORK1,SWORK2,SWORK3, &
            IBP1,JBP1,KBP1,N_TOTAL_SCALARS_IN)
            INTEGER, VALUE :: IBP1,JBP1,KBP1,N_TOTAL_SCALARS_IN
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:KBP1,1:*) :: DEL_RHO_D_DEL_Z,SWORK1,SWORK2,SWORK3
            INTEGER :: I,J,K,N

            I = (BLOCKIDX%X - 1)*BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1)*BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1)*BLOCKDIM%Z + THREADIDX%Z - 1
            IF (I<0 .OR. I>IBP1 .OR. J<0 .OR. J>JBP1 .OR. K<0 .OR. K>KBP1) RETURN

            DO N=1,N_TOTAL_SCALARS_IN
                DEL_RHO_D_DEL_Z(I,J,K,N) = 0._EB
                SWORK1(I,J,K,N) = 0._EB
                SWORK2(I,J,K,N) = 0._EB
                SWORK3(I,J,K,N) = 0._EB
            ENDDO
        END SUBROUTINE ZERO_DIFFUSIVE_WORK_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE GET_SPECIFIC_HEAT_GPU(ZZPD,CP,CP_ZD,R_H_G,TMPD,N_TRACKED_SPECIES,IBAR,JBAR,KBAR)
            INTEGER, VALUE :: N_TRACKED_SPECIES, IBAR, JBAR, KBAR
            REAL(EB) , DEVICE , DIMENSION(0:5000,1:*) :: CP_ZD
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:KBAR+1,1:*) :: ZZPD
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*)    :: CP, TMPD, R_H_G
            INTEGER :: I,J,K,ITMP,N
            REAL(EB) :: CP_LOC

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            ITMP = MIN(5000, NINT(TMPD(I,J,K)))
            CP_LOC = 0._EB
            DO N=1,N_TRACKED_SPECIES
            CP_LOC = CP_LOC + CP_ZD(ITMP,N) * ZZPD(I,J,K,N)
            ENDDO

            CP(I,J,K)    = CP_LOC
            R_H_G(I,J,K) = 1._EB / (CP_LOC * TMPD(I,J,K))
        END SUBROUTINE GET_SPECIFIC_HEAT_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE D_Z_MAX_THERMAL_GPU(D_Z_MAX,KP,CP,RHOP,CELL_INDEXD,CELL_SOLID,IBAR,JBAR,KBAR)
            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: D_Z_MAX, KP, CP, RHOP
            INTEGER, DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*)  :: CELL_INDEXD
            LOGICAL, DEVICE, DIMENSION(0:*)        :: CELL_SOLID
            INTEGER :: I,J,K,IC

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            IC = CELL_INDEXD(I,J,K)
            IF (CELL_SOLID(IC)) RETURN

            D_Z_MAX(I,J,K) = MAX(D_Z_MAX(I,J,K),KP(I,J,K)/(CP(I,J,K)*RHOP(I,J,K)))
        END SUBROUTINE D_Z_MAX_THERMAL_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE KDTDN_GPU_1(MUD,CPOPRD,KP,IBAR,JBAR,KBAR)
            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), VALUE :: CPOPRD
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: MUD, KP
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBAR+1) RETURN
            IF (J < 0 .OR. J > JBAR+1) RETURN
            IF (K < 0 .OR. K > KBAR+1) RETURN

            KP(I,J,K) = MUD(I,J,K) * CPOPRD
        END SUBROUTINE KDTDN_GPU_1


        ATTRIBUTES(GLOBAL) SUBROUTINE KDTDN_GPU_2(TMPD,KP,KDTDX,KDTDY,KDTDZ,IBAR,JBAR,KBAR,RDXND,RDYND,RDZND)
            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), DEVICE, DIMENSION(0:*)       :: RDXND, RDYND, RDZND
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: TMPD, KP, KDTDX, KDTDY, KDTDZ
            REAL(EB) :: DTDX, DTDY, DTDZ
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBAR) RETURN
            IF (J < 0 .OR. J > JBAR) RETURN
            IF (K < 0 .OR. K > KBAR) RETURN

            DTDX = (TMPD(I+1,J,K)-TMPD(I,J,K))*RDXND(I)
            KDTDX(I,J,K) = .5_EB*(KP(I+1,J,K)+KP(I,J,K))*DTDX

            DTDY = (TMPD(I,J+1,K)-TMPD(I,J,K))*RDYND(J)
            KDTDY(I,J,K) = .5_EB*(KP(I,J+1,K)+KP(I,J,K))*DTDY

            DTDZ = (TMPD(I,J,K+1)-TMPD(I,J,K))*RDZND(K)
            KDTDZ(I,J,K) = .5_EB*(KP(I,J,K+1)+KP(I,J,K))*DTDZ
        END SUBROUTINE KDTDN_GPU_2


        ATTRIBUTES(GLOBAL) SUBROUTINE KDTDN_GPU_BOUNDARY(DP, WALL_Q_CON_F, WALL_AREA_ADJUST, WALL_TOPO, WALL_RDN, KP, &
            WALL_K_G,KDTDX,KDTDY,KDTDZ,N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,IBP1,JBP1)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,IBP1,JBP1
            REAL(EB), DEVICE, DIMENSION(0:*)       :: WALL_Q_CON_F, WALL_AREA_ADJUST, WALL_RDN, WALL_K_G
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: DP,KP,KDTDX,KDTDY,KDTDZ
            INTEGER :: IW, II, JJ, KK, IIG, JJG, KKG, IERR
            REAL(EB) :: QTERM

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X

            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS) RETURN
            IF (WALL_TOPO(IW)%BOUNDARY_TYPE == NULL_BOUNDARY) RETURN

            II  = WALL_TOPO(IW)%II
            JJ  = WALL_TOPO(IW)%JJ
            KK  = WALL_TOPO(IW)%KK
            IIG = WALL_TOPO(IW)%IIG
            JJG = WALL_TOPO(IW)%JJG
            KKG = WALL_TOPO(IW)%KKG

            IF (WALL_TOPO(IW)%BOUNDARY_TYPE == INTERPOLATED_BOUNDARY) RETURN

            IF (WALL_TOPO(IW)%BOUNDARY_TYPE == OPEN_BOUNDARY) THEN
                WALL_K_G(IW) = 0.5_EB*(KP(IIG,JJG,KKG)+KP(II,JJ,KK))
                RETURN
            ENDIF

            WALL_K_G(IW) = KP(IIG,JJG,KKG)
            QTERM = -WALL_AREA_ADJUST(IW) * WALL_Q_CON_F(IW) * WALL_RDN(IW)
            IERR = atomicAdd(DP(IIG,JJG,KKG), QTERM)
            SELECT CASE(WALL_TOPO(IW)%IOR)
            CASE( 1) ; KDTDX(II  ,JJ  ,KK  ) = 0._EB
            CASE(-1) ; KDTDX(II-1,JJ  ,KK  ) = 0._EB
            CASE( 2) ; KDTDY(II  ,JJ  ,KK  ) = 0._EB
            CASE(-2) ; KDTDY(II  ,JJ-1,KK  ) = 0._EB
            CASE( 3) ; KDTDZ(II  ,JJ  ,KK  ) = 0._EB
            CASE(-3) ; KDTDZ(II  ,JJ  ,KK-1) = 0._EB
            END SELECT
        END SUBROUTINE KDTDN_GPU_BOUNDARY


        ATTRIBUTES(GLOBAL) SUBROUTINE ENTHALPY_ADVECTION_GPU_INIT( &
            H_SENS_ZD, U_DOT_DEL_RHO_H_S, FX_H_S, FY_H_S, FZ_H_S, RHO_H_S_P, &
            RHOPD, ZZPD, TMPD, N_TRACKED_SPECIES, IBP1, JBP1, KBP1)

            IMPLICIT NONE

            INTEGER, VALUE :: N_TRACKED_SPECIES, IBP1, JBP1, KBP1

            REAL(EB), DEVICE, DIMENSION(0:5000,1:*)        :: H_SENS_ZD
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*)    :: U_DOT_DEL_RHO_H_S, FX_H_S, FY_H_S, FZ_H_S
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*)    :: RHO_H_S_P, RHOPD, TMPD
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:KBP1,1:*) :: ZZPD

            INTEGER :: I, J, K
            INTEGER :: ITMP, N
            REAL(EB) :: HS_LOC, DOT_H0, DOT_DH, ZZ_N
            REAL(EB) :: TMP_LOC
            REAL(EB) :: FRAC
            REAL(EB) :: H0, H1

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBP1) RETURN
            IF (J < 0 .OR. J > JBP1) RETURN
            IF (K < 0 .OR. K > KBP1) RETURN

            ! 这三个先保留，避免 face 上未覆盖位置读到脏值
            FX_H_S(I,J,K) = 0._EB
            FY_H_S(I,J,K) = 0._EB
            FZ_H_S(I,J,K) = 0._EB
            U_DOT_DEL_RHO_H_S(I,J,K) = 0._EB

            TMP_LOC = TMPD(I,J,K)

            ITMP = MIN(4999, INT(TMP_LOC))
            FRAC = TMP_LOC - REAL(ITMP, EB)

            DOT_H0 = 0._EB
            DOT_DH = 0._EB
            DO N = 1, N_TRACKED_SPECIES
                ZZ_N = ZZPD(I,J,K,N)
                H0 = H_SENS_ZD(ITMP  ,N)
                H1 = H_SENS_ZD(ITMP+1,N)
                DOT_H0 = DOT_H0 + H0*ZZ_N
                DOT_DH = DOT_DH + (H1-H0)*ZZ_N
            ENDDO
            HS_LOC = DOT_H0 + FRAC*DOT_DH

            RHO_H_S_P(I,J,K) = RHOPD(I,J,K) * HS_LOC

        END SUBROUTINE ENTHALPY_ADVECTION_GPU_INIT


        ATTRIBUTES(GLOBAL) SUBROUTINE ENTHALPY_ADVECTION_GPU_END(U_DOT_DEL_RHO_H_S,FX_H_S,FY_H_S,FZ_H_S,RHO_H_S_P,UU,VV,WW, &
            CELL_INDEXD,CELL_SOLID,CELL_WALL_INDEX,RDX,RDY,RDZ,IBAR,JBAR,KBAR,CELL_LB,CELL_UB)
            INTEGER, VALUE :: IBAR, JBAR, KBAR,CELL_LB,CELL_UB
            INTEGER, DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: CELL_INDEXD
            INTEGER, DEVICE, DIMENSION(CELL_LB:CELL_UB,-3:*) :: CELL_WALL_INDEX
            LOGICAL, DEVICE, DIMENSION(0:*)       :: CELL_SOLID
            REAL(EB), DEVICE, DIMENSION(0:*)      :: RDX, RDY, RDZ
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: U_DOT_DEL_RHO_H_S, FX_H_S, FY_H_S, FZ_H_S, RHO_H_S_P, UU, VV, WW

            INTEGER :: I,J,K,IC
            REAL(EB) :: DU_P,DU_M,DV_P,DV_M,DW_P,DW_M
            REAL(EB) :: MXP,MXM,MYP,MYM,MZP,MZM
            REAL(EB) :: RHOHS

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            IC = CELL_INDEXD(I,J,K)
            IF (CELL_SOLID(IC)) RETURN

            RHOHS = RHO_H_S_P(I,J,K)

            MXP = MERGE(1._EB,0._EB,CELL_WALL_INDEX(IC, 1)==0)
            MXM = MERGE(1._EB,0._EB,CELL_WALL_INDEX(IC,-1)==0)
            MYP = MERGE(1._EB,0._EB,CELL_WALL_INDEX(IC, 2)==0)
            MYM = MERGE(1._EB,0._EB,CELL_WALL_INDEX(IC,-2)==0)
            MZP = MERGE(1._EB,0._EB,CELL_WALL_INDEX(IC, 3)==0)
            MZM = MERGE(1._EB,0._EB,CELL_WALL_INDEX(IC,-3)==0)

            DU_P = MXP * (FX_H_S(I  ,J  ,K  ) - RHOHS) * UU(I  ,J  ,K  )
            DU_M = MXM * (FX_H_S(I-1,J  ,K  ) - RHOHS) * UU(I-1,J  ,K  )

            DV_P = MYP * (FY_H_S(I  ,J  ,K  ) - RHOHS) * VV(I  ,J  ,K  )
            DV_M = MYM * (FY_H_S(I  ,J-1,K  ) - RHOHS) * VV(I  ,J-1,K  )

            DW_P = MZP * (FZ_H_S(I  ,J  ,K  ) - RHOHS) * WW(I  ,J  ,K  )
            DW_M = MZM * (FZ_H_S(I  ,J  ,K-1) - RHOHS) * WW(I  ,J  ,K-1)

            U_DOT_DEL_RHO_H_S(I,J,K) = U_DOT_DEL_RHO_H_S(I,J,K) + &
                (DU_P-DU_M)*RDX(I) + (DV_P-DV_M)*RDY(J) + (DW_P-DW_M)*RDZ(K)

        END SUBROUTINE ENTHALPY_ADVECTION_GPU_END


        ATTRIBUTES(GLOBAL) SUBROUTINE ENTHALPY_ADVECTION_BOUNDARY_GPU(U_DOT_DEL_RHO_H_S, RHO_H_S_P, UU, VV, WW, TMPD, &
            H_SENS_ZD, WALL_TMP_F, WALL_RHO_F, WALL_U_NORMAL, WALL_U_NORMAL_S, WALL_ZZ_F, WALL_RDN, UVW_SAVE_D, WALL_TOPO, &
            N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,N_TRACKED_SPECIES,PREDICTOR_FLAG,IBP1,JBP1,WALL_LB,WALL_UB)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,N_TRACKED_SPECIES,PREDICTOR_FLAG
            INTEGER, VALUE :: IBP1,JBP1,WALL_LB,WALL_UB
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: U_DOT_DEL_RHO_H_S,RHO_H_S_P,UU,VV,WW,TMPD
            REAL(EB), DEVICE, DIMENSION(0:5000,1:*)     :: H_SENS_ZD
            REAL(EB), DEVICE, DIMENSION(WALL_LB:WALL_UB,1:*) :: WALL_ZZ_F
            REAL(EB), DEVICE, DIMENSION(0:*)       :: WALL_TMP_F, WALL_RHO_F, WALL_U_NORMAL, WALL_U_NORMAL_S, WALL_RDN
            REAL(EB), DEVICE, DIMENSION(1:*)       :: UVW_SAVE_D
            INTEGER :: IW, IOR, II, JJ, KK, IIG, JJG, KKG, ITMP, IERR, N
            REAL(EB) :: UN, UN_P, TMP_F_GAS, FRAC, H_S, H0, H1, DOT_H0, DOT_DH, ZZ_N, DU, SGN

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS) RETURN
            IF (WALL_TOPO(IW)%BOUNDARY_TYPE == NULL_BOUNDARY) RETURN

            IOR = WALL_TOPO(IW)%IOR
            II  = WALL_TOPO(IW)%II
            JJ  = WALL_TOPO(IW)%JJ
            KK  = WALL_TOPO(IW)%KK
            IIG = WALL_TOPO(IW)%IIG
            JJG = WALL_TOPO(IW)%JJG
            KKG = WALL_TOPO(IW)%KKG
            SGN = SIGN(1._EB,REAL(IOR,EB))

            IF (WALL_TOPO(IW)%BOUNDARY_TYPE == SOLID_BOUNDARY) THEN
                IF (PREDICTOR_FLAG == 1) THEN
                    UN_P = WALL_U_NORMAL_S(IW)
                ELSE
                    UN_P = WALL_U_NORMAL(IW)
                ENDIF
                UN = -SGN*UN_P
                IF (UN_P > 0._EB) THEN
                    TMP_F_GAS = TMPD(IIG,JJG,KKG)
                ELSE
                    TMP_F_GAS = WALL_TMP_F(IW)
                ENDIF
            ELSEIF (WALL_TOPO(IW)%BOUNDARY_TYPE == INTERPOLATED_BOUNDARY) THEN
                ! Periodic/interpolated walls use the saved normal velocity from the CPU path.
                UN = UVW_SAVE_D(IW)
                TMP_F_GAS = WALL_TMP_F(IW)
            ELSE
                SELECT CASE(IOR)
                    CASE( 1); UN = UU(II  ,JJ  ,KK  )
                    CASE(-1); UN = UU(II-1,JJ  ,KK  )
                    CASE( 2); UN = VV(II  ,JJ  ,KK  )
                    CASE(-2); UN = VV(II  ,JJ-1,KK  )
                    CASE( 3); UN = WW(II  ,JJ  ,KK  )
                    CASE(-3); UN = WW(II  ,JJ  ,KK-1)
                END SELECT
                TMP_F_GAS = WALL_TMP_F(IW)
            ENDIF

            ITMP = MAX(0,MIN(4999,INT(TMP_F_GAS)))
            FRAC = TMP_F_GAS - REAL(ITMP,EB)
            DOT_H0 = 0._EB
            DOT_DH = 0._EB
            DO N=1,N_TRACKED_SPECIES
                ZZ_N = WALL_ZZ_F(IW,N)
                H0 = H_SENS_ZD(ITMP,N)
                H1 = H_SENS_ZD(ITMP+1,N)
                DOT_H0 = DOT_H0 + H0*ZZ_N
                DOT_DH = DOT_DH + (H1-H0)*ZZ_N
            ENDDO
            H_S = DOT_H0 + FRAC*DOT_DH
            DU = (WALL_RHO_F(IW)*H_S - RHO_H_S_P(IIG,JJG,KKG))*UN
            IERR = atomicAdd(U_DOT_DEL_RHO_H_S(IIG,JJG,KKG), -SGN*DU*WALL_RDN(IW))
        END SUBROUTINE ENTHALPY_ADVECTION_BOUNDARY_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE DP_ADD_DELKDELT_GPU(DP,KDTDX,KDTDY,KDTDZ,QD,QRD,RDX,RDY,RDZ,IBAR,JBAR,KBAR)
            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), DEVICE, DIMENSION(0:*)       :: RDX, RDY, RDZ
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: DP, KDTDX, KDTDY, KDTDZ, QD, QRD
            REAL(EB) :: DELKDELT
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            DELKDELT = (KDTDX(I,J,K)-KDTDX(I-1,J,K))*RDX(I) + &
                    (KDTDY(I,J,K)-KDTDY(I,J-1,K))*RDY(J) + &
                    (KDTDZ(I,J,K)-KDTDZ(I,J,K-1))*RDZ(K)

            DP(I,J,K) = DP(I,J,K) + DELKDELT + QD(I,J,K) + QRD(I,J,K)
        END SUBROUTINE DP_ADD_DELKDELT_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE DP_ADD_D_SOURCE_GPU(DP,D_SOURCE_D,IBAR,JBAR,KBAR)
            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: DP, D_SOURCE_D
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            DP(I,J,K) = DP(I,J,K) + D_SOURCE_D(I,J,K)
        END SUBROUTINE DP_ADD_D_SOURCE_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE DP_SUB_ENTHALPY_ADV_GPU(DP,U_DOT_DEL_RHO_H_S,IBAR,JBAR,KBAR)
            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: DP, U_DOT_DEL_RHO_H_S
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            DP(I,J,K) = DP(I,J,K) - U_DOT_DEL_RHO_H_S(I,J,K)
        END SUBROUTINE DP_SUB_ENTHALPY_ADV_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE DP_MUL_RTRM_GPU(DP,RTRM,R_H_G,RHOP,IBAR,JBAR,KBAR)
            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: DP, RTRM, R_H_G, RHOP
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            RTRM(I,J,K) = R_H_G(I,J,K) / RHOP(I,J,K)
            DP(I,J,K)   = RTRM(I,J,K) * DP(I,J,K)
        END SUBROUTINE DP_MUL_RTRM_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE SINGLE_SPECIES_ADVECTION_INIT_GPU( &
            RHOPD,ZZPD,RHO_Z_P,U_DOT_DEL_RHO_Z,FX_ZZ,FY_ZZ,FZ_ZZ,IBP1,JBP1,KBP1,N)
            INTEGER, VALUE :: IBP1, JBP1, KBP1, N
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: RHOPD,RHO_Z_P,U_DOT_DEL_RHO_Z
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:KBP1,1:*) :: FX_ZZ,FY_ZZ,FZ_ZZ
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:KBP1,1:*) :: ZZPD

            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBP1) RETURN
            IF (J < 0 .OR. J > JBP1) RETURN
            IF (K < 0 .OR. K > KBP1) RETURN

            RHO_Z_P(I,J,K) = RHOPD(I,J,K)*ZZPD(I,J,K,N)
            U_DOT_DEL_RHO_Z(I,J,K) = 0._EB
            FX_ZZ(I,J,K,N) = 0._EB
            FY_ZZ(I,J,K,N) = 0._EB
            FZ_ZZ(I,J,K,N) = 0._EB

        END SUBROUTINE SINGLE_SPECIES_ADVECTION_INIT_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE SINGLE_SPECIES_ADVECTION_CORRECTION_GPU( &
            DP,U_DOT_DEL_RHO_Z,DEL_RHO_D_DEL_Z,FX_ZZ,FY_ZZ,FZ_ZZ,RHO_Z_P, &
            UU,VV,WW,RHOPD,TMPD,RSUMD,R_H_GD,H_SENS_ZD_IN,SPECIES_RCON_IN,CELL_INDEXD,CELL_SOLID,CELL_WALL_INDEX, &
            RDX,RDY,RDZ,R0_LOCAL,IBAR,JBAR,KBAR,N,CELL_LB,CELL_UB)

            INTEGER, VALUE :: IBAR, JBAR, KBAR, N,CELL_LB,CELL_UB
            REAL(EB), VALUE :: R0_LOCAL
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: DP,U_DOT_DEL_RHO_Z,RHO_Z_P
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:KBAR+1,1:*) :: FX_ZZ,FY_ZZ,FZ_ZZ
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:KBAR+1,1:*) :: DEL_RHO_D_DEL_Z
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: UU,VV,WW,RHOPD,TMPD,RSUMD,R_H_GD
            REAL(EB), DEVICE, DIMENSION(0:5000,1:*) :: H_SENS_ZD_IN
            REAL(EB), DEVICE, DIMENSION(*) :: SPECIES_RCON_IN
            INTEGER, DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: CELL_INDEXD
            INTEGER, DEVICE, DIMENSION(CELL_LB:CELL_UB,-3:*) :: CELL_WALL_INDEX
            LOGICAL, DEVICE, DIMENSION(0:*) :: CELL_SOLID
            REAL(EB), DEVICE, DIMENSION(0:*) :: RDX,RDY,RDZ

            INTEGER :: I,J,K,IC,ITMP
            REAL(EB) :: DU_P,DU_M,DV_P,DV_M,DW_P,DW_M
            REAL(EB) :: MXP,MXM,MYP,MYM,MZP,MZM
            REAL(EB) :: RHOZ, UDOT, TMP_LOC, FRAC, H_S

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            IC = CELL_INDEXD(I,J,K)
            IF (CELL_SOLID(IC)) RETURN

            RHOZ = RHO_Z_P(I,J,K)

            MXP = MERGE(1._EB,0._EB,CELL_WALL_INDEX(IC, 1)==0)
            MXM = MERGE(1._EB,0._EB,CELL_WALL_INDEX(IC,-1)==0)
            MYP = MERGE(1._EB,0._EB,CELL_WALL_INDEX(IC, 2)==0)
            MYM = MERGE(1._EB,0._EB,CELL_WALL_INDEX(IC,-2)==0)
            MZP = MERGE(1._EB,0._EB,CELL_WALL_INDEX(IC, 3)==0)
            MZM = MERGE(1._EB,0._EB,CELL_WALL_INDEX(IC,-3)==0)

            DU_P = MXP * (FX_ZZ(I  ,J  ,K  ,N) - RHOZ) * UU(I  ,J  ,K)
            DU_M = MXM * (FX_ZZ(I-1,J  ,K  ,N) - RHOZ) * UU(I-1,J  ,K)
            DV_P = MYP * (FY_ZZ(I  ,J  ,K  ,N) - RHOZ) * VV(I  ,J  ,K)
            DV_M = MYM * (FY_ZZ(I  ,J-1,K  ,N) - RHOZ) * VV(I  ,J-1,K)
            DW_P = MZP * (FZ_ZZ(I  ,J  ,K  ,N) - RHOZ) * WW(I  ,J  ,K)
            DW_M = MZM * (FZ_ZZ(I  ,J  ,K-1,N) - RHOZ) * WW(I  ,J  ,K-1)

            UDOT = U_DOT_DEL_RHO_Z(I,J,K) + (DU_P-DU_M)*RDX(I) + (DV_P-DV_M)*RDY(J) + (DW_P-DW_M)*RDZ(K)
            U_DOT_DEL_RHO_Z(I,J,K) = UDOT

            TMP_LOC = TMPD(I,J,K)
            ITMP = MIN(4999,INT(TMP_LOC))
            FRAC = TMP_LOC - REAL(ITMP,EB)
            H_S = H_SENS_ZD_IN(ITMP,N) + FRAC*(H_SENS_ZD_IN(ITMP+1,N)-H_SENS_ZD_IN(ITMP,N))
            DP(I,J,K) = DP(I,J,K) + (SPECIES_RCON_IN(N)/RSUMD(I,J,K) - H_S*R_H_GD(I,J,K))* &
                (DEL_RHO_D_DEL_Z(I,J,K,N) - UDOT)/RHOPD(I,J,K)

        END SUBROUTINE SINGLE_SPECIES_ADVECTION_CORRECTION_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE SINGLE_SPECIES_ADVECTION_BOUNDARY_GPU(U_DOT_DEL_RHO_Z, RHOPD, ZZPD, UU, VV, WW, &
            WALL_RHO_F, WALL_ZZ_F, WALL_U_NORMAL, WALL_U_NORMAL_S, WALL_RDN, UVW_SAVE_D, WALL_TOPO, N_EXTERNAL_WALL_CELLS, &
            N_INTERNAL_WALL_CELLS,PREDICTOR_FLAG,N,IBP1,JBP1,KBP1,WALL_LB,WALL_UB)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,PREDICTOR_FLAG,N
            INTEGER, VALUE :: IBP1,JBP1,KBP1,WALL_LB,WALL_UB
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: U_DOT_DEL_RHO_Z,RHOPD,UU,VV,WW
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:KBP1,1:*) :: ZZPD
            REAL(EB), DEVICE, DIMENSION(0:*)          :: WALL_RHO_F, WALL_U_NORMAL, WALL_U_NORMAL_S, WALL_RDN
            REAL(EB), DEVICE, DIMENSION(WALL_LB:WALL_UB,1:*) :: WALL_ZZ_F
            REAL(EB), DEVICE, DIMENSION(1:*)          :: UVW_SAVE_D
            INTEGER :: IW, IOR, II, JJ, KK, IIG, JJG, KKG, IERR
            REAL(EB) :: UN, DU, SGN, ZZ_F_LOCAL

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS) RETURN
            IF (WALL_TOPO(IW)%BOUNDARY_TYPE == NULL_BOUNDARY) RETURN

            IOR = WALL_TOPO(IW)%IOR
            II  = WALL_TOPO(IW)%II
            JJ  = WALL_TOPO(IW)%JJ
            KK  = WALL_TOPO(IW)%KK
            IIG = WALL_TOPO(IW)%IIG
            JJG = WALL_TOPO(IW)%JJG
            KKG = WALL_TOPO(IW)%KKG
            SGN = SIGN(1._EB,REAL(IOR,EB))

            IF (WALL_TOPO(IW)%BOUNDARY_TYPE == SOLID_BOUNDARY) THEN
                IF (PREDICTOR_FLAG == 1) THEN
                    UN = -SGN*WALL_U_NORMAL_S(IW)
                ELSE
                    UN = -SGN*WALL_U_NORMAL(IW)
                ENDIF
            ELSEIF (WALL_TOPO(IW)%BOUNDARY_TYPE == INTERPOLATED_BOUNDARY) THEN
                UN = UVW_SAVE_D(IW)
            ELSE
                SELECT CASE(IOR)
                    CASE( 1); UN = UU(II  ,JJ  ,KK  )
                    CASE(-1); UN = UU(II-1,JJ  ,KK  )
                    CASE( 2); UN = VV(II  ,JJ  ,KK  )
                    CASE(-2); UN = VV(II  ,JJ-1,KK  )
                    CASE( 3); UN = WW(II  ,JJ  ,KK  )
                    CASE(-3); UN = WW(II  ,JJ  ,KK-1)
                END SELECT
            ENDIF

            ZZ_F_LOCAL = WALL_ZZ_F(IW,N)
            DU = (WALL_RHO_F(IW)*ZZ_F_LOCAL - RHOPD(IIG,JJG,KKG)*ZZPD(IIG,JJG,KKG,N))*UN
            IERR = atomicAdd(U_DOT_DEL_RHO_Z(IIG,JJG,KKG), -SGN*DU*WALL_RDN(IW))
        END SUBROUTINE SINGLE_SPECIES_ADVECTION_BOUNDARY_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE STRATIFICATION_GPU(DP,RTRM,WW,RHO_0D,GVEC3,IBAR,JBAR,KBAR)
            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), VALUE :: GVEC3
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: DP, RTRM, WW
            REAL(EB), DEVICE, DIMENSION(0:*)       :: RHO_0D
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            DP(I,J,K) = DP(I,J,K) + RTRM(I,J,K) * 0.5_EB * (WW(I,J,K)+WW(I,J,K-1)) * RHO_0D(K) * GVEC3
        END SUBROUTINE STRATIFICATION_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE D_SUM_P_SUM_ZONES_GPU( &
            DSUM_PARTIAL,PSUM_PARTIAL,DP,PBAR_P,RTRM,DXD,DYD,DZD,IBAR,JBAR,KBAR,CELL_INDEXD,CELL_SOLID,PRESSURE_ZONED, &
            N_ZONE_LOCAL,N_BLOCKS)
            INTEGER, VALUE :: IBAR, JBAR, KBAR, N_ZONE_LOCAL,N_BLOCKS
            REAL(EB), DEVICE, DIMENSION(1:N_BLOCKS,1:*) :: DSUM_PARTIAL,PSUM_PARTIAL
            REAL(EB), DEVICE, DIMENSION(0:*)        :: DXD, DYD, DZD
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*)  :: DP, RTRM
            REAL(EB), DEVICE, DIMENSION(0:KBAR+1,0:*)    :: PBAR_P
            INTEGER, DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*)   :: CELL_INDEXD, PRESSURE_ZONED
            LOGICAL, DEVICE, DIMENSION(0:*)         :: CELL_SOLID
            REAL(EB), SHARED :: DSUM_SHARED(256), PSUM_SHARED(256)
            REAL(EB) :: VC, DSUM_LOCAL, PSUM_LOCAL
            INTEGER :: TID, GID, NCELL, I, J, K, IC, IPZ, IZ, TMP, REM, OFFSET

            TID   = THREADIDX%X
            GID   = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            NCELL = IBAR*JBAR*KBAR
            IPZ = 0
            DSUM_LOCAL = 0._EB
            PSUM_LOCAL = 0._EB

            IF (GID <= NCELL) THEN
                TMP = GID - 1
                K   = TMP / (IBAR*JBAR) + 1
                REM = MOD(TMP, IBAR*JBAR)
                J   = REM / IBAR + 1
                I   = MOD(REM, IBAR) + 1

                IC = CELL_INDEXD(I,J,K)
                IPZ = PRESSURE_ZONED(I,J,K)
                IF (.NOT. CELL_SOLID(IC) .AND. IPZ >= 1 .AND. IPZ <= N_ZONE_LOCAL) THEN
                    VC = DXD(I) * DYD(J) * DZD(K)
                    DSUM_LOCAL = VC * DP(I,J,K)
                    PSUM_LOCAL = VC * (1._EB/PBAR_P(K,IPZ) - RTRM(I,J,K))
                ELSE
                    IPZ = 0
                ENDIF
            ENDIF

            DO IZ=1,N_ZONE_LOCAL
                DSUM_SHARED(TID) = MERGE(DSUM_LOCAL,0._EB,IPZ==IZ)
                PSUM_SHARED(TID) = MERGE(PSUM_LOCAL,0._EB,IPZ==IZ)
                CALL SYNCTHREADS()

                OFFSET = BLOCKDIM%X/2
                DO WHILE (OFFSET >= 1)
                    IF (TID <= OFFSET) THEN
                        DSUM_SHARED(TID) = DSUM_SHARED(TID) + DSUM_SHARED(TID+OFFSET)
                        PSUM_SHARED(TID) = PSUM_SHARED(TID) + PSUM_SHARED(TID+OFFSET)
                    ENDIF
                    CALL SYNCTHREADS()
                    OFFSET = OFFSET/2
                ENDDO

                IF (TID == 1) THEN
                    DSUM_PARTIAL(BLOCKIDX%X,IZ) = DSUM_SHARED(1)
                    PSUM_PARTIAL(BLOCKIDX%X,IZ) = PSUM_SHARED(1)
                ENDIF
                CALL SYNCTHREADS()
            ENDDO
        END SUBROUTINE D_SUM_P_SUM_ZONES_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE FINALIZE_D_SUM_P_SUM_ZONES_GPU( &
            DSUM_BLOCK,PSUM_BLOCK,DSUM_PARTIAL,PSUM_PARTIAL,N_BLOCKS,N_ZONE_LOCAL)
            INTEGER, VALUE :: N_BLOCKS, N_ZONE_LOCAL
            REAL(EB), DEVICE, DIMENSION(1:*)    :: DSUM_BLOCK, PSUM_BLOCK
            REAL(EB), DEVICE, DIMENSION(1:N_BLOCKS,1:*) :: DSUM_PARTIAL, PSUM_PARTIAL
            REAL(EB), SHARED :: DSUM_SHARED(256), PSUM_SHARED(256)
            REAL(EB) :: DSUM_LOCAL, PSUM_LOCAL
            INTEGER :: TID, IZ, IB, OFFSET

            TID = THREADIDX%X
            IZ = BLOCKIDX%X
            IF (IZ > N_ZONE_LOCAL) RETURN

            DSUM_LOCAL = 0._EB
            PSUM_LOCAL = 0._EB
            IB = TID
            DO WHILE (IB <= N_BLOCKS)
                DSUM_LOCAL = DSUM_LOCAL + DSUM_PARTIAL(IB,IZ)
                PSUM_LOCAL = PSUM_LOCAL + PSUM_PARTIAL(IB,IZ)
                IB = IB + BLOCKDIM%X
            ENDDO

            DSUM_SHARED(TID) = DSUM_LOCAL
            PSUM_SHARED(TID) = PSUM_LOCAL
            CALL SYNCTHREADS()

            OFFSET = BLOCKDIM%X/2
            DO WHILE (OFFSET >= 1)
                IF (TID <= OFFSET) THEN
                    DSUM_SHARED(TID) = DSUM_SHARED(TID) + DSUM_SHARED(TID+OFFSET)
                    PSUM_SHARED(TID) = PSUM_SHARED(TID) + PSUM_SHARED(TID+OFFSET)
                ENDIF
                CALL SYNCTHREADS()
                OFFSET = OFFSET/2
            ENDDO

            IF (TID == 1) THEN
                DSUM_BLOCK(IZ) = DSUM_SHARED(1)
                PSUM_BLOCK(IZ) = PSUM_SHARED(1)
            ENDIF
        END SUBROUTINE FINALIZE_D_SUM_P_SUM_ZONES_GPU


        !--------------------------------------------------------------------------
        ! Divergence part 2 host orchestration and density time derivative
        !--------------------------------------------------------------------------
        SUBROUTINE DIVERGENCE_PART_2_GPU(DT,NM)

            USE COMP_FUNCTIONS, ONLY: CURRENT_TIME
            IMPLICIT NONE

            INTEGER,  INTENT(IN) :: NM
            REAL(EB), INTENT(IN) :: DT

            REAL(EB), POINTER, DIMENSION(:,:,:) :: DPH, RTRMH, DIVH
            REAL(EB), POINTER, DIMENSION(:,:)   :: PBAR_PH
            REAL(EB), POINTER, DIMENSION(:)     :: D_PBAR_DT_P
            REAL(EB), DEVICE, POINTER, DIMENSION(:,:,:) :: DPD, RTRMD
            REAL(EB), DEVICE, POINTER, DIMENSION(:,:)   :: PBAR_PD
            REAL(EB), DEVICE, POINTER, DIMENSION(:)     :: D_PBAR_DT_PD
            TYPE(DIM3) :: BLKD, GRDD
            INTEGER, PARAMETER :: BS = 8
            REAL(EB) :: TNOW, RDT
            REAL(EB) :: P_EQ, SUM_P_PSUM, SUM_USUM, SUM_DSUM, SUM_PSUM
            LOGICAL :: OPEN_ZONE
            INTEGER :: IPZ,IOPZ,IW
            REAL(EB), ALLOCATABLE, DIMENSION(:) :: USUM_ADD
            REAL(EB), ALLOCATABLE, DIMENSION(:) :: WALL_U_NORMAL_H_TMP, WALL_U_NORMAL_S_H_TMP
            TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1

            IF (SOLID_PHASE_ONLY) RETURN
            IF (PERIODIC_TEST==3) RETURN
            IF (PERIODIC_TEST==4) RETURN

            TNOW = CURRENT_TIME()
            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)

            BLKD = GM%BLKD_3D_BP1
            GRDD = GM%GRDD_3D_BP1

            RDT = 1._EB / DT

            IF (PREDICTOR) THEN
                DPH      => DS
                DPD      => GM%DS_D
                PBAR_PH  => PBAR_S
                PBAR_PD  => GM%PBAR_SD
                D_PBAR_DT_P => D_PBAR_DT_S
                D_PBAR_DT_PD => GM%D_PBAR_DT_SD
            ELSE
                DPH      => D
                DPD      => GM%D_D
                PBAR_PH  => PBAR
                PBAR_PD  => GM%PBARD
                D_PBAR_DT_P => D_PBAR_DT
                D_PBAR_DT_PD => GM%D_PBAR_DTD
            ENDIF

            RTRMH => WORK1
            RTRMD => GM%WORK1D
            DIVH  => WORK1

            IF (GPU_DEBUG_H2D_D2H) THEN
                DPD     = DPH
                PBAR_PD = PBAR_PH
                RTRMD   = RTRMH
                GM%UD      = U
                GM%VD      = V
                GM%WD      = W
                IF (CORRECTOR) THEN
                    GM%USD = US
                    GM%VSD = VS
                    GM%WSD = WS
                ENDIF

                ALLOCATE(WALL_U_NORMAL_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)))
                ALLOCATE(WALL_U_NORMAL_S_H_TMP(LBOUND(WALL,1):UBOUND(WALL,1)))
                WALL_U_NORMAL_H_TMP   = 0._EB
                WALL_U_NORMAL_S_H_TMP = 0._EB
                DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
                    IF (WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE
                    B1 => BOUNDARY_PROP1(WALL(IW)%B1_INDEX)
                    WALL_U_NORMAL_H_TMP(IW)   = B1%U_NORMAL
                    WALL_U_NORMAL_S_H_TMP(IW) = B1%U_NORMAL_S
                ENDDO
                GM%WALL_U_NORMAL   = WALL_U_NORMAL_H_TMP
                GM%WALL_U_NORMAL_S = WALL_U_NORMAL_S_H_TMP
            ENDIF

            IF (N_ZONE > 0) THEN
                ALLOCATE(USUM_ADD(1:N_ZONE))
                USUM_ADD = 0._EB

                DO IPZ=1,N_ZONE
                    SUM_P_PSUM = PBAR_PH(1,IPZ)*PSUM(IPZ)
                    OPEN_ZONE  = .FALSE.
                    SUM_USUM = USUM(IPZ)
                    SUM_DSUM = DSUM(IPZ)
                    SUM_PSUM = PSUM(IPZ)
                    DO IOPZ=N_ZONE,0,-1
                        IF (IOPZ==IPZ) CYCLE
                        IF (CONNECTED_ZONES(IPZ,IOPZ)>0) THEN
                            IF (IOPZ==0) THEN
                                OPEN_ZONE = .TRUE.
                            ELSE
                                SUM_P_PSUM = SUM_P_PSUM + PBAR_PH(1,IOPZ)*PSUM(IOPZ)
                                SUM_USUM = SUM_USUM + USUM(IOPZ)
                                SUM_DSUM = SUM_DSUM + DSUM(IOPZ)
                                SUM_PSUM = SUM_PSUM + PSUM(IOPZ)
                            ENDIF
                        ENDIF
                    ENDDO
                    IF (OPEN_ZONE) THEN
                        P_EQ          = P_0(1)
                        USUM_ADD(IPZ) = PSUM(IPZ)*(PBAR_PH(1,IPZ)-P_EQ)/PRESSURE_RELAX_TIME + DSUM(IPZ) - USUM(IPZ)
                    ELSE
                        P_EQ          = SUM_P_PSUM/SUM_PSUM
                        USUM_ADD(IPZ) = PSUM(IPZ)*(PBAR_PH(1,IPZ)-P_EQ)/PRESSURE_RELAX_TIME + DSUM(IPZ) - USUM(IPZ) - &
                                        PSUM(IPZ)*(SUM_DSUM-SUM_USUM)/SUM_PSUM
                    ENDIF
                ENDDO

                DO IPZ=1,N_ZONE
                    USUM(IPZ) = USUM(IPZ) + USUM_ADD(IPZ)
                ENDDO

                DO IPZ=1,N_ZONE
                    IF (ABS(PSUM(IPZ)) > TWO_EPSILON_EB) THEN
                        D_PBAR_DT_P(IPZ) = (DSUM(IPZ) - USUM(IPZ))/PSUM(IPZ)
                    ELSE
                        D_PBAR_DT_P(IPZ) = 0._EB
                    ENDIF
                    IF (CORRECTOR) P_ZONE(IPZ)%DPSTAR = D_PBAR_DT_P(IPZ)
                ENDDO

                D_PBAR_DT_PD(0:N_ZONE) = D_PBAR_DT_P(0:N_ZONE)
                DEALLOCATE(USUM_ADD)
            ENDIF

            IF (PREDICTOR) THEN
                CALL DP_GPU_2<<<GRDD,BLKD,0,GM%STRM>>>( &
                    DPD,RTRMD,PBAR_PD,D_PBAR_DT_PD,IBAR,JBAR,KBAR,GM%CELL_INDEXD,GM%CELL_SOLID,GM%PRESSURE_ZONED)
                CALL DP_PART2_BOUNDARY_GPU<<<(N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS+256-1)/256,256,0,GM%STRM>>>(DPD, GM%WALL_TOPO_AOS, &
                    GM%WALL_U_NORMAL, GM%WALL_U_NORMAL_S, GM%RDXD, GM%RDYD, GM%RDZD, GM%RD, GM%RRND, GM%CELL_INDEXD, GM%CELL_SOLID, &
                    N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,1,IBP1,JBP1)
                CALL DDDT_GPU_PREDICTOR<<<GRDD,BLKD,0,GM%STRM>>>(RTRMD,DPD,GM%DDDTD,GM%UD,GM%VD,GM%WD,GM%RDXD,GM%RDYD,GM%RDZD,GM%RD,GM%RRND,IBAR,JBAR,KBAR,RDT)
            ELSE
                CALL DP_GPU_2<<<GRDD,BLKD,0,GM%STRM>>>( &
                    DPD,RTRMD,PBAR_PD,D_PBAR_DT_PD,IBAR,JBAR,KBAR,GM%CELL_INDEXD,GM%CELL_SOLID,GM%PRESSURE_ZONED)
                CALL DP_PART2_BOUNDARY_GPU<<<(N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS+256-1)/256,256,0,GM%STRM>>>(DPD, GM%WALL_TOPO_AOS, &
                    GM%WALL_U_NORMAL, GM%WALL_U_NORMAL_S, GM%RDXD, GM%RDYD, GM%RDZD, GM%RD, GM%RRND, GM%CELL_INDEXD, GM%CELL_SOLID, &
                    N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,0,IBP1,JBP1)
                CALL DDDT_GPU_CORRECTOR<<<GRDD,BLKD,0,GM%STRM>>>(RTRMD,DPD,GM%DDDTD,GM%UD,GM%VD,GM%WD,GM%USD,GM%VSD,GM%WSD,GM%RDXD,GM%RDYD,GM%RDZD,GM%RD,GM%RRND,IBAR,JBAR,KBAR,RDT)
            ENDIF
            IF (GPU_DEBUG_H2D_D2H) THEN
                CALL GPU_WAIT(GM%STRM)
                DDDT = GM%DDDTD
                DIVH = RTRMD
                DPH  = DPD
            ENDIF

            IF (GPU_DEBUG_H2D_D2H) DEALLOCATE(WALL_U_NORMAL_H_TMP,WALL_U_NORMAL_S_H_TMP)

            T_USED(2) = T_USED(2) + CURRENT_TIME() - TNOW

        END SUBROUTINE DIVERGENCE_PART_2_GPU



        ATTRIBUTES(GLOBAL) SUBROUTINE DP_GPU_2(DP,RTRM,PBAR_P,D_PBAR_DT_D,IBAR,JBAR,KBAR,CELL_INDEXD,CELL_SOLID,PRESSURE_ZONED)
            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: DP, RTRM
            REAL(EB), DEVICE, DIMENSION(0:KBAR+1,0:*)    :: PBAR_P
            REAL(EB), DEVICE, DIMENSION(0:*)        :: D_PBAR_DT_D
            INTEGER, DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*)  :: CELL_INDEXD, PRESSURE_ZONED
            LOGICAL, DEVICE, DIMENSION(0:*)        :: CELL_SOLID
            INTEGER :: I,J,K,IC,IPZ

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            IC = CELL_INDEXD(I,J,K)
            IF (CELL_SOLID(IC)) THEN
                DP(I,J,K) = 0._EB
                RETURN
            ENDIF

            IPZ = PRESSURE_ZONED(I,J,K)
            IF (IPZ < 1) RETURN

            DP(I,J,K) = DP(I,J,K) - (1._EB/PBAR_P(K,IPZ) - RTRM(I,J,K)) * D_PBAR_DT_D(IPZ)
        END SUBROUTINE DP_GPU_2


        ATTRIBUTES(GLOBAL) SUBROUTINE DP_PART2_BOUNDARY_GPU(DP, WALL_TOPO, WALL_U_NORMAL, WALL_U_NORMAL_S, RDXD, RDYD, RDZD, &
            RD,RRND,CELL_INDEXD,CELL_SOLID,N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,PREDICTOR_FLAG,IBP1,JBP1)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,PREDICTOR_FLAG,IBP1,JBP1
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: DP
            REAL(EB), DEVICE, DIMENSION(0:*) :: WALL_U_NORMAL, WALL_U_NORMAL_S, RDXD, RDYD, RDZD, RD, RRND
            INTEGER, DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: CELL_INDEXD
            LOGICAL, DEVICE, DIMENSION(0:*) :: CELL_SOLID

            INTEGER :: IW, II, JJ, KK, IIG, JJG, KKG, BTYPE, IOR, IC, IERR
            REAL(EB) :: UN_P, DP_ADD

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS) RETURN

            BTYPE = WALL_TOPO(IW)%BOUNDARY_TYPE
            IF (BTYPE == NULL_BOUNDARY) RETURN

            II  = WALL_TOPO(IW)%II
            JJ  = WALL_TOPO(IW)%JJ
            KK  = WALL_TOPO(IW)%KK
            IIG = WALL_TOPO(IW)%IIG
            JJG = WALL_TOPO(IW)%JJG
            KKG = WALL_TOPO(IW)%KKG

            SELECT CASE(BTYPE)
                CASE(SOLID_BOUNDARY)
                    IC = CELL_INDEXD(II,JJ,KK)
                    IF (.NOT. CELL_SOLID(IC)) RETURN

                    IF (PREDICTOR_FLAG == 1) THEN
                        UN_P = WALL_U_NORMAL_S(IW)
                    ELSE
                        UN_P = WALL_U_NORMAL(IW)
                    ENDIF

                    IOR = WALL_TOPO(IW)%IOR
                    SELECT CASE(IOR)
                        CASE( 1,-1)
                            IF (IOR == 1) THEN
                                DP_ADD = -UN_P*RDXD(II)*RRND(II)*RD(II)
                            ELSE
                                DP_ADD = -UN_P*RDXD(II)*RRND(II)*RD(II-1)
                            ENDIF
                        CASE( 2,-2)
                            DP_ADD = -UN_P*RDYD(JJ)
                        CASE( 3,-3)
                            DP_ADD = -UN_P*RDZD(KK)
                        CASE DEFAULT
                            RETURN
                    END SELECT
                    IERR = atomicAdd(DP(II,JJ,KK), DP_ADD)
                CASE(OPEN_BOUNDARY,MIRROR_BOUNDARY,INTERPOLATED_BOUNDARY)
                    DP(II,JJ,KK) = DP(IIG,JJG,KKG)
            END SELECT

        END SUBROUTINE DP_PART2_BOUNDARY_GPU


        ATTRIBUTES(GLOBAL) SUBROUTINE DDDT_GPU_PREDICTOR(DIV,DP,DDDT,UD,VD,WD,RDXD,RDYD,RDZD,RD,RRND,IBAR,JBAR,KBAR,RDT_LOCAL)
            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: DIV, DP, DDDT, UD, VD, WD
            REAL(EB), DEVICE, DIMENSION(0:*)       :: RDXD, RDYD, RDZD, RD, RRND
            REAL(EB), VALUE :: RDT_LOCAL
            INTEGER :: I,J,K
            REAL(EB) :: DIVU, XTERM, YTERM, ZTERM

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBAR+1) RETURN
            IF (J < 0 .OR. J > JBAR+1) RETURN
            IF (K < 0 .OR. K > KBAR+1) RETURN

            IF (I >= 1 .AND. I <= IBAR .AND. J >= 1 .AND. J <= JBAR .AND. K >= 1 .AND. K <= KBAR) THEN
                XTERM = RD(I)*UD(I,J,K) - RD(I-1)*UD(I-1,J,K)
                XTERM = XTERM * RDXD(I)
                XTERM = XTERM * RRND(I)
                YTERM = VD(I,J,K) - VD(I,J-1,K)
                YTERM = YTERM * RDYD(J)
                ZTERM = WD(I,J,K) - WD(I,J,K-1)
                ZTERM = ZTERM * RDZD(K)

                DIVU = XTERM
                DIVU = DIVU + YTERM
                DIVU = DIVU + ZTERM

                DIV(I,J,K) = DIVU
            ENDIF

            DDDT(I,J,K) = (DP(I,J,K)-DIV(I,J,K)) * RDT_LOCAL
        END SUBROUTINE DDDT_GPU_PREDICTOR


        ATTRIBUTES(GLOBAL) SUBROUTINE DDDT_GPU_CORRECTOR(DIV,DP,DDDT,UD,VD,WD,USD,VSD,WSD,RDXD,RDYD,RDZD,RD,RRND,IBAR,JBAR,KBAR,RDT_LOCAL)
            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), DEVICE, VOLATILE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: DIV
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: DP, DDDT, UD, VD, WD, USD, VSD, WSD
            REAL(EB), DEVICE, DIMENSION(0:*)       :: RDXD, RDYD, RDZD, RD, RRND
            REAL(EB), VALUE :: RDT_LOCAL
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBAR+1) RETURN
            IF (J < 0 .OR. J > JBAR+1) RETURN
            IF (K < 0 .OR. K > KBAR+1) RETURN

            IF (I >= 1 .AND. I <= IBAR .AND. J >= 1 .AND. J <= JBAR .AND. K >= 1 .AND. K <= KBAR) THEN
                DIV(I,J,K) = (RD(I)*UD(I,J,K)-RD(I-1)*UD(I-1,J,K))*RDXD(I)*RRND(I)
                DIV(I,J,K) = DIV(I,J,K) + (VD(I,J,K)-VD(I,J-1,K))*RDYD(J)
                DIV(I,J,K) = DIV(I,J,K) + (WD(I,J,K)-WD(I,J,K-1))*RDZD(K)
                DIV(I,J,K) = DIV(I,J,K) + (RD(I)*USD(I,J,K)-RD(I-1)*USD(I-1,J,K))*RDXD(I)*RRND(I)
                DIV(I,J,K) = DIV(I,J,K) + (VSD(I,J,K)-VSD(I,J-1,K))*RDYD(J)
                DIV(I,J,K) = DIV(I,J,K) + (WSD(I,J,K)-WSD(I,J,K-1))*RDZD(K)
            ENDIF

            DDDT(I,J,K) = (2._EB*DP(I,J,K) - DIV(I,J,K)) * RDT_LOCAL
        END SUBROUTINE DDDT_GPU_CORRECTOR

END MODULE G_DIVG
