!------------------------------------------------------------------------------
! Module: G_VELO
! Purpose: GPU viscosity, momentum flux, velocity update, and velocity boundaries.
!------------------------------------------------------------------------------
MODULE G_VELO
    USE G_MESH_POINTERS
    USE G_FUNC
    USE G_ERRO, ONLY: GPU_WAIT

    CONTAINS
        SUBROUTINE COMPUTE_VISCOSITY_GPU(NM,DT,APPLY_TO_ESTIMATED_VARIABLES)
            IMPLICIT NONE

            LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
            REAL(EB), INTENT(IN) :: DT
        
            REAL(EB), POINTER, DEVICE , DIMENSION(:,:,:) :: RHOP,UP,VP,WP, &
                                                UP_HAT,VP_HAT,WP_HAT, &
                                                UU,VV,WW

            INTEGER , INTENT(IN) :: NM
        
            REAL(EB) , DEVICE , POINTER , DIMENSION(:,:,:,:) :: ZZP

            INTEGER :: BLKD_1D , GRDD_1D
            TYPE(DIM3) :: BLKD , GRDD , BLKD_2D , GRDD_2D_X , GRDD_2D_Y , GRDD_2D_Z
            INTEGER , PARAMETER :: BS = 8
            REAL(EB) :: T_NOW

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

            T_NOW = CURRENT_TIME()

            IF(APPLY_TO_ESTIMATED_VARIABLES) THEN

            
                RHOP => GM%RHOSD
                UU   => GM%USD
                VV   => GM%VSD
                WW   => GM%WSD
                ZZP  => GM%ZZSD

                IF (GPU_DEBUG_H2D_D2H) THEN
                    GM%RHOSD = RHOS
                    GM%USD   = US
                    GM%VSD   = VS
                    GM%WSD   = WS
                    GM%ZZSD  = ZZS
                    GM%TMPD  = TMP
                    IF (N_EXTERNAL_WALL_CELLS>0) THEN
                        GM%U_GHOSTD = U_GHOST
                        GM%V_GHOSTD = V_GHOST
                        GM%W_GHOSTD = W_GHOST
                    ENDIF
                ENDIF

            ELSE

                RHOP => GM%RHOD
                UU   => GM%UD
                VV   => GM%VD
                WW   => GM%WD
                ZZP  => GM%ZZD

                IF (GPU_DEBUG_H2D_D2H) THEN
                    GM%RHOD = RHO
                    GM%UD   = U
                    GM%VD   = V
                    GM%WD   = W
                    GM%ZZD  = ZZ
                    GM%TMPD = TMP
                    IF (N_EXTERNAL_WALL_CELLS>0) THEN
                        GM%U_GHOSTD = U_GHOST
                        GM%V_GHOSTD = V_GHOST
                        GM%W_GHOSTD = W_GHOST
                    ENDIF
                ENDIF

            ENDIF


            CALL GET_VISCOSITY_GPU<<<GRDD,BLKD,0,GM%STRM>>>(ZZP,GM%MU_DNSD,GM%TMPD,GM%CELL_INDEXD,GM%CELL_SOLID,RSQ_MW_ZD,MU_RSQMW_ZD, &
                                        N_TRACKED_SPECIES,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,I_MAX_TEMP)

            UP => GM%WORK1D
            VP => GM%WORK2D
            WP => GM%WORK3D
            UP_HAT => GM%WORK4D
            VP_HAT => GM%WORK5D
            WP_HAT => GM%WORK6D

            ! Match the CPU Deardorff path: solid-wall ghost faces remain zero
            ! unless an external boundary handler overwrites them explicitly.
            CALL ZERO_UVW_FILTER_WORK_GPU<<<GRDD,BLKD,0,GM%STRM>>>( &
                UP,VP,WP,UP_HAT,VP_HAT,WP_HAT,IBP1,JBP1,KBP1)

            CALL COMPUTE_UVWP_1<<<GRDD,BLKD,0,GM%STRM>>>(UP,VP,WP,UU,VV,WW,IBAR,JBAR,KBAR)

            CALL COMPUTE_UVWPD_2_UVW<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>( &
                UP,VP,WP,GM%U_GHOSTD,GM%V_GHOSTD,GM%W_GHOSTD,GM%WALL_TOPO_AOS,N_EXTERNAL_WALL_CELLS, &
                OPEN_BOUNDARY,MIRROR_BOUNDARY,INTERPOLATED_BOUNDARY,IBP1,JBP1)

            CALL GPU_FILL_EDGES<<<(MAX(IBAR,JBAR,KBAR)+128)/128,128,0,GM%STRM>>>(UP,VP,WP,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1)
            CALL GPU_FILL_EDGES_CORNER<<<1,1,0,GM%STRM>>>(UP,VP,WP,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1)

            CALL GPU_TEST_FILTER_1_UVW<<<GRDD,BLKD,0,GM%STRM>>>( &
                K3DTD,UP_HAT,VP_HAT,WP_HAT,UP,VP,WP,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1)

            CALL GPU_TEST_FILTER_2_I_UVW<<<GRDD_2D_X,BLKD_2D,0,GM%STRM>>>( &
                UP_HAT,VP_HAT,WP_HAT,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1)
            CALL GPU_TEST_FILTER_2_J_UVW<<<GRDD_2D_Y,BLKD_2D,0,GM%STRM>>>( &
                UP_HAT,VP_HAT,WP_HAT,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1)
            CALL GPU_TEST_FILTER_2_K_UVW<<<GRDD_2D_Z,BLKD_2D,0,GM%STRM>>>( &
                UP_HAT,VP_HAT,WP_HAT,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1)

            CALL COMPUTE_MU_KRES<<<GRDD,BLKD,0,GM%STRM>>>(UP , VP , WP , UP_HAT,VP_HAT,WP_HAT,RHOP , GM%LES_FILTER_WIDTHD , GM%MU_DNSD , &
                                                GM%MUD , UU,VV,WW ,GM%KRESD,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,C_DEARDORFF)

            CALL COPY_GHOST_MU_KRES<<<(N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS+256-1)/256,256,0,GM%STRM>>> (GM%WALL_TOPO_AOS, &
                    GM%CELL_INDEXD, GM%CELL_SOLID, GM%CELL_EXTERIOR, N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS, &
                    GM%MUD, GM%KRESD, NULL_BOUNDARY, SOLID_BOUNDARY, OPEN_BOUNDARY, MIRROR_BOUNDARY, &
                    IBAR, JBAR, KBAR, IBP1, JBP1, KBP1, GM%CELL_WALL_INDEX,GM%CELL_LB,GM%CELL_UB,.FALSE.)
                          
            CALL ZERO_BOUNDARY_VISCOSITY_ACCUM_GPU<<<GRDD,BLKD,0,GM%STRM>>>( &
                GM%IWORK1D,GM%WORK7D,IBP1,JBP1,KBP1)

            CALL COMPUTE_BOUNDARY_VISCOSITY_ACCUM<<<(N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS+256-1)/256,256,0,GM%STRM>>> (GM%WALL_TOPO_AOS, &
                    GM%LES_FILTER_WIDTHD, GM%RDXD, GM%RDYD, GM%RDZD, UU, VV, WW, N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS, GM%MU_DNSD, &
                    GM%WORK7D, RHOP, GM%IWORK1D, IBAR, JBAR, KBAR)

            CALL FINALIZE_BOUNDARY_VISCOSITY<<<GRDD,BLKD,0,GM%STRM>>>(GM%WORK7D,GM%IWORK1D,GM%MUD,IBAR,JBAR,KBAR)

            CALL COPY_GHOST_MU_KRES<<<(N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS+256-1)/256,256,0,GM%STRM>>> (GM%WALL_TOPO_AOS, &
                    GM%CELL_INDEXD, GM%CELL_SOLID, GM%CELL_EXTERIOR, N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS, &
                    GM%MUD, GM%KRESD, NULL_BOUNDARY, SOLID_BOUNDARY, OPEN_BOUNDARY, MIRROR_BOUNDARY, &
                    IBAR, JBAR, KBAR, IBP1, JBP1, KBP1, GM%CELL_WALL_INDEX,GM%CELL_LB,GM%CELL_UB,.TRUE.)

            CALL COMPUTE_BOUNDARY_VISCOSITY_ADD_J<<<(JBP1+128)/128,128,0,GM%STRM>>>(IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,GM%MUD,GM%KRESD)
            CALL COMPUTE_BOUNDARY_VISCOSITY_ADD_I<<<(IBP1+128)/128,128,0,GM%STRM>>>(IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,GM%MUD,GM%KRESD)
            CALL COMPUTE_BOUNDARY_VISCOSITY_ADD_K<<<(KBP1+128)/128,128,0,GM%STRM>>>(IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,GM%MUD,GM%KRESD)


            T_USED(4) = T_USED(4) + CURRENT_TIME() - T_NOW
            IF (GPU_DEBUG_H2D_D2H) THEN
                MU     = GM%MUD
                MU_DNS = GM%MU_DNSD
                KRES   = GM%KRESD
            ENDIF

        END SUBROUTINE COMPUTE_VISCOSITY_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE GET_VISCOSITY_GPU(ZZPD,MU_DNSD,TMPD,CELL_INDEXD,CELL_SOLID, RSQ_MW_ZD,MU_RSQMW_ZD,&
                                                        N_TRACKED_SPECIES,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,I_MAX_TEMP_IN)

            INTEGER, VALUE :: N_TRACKED_SPECIES,I_MAX_TEMP_IN
            INTEGER, VALUE :: IBAR,JBAR,KBAR,IBP1,JBP1,KBP1
            REAL(EB) , DEVICE ,DIMENSION(*) :: RSQ_MW_ZD
            REAL(EB) , DEVICE ,DIMENSION(0:I_MAX_TEMP_IN,1:*) :: MU_RSQMW_ZD
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:KBP1,1:*) :: ZZPD
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*)    :: MU_DNSD, TMPD
            INTEGER , DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*)    :: CELL_INDEXD
            LOGICAL , DEVICE, DIMENSION(0:*)           :: CELL_SOLID

            INTEGER :: I,J,K,IC,N,ITMP
            REAL(EB) :: NUMER,DENOM

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            IC = CELL_INDEXD(I,J,K)
            IF (CELL_SOLID(IC)) RETURN

            ITMP = MAX(0, MIN(I_MAX_TEMP_IN, NINT(TMPD(I,J,K))))

            NUMER = 0._EB
            DENOM = 0._EB
            DO N = 1, N_TRACKED_SPECIES
                NUMER = NUMER + MU_RSQMW_ZD(ITMP,N) * ZZPD(I,J,K,N)
                DENOM = DENOM + ZZPD(I,J,K,N) * RSQ_MW_ZD(N)
            ENDDO

            MU_DNSD(I,J,K) = NUMER / DENOM

        END SUBROUTINE GET_VISCOSITY_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE ZERO_UVW_FILTER_WORK_GPU(UP,VP,WP,U_HAT,V_HAT,W_HAT,IBP1,JBP1,KBP1)
            INTEGER, VALUE :: IBP1,JBP1,KBP1
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: UP,VP,WP,U_HAT,V_HAT,W_HAT
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1)*BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1)*BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1)*BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I<0 .OR. I>IBP1 .OR. J<0 .OR. J>JBP1 .OR. K<0 .OR. K>KBP1) RETURN
            UP(I,J,K) = 0._EB
            VP(I,J,K) = 0._EB
            WP(I,J,K) = 0._EB
            U_HAT(I,J,K) = 0._EB
            V_HAT(I,J,K) = 0._EB
            W_HAT(I,J,K) = 0._EB
        END SUBROUTINE ZERO_UVW_FILTER_WORK_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE ZERO_BOUNDARY_VISCOSITY_ACCUM_GPU(IWORK,MU_SUM,IBP1,JBP1,KBP1)
            INTEGER, VALUE :: IBP1,JBP1,KBP1
            INTEGER, DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: IWORK
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: MU_SUM
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1)*BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1)*BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1)*BLOCKDIM%Z + THREADIDX%Z - 1
            IF (I<0 .OR. I>IBP1 .OR. J<0 .OR. J>JBP1 .OR. K<0 .OR. K>KBP1) RETURN
            IWORK(I,J,K) = 0
            MU_SUM(I,J,K) = 0._EB
        END SUBROUTINE ZERO_BOUNDARY_VISCOSITY_ACCUM_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE COMPUTE_UVWP_1(UPD,VPD,WPD,UUD,VVD,WWD,IBAR,JBAR,KBAR)

            INTEGER, VALUE :: IBAR,JBAR,KBAR
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: UPD,VPD,WPD
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: UUD,VVD,WWD
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            UPD(I,J,K) = 0.5_EB * ( UUD(I,J,K) + UUD(I-1,J,K) )
            VPD(I,J,K) = 0.5_EB * ( VVD(I,J,K) + VVD(I,J-1,K) )
            WPD(I,J,K) = 0.5_EB * ( WWD(I,J,K) + WWD(I,J,K-1) )

        END SUBROUTINE COMPUTE_UVWP_1
    
        ATTRIBUTES(GLOBAL) SUBROUTINE COMPUTE_UVWPD_2_UVW(UPD,VPD,WPD,U_GHOST,V_GHOST,W_GHOST,WALL_TOPO, &
            N_EXTERNAL_WALL_CELLS,OPEN_BOUNDARY,MIRROR_BOUNDARY,INTERPOLATED_BOUNDARY,IBP1,JBP1)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)

            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,OPEN_BOUNDARY,MIRROR_BOUNDARY,INTERPOLATED_BOUNDARY,IBP1,JBP1
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: UPD,VPD,WPD
            REAL(EB), DEVICE, DIMENSION(*) :: U_GHOST,V_GHOST,W_GHOST

            INTEGER :: IW, II, JJ, KK, IIG, JJG, KKG

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X

            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS) RETURN

            II  = WALL_TOPO(IW)%II
            JJ  = WALL_TOPO(IW)%JJ
            KK  = WALL_TOPO(IW)%KK
            IIG = WALL_TOPO(IW)%IIG
            JJG = WALL_TOPO(IW)%JJG
            KKG = WALL_TOPO(IW)%KKG

            IF (WALL_TOPO(IW)%BOUNDARY_TYPE == OPEN_BOUNDARY .OR. WALL_TOPO(IW)%BOUNDARY_TYPE == MIRROR_BOUNDARY) THEN
                UPD(II,JJ,KK) = UPD(IIG,JJG,KKG)
                VPD(II,JJ,KK) = VPD(IIG,JJG,KKG)
                WPD(II,JJ,KK) = WPD(IIG,JJG,KKG)
            ELSEIF (WALL_TOPO(IW)%BOUNDARY_TYPE == INTERPOLATED_BOUNDARY) THEN
                UPD(II,JJ,KK) = U_GHOST(IW)
                VPD(II,JJ,KK) = V_GHOST(IW)
                WPD(II,JJ,KK) = W_GHOST(IW)
            ENDIF

        END SUBROUTINE COMPUTE_UVWPD_2_UVW

        ATTRIBUTES(GLOBAL) SUBROUTINE GPU_FILL_EDGES(UP,VP,WP,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1)
            INTEGER , VALUE :: IBAR ,JBAR,KBAR , IBP1 ,JBP1, KBP1
            REAL(EB) ,DEVICE,DIMENSION(0:IBP1,0:JBP1,0:*):: UP,VP,WP
            INTEGER :: I,J,K,N

            N = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1

            IF (N >= 1 .AND. N <= IBAR) THEN
                I=N; J=0; K=0
                UP(I,J,K) = (UP(I,J+1,K)+UP(I,J,K+1))-UP(I,J+1,K+1)
                VP(I,J,K) = (VP(I,J+1,K)+VP(I,J,K+1))-VP(I,J+1,K+1)
                WP(I,J,K) = (WP(I,J+1,K)+WP(I,J,K+1))-WP(I,J+1,K+1)

                J=0; K=KBP1
                UP(I,J,K) = (UP(I,J+1,K)+UP(I,J,K-1))-UP(I,J+1,K-1)
                VP(I,J,K) = (VP(I,J+1,K)+VP(I,J,K-1))-VP(I,J+1,K-1)
                WP(I,J,K) = (WP(I,J+1,K)+WP(I,J,K-1))-WP(I,J+1,K-1)

                J=JBP1; K=0
                UP(I,J,K) = (UP(I,J-1,K)+UP(I,J,K+1))-UP(I,J-1,K+1)
                VP(I,J,K) = (VP(I,J-1,K)+VP(I,J,K+1))-VP(I,J-1,K+1)
                WP(I,J,K) = (WP(I,J-1,K)+WP(I,J,K+1))-WP(I,J-1,K+1)

                J=JBP1; K=KBP1
                UP(I,J,K) = (UP(I,J-1,K)+UP(I,J,K-1))-UP(I,J-1,K-1)
                VP(I,J,K) = (VP(I,J-1,K)+VP(I,J,K-1))-VP(I,J-1,K-1)
                WP(I,J,K) = (WP(I,J-1,K)+WP(I,J,K-1))-WP(I,J-1,K-1)
            ENDIF

            IF (N >= 1 .AND. N <= JBAR) THEN
                J=N; I=0; K=0
                UP(I,J,K) = (UP(I+1,J,K)+UP(I,J,K+1))-UP(I+1,J,K+1)
                VP(I,J,K) = (VP(I+1,J,K)+VP(I,J,K+1))-VP(I+1,J,K+1)
                WP(I,J,K) = (WP(I+1,J,K)+WP(I,J,K+1))-WP(I+1,J,K+1)

                I=0; K=KBP1
                UP(I,J,K) = (UP(I+1,J,K)+UP(I,J,K-1))-UP(I+1,J,K-1)
                VP(I,J,K) = (VP(I+1,J,K)+VP(I,J,K-1))-VP(I+1,J,K-1)
                WP(I,J,K) = (WP(I+1,J,K)+WP(I,J,K-1))-WP(I+1,J,K-1)

                I=IBP1; K=0
                UP(I,J,K) = (UP(I-1,J,K)+UP(I,J,K+1))-UP(I-1,J,K+1)
                VP(I,J,K) = (VP(I-1,J,K)+VP(I,J,K+1))-VP(I-1,J,K+1)
                WP(I,J,K) = (WP(I-1,J,K)+WP(I,J,K+1))-WP(I-1,J,K+1)

                I=IBP1; K=KBP1
                UP(I,J,K) = (UP(I-1,J,K)+UP(I,J,K-1))-UP(I-1,J,K-1)
                VP(I,J,K) = (VP(I-1,J,K)+VP(I,J,K-1))-VP(I-1,J,K-1)
                WP(I,J,K) = (WP(I-1,J,K)+WP(I,J,K-1))-WP(I-1,J,K-1)
            ENDIF

            IF (N >= 1 .AND. N <= KBAR) THEN
                K=N; I=0; J=0
                UP(I,J,K) = (UP(I+1,J,K)+UP(I,J+1,K))-UP(I+1,J+1,K)
                VP(I,J,K) = (VP(I+1,J,K)+VP(I,J+1,K))-VP(I+1,J+1,K)
                WP(I,J,K) = (WP(I+1,J,K)+WP(I,J+1,K))-WP(I+1,J+1,K)

                I=0; J=JBP1
                UP(I,J,K) = (UP(I+1,J,K)+UP(I,J-1,K))-UP(I+1,J-1,K)
                VP(I,J,K) = (VP(I+1,J,K)+VP(I,J-1,K))-VP(I+1,J-1,K)
                WP(I,J,K) = (WP(I+1,J,K)+WP(I,J-1,K))-WP(I+1,J-1,K)

                I=IBP1; J=0
                UP(I,J,K) = (UP(I-1,J,K)+UP(I,J+1,K))-UP(I-1,J+1,K)
                VP(I,J,K) = (VP(I-1,J,K)+VP(I,J+1,K))-VP(I-1,J+1,K)
                WP(I,J,K) = (WP(I-1,J,K)+WP(I,J+1,K))-WP(I-1,J+1,K)

                I=IBP1; J=JBP1
                UP(I,J,K) = (UP(I-1,J,K)+UP(I,J-1,K))-UP(I-1,J-1,K)
                VP(I,J,K) = (VP(I-1,J,K)+VP(I,J-1,K))-VP(I-1,J-1,K)
                WP(I,J,K) = (WP(I-1,J,K)+WP(I,J-1,K))-WP(I-1,J-1,K)
            ENDIF

        END SUBROUTINE GPU_FILL_EDGES

        ATTRIBUTES(GLOBAL) SUBROUTINE GPU_FILL_EDGES_CORNER(UP,VP,WP,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1)
            INTEGER , VALUE :: IBAR ,JBAR,KBAR , IBP1 ,JBP1, KBP1
            REAL(EB) ,DEVICE,DIMENSION(0:IBP1,0:JBP1,0:*):: UP,VP,WP

            UP(0,0,0)=2._EB*UP(1,1,1)-UP(2,2,2)
            VP(0,0,0)=2._EB*VP(1,1,1)-VP(2,2,2)
            WP(0,0,0)=2._EB*WP(1,1,1)-WP(2,2,2)
            UP(IBP1,0,0)=2._EB*UP(IBP1-1,1,1)-UP(IBP1-2,2,2)
            VP(IBP1,0,0)=2._EB*VP(IBP1-1,1,1)-VP(IBP1-2,2,2)
            WP(IBP1,0,0)=2._EB*WP(IBP1-1,1,1)-WP(IBP1-2,2,2)
            UP(0,JBP1,0)=2._EB*UP(1,JBP1-1,1)-UP(2,JBP1-2,2)
            VP(0,JBP1,0)=2._EB*VP(1,JBP1-1,1)-VP(2,JBP1-2,2)
            WP(0,JBP1,0)=2._EB*WP(1,JBP1-1,1)-WP(2,JBP1-2,2)
            UP(0,0,KBP1)=2._EB*UP(1,1,KBP1-1)-UP(2,2,KBP1-2)
            VP(0,0,KBP1)=2._EB*VP(1,1,KBP1-1)-VP(2,2,KBP1-2)
            WP(0,0,KBP1)=2._EB*WP(1,1,KBP1-1)-WP(2,2,KBP1-2)
            UP(IBP1,JBP1,0)=2._EB*UP(IBP1-1,JBP1-1,1)-UP(IBP1-2,JBP1-2,2)
            VP(IBP1,JBP1,0)=2._EB*VP(IBP1-1,JBP1-1,1)-VP(IBP1-2,JBP1-2,2)
            WP(IBP1,JBP1,0)=2._EB*WP(IBP1-1,JBP1-1,1)-WP(IBP1-2,JBP1-2,2)
            UP(IBP1,0,KBP1)=2._EB*UP(IBP1-1,1,KBP1-1)-UP(IBP1-2,2,KBP1-2)
            VP(IBP1,0,KBP1)=2._EB*VP(IBP1-1,1,KBP1-1)-VP(IBP1-2,2,KBP1-2)
            WP(IBP1,0,KBP1)=2._EB*WP(IBP1-1,1,KBP1-1)-WP(IBP1-2,2,KBP1-2)
            UP(0,JBP1,KBP1)=2._EB*UP(1,JBP1-1,KBP1-1)-UP(2,JBP1-2,KBP1-2)
            VP(0,JBP1,KBP1)=2._EB*VP(1,JBP1-1,KBP1-1)-VP(2,JBP1-2,KBP1-2)
            WP(0,JBP1,KBP1)=2._EB*WP(1,JBP1-1,KBP1-1)-WP(2,JBP1-2,KBP1-2)
            UP(IBP1,JBP1,KBP1)=2._EB*UP(IBP1-1,JBP1-1,KBP1-1)-UP(IBP1-2,JBP1-2,KBP1-2)
            VP(IBP1,JBP1,KBP1)=2._EB*VP(IBP1-1,JBP1-1,KBP1-1)-VP(IBP1-2,JBP1-2,KBP1-2)
            WP(IBP1,JBP1,KBP1)=2._EB*WP(IBP1-1,JBP1-1,KBP1-1)-WP(IBP1-2,JBP1-2,KBP1-2)

        END SUBROUTINE GPU_FILL_EDGES_CORNER

        ATTRIBUTES(GLOBAL) SUBROUTINE GPU_TEST_FILTER_1_UVW(K3DTD,U_HAT,V_HAT,W_HAT,UP,VP,WP, &
            IBAR,JBAR,KBAR,IBP1,JBP1,KBP1)
            INTEGER , VALUE :: IBAR ,JBAR,KBAR , IBP1 ,JBP1 ,KBP1
            REAL(EB) ,DEVICE:: K3DTD(-1:1,-1:1,-1:1)
            REAL(EB) ,DEVICE,DIMENSION(0:IBP1,0:JBP1,0:*):: U_HAT,V_HAT,W_HAT,UP,VP,WP
            INTEGER :: I,J,K,N,M,L
            REAL(EB) :: K3D_LOCAL,U_SUM,V_SUM,W_SUM

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1 
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            U_SUM = 0._EB
            V_SUM = 0._EB
            W_SUM = 0._EB
            DO N = -1,1
                DO M = -1,1
                    DO L = -1,1
                        K3D_LOCAL = K3DTD(L,M,N)
                        U_SUM = U_SUM + UP(I+L,J+M,K+N)*K3D_LOCAL
                        V_SUM = V_SUM + VP(I+L,J+M,K+N)*K3D_LOCAL
                        W_SUM = W_SUM + WP(I+L,J+M,K+N)*K3D_LOCAL
                    ENDDO
                ENDDO
            ENDDO
            U_HAT(I,J,K) = U_SUM
            V_HAT(I,J,K) = V_SUM
            W_HAT(I,J,K) = W_SUM

        END SUBROUTINE GPU_TEST_FILTER_1_UVW

        ATTRIBUTES(GLOBAL) SUBROUTINE GPU_TEST_FILTER_2_I_UVW(U_HAT,V_HAT,W_HAT,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1)
            INTEGER , VALUE :: IBAR ,JBAR,KBAR , IBP1 , JBP1 ,KBP1
            REAL(EB) ,DEVICE,DIMENSION(0:IBP1,0:JBP1,0:*):: U_HAT,V_HAT,W_HAT
            INTEGER :: J,K
            J = (blockIdx%x - 1) * blockDim%x + threadIdx%x - 1 
            K = (blockIdx%y - 1) * blockDim%y + threadIdx%y - 1

            if (J < 0 .or. J > JBP1) return
            if (K < 0 .or. K > KBP1) return

            U_HAT(0,J,K) = 2._EB*U_HAT(1,J,K) - U_HAT(2,J,K)
            V_HAT(0,J,K) = 2._EB*V_HAT(1,J,K) - V_HAT(2,J,K)
            W_HAT(0,J,K) = 2._EB*W_HAT(1,J,K) - W_HAT(2,J,K)
            U_HAT(IBP1,J,K) = 2._EB*U_HAT(IBP1-1,J,K) - U_HAT(IBP1-2,J,K)
            V_HAT(IBP1,J,K) = 2._EB*V_HAT(IBP1-1,J,K) - V_HAT(IBP1-2,J,K)
            W_HAT(IBP1,J,K) = 2._EB*W_HAT(IBP1-1,J,K) - W_HAT(IBP1-2,J,K)
        END SUBROUTINE GPU_TEST_FILTER_2_I_UVW

        ATTRIBUTES(GLOBAL) SUBROUTINE GPU_TEST_FILTER_2_J_UVW(U_HAT,V_HAT,W_HAT,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1)
            INTEGER , VALUE :: IBAR ,JBAR,KBAR , IBP1 , JBP1 ,KBP1
            REAL(EB) ,DEVICE,DIMENSION(0:IBP1,0:JBP1,0:*):: U_HAT,V_HAT,W_HAT
            INTEGER :: I,K
            I = (blockIdx%x - 1) * blockDim%x + threadIdx%x - 1 
            K = (blockIdx%y - 1) * blockDim%y + threadIdx%y - 1

            if (I < 0 .or. I > IBP1) return
            if (K < 0 .or. K > KBP1) return

            U_HAT(I,0,K) = 2._EB*U_HAT(I,1,K) - U_HAT(I,2,K)
            V_HAT(I,0,K) = 2._EB*V_HAT(I,1,K) - V_HAT(I,2,K)
            W_HAT(I,0,K) = 2._EB*W_HAT(I,1,K) - W_HAT(I,2,K)
            U_HAT(I,JBP1,K) = 2._EB*U_HAT(I,JBP1-1,K) - U_HAT(I,JBP1-2,K)
            V_HAT(I,JBP1,K) = 2._EB*V_HAT(I,JBP1-1,K) - V_HAT(I,JBP1-2,K)
            W_HAT(I,JBP1,K) = 2._EB*W_HAT(I,JBP1-1,K) - W_HAT(I,JBP1-2,K)
        END SUBROUTINE GPU_TEST_FILTER_2_J_UVW

        ATTRIBUTES(GLOBAL) SUBROUTINE GPU_TEST_FILTER_2_K_UVW(U_HAT,V_HAT,W_HAT,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1)
            INTEGER , VALUE :: IBAR ,JBAR,KBAR , IBP1 , JBP1 ,KBP1
            REAL(EB) ,DEVICE,DIMENSION(0:IBP1,0:JBP1,0:*):: U_HAT,V_HAT,W_HAT
            INTEGER :: I,J
            I = (blockIdx%x - 1) * blockDim%x + threadIdx%x - 1 
            J = (blockIdx%y - 1) * blockDim%y + threadIdx%y - 1

            if (I < 0 .or. I > IBP1) return
            if (J < 0 .or. J > JBP1) return

            U_HAT(I,J,0) = 2._EB*U_HAT(I,J,1) - U_HAT(I,J,2)
            V_HAT(I,J,0) = 2._EB*V_HAT(I,J,1) - V_HAT(I,J,2)
            W_HAT(I,J,0) = 2._EB*W_HAT(I,J,1) - W_HAT(I,J,2)
            U_HAT(I,J,KBP1) = 2._EB*U_HAT(I,J,KBP1-1) - U_HAT(I,J,KBP1-2)
            V_HAT(I,J,KBP1) = 2._EB*V_HAT(I,J,KBP1-1) - V_HAT(I,J,KBP1-2)
            W_HAT(I,J,KBP1) = 2._EB*W_HAT(I,J,KBP1-1) - W_HAT(I,J,KBP1-2)
        END SUBROUTINE GPU_TEST_FILTER_2_K_UVW

        ATTRIBUTES(GLOBAL) SUBROUTINE COMPUTE_MU_KRES(UPD , VPD , WPD , UP_HATD,VP_HATD,WP_HATD,RHOPD , LES_FILTER_WIDTHD , &
                                                       MU_DNSD , MUD , UUD,VVD,WWD ,KRESD,IBAR,JBAR,KBAR,IBP1,JBP1,KBP1, &
                                                       C_DEARDORFF_IN)
            INTEGER , VALUE :: IBAR ,JBAR,KBAR , IBP1 , JBP1, KBP1
            REAL(EB), VALUE :: C_DEARDORFF_IN
            REAL(EB) ,DEVICE,DIMENSION(0:IBP1,0:JBP1,0:*):: UPD,VPD,WPD,UP_HATD,VP_HATD,WP_HATD , RHOPD, LES_FILTER_WIDTHD , MU_DNSD , MUD , UUD,VVD,WWD , KRESD
            INTEGER :: I,J,K
            REAL(EB) :: DELTA,KSGS,NU_EDDY,U2,V2,W2

            i = (blockIdx%x - 1) * blockDim%x + threadIdx%x - 1
            j = (blockIdx%y - 1) * blockDim%y + threadIdx%y - 1
            k = (blockIdx%z - 1) * blockDim%z + threadIdx%z - 1

            if (i < 1 .or. i > IBAR) return
            if (j < 1 .or. j > JBAR) return
            if (k < 1 .or. k > KBAR) return

            DELTA = LES_FILTER_WIDTHD(I,J,K)
            KSGS = 0.5_EB*((UPD(I,J,K)-UP_HATD(I,J,K))**2 + &
                            (VPD(I,J,K)-VP_HATD(I,J,K))**2 + &
                            (WPD(I,J,K)-WP_HATD(I,J,K))**2)
            NU_EDDY = C_DEARDORFF_IN * DELTA * SQRT(KSGS)

            MUD(I,J,K) = MU_DNSD(I,J,K) + RHOPD(I,J,K)*NU_EDDY

            U2 = 0.25_EB*(UUD(I-1,J,K)+UUD(I,J,K))**2
            V2 = 0.25_EB*(VVD(I,J-1,K)+VVD(I,J,K))**2
            W2 = 0.25_EB*(WWD(I,J,K-1)+WWD(I,J,K))**2
            KRESD(I,J,K) = 0.5_EB*(U2+V2+W2)

        END SUBROUTINE COMPUTE_MU_KRES

        ATTRIBUTES(GLOBAL) SUBROUTINE COMPUTE_BOUNDARY_VISCOSITY_ADD_J(IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,MUD,KRESD)
            INTEGER , VALUE :: IBAR ,JBAR,KBAR , IBP1,JBP1,KBP1
            REAL(EB) ,DEVICE , DIMENSION(0:IBP1,0:JBP1,0:*):: MUD,KRESD
            INTEGER :: J

            J = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            IF (J < 0 .OR. J > JBP1) RETURN

            MUD(   0,J,   0) = MUD(   1,J,1)
            MUD(IBP1,J,   0) = MUD(IBAR,J,1)
            MUD(IBP1,J,KBP1) = MUD(IBAR,J,KBAR)
            MUD(   0,J,KBP1) = MUD(   1,J,KBAR)

            KRESD(   0,J,   0) = KRESD(   1,J,1)
            KRESD(IBP1,J,   0) = KRESD(IBAR,J,1)
            KRESD(IBP1,J,KBP1) = KRESD(IBAR,J,KBAR)
            KRESD(   0,J,KBP1) = KRESD(   1,J,KBAR)

        END SUBROUTINE COMPUTE_BOUNDARY_VISCOSITY_ADD_J

        ATTRIBUTES(GLOBAL) SUBROUTINE COMPUTE_BOUNDARY_VISCOSITY_ADD_I(IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,MUD,KRESD)
            INTEGER , VALUE :: IBAR ,JBAR,KBAR , IBP1,JBP1,KBP1
            REAL(EB) ,DEVICE , DIMENSION(0:IBP1,0:JBP1,0:*):: MUD,KRESD
            INTEGER :: I

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            IF (I < 0 .OR. I > IBP1) RETURN

            MUD(I,   0,   0) = MUD(I,   1,1)
            MUD(I,JBP1,   0) = MUD(I,JBAR,1)
            MUD(I,JBP1,KBP1) = MUD(I,JBAR,KBAR)
            MUD(I,   0,KBP1) = MUD(I,   1,KBAR)

            KRESD(I,   0,   0) = KRESD(I,   1,1)
            KRESD(I,JBP1,   0) = KRESD(I,JBAR,1)
            KRESD(I,JBP1,KBP1) = KRESD(I,JBAR,KBAR)
            KRESD(I,   0,KBP1) = KRESD(I,   1,KBAR)

        END SUBROUTINE COMPUTE_BOUNDARY_VISCOSITY_ADD_I

        ATTRIBUTES(GLOBAL) SUBROUTINE COMPUTE_BOUNDARY_VISCOSITY_ADD_K(IBAR,JBAR,KBAR,IBP1,JBP1,KBP1,MUD,KRESD)
            INTEGER , VALUE :: IBAR ,JBAR,KBAR , IBP1,JBP1,KBP1
            REAL(EB) ,DEVICE , DIMENSION(0:IBP1,0:JBP1,0:*):: MUD,KRESD
            INTEGER :: K

            K = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            IF (K < 0 .OR. K > KBP1) RETURN

            MUD(   0,   0,K) = MUD(   1,   1,K)
            MUD(IBP1,   0,K) = MUD(IBAR,   1,K)
            MUD(IBP1,JBP1,K) = MUD(IBAR,JBAR,K)
            MUD(   0,JBP1,K) = MUD(   1,JBAR,K)

            KRESD(   0,   0,K) = KRESD(   1,   1,K)
            KRESD(IBP1,   0,K) = KRESD(IBAR,   1,K)
            KRESD(IBP1,JBP1,K) = KRESD(IBAR,JBAR,K)
            KRESD(   0,JBP1,K) = KRESD(   1,JBAR,K)

        END SUBROUTINE COMPUTE_BOUNDARY_VISCOSITY_ADD_K

        ATTRIBUTES(GLOBAL) SUBROUTINE COMPUTE_BOUNDARY_VISCOSITY_ACCUM(WALL_TOPO, LES_FILTER_WIDTHD, RDXD, RDYD, RDZD, UUD, &
            VVD, WWD, N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS, MU_DNSD, MU_SUMD, RHOPD, IWORK1D, IBAR, JBAR, KBAR)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER , VALUE :: N_EXTERNAL_WALL_CELLS , N_INTERNAL_WALL_CELLS,IBAR,JBAR,KBAR
            REAL(EB) , DEVICE ,DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: LES_FILTER_WIDTHD
            REAL(EB) , DEVICE , DIMENSION(0:*) :: RDXD,RDYD,RDZD
            REAL(EB) ,DEVICE , DIMENSION(0:IBAR+1,0:JBAR+1,0:*):: UUD,VVD,WWD,RHOPD,MU_DNSD,MU_SUMD
            INTEGER , DEVICE , DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: IWORK1D

            INTEGER :: X , IIG ,JJG ,KKG , IERR
            REAL(EB) :: TEMPER
            INTEGER :: I,J,K,L
            REAL(EB) :: DUDX , DVDX , DWDX , DUDY , DVDY , DWDY , DUDZ , DVDZ , DWDZ , DELTA , NU_EDDY
            REAL(EB) :: A_IJ(3,3) , S_IJ(3,3),O_IJ(3,3),S2,O2,IV_SO,SD2,DENOM

            X = (blockIdx%x - 1) * blockDim%x + threadIdx%x

            IF (X < 1 .OR. X > N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS) RETURN
            IF (WALL_TOPO(X)%BOUNDARY_TYPE /= SOLID_BOUNDARY) RETURN

            IIG = WALL_TOPO(X)%IIG
            JJG = WALL_TOPO(X)%JJG
            KKG = WALL_TOPO(X)%KKG
            IF (IIG < 1 .OR. IIG > IBAR) RETURN
            IF (JJG < 1 .OR. JJG > JBAR) RETURN
            IF (KKG < 1 .OR. KKG > KBAR) RETURN

            DELTA = LES_FILTER_WIDTHD(IIG,JJG,KKG)

            DUDX = RDXD(IIG)*(UUD(IIG,JJG,KKG)-UUD(IIG-1,JJG,KKG))
            DVDY = RDYD(JJG)*(VVD(IIG,JJG,KKG)-VVD(IIG,JJG-1,KKG))
            DWDZ = RDZD(KKG)*(WWD(IIG,JJG,KKG)-WWD(IIG,JJG,KKG-1))
            DUDY = 0.25_EB*RDYD(JJG)*(UUD(IIG,JJG+1,KKG)-UUD(IIG,JJG-1,KKG)+UUD(IIG-1,JJG+1,KKG)-UUD(IIG-1,JJG-1,KKG))
            DUDZ = 0.25_EB*RDZD(KKG)*(UUD(IIG,JJG,KKG+1)-UUD(IIG,JJG,KKG-1)+UUD(IIG-1,JJG,KKG+1)-UUD(IIG-1,JJG,KKG-1))
            DVDX = 0.25_EB*RDXD(IIG)*(VVD(IIG+1,JJG,KKG)-VVD(IIG-1,JJG,KKG)+VVD(IIG+1,JJG-1,KKG)-VVD(IIG-1,JJG-1,KKG))
            DVDZ = 0.25_EB*RDZD(KKG)*(VVD(IIG,JJG,KKG+1)-VVD(IIG,JJG,KKG-1)+VVD(IIG,JJG-1,KKG+1)-VVD(IIG,JJG-1,KKG-1))
            DWDX = 0.25_EB*RDXD(IIG)*(WWD(IIG+1,JJG,KKG)-WWD(IIG-1,JJG,KKG)+WWD(IIG+1,JJG,KKG-1)-WWD(IIG-1,JJG,KKG-1))
            DWDY = 0.25_EB*RDYD(JJG)*(WWD(IIG,JJG+1,KKG)-WWD(IIG,JJG-1,KKG)+WWD(IIG,JJG+1,KKG-1)-WWD(IIG,JJG-1,KKG-1))

            A_IJ(1,1)=DUDX ; A_IJ(1,2)=DUDY ; A_IJ(1,3)=DUDZ
            A_IJ(2,1)=DVDX ; A_IJ(2,2)=DVDY ; A_IJ(2,3)=DVDZ
            A_IJ(3,1)=DWDX ; A_IJ(3,2)=DWDY ; A_IJ(3,3)=DWDZ

            DO J=1,3
                DO I=1,3
                    S_IJ(I,J) = 0.5_EB * ( A_IJ(I,J) + A_IJ(J,I) )
                    O_IJ(I,J) = 0.5_EB * ( A_IJ(I,J) - A_IJ(J,I) )
                ENDDO
            ENDDO

            S2 = 0._EB
            O2 = 0._EB
            DO J=1,3
                DO I=1,3
                    S2 = S2 + S_IJ(I,J)*S_IJ(I,J)
                    O2 = O2 + O_IJ(I,J)*O_IJ(I,J)
                ENDDO
            ENDDO

            IV_SO = 0._EB
            DO L=1,3
                DO K=1,3
                    DO J=1,3
                        DO I=1,3
                            IV_SO = IV_SO + S_IJ(I,K)*S_IJ(K,J)*O_IJ(J,L)*O_IJ(L,I)
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO

            SD2 = (1._EB/6._EB)*(S2*S2 + O2*O2) + (2._EB/3._EB)*S2*O2 + 2._EB*IV_SO
            IF (SD2 < 0._EB) SD2 = 0._EB

            DENOM = S2**2.5_EB + SD2**1.25_EB
            IF (DENOM > TWO_EPSILON_EB) THEN
                NU_EDDY = (0.60_EB*DELTA)**2 * SD2**1.5_EB / DENOM
            ELSE
                NU_EDDY = 0._EB
            ENDIF

            TEMPER = MU_DNSD(IIG,JJG,KKG) + RHOPD(IIG,JJG,KKG)*NU_EDDY

            IERR = atomicAdd(IWORK1D(IIG,JJG,KKG), 1)
            IERR = atomicAdd(MU_SUMD(IIG,JJG,KKG), TEMPER)

        END SUBROUTINE COMPUTE_BOUNDARY_VISCOSITY_ACCUM

        ATTRIBUTES(GLOBAL) SUBROUTINE FINALIZE_BOUNDARY_VISCOSITY(MU_SUMD,IWORK1D,MUD,IBAR,JBAR,KBAR)
            INTEGER, VALUE :: IBAR,JBAR,KBAR
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: MU_SUMD, MUD
            INTEGER , DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: IWORK1D

            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            IF (IWORK1D(I,J,K) > 0) THEN
                MUD(I,J,K) = MU_SUMD(I,J,K) / REAL(IWORK1D(I,J,K),EB)
            ENDIF
        END SUBROUTINE FINALIZE_BOUNDARY_VISCOSITY

        ATTRIBUTES(GLOBAL) SUBROUTINE COPY_GHOST_MU_KRES(WALL_TOPO, CELL_INDEXD, CELL_SOLID, CELL_EXTERIOR, &
            N_EXTERNAL_WALL_CELLS, N_INTERNAL_WALL_CELLS, MUD, KRESD, NULL_BOUNDARY, SOLID_BOUNDARY, OPEN_BOUNDARY, &
            MIRROR_BOUNDARY, IBAR, JBAR, KBAR, IBP1, JBP1, KBP1, CELL_WALL_INDEX,CELL_LB,CELL_UB,SOLID_GHOSTS)
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER , VALUE :: N_EXTERNAL_WALL_CELLS,N_INTERNAL_WALL_CELLS,NULL_BOUNDARY,SOLID_BOUNDARY,OPEN_BOUNDARY,MIRROR_BOUNDARY
            LOGICAL, VALUE :: SOLID_GHOSTS
            INTEGER , VALUE :: IBAR,JBAR,KBAR,IBP1,JBP1,KBP1
            INTEGER, VALUE :: CELL_LB,CELL_UB
            INTEGER, DEVICE :: CELL_WALL_INDEX(CELL_LB:CELL_UB,-3:3)
            INTEGER , DEVICE , DIMENSION(0:IBP1,0:JBP1,0:*) :: CELL_INDEXD
            LOGICAL , DEVICE , DIMENSION(0:*) :: CELL_SOLID,CELL_EXTERIOR
            REAL(EB) , DEVICE , DIMENSION(0:IBP1,0:JBP1,0:*) :: MUD,KRESD

            INTEGER :: X,II,JJ,KK,IIG,JJG,KKG,IC,BTYPE,ICG,IOR,IW
            LOGICAL :: AFTER_SOLID

            X = (blockIdx%x - 1) * blockDim%x + threadIdx%x

            IF (X < 1 .OR. X > N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS) RETURN
            BTYPE = WALL_TOPO(X)%BOUNDARY_TYPE
            IF (BTYPE == NULL_BOUNDARY) RETURN

            II  = WALL_TOPO(X)%II
            JJ  = WALL_TOPO(X)%JJ
            KK  = WALL_TOPO(X)%KK
            IIG = WALL_TOPO(X)%IIG
            JJG = WALL_TOPO(X)%JJG
            KKG = WALL_TOPO(X)%KKG

            IF (II < 0 .OR. II > IBP1) RETURN
            IF (JJ < 0 .OR. JJ > JBP1) RETURN
            IF (KK < 0 .OR. KK > KBP1) RETURN
            IF (IIG < 1 .OR. IIG > IBAR) RETURN
            IF (JJG < 1 .OR. JJG > JBAR) RETURN
            IF (KKG < 1 .OR. KKG > KBAR) RETURN

            IC = CELL_INDEXD(II,JJ,KK)
            IF (SOLID_GHOSTS) THEN
                IF (BTYPE == SOLID_BOUNDARY .AND. CELL_SOLID(IC)) MUD(II,JJ,KK) = MUD(IIG,JJG,KKG)
            ELSE
                IF (CELL_SOLID(IC) .OR. CELL_EXTERIOR(IC)) KRESD(II,JJ,KK) = KRESD(IIG,JJG,KKG)
            ENDIF
            IF (BTYPE == OPEN_BOUNDARY .OR. BTYPE == MIRROR_BOUNDARY) THEN
                ! Preserve whether the CPU WALL loop has already updated the adjacent gas cell.
                AFTER_SOLID = .FALSE.
                ICG = CELL_INDEXD(IIG,JJG,KKG)
                DO IOR=-3,3
                    IW = CELL_WALL_INDEX(ICG,IOR)
                    IF (IW<1 .OR. IW>=X) CYCLE
                    IF (WALL_TOPO(IW)%BOUNDARY_TYPE==SOLID_BOUNDARY) AFTER_SOLID = .TRUE.
                ENDDO
                IF (SOLID_GHOSTS .EQV. AFTER_SOLID) MUD(II,JJ,KK) = MUD(IIG,JJG,KKG)
            ENDIF

        END SUBROUTINE COPY_GHOST_MU_KRES
        !--------------------------------------------------------------------------
        ! Momentum-flux host orchestration and kernels
        !--------------------------------------------------------------------------
        SUBROUTINE VELOCITY_FLUX_GPU(T,DT,NM,APPLY_TO_ESTIMATED_VARIABLES)
            USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
            USE PHYSICAL_FUNCTIONS, ONLY: COMPUTE_WIND_COMPONENTS
            REAL(EB) :: T_NOW
            INTEGER, INTENT(IN) :: NM
            REAL(EB), INTENT(IN) :: T,DT
            LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
            REAL(EB), DEVICE , POINTER, DIMENSION(:,:,:) :: TXY,TXZ,TYZ,OMX,OMY,OMZ,UU,VV,WW,RHOP,DP
            TYPE(DIM3) :: BLKD,GRDD
            INTEGER :: I
            REAL(EB) :: GX_H(0:IBAR_MAX),GY_H(0:IBAR_MAX),GZ_H(0:IBAR_MAX)

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)

            BLKD = GM%BLKD_3D_BP1
            GRDD = GM%GRDD_3D_BP1

            T_NOW=CURRENT_TIME()

            IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
                UU => GM%USD
                VV => GM%VSD
                WW => GM%WSD
                DP => GM%DS_D
                RHOP => GM%RHOSD

                IF (GPU_DEBUG_H2D_D2H) THEN
                    GM%USD   = US
                    GM%VSD   = VS
                    GM%WSD   = WS
                    GM%RHOSD = RHOS
                    GM%DS_D  = DS
                    GM%MUD   = MU
                    GM%FVXD  = FVX
                    GM%FVYD  = FVY
                    GM%FVZD  = FVZ
                    GM%RHO_0D = RHO_0
                ENDIF
            ELSE
                UU => GM%UD
                VV => GM%VD
                WW => GM%WD
                DP => GM%D_D
                RHOP => GM%RHOD

                IF (GPU_DEBUG_H2D_D2H) THEN
                    GM%UD    = U
                    GM%VD    = V
                    GM%WD    = W
                    GM%RHOD  = RHO
                    GM%D_D   = D
                    GM%MUD   = MU
                    GM%FVXD  = FVX
                    GM%FVYD  = FVY
                    GM%FVZD  = FVZ
                    GM%RHO_0D = RHO_0
                ENDIF
            ENDIF

            IF (GPU_DEBUG_H2D_D2H) THEN
                IF (.NOT.SPATIAL_GRAVITY_VARIATION) THEN
                    GX_H(0:IBAR) = EVALUATE_RAMP(T,I_RAMP_GX)*GVEC(1)
                    GY_H(0:IBAR) = EVALUATE_RAMP(T,I_RAMP_GY)*GVEC(2)
                    GZ_H(0:IBAR) = EVALUATE_RAMP(T,I_RAMP_GZ)*GVEC(3)
                ELSE
                    DO I=0,IBAR
                        GX_H(I) = EVALUATE_RAMP(X(I),I_RAMP_GX)*GVEC(1)
                        GY_H(I) = EVALUATE_RAMP(X(I),I_RAMP_GY)*GVEC(2)
                        GZ_H(I) = EVALUATE_RAMP(X(I),I_RAMP_GZ)*GVEC(3)
                    ENDDO
                ENDIF
                GM%GXD(0:IBAR) = GX_H(0:IBAR)
                GM%GYD(0:IBAR) = GY_H(0:IBAR)
                GM%GZD(0:IBAR) = GZ_H(0:IBAR)
            ENDIF

            TXY => GM%WORK1D
            TXZ => GM%WORK2D
            TYZ => GM%WORK3D
            OMX => GM%WORK4D
            OMY => GM%WORK5D
            OMZ => GM%WORK6D

            IF (GPU_DEBUG_H2D_D2H) CALL SYNC_EDGE_OMEGA_TAU_TO_GPU(NM)

            CALL VELOCITY_FLUX_KERNEL_1<<<GRDD,BLKD,0,GM%STRM>>>(OMX,OMY,OMZ,TXY,TXZ,TYZ,UU,VV,WW,GM%RDXND,GM%RDYND,GM%RDZND,GM%MUD,IBAR,JBAR,KBAR)
            CALL VELOCITY_FLUX_KERNEL_2X<<<GRDD,BLKD,0,GM%STRM>>>(OMX,OMY,OMZ,TXY,TXZ,TYZ,UU,VV,WW,RHOP,DP,GM%CELL_INDEXD,GM%CELL_EDGE_INDEX,GM%EDGE_OMEGA,GM%EDGE_TAU, &
                                                        GM%RDXND,GM%RDYD,GM%RDZD,GM%MUD,GM%GXD,IBAR,JBAR,KBAR,GM%RHO_0D,GM%FVXD,GM%WORK10D, &
                                                        GM%CELL_LB,GM%CELL_UB,GM%EDGE_LB,GM%EDGE_UB)
            CALL VELOCITY_FLUX_KERNEL_2Y<<<GRDD,BLKD,0,GM%STRM>>>(OMX,OMY,OMZ,TXY,TXZ,TYZ,UU,VV,WW,RHOP,DP,GM%CELL_INDEXD,GM%CELL_EDGE_INDEX,GM%EDGE_OMEGA,GM%EDGE_TAU, &
                                                        GM%RDXD,GM%RDYND,GM%RDZD,GM%MUD,GM%GYD,IBAR,JBAR,KBAR,GM%RHO_0D,GM%FVYD,GM%WORK11D, &
                                                        GM%CELL_LB,GM%CELL_UB,GM%EDGE_LB,GM%EDGE_UB)
            CALL VELOCITY_FLUX_KERNEL_2Z<<<GRDD,BLKD,0,GM%STRM>>>(OMX,OMY,OMZ,TXY,TXZ,TYZ,UU,VV,WW,RHOP,DP,GM%CELL_INDEXD,GM%CELL_EDGE_INDEX,GM%EDGE_OMEGA,GM%EDGE_TAU, &
                                                        GM%RDXD,GM%RDYD,GM%RDZND,GM%MUD,GM%GZD,IBAR,JBAR,KBAR,GM%RHO_0D,GM%FVZD,GM%WORK12D, &
                                                        GM%CELL_LB,GM%CELL_UB,GM%EDGE_LB,GM%EDGE_UB)

            IF (OPEN_WIND_BOUNDARY) THEN
                CALL COMPUTE_WIND_COMPONENTS(T,NM)
                GM%U_WINDD = U_WIND
                GM%V_WINDD = V_WIND
                GM%W_WINDD = W_WIND
            ENDIF
            IF (ANY(ABS(FVEC)>TWO_EPSILON_EB) .OR. CTRL_DIRECT_FORCE .OR. ANY(ABS(OVEC)>TWO_EPSILON_EB)) &
                CALL WIND_FORCES_GPU(T,UU,VV,WW,RHOP)

            IF (GPU_DEBUG_H2D_D2H) THEN
                CALL GPU_WAIT(GM%STRM)
                FVX = GM%FVXD
                FVY = GM%FVYD
                FVZ = GM%FVZD
            ENDIF
            T_USED(4) = T_USED(4) + CURRENT_TIME() - T_NOW

        END SUBROUTINE VELOCITY_FLUX_GPU

        SUBROUTINE WIND_FORCES_GPU(T,UU,VV,WW,RHOP)
            USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP
            USE CONTROL_VARIABLES, ONLY: CONTROL,N_CTRL
            REAL(EB), INTENT(IN) :: T
            REAL(EB), DEVICE :: UU(0:IBP1,0:JBP1,0:KBP1),VV(0:IBP1,0:JBP1,0:KBP1), &
                WW(0:IBP1,0:JBP1,0:KBP1),RHOP(0:IBP1,0:JBP1,0:KBP1)
            REAL(EB) :: FORCE(3),DIRECTION_FACTOR(3),THETA
            INTEGER :: N,RAMP_INDEX(3)
            LOGICAL :: ROTATING

            IF (CTRL_DIRECT_FORCE) THEN
                DO N=1,N_CTRL
                    WHERE (CONTROL(N)%CONTROL_FORCE) FVEC = FVEC-CONTROL(N)%INSTANT_VALUE
                ENDDO
            ENDIF
            DIRECTION_FACTOR = 1._EB
            IF (I_RAMP_DIRECTION_T/=0) THEN
                THETA = EVALUATE_RAMP(T,I_RAMP_DIRECTION_T)*DEG2RAD
                DIRECTION_FACTOR(1) = -SIN(THETA)
                DIRECTION_FACTOR(2) = -COS(THETA)
            ENDIF
            RAMP_INDEX = (/I_RAMP_FVX_T,I_RAMP_FVY_T,I_RAMP_FVZ_T/)
            FORCE = 0._EB
            DO N=1,3
                IF (ABS(FVEC(N))<=TWO_EPSILON_EB) CYCLE
                IF (RAMP_INDEX(N)<=0) RAMP_INDEX(N) = MAX(0,I_RAMP_PGF_T)
                FORCE(N) = FVEC(N)*EVALUATE_RAMP(T,RAMP_INDEX(N))*DIRECTION_FACTOR(N)
            ENDDO

            ROTATING = ANY(ABS(OVEC)>TWO_EPSILON_EB)
            IF (ROTATING) THEN
                IF (GPU_DEBUG_H2D_D2H) THEN
                    GM%U_GHOSTD = U_GHOST
                    GM%V_GHOSTD = V_GHOST
                    GM%W_GHOSTD = W_GHOST
                ENDIF
                GM%WORK7D = 0._EB
                GM%WORK8D = 0._EB
                GM%WORK9D = 0._EB
                CALL COMPUTE_UVWP_1<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1,0,GM%STRM>>>( &
                    GM%WORK7D,GM%WORK8D,GM%WORK9D,UU,VV,WW,IBAR,JBAR,KBAR)
                IF (N_EXTERNAL_WALL_CELLS>0) CALL WIND_FORCE_GHOST_GPU<<<GM%GRDD_EXTERNAL_WALL_1D,GPU_BLKD_1D,0,GM%STRM>>>( &
                    GM%WORK7D,GM%WORK8D,GM%WORK9D,GM%U_GHOSTD,GM%V_GHOSTD,GM%W_GHOSTD,GM%WALL_TOPO_AOS, &
                    N_EXTERNAL_WALL_CELLS,IBP1,JBP1)
            ENDIF
            CALL WIND_FORCE_KERNEL_GPU<<<GM%GRDD_3D_BP1,GM%BLKD_3D_BP1,0,GM%STRM>>>( &
                GM%FVXD,GM%FVYD,GM%FVZD,RHOP,GM%WORK7D,GM%WORK8D,GM%WORK9D, &
                FORCE(1),FORCE(2),FORCE(3),OVEC(1),OVEC(2),OVEC(3),ROTATING,IBAR,JBAR,KBAR)
        END SUBROUTINE WIND_FORCES_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE WIND_FORCE_GHOST_GPU(UP,VP,WP,U_GHOST,V_GHOST,W_GHOST,WALL_TOPO,NW,IBP1,JBP1)
            INTEGER, VALUE :: NW,IBP1,JBP1
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            REAL(EB), DEVICE :: UP(0:IBP1,0:JBP1,0:*),VP(0:IBP1,0:JBP1,0:*),WP(0:IBP1,0:JBP1,0:*)
            REAL(EB), DEVICE :: U_GHOST(1:*),V_GHOST(1:*),W_GHOST(1:*)
            INTEGER :: IW,I,J,K
            IW = (BLOCKIDX%X-1)*BLOCKDIM%X+THREADIDX%X
            IF (IW>NW) RETURN
            I = WALL_TOPO(IW)%II
            J = WALL_TOPO(IW)%JJ
            K = WALL_TOPO(IW)%KK
            UP(I,J,K) = U_GHOST(IW)
            VP(I,J,K) = V_GHOST(IW)
            WP(I,J,K) = W_GHOST(IW)
        END SUBROUTINE WIND_FORCE_GHOST_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE WIND_FORCE_KERNEL_GPU(FVX,FVY,FVZ,RHOP,UP,VP,WP,FX,FY,FZ,OX,OY,OZ,ROTATING,IBAR,JBAR,KBAR)
            INTEGER, VALUE :: IBAR,JBAR,KBAR
            LOGICAL, VALUE :: ROTATING
            REAL(EB), VALUE :: FX,FY,FZ,OX,OY,OZ
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: FVX,FVY,FVZ,RHOP,UP,VP,WP
            INTEGER :: I,J,K
            I = (BLOCKIDX%X-1)*BLOCKDIM%X+THREADIDX%X-1
            J = (BLOCKIDX%Y-1)*BLOCKDIM%Y+THREADIDX%Y-1
            K = (BLOCKIDX%Z-1)*BLOCKDIM%Z+THREADIDX%Z-1
            IF (I>IBAR .OR. J>JBAR .OR. K>KBAR) RETURN
            IF (J>=1 .AND. K>=1) THEN
                FVX(I,J,K) = FVX(I,J,K)-2._EB*FX/(RHOP(I,J,K)+RHOP(I+1,J,K))
                IF (ROTATING) FVX(I,J,K) = FVX(I,J,K)+OY*(WP(I,J,K)+WP(I+1,J,K))-OZ*(VP(I,J,K)+VP(I+1,J,K))
            ENDIF
            IF (I>=1 .AND. K>=1) THEN
                FVY(I,J,K) = FVY(I,J,K)-2._EB*FY/(RHOP(I,J,K)+RHOP(I,J+1,K))
                IF (ROTATING) FVY(I,J,K) = FVY(I,J,K)+OZ*(UP(I,J,K)+UP(I,J+1,K))-OX*(WP(I,J,K)+WP(I,J+1,K))
            ENDIF
            IF (I>=1 .AND. J>=1) THEN
                FVZ(I,J,K) = FVZ(I,J,K)-2._EB*FZ/(RHOP(I,J,K)+RHOP(I,J,K+1))
                IF (ROTATING) FVZ(I,J,K) = FVZ(I,J,K)+OX*(VP(I,J,K)+VP(I,J,K+1))-OY*(UP(I,J,K)+UP(I,J,K+1))
            ENDIF
        END SUBROUTINE WIND_FORCE_KERNEL_GPU

        SUBROUTINE SYNC_EDGE_OMEGA_TAU_TO_GPU(NM)
            INTEGER, INTENT(IN) :: NM
            INTEGER :: EDGE_LB, EDGE_UB, IE
            REAL(EB), ALLOCATABLE :: EDGE_OMEGA_H(:,:), EDGE_TAU_H(:,:)

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)
            IF (.NOT.ALLOCATED(GM%EDGE_OMEGA)) RETURN
            IF (.NOT.ALLOCATED(GM%EDGE_TAU))   RETURN

            EDGE_LB = LBOUND(EDGE,1)
            EDGE_UB = EDGE_COUNT(NM)
            IF (EDGE_UB < EDGE_LB) RETURN

            ALLOCATE(EDGE_OMEGA_H(EDGE_LB:EDGE_UB,-2:2))
            ALLOCATE(EDGE_TAU_H  (EDGE_LB:EDGE_UB,-2:2))

            DO IE = EDGE_LB, EDGE_UB
                EDGE_OMEGA_H(IE,-2:2) = EDGE(IE)%OMEGA(-2:2)
                EDGE_TAU_H  (IE,-2:2) = EDGE(IE)%TAU  (-2:2)
            ENDDO

            DO IE=-2,2
                GM%EDGE_OMEGA(EDGE_LB:EDGE_UB,IE) = EDGE_OMEGA_H(:,IE)
                GM%EDGE_TAU  (EDGE_LB:EDGE_UB,IE) = EDGE_TAU_H(:,IE)
            ENDDO

            DEALLOCATE(EDGE_OMEGA_H,EDGE_TAU_H)

        END SUBROUTINE SYNC_EDGE_OMEGA_TAU_TO_GPU

        SUBROUTINE SYNC_EDGE_VELOCITY_AVG_D2H(NM)
            INTEGER, INTENT(IN) :: NM
            INTEGER :: EDGE_LB,EDGE_UB,IE
            REAL(EB), ALLOCATABLE :: EDGE_U_AVG_H(:),EDGE_V_AVG_H(:),EDGE_W_AVG_H(:)

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)
            IF (.NOT.ALLOCATED(GM%EDGE_U_AVG)) RETURN
            IF (.NOT.ALLOCATED(GM%EDGE_V_AVG)) RETURN
            IF (.NOT.ALLOCATED(GM%EDGE_W_AVG)) RETURN

            EDGE_LB = LBOUND(EDGE,1)
            EDGE_UB = EDGE_COUNT(NM)
            IF (EDGE_UB < EDGE_LB) RETURN

            CALL GPU_WAIT(GM%STRM)

            ALLOCATE(EDGE_U_AVG_H(EDGE_LB:EDGE_UB))
            ALLOCATE(EDGE_V_AVG_H(EDGE_LB:EDGE_UB))
            ALLOCATE(EDGE_W_AVG_H(EDGE_LB:EDGE_UB))
            EDGE_U_AVG_H = GM%EDGE_U_AVG(EDGE_LB:EDGE_UB)
            EDGE_V_AVG_H = GM%EDGE_V_AVG(EDGE_LB:EDGE_UB)
            EDGE_W_AVG_H = GM%EDGE_W_AVG(EDGE_LB:EDGE_UB)

            DO IE=EDGE_LB,EDGE_UB
                EDGE(IE)%U_AVG = EDGE_U_AVG_H(IE)
                EDGE(IE)%V_AVG = EDGE_V_AVG_H(IE)
                EDGE(IE)%W_AVG = EDGE_W_AVG_H(IE)
            ENDDO

            DEALLOCATE(EDGE_U_AVG_H,EDGE_V_AVG_H,EDGE_W_AVG_H)

        END SUBROUTINE SYNC_EDGE_VELOCITY_AVG_D2H

        ATTRIBUTES(GLOBAL) SUBROUTINE VELOCITY_FLUX_KERNEL_1(OMX,OMY,OMZ,TXY,TXZ,TYZ,UU,VV,WW,RDXN,RDYN,RDZN,MU ,IBAR,JBAR,KBAR)
            IMPLICIT NONE
            INTEGER , VALUE :: IBAR , JBAR ,KBAR
            REAL(EB) , DEVICE , DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: OMX , OMY , OMZ , TXY , TXZ , TYZ , UU , VV , WW , MU
            REAL(EB) , DEVICE , DIMENSION(0:*) :: RDXN , RDYN , RDZN
            REAL(EB) ::DWDX, DUDY , DVDX , DUDZ , DVDZ , DWDY , MUX ,MUY ,MUZ
            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBAR) RETURN
            IF (J < 0 .OR. J > JBAR) RETURN
            IF (K < 0 .OR. K > KBAR) RETURN

            DUDY = RDYN(J)*(UU(I,J+1,K)-UU(I,J,K))
            DVDX = RDXN(I)*(VV(I+1,J,K)-VV(I,J,K))
            DUDZ = RDZN(K)*(UU(I,J,K+1)-UU(I,J,K))
            DWDX = RDXN(I)*(WW(I+1,J,K)-WW(I,J,K))
            DVDZ = RDZN(K)*(VV(I,J,K+1)-VV(I,J,K))
            DWDY = RDYN(J)*(WW(I,J+1,K)-WW(I,J,K))
            OMX(I,J,K) = DWDY - DVDZ
            OMY(I,J,K) = DUDZ - DWDX
            OMZ(I,J,K) = DVDX - DUDY
            MUX = 0.25_EB*(MU(I,J+1,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I,J+1,K+1))
            MUY = 0.25_EB*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I+1,J,K+1))
            MUZ = 0.25_EB*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J+1,K)+MU(I+1,J+1,K))
            TXY(I,J,K) = MUZ*(DVDX + DUDY)
            TXZ(I,J,K) = MUY*(DUDZ + DWDX)
            TYZ(I,J,K) = MUX*(DVDZ + DWDY)

        END SUBROUTINE VELOCITY_FLUX_KERNEL_1

        ATTRIBUTES(GLOBAL) SUBROUTINE VELOCITY_FLUX_KERNEL_2X(OMX,OMY,OMZ,TXY,TXZ,TYZ,UU,VV,WW,RHOP,DP,CELL_INDEXD,CELL_EDGE_INDEX,EDGE_OMEGA,EDGE_TAU, &
                                                              RDXN,RDY,RDZ,MU,GX,IBAR,JBAR,KBAR,RHO_0,FVX,RRHO_X, &
                                                              CELL_LB,CELL_UB,EDGE_LB,EDGE_UB)
            INTEGER , VALUE :: IBAR,JBAR,KBAR,CELL_LB,CELL_UB,EDGE_LB,EDGE_UB
            REAL(EB) , DEVICE , DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: OMX , OMY , OMZ , TXY , TXZ , TYZ , UU , VV , WW , MU , RHOP , DP ,FVX,RRHO_X
            INTEGER , DEVICE , DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: CELL_INDEXD
            INTEGER , DEVICE , DIMENSION(CELL_LB:CELL_UB,1:*) :: CELL_EDGE_INDEX
            REAL(EB), DEVICE, DIMENSION(EDGE_LB:EDGE_UB,-2:*) :: EDGE_OMEGA,EDGE_TAU
            REAL(EB) , DEVICE , DIMENSION(0:*) :: GX , RHO_0 , RDXN , RDY , RDZ
            REAL(EB) :: WP , WM ,VP ,VM ,OMYP , OMYM , OMZP , OMZM , TXZP ,TXZM , TXYP , TXYM,WOMY , VOMZ , RRHO , DVDY , DWDZ , TXXP , TXXM , DTXXDX , DTXYDY , DTXZDZ , VTRM , FOTH
            REAL(EB) :: EDGE_VAL
            INTEGER :: I,J,K,IC,IEYP,IEYM,IEZP,IEZM

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 0 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN

            FOTH = 4._EB / 3._EB

            WP    = WW(I,J,K)   + WW(I+1,J,K)
            WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
            VP    = VV(I,J,K)   + VV(I+1,J,K)
            VM    = VV(I,J-1,K) + VV(I+1,J-1,K)
            OMYP  = OMY(I,J,K)
            OMYM  = OMY(I,J,K-1)
            OMZP  = OMZ(I,J,K)
            OMZM  = OMZ(I,J-1,K)
            TXZP  = TXZ(I,J,K)
            TXZM  = TXZ(I,J,K-1)
            TXYP  = TXY(I,J,K)
            TXYM  = TXY(I,J-1,K)
            IC    = CELL_INDEXD(I,J,K)
            IEYP  = CELL_EDGE_INDEX(IC,8)
            IEYM  = CELL_EDGE_INDEX(IC,6)
            IEZP  = CELL_EDGE_INDEX(IC,12)
            IEZM  = CELL_EDGE_INDEX(IC,10)

            EDGE_VAL = EDGE_OMEGA(IEYP,-1)
            IF (EDGE_VAL > -1.E5_EB) THEN
                OMYP = EDGE_VAL
                TXZP = EDGE_TAU(IEYP,-1)
            ENDIF

            EDGE_VAL = EDGE_OMEGA(IEYM,1)
            IF (EDGE_VAL > -1.E5_EB) THEN
                OMYM = EDGE_VAL
                TXZM = EDGE_TAU(IEYM,1)
            ENDIF

            EDGE_VAL = EDGE_OMEGA(IEZP,-2)
            IF (EDGE_VAL > -1.E5_EB) THEN
                OMZP = EDGE_VAL
                TXYP = EDGE_TAU(IEZP,-2)
            ENDIF

            EDGE_VAL = EDGE_OMEGA(IEZM,2)
            IF (EDGE_VAL > -1.E5_EB) THEN
                OMZM = EDGE_VAL
                TXYM = EDGE_TAU(IEZM,2)
            ENDIF

            WOMY  = WP*OMYP + WM*OMYM
            VOMZ  = VP*OMZP + VM*OMZM
            RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
            RRHO_X(I,J,K) = RRHO
            DVDY  = (VV(I+1,J,K)-VV(I+1,J-1,K))*RDY(J)
            DWDZ  = (WW(I+1,J,K)-WW(I+1,J,K-1))*RDZ(K)
            TXXP  = MU(I+1,J,K)*( FOTH*DP(I+1,J,K) - 2._EB*(DVDY+DWDZ) )
            DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
            DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
            TXXM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DVDY+DWDZ) )
            DTXXDX= RDXN(I)*(TXXP-TXXM)
            DTXYDY= RDY(J) *(TXYP-TXYM)
            DTXZDZ= RDZ(K) *(TXZP-TXZM)
            VTRM  = DTXXDX + DTXYDY + DTXZDZ
            FVX(I,J,K) = 0.25_EB*(WOMY - VOMZ) - GX(I) + RRHO*(GX(I)*RHO_0(K) - VTRM)
        END SUBROUTINE VELOCITY_FLUX_KERNEL_2X

        ATTRIBUTES(GLOBAL) SUBROUTINE VELOCITY_FLUX_KERNEL_2Y(OMX,OMY,OMZ,TXY,TXZ,TYZ,UU,VV,WW,RHOP,DP,CELL_INDEXD,CELL_EDGE_INDEX,EDGE_OMEGA,EDGE_TAU, &
                                                              RDX,RDYN,RDZ,MU,GY,IBAR,JBAR,KBAR,RHO_0,FVY,RRHO_Y, &
                                                              CELL_LB,CELL_UB,EDGE_LB,EDGE_UB)
                IMPLICIT NONE
            INTEGER , VALUE :: IBAR,JBAR,KBAR,CELL_LB,CELL_UB,EDGE_LB,EDGE_UB
            REAL(EB) , DEVICE , DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: OMX , OMY , OMZ , TXY , TXZ , TYZ , UU , VV , WW , MU , RHOP , DP ,FVY,RRHO_Y
            INTEGER , DEVICE , DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: CELL_INDEXD
            INTEGER , DEVICE , DIMENSION(CELL_LB:CELL_UB,1:*) :: CELL_EDGE_INDEX
            REAL(EB), DEVICE, DIMENSION(EDGE_LB:EDGE_UB,-2:*) :: EDGE_OMEGA,EDGE_TAU
            REAL(EB) , DEVICE , DIMENSION(0:*) :: GY , RHO_0 , RDX , RDYN , RDZ
            REAL(EB) :: UP , UM ,WP ,WM ,OMXP , OMXM , OMZP , OMZM , TYZP ,TYZM , TXYP , TXYM,WOMX , UOMZ , RRHO , DUDX , DWDZ , TYYP , TYYM , DTXYDX , DTYYDY , DTYZDZ , VTRM,FOTH
            REAL(EB) :: EDGE_VAL
            INTEGER :: I,J,K,IC,IEXP,IEXM,IEZP,IEZM

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 0 .OR. J > JBAR) RETURN
            IF (K < 1 .OR. K > KBAR) RETURN
            FOTH = 4._EB / 3._EB

            UP    = UU(I,J,K)   + UU(I,J+1,K)
            UM    = UU(I-1,J,K) + UU(I-1,J+1,K)
            WP    = WW(I,J,K)   + WW(I,J+1,K)
            WM    = WW(I,J,K-1) + WW(I,J+1,K-1)
            OMXP  = OMX(I,J,K)
            OMXM  = OMX(I,J,K-1)
            OMZP  = OMZ(I,J,K)
            OMZM  = OMZ(I-1,J,K)
            TYZP  = TYZ(I,J,K)
            TYZM  = TYZ(I,J,K-1)
            TXYP  = TXY(I,J,K)
            TXYM  = TXY(I-1,J,K)
            IC    = CELL_INDEXD(I,J,K)
            IEXP  = CELL_EDGE_INDEX(IC,4)
            IEXM  = CELL_EDGE_INDEX(IC,2)
            IEZP  = CELL_EDGE_INDEX(IC,12)
            IEZM  = CELL_EDGE_INDEX(IC,11)

            EDGE_VAL = EDGE_OMEGA(IEXP,-2)
            IF (EDGE_VAL > -1.E5_EB) THEN
                OMXP = EDGE_VAL
                TYZP = EDGE_TAU(IEXP,-2)
            ENDIF

            EDGE_VAL = EDGE_OMEGA(IEXM,2)
            IF (EDGE_VAL > -1.E5_EB) THEN
                OMXM = EDGE_VAL
                TYZM = EDGE_TAU(IEXM,2)
            ENDIF

            EDGE_VAL = EDGE_OMEGA(IEZP,-1)
            IF (EDGE_VAL > -1.E5_EB) THEN
                OMZP = EDGE_VAL
                TXYP = EDGE_TAU(IEZP,-1)
            ENDIF

            EDGE_VAL = EDGE_OMEGA(IEZM,1)
            IF (EDGE_VAL > -1.E5_EB) THEN
                OMZM = EDGE_VAL
                TXYM = EDGE_TAU(IEZM,1)
            ENDIF

            WOMX  = WP*OMXP + WM*OMXM
            UOMZ  = UP*OMZP + UM*OMZM
            RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J+1,K))
            RRHO_Y(I,J,K) = RRHO
            DUDX  = (UU(I,J+1,K)-UU(I-1,J+1,K))*RDX(I)
            DWDZ  = (WW(I,J+1,K)-WW(I,J+1,K-1))*RDZ(K)
            TYYP  = MU(I,J+1,K)*( FOTH*DP(I,J+1,K) - 2._EB*(DUDX+DWDZ) )
            DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
            DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
            TYYM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DWDZ) )
            DTXYDX= RDX(I) *(TXYP-TXYM)
            DTYYDY= RDYN(J)*(TYYP-TYYM)
            DTYZDZ= RDZ(K) *(TYZP-TYZM)
            VTRM  = DTXYDX + DTYYDY + DTYZDZ
            FVY(I,J,K) = 0.25_EB*(UOMZ - WOMX) - GY(I) + RRHO*(GY(I)*RHO_0(K) - VTRM)
        END SUBROUTINE VELOCITY_FLUX_KERNEL_2Y

        ATTRIBUTES(GLOBAL) SUBROUTINE VELOCITY_FLUX_KERNEL_2Z(OMX,OMY,OMZ,TXY,TXZ,TYZ,UU,VV,WW,RHOP,DP,CELL_INDEXD,CELL_EDGE_INDEX,EDGE_OMEGA,EDGE_TAU, &
                                                              RDX,RDY,RDZN,MU,GZ,IBAR,JBAR,KBAR,RHO_0,FVZ,RRHO_Z, &
                                                              CELL_LB,CELL_UB,EDGE_LB,EDGE_UB)
                IMPLICIT NONE
            INTEGER , VALUE :: IBAR,JBAR,KBAR,CELL_LB,CELL_UB,EDGE_LB,EDGE_UB
            REAL(EB) , DEVICE , DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: OMX , OMY , OMZ , TXY , TXZ , TYZ , UU , VV , WW , MU , RHOP , DP ,FVZ,RRHO_Z
            INTEGER , DEVICE , DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: CELL_INDEXD
            INTEGER , DEVICE , DIMENSION(CELL_LB:CELL_UB,1:*) :: CELL_EDGE_INDEX
            REAL(EB), DEVICE, DIMENSION(EDGE_LB:EDGE_UB,-2:*) :: EDGE_OMEGA,EDGE_TAU
            REAL(EB) , DEVICE , DIMENSION(0:*) :: GZ , RHO_0 , RDX , RDY , RDZN
            REAL(EB) :: UP , UM ,VP ,VM ,OMYP , OMYM , OMXP , OMXM , TXZP ,TXZM , TYZP , TYZM , UOMY , VOMX , RRHO , DUDX , DVDY , TZZP , TZZM , DTXZDX , DTYZDY , DTZZDZ , VTRM , FOTH
            REAL(EB) :: EDGE_VAL
            INTEGER :: I,J,K,IC,IEXP,IEXM,IEYP,IEYM

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I < 1 .OR. I > IBAR) RETURN
            IF (J < 1 .OR. J > JBAR) RETURN
            IF (K < 0 .OR. K > KBAR) RETURN
            FOTH = 4._EB / 3._EB

            UP    = UU(I,J,K)   + UU(I,J,K+1)
            UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
            VP    = VV(I,J,K)   + VV(I,J,K+1)
            VM    = VV(I,J-1,K) + VV(I,J-1,K+1)
            OMYP  = OMY(I,J,K)
            OMYM  = OMY(I-1,J,K)
            OMXP  = OMX(I,J,K)
            OMXM  = OMX(I,J-1,K)
            TXZP  = TXZ(I,J,K)
            TXZM  = TXZ(I-1,J,K)
            TYZP  = TYZ(I,J,K)
            TYZM  = TYZ(I,J-1,K)
            IC    = CELL_INDEXD(I,J,K)
            IEXP  = CELL_EDGE_INDEX(IC,4)
            IEXM  = CELL_EDGE_INDEX(IC,3)
            IEYP  = CELL_EDGE_INDEX(IC,8)
            IEYM  = CELL_EDGE_INDEX(IC,7)

            EDGE_VAL = EDGE_OMEGA(IEXP,-1)
            IF (EDGE_VAL > -1.E5_EB) THEN
                OMXP = EDGE_VAL
                TYZP = EDGE_TAU(IEXP,-1)
            ENDIF

            EDGE_VAL = EDGE_OMEGA(IEXM,1)
            IF (EDGE_VAL > -1.E5_EB) THEN
                OMXM = EDGE_VAL
                TYZM = EDGE_TAU(IEXM,1)
            ENDIF

            EDGE_VAL = EDGE_OMEGA(IEYP,-2)
            IF (EDGE_VAL > -1.E5_EB) THEN
                OMYP = EDGE_VAL
                TXZP = EDGE_TAU(IEYP,-2)
            ENDIF

            EDGE_VAL = EDGE_OMEGA(IEYM,2)
            IF (EDGE_VAL > -1.E5_EB) THEN
                OMYM = EDGE_VAL
                TXZM = EDGE_TAU(IEYM,2)
            ENDIF

            UOMY  = UP*OMYP + UM*OMYM
            VOMX  = VP*OMXP + VM*OMXM
            RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
            RRHO_Z(I,J,K) = RRHO
            DUDX  = (UU(I,J,K+1)-UU(I-1,J,K+1))*RDX(I)
            DVDY  = (VV(I,J,K+1)-VV(I,J-1,K+1))*RDY(J)
            TZZP  = MU(I,J,K+1)*( FOTH*DP(I,J,K+1) - 2._EB*(DUDX+DVDY) )
            DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
            DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
            TZZM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DVDY) )
            DTXZDX= RDX(I) *(TXZP-TXZM)
            DTYZDY= RDY(J) *(TYZP-TYZM)
            DTZZDZ= RDZN(K)*(TZZP-TZZM)
            VTRM  = DTXZDX + DTYZDY + DTZZDZ
            FVZ(I,J,K) = 0.25_EB*(VOMX - UOMY) - GZ(I) + RRHO*(GZ(I)*0.5_EB*(RHO_0(K)+RHO_0(K+1)) - VTRM)
        END SUBROUTINE VELOCITY_FLUX_KERNEL_2Z
        !--------------------------------------------------------------------------
        ! Velocity predictor, stability reduction, and corrector
        !--------------------------------------------------------------------------
        SUBROUTINE VELOCITY_PREDICTOR_GPU(T,DT,DT_NEW,NM)
            USE TURBULENCE, ONLY: COMPRESSION_WAVE
            USE MANUFACTURED_SOLUTIONS, ONLY: UF_MMS,WF_MMS,VD2D_MMS_U,VD2D_MMS_V
            USE CC_SCALARS, ONLY : CC_PROJECT_VELOCITY

            REAL(EB) :: T_NOW
            INTEGER, INTENT(IN) :: NM
            REAL(EB), INTENT(IN) :: T,DT
            REAL(EB) :: DT_NEW(NMESHES)

            INTEGER , PARAMETER :: BS = 8
            TYPE(DIM3) :: BLKD , GRDD

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)
            BLKD = GM%BLKD_3D_BP1
            GRDD = GM%GRDD_3D_BP1

            IF (SOLID_PHASE_ONLY) RETURN

            T_NOW=CURRENT_TIME()
            IF (GPU_DEBUG_H2D_D2H) THEN
                GM%USD   = US
                GM%VSD   = VS
                GM%WSD   = WS
                GM%UD    = U
                GM%VD    = V
                GM%WD    = W
                GM%FVXD  = FVX
                GM%FVYD  = FVY
                GM%FVZD  = FVZ
                GM%HD    = H
                GM%DS_D  = DS
                GM%MUD   = MU
                GM%RHOSD = RHOS
                IF (CHECK_VN) GM%D_Z_MAXD = D_Z_MAX
            ENDIF
            IF (FREEZE_VELOCITY) THEN
                GM%USD = GM%UD
                GM%VSD = GM%VD
                GM%WSD = GM%WD
            ELSE
                CALL UVW_PREDICTOR_KERNEL<<<GRDD,BLKD,0,GM%STRM>>>( &
                    GM%USD,GM%VSD,GM%WSD,GM%UD,GM%VD,GM%WD,GM%FVXD,GM%FVYD,GM%FVZD, &
                    GM%RDXND,GM%RDYND,GM%RDZND,GM%HD,DT,IBAR,JBAR,KBAR)
            ENDIF

            IF (GPU_DEBUG_H2D_D2H) THEN
                CALL GPU_WAIT(GM%STRM)
                US = GM%USD
                VS = GM%VSD
                WS = GM%WSD
            ENDIF


            T_USED(4)=T_USED(4)+CURRENT_TIME()-T_NOW

            ! Check the stability criteria, and if the time step is too small, send back a signal to kill the job

            CALL CHECK_STABILITY_GPU(DT,DT_NEW,T,NM)

            IF (DT_NEW(NM)<DT_INITIAL*LIMITING_DT_RATIO .AND. (T+DT_NEW(NM)<(T_END-TWO_EPSILON_EB))) STOP_STATUS = INSTABILITY_STOP

            CONTAINS

            SUBROUTINE CHECK_STABILITY_GPU(DT,DT_NEW,T,NM)

                USE CC_SCALARS, ONLY : CHECK_CFLVN_LINKED_CELLS
                USE OUTPUT_CLOCKS, ONLY: RAMP_TIME_INDEX,RAMP_DT_INDEX
                USE MATH_FUNCTIONS, ONLY: EVALUATE_RAMP

                INTEGER, INTENT(IN) :: NM
                REAL(EB), INTENT(IN) :: DT,T
                REAL(EB) :: UODX,VODY,WODZ,UVW,UVWMAX,R_DX2,MU_MAX,MUTRM,PART_CFL,MU_TMP, UVWMAX_TMP, DT_CLIP, T_NOW
                REAL(EB) :: DT_NEW(NMESHES)
                INTEGER  :: I,J,K,IW,IIG,JJG,KKG, ICFL_TMP, JCFL_TMP, KCFL_TMP
                REAL(EB), PARAMETER :: DT_EPS = 1.E-10_EB
                TYPE(WALL_TYPE), POINTER :: WC
                TYPE(BOUNDARY_COORD_TYPE), POINTER :: BC
                TYPE(BOUNDARY_PROP1_TYPE), POINTER :: B1
                TYPE(RAMPS_TYPE), POINTER :: RP

                ! GPU reduction workspace.
                INTEGER, PARAMETER :: BS1 = 256
                INTEGER :: NCELL,NBLK
                TYPE(DIM3) :: BLOCK1, GRID1

                T_NOW = CURRENT_TIME()

                UVWMAX = 0._EB
                VN     = 0._EB
                MUTRM  = 1.E-9_EB
                R_DX2  = 1.E-9_EB
                ICFL   = 0; JCFL   = 0; KCFL   = 0
                I_VN   = 0; J_VN   = 0; K_VN   = 0

                !--------------------------------------------------------------
                ! CFL maximum reduction
                !--------------------------------------------------------------
                NCELL  = IBAR*JBAR*KBAR
                NBLK   = (NCELL + BS1 - 1) / BS1
                BLOCK1 = DIM3(BS1,1,1)
                GRID1  = DIM3(NBLK,1,1)

                IF (.NOT.ALLOCATED(GM%CFL_BLOCKMAX_D)) THEN
                    ALLOCATE(GM%CFL_BLOCKMAX_D(NBLK),GM%CFL_BLOCKMAX_H(NBLK))
                    ALLOCATE(GM%CFL_BLOCKI_D(NBLK),GM%CFL_BLOCKI_H(NBLK))
                    ALLOCATE(GM%CFL_BLOCKJ_D(NBLK),GM%CFL_BLOCKJ_H(NBLK))
                    ALLOCATE(GM%CFL_BLOCKK_D(NBLK),GM%CFL_BLOCKK_H(NBLK))
                ELSEIF (SIZE(GM%CFL_BLOCKMAX_D) < NBLK) THEN
                    DEALLOCATE(GM%CFL_BLOCKMAX_D,GM%CFL_BLOCKMAX_H, &
                               GM%CFL_BLOCKI_D,GM%CFL_BLOCKI_H, &
                               GM%CFL_BLOCKJ_D,GM%CFL_BLOCKJ_H, &
                               GM%CFL_BLOCKK_D,GM%CFL_BLOCKK_H)
                    ALLOCATE(GM%CFL_BLOCKMAX_D(NBLK),GM%CFL_BLOCKMAX_H(NBLK))
                    ALLOCATE(GM%CFL_BLOCKI_D(NBLK),GM%CFL_BLOCKI_H(NBLK))
                    ALLOCATE(GM%CFL_BLOCKJ_D(NBLK),GM%CFL_BLOCKJ_H(NBLK))
                    ALLOCATE(GM%CFL_BLOCKK_D(NBLK),GM%CFL_BLOCKK_H(NBLK))
                ENDIF

                CALL CFL_SCAN_BLOCKMAX_KERNEL<<<GRID1,BLOCK1,0,GM%STRM>>>( &
                    GM%CFL_BLOCKMAX_D, GM%CFL_BLOCKI_D, GM%CFL_BLOCKJ_D, GM%CFL_BLOCKK_D, &
                    GM%USD, GM%VSD, GM%WSD, GM%DS_D, GM%RDXD, GM%RDYD, GM%RDZD, GM%CELL_INDEXD, GM%CELL_SOLID, IBAR, JBAR, KBAR, CFL_VELOCITY_NORM)

                IF (.NOT.ALLOCATED(GM%FINAL_BLOCK_RESULT_D)) ALLOCATE(GM%FINAL_BLOCK_RESULT_D(4))

                CALL FINAL_BLOCKMAX_REDUCE_KERNEL<<<1,BS1,0,GM%STRM>>>( &
                    GM%FINAL_BLOCK_RESULT_D, &
                    GM%CFL_BLOCKMAX_D,GM%CFL_BLOCKI_D,GM%CFL_BLOCKJ_D,GM%CFL_BLOCKK_D,NBLK)

                CALL GPU_WAIT(GM%STRM)
                GM%FINAL_BLOCK_RESULT_H = GM%FINAL_BLOCK_RESULT_D
                UVWMAX = GM%FINAL_BLOCK_RESULT_H(1)
                ICFL = NINT(GM%FINAL_BLOCK_RESULT_H(2))
                JCFL = NINT(GM%FINAL_BLOCK_RESULT_H(3))
                KCFL = NINT(GM%FINAL_BLOCK_RESULT_H(4))

                CFL = DT*UVWMAX

                ! Include surface vegetation drag if necessary
                PART_CFL = 0._EB

                !--------------------------------------------------------------
                ! VN maximum reduction
                !--------------------------------------------------------------
                PARABOLIC_IF: IF (CHECK_VN) THEN

                    IF (.NOT.ALLOCATED(GM%VN_BLOCKMAX_D)) THEN
                        ALLOCATE(GM%VN_BLOCKMAX_D(NBLK),GM%VN_BLOCKMAX_H(NBLK))
                        ALLOCATE(GM%VN_BLOCKI_D(NBLK),GM%VN_BLOCKI_H(NBLK))
                        ALLOCATE(GM%VN_BLOCKJ_D(NBLK),GM%VN_BLOCKJ_H(NBLK))
                        ALLOCATE(GM%VN_BLOCKK_D(NBLK),GM%VN_BLOCKK_H(NBLK))
                    ELSEIF (SIZE(GM%VN_BLOCKMAX_D) < NBLK) THEN
                        DEALLOCATE(GM%VN_BLOCKMAX_D,GM%VN_BLOCKMAX_H, &
                                   GM%VN_BLOCKI_D,GM%VN_BLOCKI_H, &
                                   GM%VN_BLOCKJ_D,GM%VN_BLOCKJ_H, &
                                   GM%VN_BLOCKK_D,GM%VN_BLOCKK_H)
                        ALLOCATE(GM%VN_BLOCKMAX_D(NBLK),GM%VN_BLOCKMAX_H(NBLK))
                        ALLOCATE(GM%VN_BLOCKI_D(NBLK),GM%VN_BLOCKI_H(NBLK))
                        ALLOCATE(GM%VN_BLOCKJ_D(NBLK),GM%VN_BLOCKJ_H(NBLK))
                        ALLOCATE(GM%VN_BLOCKK_D(NBLK),GM%VN_BLOCKK_H(NBLK))
                    ENDIF

                    CALL VN_SCAN_BLOCKMAX_KERNEL<<<GRID1,BLOCK1,0,GM%STRM>>>( &
                        GM%VN_BLOCKMAX_D, GM%VN_BLOCKI_D, GM%VN_BLOCKJ_D, GM%VN_BLOCKK_D, &
                        GM%D_Z_MAXD, GM%MUD, GM%RHOSD, GM%CELL_INDEXD, GM%CELL_SOLID, IBAR, JBAR, KBAR)

                    CALL FINAL_BLOCKMAX_REDUCE_KERNEL<<<1,BS1,0,GM%STRM>>>( &
                        GM%FINAL_BLOCK_RESULT_D, &
                        GM%VN_BLOCKMAX_D,GM%VN_BLOCKI_D,GM%VN_BLOCKJ_D,GM%VN_BLOCKK_D,NBLK)

                    CALL GPU_WAIT(GM%STRM)
                    GM%FINAL_BLOCK_RESULT_H = GM%FINAL_BLOCK_RESULT_D
                    MU_MAX = GM%FINAL_BLOCK_RESULT_H(1)
                    I_VN = NINT(GM%FINAL_BLOCK_RESULT_H(2))
                    J_VN = NINT(GM%FINAL_BLOCK_RESULT_H(3))
                    K_VN = NINT(GM%FINAL_BLOCK_RESULT_H(4))

                    IF (TWO_D) THEN
                        R_DX2 = RDX(I_VN)**2 + RDZ(K_VN)**2
                    ELSE
                        R_DX2 = RDX(I_VN)**2 + RDY(J_VN)**2 + RDZ(K_VN)**2
                    ENDIF

                    MUTRM = MU_MAX
                    VN = DT*2._EB*R_DX2*MUTRM

                ENDIF PARABOLIC_IF

                ! Attempt DT restriction to avoid clippings
                DT_CLIP = HUGE(1._EB)

                RAMP_TIME_IF: IF (RAMP_TIME_INDEX>0) THEN

                    ! Keep the RAMP_TIME_INDEX-controlled time step unchanged.

                ELSE RAMP_TIME_IF

                    ! Adjust time step size if necessary
                    IF ((CFL<CFL_MAX .AND. VN<VN_MAX .AND. PART_CFL<PARTICLE_CFL_MAX) .OR. LOCK_TIME_STEP) THEN
                        DT_NEW(NM) = DT
                        IF (CFL<=CFL_MIN .AND. VN<VN_MIN .AND. PART_CFL<PARTICLE_CFL_MIN .AND. .NOT.LOCK_TIME_STEP) THEN
                            SELECT CASE (RESTRICT_TIME_STEP)
                            CASE (.TRUE.)
                                DT_NEW(NM) = MIN(1.1_EB*DT,DT_INITIAL)
                            CASE (.FALSE.)
                                DT_NEW(NM) = 1.1_EB*DT
                            END SELECT
                            CHANGE_TIME_STEP_INDEX(NM) = 1
                        ELSE
                            CHANGE_TIME_STEP_INDEX(NM) = 0
                        ENDIF
                    ELSE
                        DT_NEW(NM) = 0.9_EB*MIN( CFL_MAX/MAX(UVWMAX,DT_EPS)               , &
                                                VN_MAX/(2._EB*R_DX2*MAX(MUTRM,DT_EPS))   , &
                                                PARTICLE_CFL_MAX/MAX(PART_UVWMAX,DT_EPS) , &
                                                DT_CLIP)
                        CHANGE_TIME_STEP_INDEX(NM) = -1
                    ENDIF

                    IF (RAMP_DT_INDEX > 0) DT_NEW(NM) = MIN(DT_NEW(NM),EVALUATE_RAMP(T,RAMP_DT_INDEX))

                ENDIF RAMP_TIME_IF

                T_USED(4)=T_USED(4)+CURRENT_TIME()-T_NOW

            END SUBROUTINE CHECK_STABILITY_GPU

        END SUBROUTINE VELOCITY_PREDICTOR_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE CFL_SCAN_BLOCKMAX_KERNEL( &
            BLOCKMAX, BLOCKI, BLOCKJ, BLOCKK, &
            US, VS, WS, DS, RDX, RDY, RDZ, CELL_INDEX, CELL_SOLID, IBAR, JBAR, KBAR, CFL_NORM)

            IMPLICIT NONE

            INTEGER, VALUE :: IBAR, JBAR, KBAR, CFL_NORM
            REAL(EB), DEVICE, DIMENSION(*)         :: BLOCKMAX
            INTEGER,  DEVICE, DIMENSION(*)         :: BLOCKI, BLOCKJ, BLOCKK
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*)  :: US, VS, WS, DS
            REAL(EB), DEVICE, DIMENSION(0:*)        :: RDX, RDY, RDZ
            INTEGER,  DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*)  :: CELL_INDEX
            LOGICAL,  DEVICE, DIMENSION(0:*)         :: CELL_SOLID

            REAL(EB), SHARED :: SMAX(256)
            INTEGER,  SHARED :: SI(256), SJ(256), SK(256)

            INTEGER :: TID, GID, NCELL
            INTEGER :: I, J, K, TMP, REM, IC
            REAL(EB) :: UODX, VODY, WODZ, UVW, VAL

            TID   = THREADIDX%X
            GID   = (BLOCKIDX%X-1) * BLOCKDIM%X + THREADIDX%X
            NCELL = IBAR*JBAR*KBAR

            VAL = -1._EB
            I = 0
            J = 0
            K = 0

            IF (GID <= NCELL) THEN
                TMP = GID - 1
                K   = TMP / (IBAR*JBAR) + 1
                REM = MOD(TMP, IBAR*JBAR)
                J   = REM / IBAR + 1
                I   = MOD(REM, IBAR) + 1

                IC = CELL_INDEX(I,J,K)

                IF (.NOT. CELL_SOLID(IC)) THEN
                    UODX = MAX(ABS(US(I-1,J,K)), ABS(US(I,J,K))) * RDX(I)
                    VODY = MAX(ABS(VS(I,J-1,K)), ABS(VS(I,J,K))) * RDY(J)
                    WODZ = MAX(ABS(WS(I,J,K-1)), ABS(WS(I,J,K))) * RDZ(K)

                    SELECT CASE (CFL_NORM)
                        CASE (0)
                            UVW = MAX(UODX,VODY,WODZ) + ABS(DS(I,J,K))
                        CASE (1)
                            UVW = UODX + VODY + WODZ + ABS(DS(I,J,K))
                        CASE (2)
                            UVW = SQRT(UODX*UODX + VODY*VODY + WODZ*WODZ) + ABS(DS(I,J,K))
                        CASE (3)
                            UVW = MAX(UODX,VODY,WODZ)
                        CASE DEFAULT
                            UVW = SQRT(UODX*UODX + VODY*VODY + WODZ*WODZ) + ABS(DS(I,J,K))
                    END SELECT
                    VAL = UVW
                ENDIF
            ENDIF

            SMAX(TID) = VAL
            SI(TID)   = I
            SJ(TID)   = J
            SK(TID)   = K
            CALL SYNCTHREADS()

            IF (TID <= 128) THEN
                IF (SMAX(TID+128) > SMAX(TID)) THEN
                    SMAX(TID) = SMAX(TID+128)
                    SI(TID) = SI(TID+128)
                    SJ(TID) = SJ(TID+128)
                    SK(TID) = SK(TID+128)
                ENDIF
            ENDIF
            CALL SYNCTHREADS()

            IF (TID <= 64) THEN
                IF (SMAX(TID+64) > SMAX(TID)) THEN
                    SMAX(TID) = SMAX(TID+64)
                    SI(TID) = SI(TID+64)
                    SJ(TID) = SJ(TID+64)
                    SK(TID) = SK(TID+64)
                ENDIF
            ENDIF
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
                BLOCKI(BLOCKIDX%X)   = SI(1)
                BLOCKJ(BLOCKIDX%X)   = SJ(1)
                BLOCKK(BLOCKIDX%X)   = SK(1)
            ENDIF

        END SUBROUTINE CFL_SCAN_BLOCKMAX_KERNEL

        ATTRIBUTES(GLOBAL) SUBROUTINE VN_SCAN_BLOCKMAX_KERNEL( &
            BLOCKMAX, BLOCKI, BLOCKJ, BLOCKK, &
            D_Z_MAX, MU, RHOS, CELL_INDEX, CELL_SOLID, IBAR, JBAR, KBAR)

            IMPLICIT NONE

            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), DEVICE, DIMENSION(*)         :: BLOCKMAX
            INTEGER,  DEVICE, DIMENSION(*)         :: BLOCKI, BLOCKJ, BLOCKK
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*)  :: D_Z_MAX, MU, RHOS
            INTEGER,  DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*)  :: CELL_INDEX
            LOGICAL,  DEVICE, DIMENSION(0:*)         :: CELL_SOLID

            REAL(EB), SHARED :: SMAX(256)
            INTEGER,  SHARED :: SI(256), SJ(256), SK(256)

            INTEGER :: TID, GID, NCELL
            INTEGER :: I, J, K, TMP, REM, IC
            REAL(EB) :: MU_TMP, VAL
            REAL(EB), PARAMETER :: RHOS_EPS = 1.E-12_EB

            TID   = THREADIDX%X
            GID   = (BLOCKIDX%X-1) * BLOCKDIM%X + THREADIDX%X
            NCELL = IBAR*JBAR*KBAR

            VAL = -1._EB
            I = 0
            J = 0
            K = 0

            IF (GID <= NCELL) THEN
                TMP = GID - 1
                K   = TMP / (IBAR*JBAR) + 1
                REM = MOD(TMP, IBAR*JBAR)
                J   = REM / IBAR + 1
                I   = MOD(REM, IBAR) + 1

                IC = CELL_INDEX(I,J,K)

                IF (.NOT. CELL_SOLID(IC)) THEN
                    MU_TMP = MAX(D_Z_MAX(I,J,K), MU(I,J,K)/MAX(RHOS(I,J,K), RHOS_EPS))
                    VAL = MU_TMP
                ENDIF
            ENDIF

            SMAX(TID) = VAL
            SI(TID)   = I
            SJ(TID)   = J
            SK(TID)   = K
            CALL SYNCTHREADS()

            IF (TID <= 128) THEN
                IF (SMAX(TID+128) > SMAX(TID)) THEN
                    SMAX(TID) = SMAX(TID+128)
                    SI(TID) = SI(TID+128)
                    SJ(TID) = SJ(TID+128)
                    SK(TID) = SK(TID+128)
                ENDIF
            ENDIF
            CALL SYNCTHREADS()

            IF (TID <= 64) THEN
                IF (SMAX(TID+64) > SMAX(TID)) THEN
                    SMAX(TID) = SMAX(TID+64)
                    SI(TID) = SI(TID+64)
                    SJ(TID) = SJ(TID+64)
                    SK(TID) = SK(TID+64)
                ENDIF
            ENDIF
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
                BLOCKI(BLOCKIDX%X)   = SI(1)
                BLOCKJ(BLOCKIDX%X)   = SJ(1)
                BLOCKK(BLOCKIDX%X)   = SK(1)
            ENDIF

        END SUBROUTINE VN_SCAN_BLOCKMAX_KERNEL

        ATTRIBUTES(GLOBAL) SUBROUTINE UVW_PREDICTOR_KERNEL(US,VS,WS,U,V,W,FVX,FVY,FVZ,RDXN,RDYN,RDZN,H,DT, &
            IBAR,JBAR,KBAR)

            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: US,VS,WS,U,V,W,FVX,FVY,FVZ,H
            REAL(EB), DEVICE, DIMENSION(0:*) :: RDXN,RDYN,RDZN
            REAL(EB), VALUE :: DT

            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I>=0 .AND. I<=IBAR .AND. J>=1 .AND. J<=JBAR .AND. K>=1 .AND. K<=KBAR) &
                US(I,J,K) = U(I,J,K) - DT*(FVX(I,J,K) + RDXN(I)*(H(I+1,J,K)-H(I,J,K)))

            IF (I>=1 .AND. I<=IBAR .AND. J>=0 .AND. J<=JBAR .AND. K>=1 .AND. K<=KBAR) &
                VS(I,J,K) = V(I,J,K) - DT*(FVY(I,J,K) + RDYN(J)*(H(I,J+1,K)-H(I,J,K)))

            IF (I>=1 .AND. I<=IBAR .AND. J>=1 .AND. J<=JBAR .AND. K>=0 .AND. K<=KBAR) &
                WS(I,J,K) = W(I,J,K) - DT*(FVZ(I,J,K) + RDZN(K)*(H(I,J,K+1)-H(I,J,K)))

        END SUBROUTINE UVW_PREDICTOR_KERNEL

        SUBROUTINE VELOCITY_CORRECTOR_GPU(T,DT,NM)

            USE TURBULENCE, ONLY: COMPRESSION_WAVE
            USE MANUFACTURED_SOLUTIONS, ONLY: UF_MMS,WF_MMS,VD2D_MMS_U,VD2D_MMS_V
            USE CC_SCALARS, ONLY : CC_PROJECT_VELOCITY

            REAL(EB) :: T_NOW,XHAT,ZHAT
            INTEGER, INTENT(IN) :: NM
            REAL(EB), INTENT(IN) :: T,DT
            INTEGER , PARAMETER :: BS = 8
            TYPE(DIM3) :: BLKD , GRDD

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)
            BLKD = GM%BLKD_3D_BP1
            GRDD = GM%GRDD_3D_BP1

            IF (SOLID_PHASE_ONLY) RETURN

            T_NOW=CURRENT_TIME()

            IF (GPU_DEBUG_H2D_D2H) THEN
                GM%UD    = U
                GM%VD    = V
                GM%WD    = W
                GM%USD   = US
                GM%VSD   = VS
                GM%WSD   = WS
                GM%FVXD  = FVX
                GM%FVYD  = FVY
                GM%FVZD  = FVZ
                GM%HSD   = HS
            ENDIF
            IF (FREEZE_VELOCITY) THEN
                GM%UD = GM%USD
                GM%VD = GM%VSD
                GM%WD = GM%WSD
            ELSE
                CALL UVW_CORRECTOR_KERNEL<<<GRDD,BLKD,0,GM%STRM>>>( &
                    GM%UD,GM%VD,GM%WD,GM%USD,GM%VSD,GM%WSD,GM%FVXD,GM%FVYD,GM%FVZD, &
                    GM%RDXND,GM%RDYND,GM%RDZND,GM%HSD,DT,IBAR,JBAR,KBAR)
            ENDIF

            IF (GPU_DEBUG_H2D_D2H) THEN
                U = GM%UD
                V = GM%VD
                W = GM%WD
            ENDIF

            T_USED(4)=T_USED(4)+CURRENT_TIME()-T_NOW
        END SUBROUTINE VELOCITY_CORRECTOR_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE UVW_CORRECTOR_KERNEL(U,V,W,US,VS,WS,FVX,FVY,FVZ,RDXN,RDYN,RDZN,HS,DT, &
            IBAR,JBAR,KBAR)

            INTEGER, VALUE :: IBAR, JBAR, KBAR
            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: U,V,W,US,VS,WS,FVX,FVY,FVZ,HS
            REAL(EB), DEVICE, DIMENSION(0:*) :: RDXN,RDYN,RDZN
            REAL(EB), VALUE :: DT

            INTEGER :: I,J,K

            I = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            J = (BLOCKIDX%Y - 1) * BLOCKDIM%Y + THREADIDX%Y - 1
            K = (BLOCKIDX%Z - 1) * BLOCKDIM%Z + THREADIDX%Z - 1

            IF (I>=0 .AND. I<=IBAR .AND. J>=1 .AND. J<=JBAR .AND. K>=1 .AND. K<=KBAR) &
                U(I,J,K) = 0.5_EB*(U(I,J,K)+US(I,J,K)-DT*(FVX(I,J,K)+RDXN(I)*(HS(I+1,J,K)-HS(I,J,K))))

            IF (I>=1 .AND. I<=IBAR .AND. J>=0 .AND. J<=JBAR .AND. K>=1 .AND. K<=KBAR) &
                V(I,J,K) = 0.5_EB*(V(I,J,K)+VS(I,J,K)-DT*(FVY(I,J,K)+RDYN(J)*(HS(I,J+1,K)-HS(I,J,K))))

            IF (I>=1 .AND. I<=IBAR .AND. J>=1 .AND. J<=JBAR .AND. K>=0 .AND. K<=KBAR) &
                W(I,J,K) = 0.5_EB*(W(I,J,K)+WS(I,J,K)-DT*(FVZ(I,J,K)+RDZN(K)*(HS(I,J,K+1)-HS(I,J,K))))

        END SUBROUTINE UVW_CORRECTOR_KERNEL
        !--------------------------------------------------------------------------
        ! Velocity-boundary device helper and host orchestration
        !--------------------------------------------------------------------------
        ATTRIBUTES(DEVICE) SUBROUTINE WALL_MODEL_DEVICE( &
            SLIP_FACTOR,U_TAU,Y_PLUS, &
            NU,S,Y_EXTERNAL_POINT,U_EXTERNAL_POINT, &
            SIM_MODE_LOCAL,DNS_MODE_LOCAL, &
            VON_KARMAN_CONSTANT_LOCAL,Y_WERNER_WENGLE_LOCAL)

            IMPLICIT NONE

            REAL(EB), INTENT(OUT) :: SLIP_FACTOR,U_TAU,Y_PLUS
            REAL(EB), INTENT(IN)  :: NU,S,Y_EXTERNAL_POINT,U_EXTERNAL_POINT
            INTEGER , VALUE       :: SIM_MODE_LOCAL,DNS_MODE_LOCAL
            REAL(EB), VALUE       :: VON_KARMAN_CONSTANT_LOCAL,Y_WERNER_WENGLE_LOCAL

            REAL(EB), PARAMETER :: B=5.2_EB,BTILDE_MAX=9.5_EB,BTILDE_ROUGH=8.5_EB
            REAL(EB), PARAMETER :: S0=1._EB,S1=5.83_EB,S2=30._EB,EPS=1.E-10_EB
            REAL(EB) :: U,Y_CELL_CENTER,TAU_W,BTILDE,S_PLUS,DUDY,DY,RKAPPA,DELTA_NU
            INTEGER :: ITER

            RKAPPA = 1._EB / VON_KARMAN_CONSTANT_LOCAL
            DY = 2._EB*Y_EXTERNAL_POINT
            Y_CELL_CENTER = Y_EXTERNAL_POINT
            U = U_EXTERNAL_POINT

            DUDY = ABS(U) / MAX(Y_CELL_CENTER,EPS)
            TAU_W = NU * DUDY
            U_TAU = SQRT(ABS(TAU_W))
            DELTA_NU = NU / (U_TAU + EPS)
            Y_PLUS = Y_CELL_CENTER / (DELTA_NU + EPS)
            SLIP_FACTOR = -1._EB

            IF (SIM_MODE_LOCAL /= DNS_MODE_LOCAL) THEN
                DO ITER=1,3
                    S_PLUS = S / (DELTA_NU + EPS)

                    IF (S_PLUS < S0) THEN
                        Y_PLUS = Y_CELL_CENTER / (DELTA_NU + EPS)

                        IF (Y_PLUS < Y_WERNER_WENGLE_LOCAL) THEN
                            TAU_W = ( U / MAX(Y_PLUS,EPS) )**2
                            U_TAU = SQRT(ABS(TAU_W))
                            DUDY  = ABS(U) / MAX(Y_CELL_CENTER,EPS)
                        ELSE
                            TAU_W = ( U / (RKAPPA*LOG(MAX(Y_PLUS,1.000001_EB)) + B) )**2
                            U_TAU = SQRT(ABS(TAU_W))
                            DUDY  = U_TAU*RKAPPA / MAX(Y_CELL_CENTER,EPS)
                        ENDIF
                    ELSE
                        IF (S_PLUS < S1) THEN
                            BTILDE = B + RKAPPA*LOG(S_PLUS)
                        ELSEIF (S_PLUS < S2) THEN
                            BTILDE = BTILDE_MAX
                        ELSE
                            BTILDE = BTILDE_ROUGH
                        ENDIF
                        Y_PLUS = Y_CELL_CENTER / S
                        TAU_W = ( U / (RKAPPA*LOG(Y_PLUS) + BTILDE) )**2
                        U_TAU = SQRT(ABS(TAU_W))
                        DUDY  = U_TAU*RKAPPA / Y_CELL_CENTER
                    ENDIF

                    DELTA_NU = NU / (U_TAU + EPS)
                ENDDO

                SLIP_FACTOR = MAX(-1._EB,MIN(1._EB,1._EB - DUDY*DY/(ABS(U)+EPS)))
            ENDIF

        END SUBROUTINE WALL_MODEL_DEVICE

        ATTRIBUTES(DEVICE) SUBROUTINE WALL_MODEL_OUTPUT_DEVICE( &
            SLIP_FACTOR,U_TAU,Y_PLUS,NU,S,Y_EXTERNAL_POINT,U_EXTERNAL_POINT, &
            SIM_MODE_LOCAL,DNS_MODE_LOCAL,VON_KARMAN_CONSTANT_LOCAL,Y_WERNER_WENGLE_LOCAL)
            IMPLICIT NONE

            REAL(EB), INTENT(OUT) :: SLIP_FACTOR,U_TAU,Y_PLUS
            REAL(EB), INTENT(IN) :: NU,S,Y_EXTERNAL_POINT,U_EXTERNAL_POINT
            INTEGER, VALUE :: SIM_MODE_LOCAL,DNS_MODE_LOCAL
            REAL(EB), VALUE :: VON_KARMAN_CONSTANT_LOCAL,Y_WERNER_WENGLE_LOCAL

            REAL(EB), PARAMETER :: B=5.2_EB,BTILDE_MAX=9.5_EB,BTILDE_ROUGH=8.5_EB
            REAL(EB), PARAMETER :: S0=1._EB,S1=5.83_EB,S2=30._EB,EPS=1.E-10_EB
            REAL(EB) :: U,Y_CELL_CENTER,TAU_W,BTILDE,S_PLUS,DUDY,DY,RKAPPA,DELTA_NU
            INTEGER :: ITER

            RKAPPA = 1._EB / VON_KARMAN_CONSTANT_LOCAL
            DY = 2._EB*Y_EXTERNAL_POINT
            Y_CELL_CENTER = Y_EXTERNAL_POINT
            U = U_EXTERNAL_POINT

            DUDY = ABS(U) / MAX(Y_CELL_CENTER,EPS)
            TAU_W = NU*DUDY
            U_TAU = SQRT(ABS(TAU_W))
            DELTA_NU = NU/(U_TAU+EPS)
            Y_PLUS = Y_CELL_CENTER/(DELTA_NU+EPS)
            SLIP_FACTOR = -1._EB

            IF (SIM_MODE_LOCAL/=DNS_MODE_LOCAL) THEN
                DO ITER=1,3
                    S_PLUS = S/(DELTA_NU+EPS)
                    IF (S_PLUS<S0) THEN
                        Y_PLUS = Y_CELL_CENTER/(DELTA_NU+EPS)
                        IF (Y_PLUS<Y_WERNER_WENGLE_LOCAL) THEN
                            TAU_W = (U/MAX(Y_PLUS,EPS))**2
                            U_TAU = SQRT(ABS(TAU_W))
                            DUDY = ABS(U)/MAX(Y_CELL_CENTER,EPS)
                        ELSE
                            TAU_W = (U/(RKAPPA*LOG(MAX(Y_PLUS,1.000001_EB))+B))**2
                            U_TAU = SQRT(ABS(TAU_W))
                            DUDY = U_TAU*RKAPPA/MAX(Y_CELL_CENTER,EPS)
                        ENDIF
                    ELSE
                        IF (S_PLUS<S1) THEN
                            BTILDE = B+RKAPPA*LOG(MAX(S_PLUS,EPS))
                        ELSEIF (S_PLUS<S2) THEN
                            BTILDE = BTILDE_MAX
                        ELSE
                            BTILDE = BTILDE_ROUGH
                        ENDIF
                        Y_PLUS = Y_CELL_CENTER/MAX(S,EPS)
                        TAU_W = (U/(RKAPPA*LOG(MAX(Y_PLUS,1.000001_EB))+BTILDE))**2
                        U_TAU = SQRT(ABS(TAU_W))
                        DUDY = U_TAU*RKAPPA/MAX(Y_CELL_CENTER,EPS)
                    ENDIF
                    DELTA_NU = NU/(U_TAU+EPS)
                ENDDO
                SLIP_FACTOR = MAX(-1._EB,MIN(1._EB,1._EB-DUDY*DY/(ABS(U)+EPS)))
            ENDIF

        END SUBROUTINE WALL_MODEL_OUTPUT_DEVICE

        SUBROUTINE COMPUTE_WALL_MODEL_OUTPUT_GPU(NM)
            IMPLICIT NONE

            INTEGER, INTENT(IN) :: NM
            INTEGER :: TOTAL_WALL_CELLS
            TYPE(DIM3) :: BLK,GRD

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)
            TOTAL_WALL_CELLS = N_EXTERNAL_WALL_CELLS + N_INTERNAL_WALL_CELLS
            IF (TOTAL_WALL_CELLS<=0) RETURN

            BLK = DIM3(GPU_BLKD_1D,1,1)
            GRD = DIM3((TOTAL_WALL_CELLS+GPU_BLKD_1D-1)/GPU_BLKD_1D,1,1)
            CALL WALL_MODEL_OUTPUT_KERNEL<<<GRD,BLK,0,GM%STRM>>>(GM%WALL_TOPO_AOS,GM%WALL_SF_INDEX,GM%WALL_RDN, &
                GM%WALL_U_TANG,GM%MU_DNSD,GM%RHOD,SURF_ROUGHNESS,GM%WALL_U_TAU,GM%WALL_Y_PLUS,TOTAL_WALL_CELLS, &
                IBP1,JBP1,SOLID_BOUNDARY,SIM_MODE,DNS_MODE,VON_KARMAN_CONSTANTD,Y_WERNER_WENGLED)

        END SUBROUTINE COMPUTE_WALL_MODEL_OUTPUT_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE WALL_MODEL_OUTPUT_KERNEL(WALL_TOPO,WALL_SF_INDEX,WALL_RDN,WALL_U_TANG,MU_DNSD,RHOD, &
            SURF_ROUGHNESS,WALL_U_TAU,WALL_Y_PLUS,TOTAL_WALL_CELLS,IBP1,JBP1,SOLID_BOUNDARY_LOCAL,SIM_MODE_LOCAL,DNS_MODE_LOCAL, &
            VON_KARMAN_CONSTANT_LOCAL,Y_WERNER_WENGLE_LOCAL)
            IMPLICIT NONE

            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, DEVICE :: WALL_SF_INDEX(0:*)
            REAL(EB), DEVICE :: WALL_RDN(0:*),WALL_U_TANG(0:*),SURF_ROUGHNESS(0:*),WALL_U_TAU(0:*),WALL_Y_PLUS(0:*)
            REAL(EB), DEVICE :: MU_DNSD(0:IBP1,0:JBP1,0:*),RHOD(0:IBP1,0:JBP1,0:*)
            INTEGER, VALUE :: TOTAL_WALL_CELLS,IBP1,JBP1,SOLID_BOUNDARY_LOCAL,SIM_MODE_LOCAL,DNS_MODE_LOCAL
            REAL(EB), DEVICE :: VON_KARMAN_CONSTANT_LOCAL,Y_WERNER_WENGLE_LOCAL

            INTEGER :: IW,IIG,JJG,KKG,SF_INDEX
            REAL(EB) :: SLIP_FACTOR,NU,VK_LOCAL,YW_LOCAL

            IW = (BLOCKIDX%X-1)*BLOCKDIM%X + THREADIDX%X
            IF (IW<1 .OR. IW>TOTAL_WALL_CELLS) RETURN

            IF (WALL_TOPO(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY_LOCAL) THEN
                WALL_U_TAU(IW) = 0._EB
                WALL_Y_PLUS(IW) = 1._EB
                RETURN
            ENDIF

            IIG = WALL_TOPO(IW)%IIG
            JJG = WALL_TOPO(IW)%JJG
            KKG = WALL_TOPO(IW)%KKG
            SF_INDEX = WALL_SF_INDEX(IW)
            NU = MU_DNSD(IIG,JJG,KKG) / MAX(RHOD(IIG,JJG,KKG),1.E-20_EB)
            VK_LOCAL = VON_KARMAN_CONSTANT_LOCAL
            YW_LOCAL = Y_WERNER_WENGLE_LOCAL

            CALL WALL_MODEL_OUTPUT_DEVICE(SLIP_FACTOR,WALL_U_TAU(IW),WALL_Y_PLUS(IW),NU,SURF_ROUGHNESS(SF_INDEX), &
                0.5_EB*WALL_RDN(IW)**(-1),WALL_U_TANG(IW),SIM_MODE_LOCAL,DNS_MODE_LOCAL, &
                VK_LOCAL,YW_LOCAL)

        END SUBROUTINE WALL_MODEL_OUTPUT_KERNEL

        SUBROUTINE VELOCITY_BC_GPU(T,NM,APPLY_TO_ESTIMATED_VARIABLES)
            IMPLICIT NONE
            REAL(EB) :: T
            INTEGER :: IE,IOR,NOM,AXIS,NEDGE_PREP,STRIDE_1,STRIDE_2
            INTEGER, INTENT(IN) :: NM
            LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
            TYPE(OMESH_TYPE), POINTER :: OM
            TYPE(GPU_OMESH_TYPE), POINTER :: GOM_LOCAL

            TYPE(DIM3) :: BLK,GRD,GRD_AXIS1,GRD_AXIS2,GRD_AXIS3,GRD_PREP
            INTEGER :: EDGE_ACTIVE_LB, EDGE_ACTIVE_UB, NEDGE_ACTIVE
            REAL(EB), ALLOCATABLE :: EDGE_OMEGA_H(:,:), EDGE_TAU_H(:,:)
            REAL(EB), DEVICE, POINTER, CONTIGUOUS :: OM_VEL_1D(:,:,:),OM_VEL_2D(:,:,:)
            INTEGER, DEVICE, POINTER, CONTIGUOUS :: EDGE_AXIS_LISTD(:)

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)

            IF (OPEN_WIND_BOUNDARY) THEN
                GM%U_WINDD = U_WIND
                GM%V_WINDD = V_WIND
                GM%W_WINDD = W_WIND
            ENDIF

            DRAG_UVWMAX = 0._EB

            EDGE_ACTIVE_LB = 1
            EDGE_ACTIVE_UB = EDGE_COUNT(NM)
            NEDGE_ACTIVE   = EDGE_ACTIVE_UB - EDGE_ACTIVE_LB + 1

            IF (NEDGE_ACTIVE <= 0) RETURN

            BLK = DIM3(64,1,1)
            GRD = DIM3((NEDGE_ACTIVE + BLK%X - 1)/BLK%X,1,1)
            GRD_AXIS1 = DIM3((GM%N_EDGE_AXIS1 + BLK%X - 1)/BLK%X,1,1)
            GRD_AXIS2 = DIM3((GM%N_EDGE_AXIS2 + BLK%X - 1)/BLK%X,1,1)
            GRD_AXIS3 = DIM3((GM%N_EDGE_AXIS3 + BLK%X - 1)/BLK%X,1,1)

            GM%EDGE_GHOST_1 = 0._EB
            GM%EDGE_GHOST_2 = 0._EB
            DO NOM=1,NMESHES
                IF (.NOT.ALLOCATED(GM%OMESH)) EXIT
                GOM_LOCAL => GM%OMESH(NOM)
                IF (.NOT.GOM_LOCAL%ACTIVE) CYCLE
                OM => MESHES(NM)%OMESH(NOM)

                IF (GPU_DEBUG_H2D_D2H) THEN
                    IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
                        GOM_LOCAL%US_D = OM%US
                        GOM_LOCAL%VS_D = OM%VS
                        GOM_LOCAL%WS_D = OM%WS
                    ELSE
                        GOM_LOCAL%U_D = OM%U
                        GOM_LOCAL%V_D = OM%V
                        GOM_LOCAL%W_D = OM%W
                    ENDIF
                ENDIF

                DO AXIS=1,3
                    SELECT CASE(AXIS)
                        CASE(1)
                            IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
                                OM_VEL_1D => GOM_LOCAL%WS_D
                                OM_VEL_2D => GOM_LOCAL%VS_D
                            ELSE
                                OM_VEL_1D => GOM_LOCAL%W_D
                                OM_VEL_2D => GOM_LOCAL%V_D
                            ENDIF
                            EDGE_AXIS_LISTD => GM%EDGE_AXIS1_LIST
                            NEDGE_PREP = GM%N_EDGE_AXIS1
                            GRD_PREP = GRD_AXIS1
                            STRIDE_1 = SIZE(OM_VEL_1D,1)*SIZE(OM_VEL_1D,2)
                            STRIDE_2 = SIZE(OM_VEL_1D,1)
                        CASE(2)
                            IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
                                OM_VEL_1D => GOM_LOCAL%US_D
                                OM_VEL_2D => GOM_LOCAL%WS_D
                            ELSE
                                OM_VEL_1D => GOM_LOCAL%U_D
                                OM_VEL_2D => GOM_LOCAL%W_D
                            ENDIF
                            EDGE_AXIS_LISTD => GM%EDGE_AXIS2_LIST
                            NEDGE_PREP = GM%N_EDGE_AXIS2
                            GRD_PREP = GRD_AXIS2
                            STRIDE_1 = 1
                            STRIDE_2 = SIZE(OM_VEL_1D,1)*SIZE(OM_VEL_1D,2)
                        CASE(3)
                            IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
                                OM_VEL_1D => GOM_LOCAL%VS_D
                                OM_VEL_2D => GOM_LOCAL%US_D
                            ELSE
                                OM_VEL_1D => GOM_LOCAL%V_D
                                OM_VEL_2D => GOM_LOCAL%U_D
                            ENDIF
                            EDGE_AXIS_LISTD => GM%EDGE_AXIS3_LIST
                            NEDGE_PREP = GM%N_EDGE_AXIS3
                            GRD_PREP = GRD_AXIS3
                            STRIDE_1 = SIZE(OM_VEL_1D,1)
                            STRIDE_2 = 1
                    END SELECT
                    IF (NEDGE_PREP <= 0) CYCLE
                    CALL PREPARE_VELOCITY_BC_EDGE_GHOST_GPU<<<GRD_PREP,BLK,0,GM%STRM>>>(GM%EDGE_GHOST_1,GM%EDGE_GHOST_2,OM_VEL_1D, &
                        OM_VEL_2D,EDGE_AXIS_LISTD,NEDGE_PREP,GM%EDGE_NOM_1,GM%EDGE_NOM_2,GM%EDGE_IIO_1,GM%EDGE_JJO_1,GM%EDGE_KKO_1, &
                        GM%EDGE_IIO_2,GM%EDGE_JJO_2,GM%EDGE_KKO_2,GM%EDGE_INTERP_1,GM%EDGE_INTERP_2,NOM,STRIDE_1,STRIDE_2, &
                        SIZE(OM_VEL_1D,1),SIZE(OM_VEL_1D,2),LBOUND(OM_VEL_1D,1),LBOUND(OM_VEL_1D,2), &
                        LBOUND(OM_VEL_1D,3),GM%EDGE_LB)
                ENDDO
            ENDDO

            IF (APPLY_TO_ESTIMATED_VARIABLES) THEN

                IF (GPU_DEBUG_H2D_D2H) THEN
                    GM%USD   = US
                    GM%VSD   = VS
                    GM%WSD   = WS
                    GM%MUD   = MU
                    GM%RHOSD = RHOS
                    GM%TMPD  = TMP
                    CALL SYNC_EDGE_OMEGA_TAU_TO_GPU(NM)
                ENDIF

                IF (GM%N_EDGE_AXIS1 > 0) &
                    CALL VELOCITY_BC_GPU_GHOST_KERNEL_AXIS<<<GRD_AXIS1,BLK,0,GM%STRM>>>( &
                        GM%USD,GM%VSD,GM%WSD,GM%RHOSD,GM%TMPD, &
                        GM%EDGE_GHOST_1,GM%EDGE_GHOST_2, &
                        GM%EDGE_U_AVG,GM%EDGE_V_AVG,GM%EDGE_W_AVG, &
                        GM%EDGE_OMEGA,GM%EDGE_TAU,GM%EDGE_SIDE_ISGN,GM%EDGE_SIDE_IOR, &
                        GM%EDGE_SIDE_IS2,GM%EDGE_SIDE_ICD,GM%EDGE_SIDE_IWM,GM%EDGE_SIDE_IWMI, &
                        GM%EDGE_SIDE_IWP,GM%EDGE_SIDE_IWPI, &
                        GM%EDGE_AXIS1_LIST,GM%N_EDGE_AXIS1, &
                        GM%EDGE_I,GM%EDGE_J,GM%EDGE_K, &
                        GM%EDGE_CELL_INDEX_MM,GM%EDGE_CELL_INDEX_PM,GM%EDGE_CELL_INDEX_MP,GM%EDGE_CELL_INDEX_PP, &
                        GM%EDGE_NOM_1,GM%EDGE_NOM_2, &
                        GM%EDGE_EXTERNAL, &
                        GM%CELL_I,GM%CELL_J,GM%CELL_K,GM%CELL_SOLID,GM%CELL_EXTERIOR, &
                        GM%WALL_TOPO_AOS,GM%WALL_SF_INDEX, &
                        SURF_VELOCITY_BC_INDEX,SURF_SPECIFIED_NORMAL_VELOCITY,SURF_SPECIFIED_TANGENTIAL_VELOCITY, &
                        SURF_SPECIFIED_NORMAL_GRADIENT,SURF_ROUGHNESS,SURF_VEL_T_1,SURF_VEL_T_2, &
                        SF_RAMP_VELO_INDEX,SF_RAMP_VELO_TAU,RAMP_T_MIN_D,RAMP_SPAN_D,RAMP_RDT_D, &
                        RAMP_NPOINTS_D,RAMP_OFFSET_D,RAMP_INTERP_D,T, &
                        GM%WALL_RDN, &
                        MU_RSQMW_ZD,RSQ_MW_ZD, &
                        GM%DXD,GM%DYD,GM%DZD, &
                        VON_KARMAN_CONSTANTD,Y_WERNER_WENGLED, &
                        IBAR,JBAR,KBAR,I_MAX_TEMP, &
                        SOLID_BOUNDARY,FREE_SLIP_BC,NO_SLIP_BC,WALL_MODEL_BC, &
                        INTERPOLATED_BOUNDARY, &
                        SIM_MODE,DNS_MODE,1,GM%EDGE_LB,GM%EDGE_UB,0,OPEN_WIND_BOUNDARY,GM%U_WINDD,GM%V_WINDD,GM%W_WINDD)
                IF (GM%N_EDGE_AXIS2 > 0) &
                    CALL VELOCITY_BC_GPU_GHOST_KERNEL_AXIS<<<GRD_AXIS2,BLK,0,GM%STRM>>>( &
                        GM%USD,GM%VSD,GM%WSD,GM%RHOSD,GM%TMPD, &
                        GM%EDGE_GHOST_1,GM%EDGE_GHOST_2, &
                        GM%EDGE_U_AVG,GM%EDGE_V_AVG,GM%EDGE_W_AVG, &
                        GM%EDGE_OMEGA,GM%EDGE_TAU,GM%EDGE_SIDE_ISGN,GM%EDGE_SIDE_IOR, &
                        GM%EDGE_SIDE_IS2,GM%EDGE_SIDE_ICD,GM%EDGE_SIDE_IWM,GM%EDGE_SIDE_IWMI, &
                        GM%EDGE_SIDE_IWP,GM%EDGE_SIDE_IWPI, &
                        GM%EDGE_AXIS2_LIST,GM%N_EDGE_AXIS2, &
                        GM%EDGE_I,GM%EDGE_J,GM%EDGE_K, &
                        GM%EDGE_CELL_INDEX_MM,GM%EDGE_CELL_INDEX_PM,GM%EDGE_CELL_INDEX_MP,GM%EDGE_CELL_INDEX_PP, &
                        GM%EDGE_NOM_1,GM%EDGE_NOM_2, &
                        GM%EDGE_EXTERNAL, &
                        GM%CELL_I,GM%CELL_J,GM%CELL_K,GM%CELL_SOLID,GM%CELL_EXTERIOR, &
                        GM%WALL_TOPO_AOS,GM%WALL_SF_INDEX, &
                        SURF_VELOCITY_BC_INDEX,SURF_SPECIFIED_NORMAL_VELOCITY,SURF_SPECIFIED_TANGENTIAL_VELOCITY, &
                        SURF_SPECIFIED_NORMAL_GRADIENT,SURF_ROUGHNESS,SURF_VEL_T_1,SURF_VEL_T_2, &
                        SF_RAMP_VELO_INDEX,SF_RAMP_VELO_TAU,RAMP_T_MIN_D,RAMP_SPAN_D,RAMP_RDT_D, &
                        RAMP_NPOINTS_D,RAMP_OFFSET_D,RAMP_INTERP_D,T, &
                        GM%WALL_RDN, &
                        MU_RSQMW_ZD,RSQ_MW_ZD, &
                        GM%DXD,GM%DYD,GM%DZD, &
                        VON_KARMAN_CONSTANTD,Y_WERNER_WENGLED, &
                        IBAR,JBAR,KBAR,I_MAX_TEMP, &
                        SOLID_BOUNDARY,FREE_SLIP_BC,NO_SLIP_BC,WALL_MODEL_BC, &
                        INTERPOLATED_BOUNDARY, &
                        SIM_MODE,DNS_MODE,2,GM%EDGE_LB,GM%EDGE_UB,0,OPEN_WIND_BOUNDARY,GM%U_WINDD,GM%V_WINDD,GM%W_WINDD)
                IF (GM%N_EDGE_AXIS3 > 0) &
                    CALL VELOCITY_BC_GPU_GHOST_KERNEL_AXIS<<<GRD_AXIS3,BLK,0,GM%STRM>>>( &
                        GM%USD,GM%VSD,GM%WSD,GM%RHOSD,GM%TMPD, &
                        GM%EDGE_GHOST_1,GM%EDGE_GHOST_2, &
                        GM%EDGE_U_AVG,GM%EDGE_V_AVG,GM%EDGE_W_AVG, &
                        GM%EDGE_OMEGA,GM%EDGE_TAU,GM%EDGE_SIDE_ISGN,GM%EDGE_SIDE_IOR, &
                        GM%EDGE_SIDE_IS2,GM%EDGE_SIDE_ICD,GM%EDGE_SIDE_IWM,GM%EDGE_SIDE_IWMI, &
                        GM%EDGE_SIDE_IWP,GM%EDGE_SIDE_IWPI, &
                        GM%EDGE_AXIS3_LIST,GM%N_EDGE_AXIS3, &
                        GM%EDGE_I,GM%EDGE_J,GM%EDGE_K, &
                        GM%EDGE_CELL_INDEX_MM,GM%EDGE_CELL_INDEX_PM,GM%EDGE_CELL_INDEX_MP,GM%EDGE_CELL_INDEX_PP, &
                        GM%EDGE_NOM_1,GM%EDGE_NOM_2, &
                        GM%EDGE_EXTERNAL, &
                        GM%CELL_I,GM%CELL_J,GM%CELL_K,GM%CELL_SOLID,GM%CELL_EXTERIOR, &
                        GM%WALL_TOPO_AOS,GM%WALL_SF_INDEX, &
                        SURF_VELOCITY_BC_INDEX,SURF_SPECIFIED_NORMAL_VELOCITY,SURF_SPECIFIED_TANGENTIAL_VELOCITY, &
                        SURF_SPECIFIED_NORMAL_GRADIENT,SURF_ROUGHNESS,SURF_VEL_T_1,SURF_VEL_T_2, &
                        SF_RAMP_VELO_INDEX,SF_RAMP_VELO_TAU,RAMP_T_MIN_D,RAMP_SPAN_D,RAMP_RDT_D, &
                        RAMP_NPOINTS_D,RAMP_OFFSET_D,RAMP_INTERP_D,T, &
                        GM%WALL_RDN, &
                        MU_RSQMW_ZD,RSQ_MW_ZD, &
                        GM%DXD,GM%DYD,GM%DZD, &
                        VON_KARMAN_CONSTANTD,Y_WERNER_WENGLED, &
                        IBAR,JBAR,KBAR,I_MAX_TEMP, &
                        SOLID_BOUNDARY,FREE_SLIP_BC,NO_SLIP_BC,WALL_MODEL_BC, &
                        INTERPOLATED_BOUNDARY, &
                        SIM_MODE,DNS_MODE,3,GM%EDGE_LB,GM%EDGE_UB,0,OPEN_WIND_BOUNDARY,GM%U_WINDD,GM%V_WINDD,GM%W_WINDD)
                CALL VELOCITY_BC_GPU_TAU_OMEGA_KERNEL<<<GRD,BLK,0,GM%STRM>>>( &
                    GM%EDGE_OMEGA,GM%EDGE_TAU,GM%EDGE_LB,GM%EDGE_UB,EDGE_ACTIVE_LB,EDGE_ACTIVE_UB)
            ELSE
                IF (GPU_DEBUG_H2D_D2H) THEN
                    GM%UD   = U
                    GM%VD   = V
                    GM%WD   = W
                    GM%MUD  = MU
                    GM%RHOD = RHO
                    GM%TMPD = TMP
                    CALL SYNC_EDGE_OMEGA_TAU_TO_GPU(NM)
                ENDIF

                IF (GM%N_EDGE_AXIS1 > 0) &
                    CALL VELOCITY_BC_GPU_GHOST_KERNEL_AXIS<<<GRD_AXIS1,BLK,0,GM%STRM>>>( &
                        GM%UD,GM%VD,GM%WD,GM%RHOD,GM%TMPD, &
                        GM%EDGE_GHOST_1,GM%EDGE_GHOST_2, &
                        GM%EDGE_U_AVG,GM%EDGE_V_AVG,GM%EDGE_W_AVG, &
                        GM%EDGE_OMEGA,GM%EDGE_TAU,GM%EDGE_SIDE_ISGN,GM%EDGE_SIDE_IOR, &
                        GM%EDGE_SIDE_IS2,GM%EDGE_SIDE_ICD,GM%EDGE_SIDE_IWM,GM%EDGE_SIDE_IWMI, &
                        GM%EDGE_SIDE_IWP,GM%EDGE_SIDE_IWPI, &
                        GM%EDGE_AXIS1_LIST,GM%N_EDGE_AXIS1, &
                        GM%EDGE_I,GM%EDGE_J,GM%EDGE_K, &
                        GM%EDGE_CELL_INDEX_MM,GM%EDGE_CELL_INDEX_PM,GM%EDGE_CELL_INDEX_MP,GM%EDGE_CELL_INDEX_PP, &
                        GM%EDGE_NOM_1,GM%EDGE_NOM_2, &
                        GM%EDGE_EXTERNAL, &
                        GM%CELL_I,GM%CELL_J,GM%CELL_K,GM%CELL_SOLID,GM%CELL_EXTERIOR, &
                        GM%WALL_TOPO_AOS,GM%WALL_SF_INDEX, &
                        SURF_VELOCITY_BC_INDEX,SURF_SPECIFIED_NORMAL_VELOCITY,SURF_SPECIFIED_TANGENTIAL_VELOCITY, &
                        SURF_SPECIFIED_NORMAL_GRADIENT,SURF_ROUGHNESS,SURF_VEL_T_1,SURF_VEL_T_2, &
                        SF_RAMP_VELO_INDEX,SF_RAMP_VELO_TAU,RAMP_T_MIN_D,RAMP_SPAN_D,RAMP_RDT_D, &
                        RAMP_NPOINTS_D,RAMP_OFFSET_D,RAMP_INTERP_D,T, &
                        GM%WALL_RDN, &
                        MU_RSQMW_ZD,RSQ_MW_ZD, &
                        GM%DXD,GM%DYD,GM%DZD, &
                        VON_KARMAN_CONSTANTD,Y_WERNER_WENGLED, &
                        IBAR,JBAR,KBAR,I_MAX_TEMP, &
                        SOLID_BOUNDARY,FREE_SLIP_BC,NO_SLIP_BC,WALL_MODEL_BC, &
                        INTERPOLATED_BOUNDARY, &
                        SIM_MODE,DNS_MODE,1,GM%EDGE_LB,GM%EDGE_UB,1,OPEN_WIND_BOUNDARY,GM%U_WINDD,GM%V_WINDD,GM%W_WINDD)
                IF (GM%N_EDGE_AXIS2 > 0) &
                    CALL VELOCITY_BC_GPU_GHOST_KERNEL_AXIS<<<GRD_AXIS2,BLK,0,GM%STRM>>>( &
                        GM%UD,GM%VD,GM%WD,GM%RHOD,GM%TMPD, &
                        GM%EDGE_GHOST_1,GM%EDGE_GHOST_2, &
                        GM%EDGE_U_AVG,GM%EDGE_V_AVG,GM%EDGE_W_AVG, &
                        GM%EDGE_OMEGA,GM%EDGE_TAU,GM%EDGE_SIDE_ISGN,GM%EDGE_SIDE_IOR, &
                        GM%EDGE_SIDE_IS2,GM%EDGE_SIDE_ICD,GM%EDGE_SIDE_IWM,GM%EDGE_SIDE_IWMI, &
                        GM%EDGE_SIDE_IWP,GM%EDGE_SIDE_IWPI, &
                        GM%EDGE_AXIS2_LIST,GM%N_EDGE_AXIS2, &
                        GM%EDGE_I,GM%EDGE_J,GM%EDGE_K, &
                        GM%EDGE_CELL_INDEX_MM,GM%EDGE_CELL_INDEX_PM,GM%EDGE_CELL_INDEX_MP,GM%EDGE_CELL_INDEX_PP, &
                        GM%EDGE_NOM_1,GM%EDGE_NOM_2, &
                        GM%EDGE_EXTERNAL, &
                        GM%CELL_I,GM%CELL_J,GM%CELL_K,GM%CELL_SOLID,GM%CELL_EXTERIOR, &
                        GM%WALL_TOPO_AOS,GM%WALL_SF_INDEX, &
                        SURF_VELOCITY_BC_INDEX,SURF_SPECIFIED_NORMAL_VELOCITY,SURF_SPECIFIED_TANGENTIAL_VELOCITY, &
                        SURF_SPECIFIED_NORMAL_GRADIENT,SURF_ROUGHNESS,SURF_VEL_T_1,SURF_VEL_T_2, &
                        SF_RAMP_VELO_INDEX,SF_RAMP_VELO_TAU,RAMP_T_MIN_D,RAMP_SPAN_D,RAMP_RDT_D, &
                        RAMP_NPOINTS_D,RAMP_OFFSET_D,RAMP_INTERP_D,T, &
                        GM%WALL_RDN, &
                        MU_RSQMW_ZD,RSQ_MW_ZD, &
                        GM%DXD,GM%DYD,GM%DZD, &
                        VON_KARMAN_CONSTANTD,Y_WERNER_WENGLED, &
                        IBAR,JBAR,KBAR,I_MAX_TEMP, &
                        SOLID_BOUNDARY,FREE_SLIP_BC,NO_SLIP_BC,WALL_MODEL_BC, &
                        INTERPOLATED_BOUNDARY, &
                        SIM_MODE,DNS_MODE,2,GM%EDGE_LB,GM%EDGE_UB,1,OPEN_WIND_BOUNDARY,GM%U_WINDD,GM%V_WINDD,GM%W_WINDD)
                IF (GM%N_EDGE_AXIS3 > 0) &
                    CALL VELOCITY_BC_GPU_GHOST_KERNEL_AXIS<<<GRD_AXIS3,BLK,0,GM%STRM>>>( &
                        GM%UD,GM%VD,GM%WD,GM%RHOD,GM%TMPD, &
                        GM%EDGE_GHOST_1,GM%EDGE_GHOST_2, &
                        GM%EDGE_U_AVG,GM%EDGE_V_AVG,GM%EDGE_W_AVG, &
                        GM%EDGE_OMEGA,GM%EDGE_TAU,GM%EDGE_SIDE_ISGN,GM%EDGE_SIDE_IOR, &
                        GM%EDGE_SIDE_IS2,GM%EDGE_SIDE_ICD,GM%EDGE_SIDE_IWM,GM%EDGE_SIDE_IWMI, &
                        GM%EDGE_SIDE_IWP,GM%EDGE_SIDE_IWPI, &
                        GM%EDGE_AXIS3_LIST,GM%N_EDGE_AXIS3, &
                        GM%EDGE_I,GM%EDGE_J,GM%EDGE_K, &
                        GM%EDGE_CELL_INDEX_MM,GM%EDGE_CELL_INDEX_PM,GM%EDGE_CELL_INDEX_MP,GM%EDGE_CELL_INDEX_PP, &
                        GM%EDGE_NOM_1,GM%EDGE_NOM_2, &
                        GM%EDGE_EXTERNAL, &
                        GM%CELL_I,GM%CELL_J,GM%CELL_K,GM%CELL_SOLID,GM%CELL_EXTERIOR, &
                        GM%WALL_TOPO_AOS,GM%WALL_SF_INDEX, &
                        SURF_VELOCITY_BC_INDEX,SURF_SPECIFIED_NORMAL_VELOCITY,SURF_SPECIFIED_TANGENTIAL_VELOCITY, &
                        SURF_SPECIFIED_NORMAL_GRADIENT,SURF_ROUGHNESS,SURF_VEL_T_1,SURF_VEL_T_2, &
                        SF_RAMP_VELO_INDEX,SF_RAMP_VELO_TAU,RAMP_T_MIN_D,RAMP_SPAN_D,RAMP_RDT_D, &
                        RAMP_NPOINTS_D,RAMP_OFFSET_D,RAMP_INTERP_D,T, &
                        GM%WALL_RDN, &
                        MU_RSQMW_ZD,RSQ_MW_ZD, &
                        GM%DXD,GM%DYD,GM%DZD, &
                        VON_KARMAN_CONSTANTD,Y_WERNER_WENGLED, &
                        IBAR,JBAR,KBAR,I_MAX_TEMP, &
                        SOLID_BOUNDARY,FREE_SLIP_BC,NO_SLIP_BC,WALL_MODEL_BC, &
                        INTERPOLATED_BOUNDARY, &
                        SIM_MODE,DNS_MODE,3,GM%EDGE_LB,GM%EDGE_UB,1,OPEN_WIND_BOUNDARY,GM%U_WINDD,GM%V_WINDD,GM%W_WINDD)
                CALL VELOCITY_BC_GPU_TAU_OMEGA_KERNEL<<<GRD,BLK,0,GM%STRM>>>( &
                    GM%EDGE_OMEGA,GM%EDGE_TAU,GM%EDGE_LB,GM%EDGE_UB,EDGE_ACTIVE_LB,EDGE_ACTIVE_UB)
            ENDIF

            IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
                IF (GPU_DEBUG_H2D_D2H) THEN
                    US = GM%USD
                    VS = GM%VSD
                    WS = GM%WSD
                ENDIF
            ELSE
                IF (GPU_DEBUG_H2D_D2H) THEN
                    U = GM%UD
                    V = GM%VD
                    W = GM%WD
                ENDIF
            ENDIF

            IF (GPU_DEBUG_H2D_D2H) THEN
                CALL GPU_WAIT(GM%STRM)
                ALLOCATE(EDGE_OMEGA_H(EDGE_ACTIVE_LB:EDGE_ACTIVE_UB,-2:2))
                ALLOCATE(EDGE_TAU_H  (EDGE_ACTIVE_LB:EDGE_ACTIVE_UB,-2:2))
                EDGE_OMEGA_H = GM%EDGE_OMEGA(EDGE_ACTIVE_LB:EDGE_ACTIVE_UB,-2:2)
                EDGE_TAU_H   = GM%EDGE_TAU  (EDGE_ACTIVE_LB:EDGE_ACTIVE_UB,-2:2)
                DO IE = EDGE_ACTIVE_LB, EDGE_ACTIVE_UB
                    DO IOR = -2, 2
                        EDGE(IE)%OMEGA(IOR) = EDGE_OMEGA_H(IE,IOR)
                        EDGE(IE)%TAU  (IOR) = EDGE_TAU_H  (IE,IOR)
                    ENDDO
                ENDDO
                DEALLOCATE(EDGE_OMEGA_H, EDGE_TAU_H)
            ENDIF

        END SUBROUTINE VELOCITY_BC_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE PREPARE_VELOCITY_BC_EDGE_GHOST_GPU(EDGE_GHOST_1,EDGE_GHOST_2,OM_VEL_1,OM_VEL_2, &
            EDGE_AXIS_LIST,NEDGE_AXIS,EDGE_NOM_1,EDGE_NOM_2,EDGE_IIO_1,EDGE_JJO_1,EDGE_KKO_1,EDGE_IIO_2,EDGE_JJO_2, &
            EDGE_KKO_2,EDGE_INTERP_1,EDGE_INTERP_2,NOM_FILTER,STRIDE_1,STRIDE_2,OM_NI,OM_NJ,OM_LB_I,OM_LB_J,OM_LB_K, &
            EDGE_LB)
            IMPLICIT NONE

            INTEGER, VALUE :: NEDGE_AXIS,NOM_FILTER,STRIDE_1,STRIDE_2,OM_NI,OM_NJ,OM_LB_I,OM_LB_J,OM_LB_K,EDGE_LB
            REAL(EB), DEVICE, DIMENSION(*) :: EDGE_GHOST_1,EDGE_GHOST_2
            REAL(EB), DEVICE, DIMENSION(*) :: OM_VEL_1,OM_VEL_2
            INTEGER, DEVICE, DIMENSION(*) :: EDGE_AXIS_LIST
            INTEGER, DEVICE, DIMENSION(*) :: EDGE_NOM_1,EDGE_NOM_2
            INTEGER, DEVICE, DIMENSION(*) :: EDGE_IIO_1,EDGE_JJO_1,EDGE_KKO_1
            INTEGER, DEVICE, DIMENSION(*) :: EDGE_IIO_2,EDGE_JJO_2,EDGE_KKO_2
            REAL(EB), DEVICE, DIMENSION(*) :: EDGE_INTERP_1,EDGE_INTERP_2
            INTEGER :: IDX,IE,IE_IDX,IDX_C
            REAL(EB) :: WGT

            IDX = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IDX < 1 .OR. IDX > NEDGE_AXIS) RETURN
            IE = EDGE_AXIS_LIST(IDX)
            IE_IDX = 1 + IE - EDGE_LB

            IF (ABS(EDGE_NOM_1(IE_IDX)) == NOM_FILTER) THEN
                IDX_C = 1 + (EDGE_IIO_1(IE_IDX)-OM_LB_I) + &
                        OM_NI*((EDGE_JJO_1(IE_IDX)-OM_LB_J) + OM_NJ*(EDGE_KKO_1(IE_IDX)-OM_LB_K))
                WGT = EDGE_INTERP_1(IE_IDX)
                EDGE_GHOST_1(IE_IDX) = WGT*OM_VEL_1(IDX_C) + (1._EB-WGT)*OM_VEL_1(IDX_C-STRIDE_1)
            ENDIF

            IF (ABS(EDGE_NOM_2(IE_IDX)) == NOM_FILTER) THEN
                IDX_C = 1 + (EDGE_IIO_2(IE_IDX)-OM_LB_I) + &
                        OM_NI*((EDGE_JJO_2(IE_IDX)-OM_LB_J) + OM_NJ*(EDGE_KKO_2(IE_IDX)-OM_LB_K))
                WGT = EDGE_INTERP_2(IE_IDX)
                EDGE_GHOST_2(IE_IDX) = WGT*OM_VEL_2(IDX_C) + (1._EB-WGT)*OM_VEL_2(IDX_C-STRIDE_2)
            ENDIF

        END SUBROUTINE PREPARE_VELOCITY_BC_EDGE_GHOST_GPU

        ATTRIBUTES(DEVICE) SUBROUTINE STORE_EDGE_VELOCITY_AVG(EDGE_U_AVG,EDGE_V_AVG,EDGE_W_AVG,IE,IEC,ICD, &
            VEL_GHOST,VEL_GAS)
            IMPLICIT NONE
            REAL(EB), DEVICE, DIMENSION(0:*) :: EDGE_U_AVG,EDGE_V_AVG,EDGE_W_AVG
            INTEGER, VALUE :: IE,IEC,ICD
            REAL(EB), VALUE :: VEL_GHOST,VEL_GAS
            REAL(EB) :: VEL_AVG

            VEL_AVG = 0.5_EB*(VEL_GHOST+VEL_GAS)
            SELECT CASE(IEC)
            CASE(1)
                IF (ICD==1) THEN
                    EDGE_W_AVG(IE) = VEL_AVG
                ELSE
                    EDGE_V_AVG(IE) = VEL_AVG
                ENDIF
            CASE(2)
                IF (ICD==1) THEN
                    EDGE_U_AVG(IE) = VEL_AVG
                ELSE
                    EDGE_W_AVG(IE) = VEL_AVG
                ENDIF
            CASE(3)
                IF (ICD==1) THEN
                    EDGE_V_AVG(IE) = VEL_AVG
                ELSE
                    EDGE_U_AVG(IE) = VEL_AVG
                ENDIF
            END SELECT

        END SUBROUTINE STORE_EDGE_VELOCITY_AVG

        !--------------------------------------------------------------------------
        ! One-thread-per-edge velocity ghost kernel
        !--------------------------------------------------------------------------
        ATTRIBUTES(GLOBAL) SUBROUTINE VELOCITY_BC_GPU_GHOST_KERNEL_AXIS(UU,VV,WW,RHOPD,TMPD_IN,EDGE_GHOST_1,EDGE_GHOST_2, &
            EDGE_U_AVG,EDGE_V_AVG,EDGE_W_AVG, &
            EDGE_OMEGA_IN,EDGE_TAU_IN,EDGE_SIDE_ISGN_IN,EDGE_SIDE_IOR_IN,EDGE_SIDE_IS2_IN,EDGE_SIDE_ICD_IN, &
            EDGE_SIDE_IWM_IN,EDGE_SIDE_IWMI_IN,EDGE_SIDE_IWP_IN,EDGE_SIDE_IWPI_IN, &
            EDGE_AXIS_LIST, NEDGE_AXIS, EDGE_I, EDGE_J, EDGE_K, EDGE_CELL_INDEX_MM, EDGE_CELL_INDEX_PM, &
            EDGE_CELL_INDEX_MP,EDGE_CELL_INDEX_PP,EDGE_NOM_1,EDGE_NOM_2,EDGE_EXTERNAL,CELL_I,CELL_J,CELL_K, &
            CELL_SOLID, CELL_EXTERIOR, WALL_TOPO, WALL_SF_INDEX, SURF_VELOCITY_BC_INDEX, SURF_SPECIFIED_NORMAL_VELOCITY, &
            SURF_SPECIFIED_TANGENTIAL_VELOCITY, SURF_SPECIFIED_NORMAL_GRADIENT, SURF_ROUGHNESS, SURF_VEL_T_1, &
            SURF_VEL_T_2, SF_RAMP_VELO_INDEX, SF_RAMP_VELO_TAU, RAMP_T_MIN, RAMP_SPAN, RAMP_RDT, RAMP_NPOINTS, &
            RAMP_OFFSET, RAMP_INTERP, T_IN, WALL_RDN, MU_RSQMW_ZD, RSQ_MW_ZD, DXD, DYD, DZD, VON_KARMAN_CONSTANTD, &
            Y_WERNER_WENGLED, IBAR, &
            JBAR, KBAR, I_MAX_TEMP_LOCAL, SOLID_BOUNDARY_LOCAL, FREE_SLIP_BC_LOCAL, NO_SLIP_BC_LOCAL, WALL_MODEL_BC_LOCAL, &
            INTERPOLATED_BOUNDARY_LOCAL,SIM_MODE_LOCAL,DNS_MODE_LOCAL,IEC_TARGET_LOCAL,EDGE_ACTIVE_LB,EDGE_ACTIVE_UB, &
            CORRECTOR_LOCAL,WIND_BOUNDARY,U_WIND,V_WIND,W_WIND)
            IMPLICIT NONE
            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)

            REAL(EB), DEVICE, DIMENSION(0:IBAR+1,0:JBAR+1,0:*) :: UU,VV,WW,RHOPD,TMPD_IN
            REAL(EB), DEVICE, DIMENSION(0:*) :: EDGE_GHOST_1,EDGE_GHOST_2
            REAL(EB), DEVICE, DIMENSION(0:*) :: EDGE_U_AVG,EDGE_V_AVG,EDGE_W_AVG
            REAL(EB), DEVICE, DIMENSION(EDGE_ACTIVE_LB:EDGE_ACTIVE_UB,-2:*) :: EDGE_OMEGA_IN,EDGE_TAU_IN
            INTEGER, DEVICE, DIMENSION(EDGE_ACTIVE_LB:EDGE_ACTIVE_UB,1:*) :: EDGE_SIDE_ISGN_IN,EDGE_SIDE_IOR_IN
            INTEGER, DEVICE, DIMENSION(EDGE_ACTIVE_LB:EDGE_ACTIVE_UB,1:*) :: EDGE_SIDE_IS2_IN,EDGE_SIDE_ICD_IN
            INTEGER, DEVICE, DIMENSION(EDGE_ACTIVE_LB:EDGE_ACTIVE_UB,1:*) :: EDGE_SIDE_IWM_IN,EDGE_SIDE_IWMI_IN
            INTEGER, DEVICE, DIMENSION(EDGE_ACTIVE_LB:EDGE_ACTIVE_UB,1:*) :: EDGE_SIDE_IWP_IN,EDGE_SIDE_IWPI_IN

            INTEGER , DEVICE, DIMENSION(*)  :: EDGE_AXIS_LIST
            INTEGER , DEVICE, DIMENSION(0:*) :: EDGE_I,EDGE_J,EDGE_K
            INTEGER , DEVICE, DIMENSION(0:*) :: EDGE_CELL_INDEX_MM,EDGE_CELL_INDEX_PM,EDGE_CELL_INDEX_MP,EDGE_CELL_INDEX_PP
            INTEGER , DEVICE, DIMENSION(0:*) :: EDGE_NOM_1,EDGE_NOM_2
            LOGICAL , DEVICE, DIMENSION(0:*) :: EDGE_EXTERNAL

            INTEGER , DEVICE, DIMENSION(0:*) :: CELL_I,CELL_J,CELL_K
            LOGICAL , DEVICE, DIMENSION(0:*) :: CELL_SOLID,CELL_EXTERIOR

            INTEGER , DEVICE, DIMENSION(0:*) :: WALL_SF_INDEX

            INTEGER , DEVICE, DIMENSION(0:*) :: SURF_VELOCITY_BC_INDEX,SF_RAMP_VELO_INDEX,RAMP_NPOINTS,RAMP_OFFSET
            LOGICAL , DEVICE, DIMENSION(0:*) :: SURF_SPECIFIED_NORMAL_VELOCITY, SURF_SPECIFIED_TANGENTIAL_VELOCITY
            LOGICAL , DEVICE, DIMENSION(0:*) :: SURF_SPECIFIED_NORMAL_GRADIENT
            REAL(EB), DEVICE, DIMENSION(0:*) :: SURF_ROUGHNESS,SURF_VEL_T_1,SURF_VEL_T_2,SF_RAMP_VELO_TAU
            REAL(EB), DEVICE, DIMENSION(0:*) :: RAMP_T_MIN,RAMP_SPAN,RAMP_RDT
            REAL(EB), DEVICE, DIMENSION(1:*) :: RAMP_INTERP
            REAL(EB), DEVICE, DIMENSION(0:*) :: WALL_RDN

            REAL(EB), DEVICE, DIMENSION(0:I_MAX_TEMP_LOCAL,1:*) :: MU_RSQMW_ZD
            REAL(EB), DEVICE, DIMENSION(*)    :: RSQ_MW_ZD
            REAL(EB), DEVICE, DIMENSION(0:*)   :: DXD,DYD,DZD

            REAL(EB), DEVICE :: VON_KARMAN_CONSTANTD,Y_WERNER_WENGLED
            REAL(EB), VALUE :: T_IN
            LOGICAL, VALUE :: WIND_BOUNDARY
            REAL(EB), DEVICE, DIMENSION(0:*) :: U_WIND,V_WIND,W_WIND

            INTEGER, VALUE :: IBAR,JBAR,KBAR,I_MAX_TEMP_LOCAL
            INTEGER, VALUE :: SOLID_BOUNDARY_LOCAL,FREE_SLIP_BC_LOCAL,NO_SLIP_BC_LOCAL,WALL_MODEL_BC_LOCAL
            INTEGER, VALUE :: INTERPOLATED_BOUNDARY_LOCAL
            INTEGER, VALUE :: SIM_MODE_LOCAL,DNS_MODE_LOCAL
            INTEGER, VALUE :: NEDGE_AXIS,IEC_TARGET_LOCAL,EDGE_ACTIVE_LB,EDGE_ACTIVE_UB,CORRECTOR_LOCAL

            INTEGER :: IDX,IE,L
            INTEGER :: II,JJ,KK,IEC
            INTEGER :: ICMM,ICPM,ICMP,ICPP
            INTEGER :: I_SGN,IS,IS2,IOR,ICD,ICD_SGN,IVL,SIDE
            INTEGER :: IWM,IWP,IWMI,IWPI
            INTEGER :: IIGM,JJGM,KKGM,IIGP,JJGP,KKGP
            INTEGER :: NOM_ICD
            INTEGER :: SURF_INDEXM,SURF_INDEXP,SURF_INDEX_USE
            INTEGER :: BOUNDARY_TYPE_M,BOUNDARY_TYPE_P
            INTEGER :: VELOCITY_BC_INDEX
            INTEGER :: ITMP
            REAL(EB) :: UUP(2),UUM(2),DXX(2)
            REAL(EB) :: VEL_GAS,VEL_GHOST,VEL_T,RAMP_T
            REAL(EB) :: SLIP_COEF,U_TAU,Y_PLUS
            REAL(EB) :: MUA,MU_WALL,RHO_WALL
            REAL(EB) :: VK_LOCAL,YW_LOCAL
            REAL(EB) :: U_NORM,WIND_NORMAL

            LOGICAL  :: ALL_BLOCKED, OPEN_M, OPEN_P, INTERPOLATED_EDGE, WIND_INFLOW


            IDX = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IDX < 1 .OR. IDX > NEDGE_AXIS) RETURN
            IE = EDGE_AXIS_LIST(IDX)

            IF (CORRECTOR_LOCAL /= 0) THEN
                EDGE_U_AVG(IE) = -1.E6_EB
                EDGE_V_AVG(IE) = -1.E6_EB
                EDGE_W_AVG(IE) = -1.E6_EB
            ENDIF

            DO L=-2,2
                EDGE_OMEGA_IN(IE,L) = -1.E6_EB
                EDGE_TAU_IN(IE,L) = -1.E6_EB
            ENDDO

            ICMM = EDGE_CELL_INDEX_MM(IE)
            ICPM = EDGE_CELL_INDEX_PM(IE)
            ICMP = EDGE_CELL_INDEX_MP(IE)
            ICPP = EDGE_CELL_INDEX_PP(IE)

            ALL_BLOCKED = ( (CELL_SOLID(ICMM) .OR. CELL_EXTERIOR(ICMM)) .AND. &
                            (CELL_SOLID(ICPM) .OR. CELL_EXTERIOR(ICPM)) .AND. &
                            (CELL_SOLID(ICMP) .OR. CELL_EXTERIOR(ICMP)) .AND. &
                            (CELL_SOLID(ICPP) .OR. CELL_EXTERIOR(ICPP)) )
            IF (ALL_BLOCKED) RETURN

            II  = EDGE_I(IE)
            JJ  = EDGE_J(IE)
            KK  = EDGE_K(IE)
            IEC = IEC_TARGET_LOCAL

            SELECT CASE(IEC)
            CASE(1)
                UUP(1) = VV(II,JJ,KK+1)
                UUM(1) = VV(II,JJ,KK)
                UUP(2) = WW(II,JJ+1,KK)
                UUM(2) = WW(II,JJ,KK)
                DXX(1) = DYD(JJ)
                DXX(2) = DZD(KK)
            CASE(2)
                UUP(1) = WW(II+1,JJ,KK)
                UUM(1) = WW(II,JJ,KK)
                UUP(2) = UU(II,JJ,KK+1)
                UUM(2) = UU(II,JJ,KK)
                DXX(1) = DZD(KK)
                DXX(2) = DXD(II)
            CASE(3)
                UUP(1) = UU(II,JJ+1,KK)
                UUM(1) = UU(II,JJ,KK)
                UUP(2) = VV(II+1,JJ,KK)
                UUM(2) = VV(II,JJ,KK)
                DXX(1) = DXD(II)
                DXX(2) = DYD(JJ)
            END SELECT

            VK_LOCAL = VON_KARMAN_CONSTANTD
            YW_LOCAL = Y_WERNER_WENGLED

            DO SIDE = 1, 4
                    I_SGN = EDGE_SIDE_ISGN_IN(IE,SIDE)
                    IOR   = EDGE_SIDE_IOR_IN (IE,SIDE)
                    IS2   = EDGE_SIDE_IS2_IN (IE,SIDE)
                    ICD   = EDGE_SIDE_ICD_IN (IE,SIDE)
                    ICD_SGN = I_SGN*ICD
                    IWM  = EDGE_SIDE_IWM_IN (IE,SIDE)
                    IWMI = EDGE_SIDE_IWMI_IN(IE,SIDE)
                    IWP  = EDGE_SIDE_IWP_IN (IE,SIDE)
                    IWPI = EDGE_SIDE_IWPI_IN(IE,SIDE)

                    IF (IWM == 0 .AND. IWP == 0) CYCLE

                    IF (IWM > 0 .AND. IWMI > 0) THEN
                        IF (WALL_TOPO(IWMI)%BOUNDARY_TYPE == SOLID_BOUNDARY_LOCAL) THEN
                            IF (SURF_VELOCITY_BC_INDEX(WALL_SF_INDEX(IWM)) /= FREE_SLIP_BC_LOCAL) CYCLE
                        ENDIF
                    ENDIF

                    IF (IWP > 0 .AND. IWPI > 0) THEN
                        IF (WALL_TOPO(IWPI)%BOUNDARY_TYPE == SOLID_BOUNDARY_LOCAL) THEN
                            IF (SURF_VELOCITY_BC_INDEX(WALL_SF_INDEX(IWP)) /= FREE_SLIP_BC_LOCAL) CYCLE
                        ENDIF
                    ENDIF

                    OPEN_M = (IWM == 0)
                    OPEN_P = (IWP == 0)
                    IF (IWM > 0) OPEN_M = (WALL_TOPO(IWM)%BOUNDARY_TYPE == OPEN_BOUNDARY)
                    IF (IWP > 0) OPEN_P = (WALL_TOPO(IWP)%BOUNDARY_TYPE == OPEN_BOUNDARY)

                    IF (OPEN_M .AND. OPEN_P) THEN
                        U_NORM = 0._EB
                        WIND_NORMAL = 0._EB
                        IF (WIND_BOUNDARY) THEN
                            SELECT CASE(IEC)
                                CASE(1)
                                    IF (JJ==0    .AND. IOR== 2) U_NORM = 0.5_EB*(VV(II,0,KK)+VV(II,0,KK+1))
                                    IF (JJ==JBAR .AND. IOR==-2) U_NORM = 0.5_EB*(VV(II,JBAR,KK)+VV(II,JBAR,KK+1))
                                    IF (KK==0    .AND. IOR== 3) U_NORM = 0.5_EB*(WW(II,JJ,0)+WW(II,JJ+1,0))
                                    IF (KK==KBAR .AND. IOR==-3) U_NORM = 0.5_EB*(WW(II,JJ,KBAR)+WW(II,JJ+1,KBAR))
                                CASE(2)
                                    IF (II==0    .AND. IOR== 1) U_NORM = 0.5_EB*(UU(0,JJ,KK)+UU(0,JJ,KK+1))
                                    IF (II==IBAR .AND. IOR==-1) U_NORM = 0.5_EB*(UU(IBAR,JJ,KK)+UU(IBAR,JJ,KK+1))
                                    IF (KK==0    .AND. IOR== 3) U_NORM = 0.5_EB*(WW(II,JJ,0)+WW(II+1,JJ,0))
                                    IF (KK==KBAR .AND. IOR==-3) U_NORM = 0.5_EB*(WW(II,JJ,KBAR)+WW(II+1,JJ,KBAR))
                                CASE(3)
                                    IF (II==0    .AND. IOR== 1) U_NORM = 0.5_EB*(UU(0,JJ,KK)+UU(0,JJ+1,KK))
                                    IF (II==IBAR .AND. IOR==-1) U_NORM = 0.5_EB*(UU(IBAR,JJ,KK)+UU(IBAR,JJ+1,KK))
                                    IF (JJ==0    .AND. IOR== 2) U_NORM = 0.5_EB*(VV(II,0,KK)+VV(II+1,0,KK))
                                    IF (JJ==JBAR .AND. IOR==-2) U_NORM = 0.5_EB*(VV(II,JBAR,KK)+VV(II+1,JBAR,KK))
                            END SELECT
                            SELECT CASE(ABS(IOR))
                                CASE(1); WIND_NORMAL = U_WIND(KK)
                                CASE(2); WIND_NORMAL = V_WIND(KK)
                                CASE(3); WIND_NORMAL = W_WIND(KK)
                            END SELECT
                        ENDIF
                        WIND_INFLOW = WIND_BOUNDARY .AND. REAL(IOR,EB)*WIND_NORMAL>=0._EB .AND. REAL(IOR,EB)*U_NORM>=0._EB
                        SELECT CASE(IEC)
                            CASE(1)
                                IF (JJ==0    .AND. IOR== 2) WW(II,0,KK)      = MERGE(W_WIND(KK),WW(II,1,KK),WIND_INFLOW)
                                IF (JJ==JBAR .AND. IOR==-2) WW(II,JBAR+1,KK) = MERGE(W_WIND(KK),WW(II,JBAR,KK),WIND_INFLOW)
                                IF (KK==0    .AND. IOR== 3) VV(II,JJ,0)      = MERGE(V_WIND(KK),VV(II,JJ,1),WIND_INFLOW)
                                IF (KK==KBAR .AND. IOR==-3) VV(II,JJ,KBAR+1) = MERGE(V_WIND(KK),VV(II,JJ,KBAR),WIND_INFLOW)
                            CASE(2)
                                IF (II==0    .AND. IOR== 1) WW(0,JJ,KK)      = MERGE(W_WIND(KK),WW(1,JJ,KK),WIND_INFLOW)
                                IF (II==IBAR .AND. IOR==-1) WW(IBAR+1,JJ,KK) = MERGE(W_WIND(KK),WW(IBAR,JJ,KK),WIND_INFLOW)
                                IF (KK==0    .AND. IOR== 3) UU(II,JJ,0)      = MERGE(U_WIND(KK),UU(II,JJ,1),WIND_INFLOW)
                                IF (KK==KBAR .AND. IOR==-3) UU(II,JJ,KBAR+1) = MERGE(U_WIND(KK),UU(II,JJ,KBAR),WIND_INFLOW)
                            CASE(3)
                                IF (II==0    .AND. IOR== 1) VV(0,JJ,KK)      = MERGE(V_WIND(KK),VV(1,JJ,KK),WIND_INFLOW)
                                IF (II==IBAR .AND. IOR==-1) VV(IBAR+1,JJ,KK) = MERGE(V_WIND(KK),VV(IBAR,JJ,KK),WIND_INFLOW)
                                IF (JJ==0    .AND. IOR== 2) UU(II,0,KK)      = MERGE(U_WIND(KK),UU(II,1,KK),WIND_INFLOW)
                                IF (JJ==JBAR .AND. IOR==-2) UU(II,JBAR+1,KK) = MERGE(U_WIND(KK),UU(II,JBAR,KK),WIND_INFLOW)
                        END SELECT

                        IF (IWM /= 0 .AND. IWP /= 0) RETURN
                        CYCLE
                    ENDIF

                    IF (IWM > 0) THEN
                        BOUNDARY_TYPE_M = WALL_TOPO(IWM)%BOUNDARY_TYPE
                        SURF_INDEXM     = WALL_SF_INDEX(IWM)
                    ELSE
                        BOUNDARY_TYPE_M = WALL_TOPO(IWP)%BOUNDARY_TYPE
                        SURF_INDEXM     = WALL_SF_INDEX(IWP)
                    ENDIF

                    IF (IWP > 0) THEN
                        BOUNDARY_TYPE_P = WALL_TOPO(IWP)%BOUNDARY_TYPE
                        SURF_INDEXP     = WALL_SF_INDEX(IWP)
                    ELSE
                        BOUNDARY_TYPE_P = WALL_TOPO(IWM)%BOUNDARY_TYPE
                        SURF_INDEXP     = WALL_SF_INDEX(IWM)
                    ENDIF

                    IF (ICD == 1) THEN
                        NOM_ICD = EDGE_NOM_1(IE)
                    ELSE
                        NOM_ICD = EDGE_NOM_2(IE)
                    ENDIF

                    INTERPOLATED_EDGE = (NOM_ICD /= 0) .AND. &
                        (BOUNDARY_TYPE_M /= SOLID_BOUNDARY_LOCAL .AND. BOUNDARY_TYPE_P /= SOLID_BOUNDARY_LOCAL) .AND. &
                        (BOUNDARY_TYPE_M == INTERPOLATED_BOUNDARY_LOCAL .OR. BOUNDARY_TYPE_P == INTERPOLATED_BOUNDARY_LOCAL)

                    IF (INTERPOLATED_EDGE) THEN
                        IVL = 3 - ICD
                        IF (IOR < 0) THEN
                            VEL_GAS = UUM(IVL)
                        ELSE
                            VEL_GAS = UUP(IVL)
                        ENDIF
                        IF (ICD == 1) VEL_GHOST = EDGE_GHOST_1(IE)
                        IF (ICD == 2) VEL_GHOST = EDGE_GHOST_2(IE)

                        SELECT CASE(IEC)
                        CASE(1)
                            IF (JJ==0    .AND. IOR== 2) WW(II,JJ,KK)   = VEL_GHOST
                            IF (JJ==JBAR .AND. IOR==-2) WW(II,JJ+1,KK) = VEL_GHOST
                            IF (KK==0    .AND. IOR== 3) VV(II,JJ,KK)   = VEL_GHOST
                            IF (KK==KBAR .AND. IOR==-3) VV(II,JJ,KK+1) = VEL_GHOST
                        CASE(2)
                            IF (II==0    .AND. IOR== 1) WW(II,JJ,KK)   = VEL_GHOST
                            IF (II==IBAR .AND. IOR==-1) WW(II+1,JJ,KK) = VEL_GHOST
                            IF (KK==0    .AND. IOR== 3) UU(II,JJ,KK)   = VEL_GHOST
                            IF (KK==KBAR .AND. IOR==-3) UU(II,JJ,KK+1) = VEL_GHOST
                        CASE(3)
                            IF (II==0    .AND. IOR== 1) VV(II,JJ,KK)   = VEL_GHOST
                            IF (II==IBAR .AND. IOR==-1) VV(II+1,JJ,KK) = VEL_GHOST
                            IF (JJ==0    .AND. IOR== 2) UU(II,JJ,KK)   = VEL_GHOST
                            IF (JJ==JBAR .AND. IOR==-2) UU(II,JJ+1,KK) = VEL_GHOST
                        END SELECT

                        IF (CORRECTOR_LOCAL /= 0) CALL STORE_EDGE_VELOCITY_AVG( &
                            EDGE_U_AVG,EDGE_V_AVG,EDGE_W_AVG,IE,IEC,ICD,VEL_GHOST,VEL_GAS)

                        CYCLE
                    ENDIF

                    IF (SURF_SPECIFIED_NORMAL_VELOCITY(SURF_INDEXM) .OR. SURF_SPECIFIED_NORMAL_GRADIENT(SURF_INDEXM)) THEN
                        SURF_INDEX_USE = SURF_INDEXM
                    ELSEIF (SURF_SPECIFIED_NORMAL_VELOCITY(SURF_INDEXP) .OR. SURF_SPECIFIED_NORMAL_GRADIENT(SURF_INDEXP)) THEN
                        SURF_INDEX_USE = SURF_INDEXP
                    ELSE
                        SURF_INDEX_USE = MAX(SURF_INDEXM,SURF_INDEXP)
                    ENDIF

                    VELOCITY_BC_INDEX = SURF_VELOCITY_BC_INDEX(SURF_INDEX_USE)
                    IF (IWM > 0 .AND. SURF_ROUGHNESS(SURF_INDEX_USE) > 2._EB/WALL_RDN(IWM)) &
                        VELOCITY_BC_INDEX = NO_SLIP_BC_LOCAL

                    IF (VELOCITY_BC_INDEX /= WALL_MODEL_BC_LOCAL .AND. VELOCITY_BC_INDEX /= FREE_SLIP_BC_LOCAL .AND. &
                        VELOCITY_BC_INDEX /= NO_SLIP_BC_LOCAL) CYCLE

                    IVL = 3 - ICD

                    IF (IOR < 0) THEN
                        VEL_GAS = UUM(IVL)
                        IF (ICD == 1) THEN
                            IIGM = CELL_I(ICMM); JJGM = CELL_J(ICMM); KKGM = CELL_K(ICMM)
                            IIGP = CELL_I(ICMP); JJGP = CELL_J(ICMP); KKGP = CELL_K(ICMP)
                        ELSE
                            IIGM = CELL_I(ICMM); JJGM = CELL_J(ICMM); KKGM = CELL_K(ICMM)
                            IIGP = CELL_I(ICPM); JJGP = CELL_J(ICPM); KKGP = CELL_K(ICPM)
                        ENDIF
                    ELSE
                        VEL_GAS = UUP(IVL)
                        IF (ICD == 1) THEN
                            IIGM = CELL_I(ICPM); JJGM = CELL_J(ICPM); KKGM = CELL_K(ICPM)
                        ELSE
                            IIGM = CELL_I(ICMP); JJGM = CELL_J(ICMP); KKGM = CELL_K(ICMP)
                        ENDIF
                        IIGP = CELL_I(ICPP); JJGP = CELL_J(ICPP); KKGP = CELL_K(ICPP)
                    ENDIF

                    VEL_T = 0._EB
                    IF (SURF_SPECIFIED_TANGENTIAL_VELOCITY(SURF_INDEX_USE)) THEN
                        RAMP_T = EVALUATE_RAMP_DEVICE(T_IN,SF_RAMP_VELO_INDEX(SURF_INDEX_USE), &
                            SF_RAMP_VELO_TAU(SURF_INDEX_USE),RAMP_T_MIN,RAMP_SPAN,RAMP_RDT, &
                            RAMP_NPOINTS,RAMP_OFFSET,RAMP_INTERP)
                        IF (IEC==1 .OR. (IEC==2 .AND. ICD==2)) VEL_T = RAMP_T*SURF_VEL_T_2(SURF_INDEX_USE)
                        IF (IEC==3 .OR. (IEC==2 .AND. ICD==1)) VEL_T = RAMP_T*SURF_VEL_T_1(SURF_INDEX_USE)
                    ENDIF
                    MUA = 0.5_EB*(MUD(IIGM,JJGM,KKGM) + MUD(IIGP,JJGP,KKGP))

                    IF (VELOCITY_BC_INDEX == FREE_SLIP_BC_LOCAL) THEN
                        VEL_GHOST = VEL_GAS

                        EDGE_OMEGA_IN(IE,ICD_SGN) = 0._EB
                        EDGE_TAU_IN  (IE,ICD_SGN) = 0._EB

                        SELECT CASE(IEC)
                        CASE(1)
                            IF (JJ==0    .AND. IOR== 2) WW(II,JJ,KK)   = VEL_GHOST
                            IF (JJ==JBAR .AND. IOR==-2) WW(II,JJ+1,KK) = VEL_GHOST
                            IF (KK==0    .AND. IOR== 3) VV(II,JJ,KK)   = VEL_GHOST
                            IF (KK==KBAR .AND. IOR==-3) VV(II,JJ,KK+1) = VEL_GHOST
                        CASE(2)
                            IF (II==0    .AND. IOR== 1) WW(II,JJ,KK)   = VEL_GHOST
                            IF (II==IBAR .AND. IOR==-1) WW(II+1,JJ,KK) = VEL_GHOST
                            IF (KK==0    .AND. IOR== 3) UU(II,JJ,KK)   = VEL_GHOST
                            IF (KK==KBAR .AND. IOR==-3) UU(II,JJ,KK+1) = VEL_GHOST
                        CASE(3)
                            IF (II==0    .AND. IOR== 1) VV(II,JJ,KK)   = VEL_GHOST
                            IF (II==IBAR .AND. IOR==-1) VV(II+1,JJ,KK) = VEL_GHOST
                            IF (JJ==0    .AND. IOR== 2) UU(II,JJ,KK)   = VEL_GHOST
                            IF (JJ==JBAR .AND. IOR==-2) UU(II,JJ+1,KK) = VEL_GHOST
                        END SELECT

                        IF (CORRECTOR_LOCAL /= 0) CALL STORE_EDGE_VELOCITY_AVG( &
                            EDGE_U_AVG,EDGE_V_AVG,EDGE_W_AVG,IE,IEC,ICD,VEL_GHOST,VEL_GAS)

                        CYCLE
                    ENDIF

                    IF (VELOCITY_BC_INDEX == NO_SLIP_BC_LOCAL) THEN
                        VEL_GHOST = 2._EB*VEL_T - VEL_GAS

                        EDGE_OMEGA_IN(IE,ICD_SGN) = I_SGN*(VEL_GAS - VEL_GHOST)/DXX(ICD)
                        EDGE_TAU_IN  (IE,ICD_SGN) = MUA*EDGE_OMEGA_IN(IE,ICD_SGN)

                        SELECT CASE(IEC)
                        CASE(1)
                            IF (JJ==0    .AND. IOR== 2) WW(II,JJ,KK)   = VEL_GHOST
                            IF (JJ==JBAR .AND. IOR==-2) WW(II,JJ+1,KK) = VEL_GHOST
                            IF (KK==0    .AND. IOR== 3) VV(II,JJ,KK)   = VEL_GHOST
                            IF (KK==KBAR .AND. IOR==-3) VV(II,JJ,KK+1) = VEL_GHOST
                        CASE(2)
                            IF (II==0    .AND. IOR== 1) WW(II,JJ,KK)   = VEL_GHOST
                            IF (II==IBAR .AND. IOR==-1) WW(II+1,JJ,KK) = VEL_GHOST
                            IF (KK==0    .AND. IOR== 3) UU(II,JJ,KK)   = VEL_GHOST
                            IF (KK==KBAR .AND. IOR==-3) UU(II,JJ,KK+1) = VEL_GHOST
                        CASE(3)
                            IF (II==0    .AND. IOR== 1) VV(II,JJ,KK)   = VEL_GHOST
                            IF (II==IBAR .AND. IOR==-1) VV(II+1,JJ,KK) = VEL_GHOST
                            IF (JJ==0    .AND. IOR== 2) UU(II,JJ,KK)   = VEL_GHOST
                            IF (JJ==JBAR .AND. IOR==-2) UU(II,JJ+1,KK) = VEL_GHOST
                        END SELECT

                        IF (CORRECTOR_LOCAL /= 0) CALL STORE_EDGE_VELOCITY_AVG( &
                            EDGE_U_AVG,EDGE_V_AVG,EDGE_W_AVG,IE,IEC,ICD,VEL_GHOST,VEL_GAS)

                        CYCLE
                    ENDIF

                    ITMP = MIN(I_MAX_TEMP_LOCAL, NINT(0.5_EB*(TMPD_IN(IIGM,JJGM,KKGM)+TMPD_IN(IIGP,JJGP,KKGP))))
                    MU_WALL  = MU_RSQMW_ZD(ITMP,1) / RSQ_MW_ZD(1)
                    RHO_WALL = 0.5_EB*(RHOPD(IIGM,JJGM,KKGM) + RHOPD(IIGP,JJGP,KKGP))

                    CALL WALL_MODEL_DEVICE( &
                        SLIP_COEF,U_TAU,Y_PLUS, &
                        MU_WALL/RHO_WALL, &
                        SURF_ROUGHNESS(SURF_INDEX_USE), &
                        0.5_EB*DXX(ICD), &
                        VEL_GAS-VEL_T, &
                        SIM_MODE_LOCAL,DNS_MODE_LOCAL, &
                        VK_LOCAL,YW_LOCAL )

                    IF ((IWM==0 .OR. IWP==0) .AND. (.NOT. EDGE_EXTERNAL(IE))) THEN
                        SLIP_COEF = 0._EB
                    ENDIF

                    VEL_GHOST = VEL_T + SLIP_COEF*(VEL_GAS-VEL_T)

                    EDGE_OMEGA_IN(IE,ICD_SGN) = I_SGN*(VEL_GAS - VEL_GHOST)/DXX(ICD)
                    EDGE_TAU_IN  (IE,ICD_SGN) = RHO_WALL*U_TAU**2 * SIGN(1._EB,EDGE_OMEGA_IN(IE,ICD_SGN))

                    SELECT CASE(IEC)
                    CASE(1)
                        IF (JJ==0    .AND. IOR== 2) WW(II,JJ,KK)   = VEL_GHOST
                        IF (JJ==JBAR .AND. IOR==-2) WW(II,JJ+1,KK) = VEL_GHOST
                        IF (KK==0    .AND. IOR== 3) VV(II,JJ,KK)   = VEL_GHOST
                        IF (KK==KBAR .AND. IOR==-3) VV(II,JJ,KK+1) = VEL_GHOST
                    CASE(2)
                        IF (II==0    .AND. IOR== 1) WW(II,JJ,KK)   = VEL_GHOST
                        IF (II==IBAR .AND. IOR==-1) WW(II+1,JJ,KK) = VEL_GHOST
                        IF (KK==0    .AND. IOR== 3) UU(II,JJ,KK)   = VEL_GHOST
                        IF (KK==KBAR .AND. IOR==-3) UU(II,JJ,KK+1) = VEL_GHOST
                    CASE(3)
                        IF (II==0    .AND. IOR== 1) VV(II,JJ,KK)   = VEL_GHOST
                        IF (II==IBAR .AND. IOR==-1) VV(II+1,JJ,KK) = VEL_GHOST
                        IF (JJ==0    .AND. IOR== 2) UU(II,JJ,KK)   = VEL_GHOST
                        IF (JJ==JBAR .AND. IOR==-2) UU(II,JJ+1,KK) = VEL_GHOST
                    END SELECT

                    IF (CORRECTOR_LOCAL /= 0) CALL STORE_EDGE_VELOCITY_AVG( &
                        EDGE_U_AVG,EDGE_V_AVG,EDGE_W_AVG,IE,IEC,ICD,VEL_GHOST,VEL_GAS)


            ENDDO

        END SUBROUTINE VELOCITY_BC_GPU_GHOST_KERNEL_AXIS

        ATTRIBUTES(GLOBAL) SUBROUTINE VELOCITY_BC_GPU_TAU_OMEGA_KERNEL( &
            EDGE_OMEGA_IN,EDGE_TAU_IN,EDGE_STORAGE_LB,EDGE_STORAGE_UB,EDGE_ACTIVE_LB,EDGE_ACTIVE_UB)
            IMPLICIT NONE

            INTEGER, VALUE :: EDGE_STORAGE_LB,EDGE_STORAGE_UB,EDGE_ACTIVE_LB,EDGE_ACTIVE_UB
            REAL(EB), DEVICE, DIMENSION(EDGE_STORAGE_LB:EDGE_STORAGE_UB,-2:*) :: EDGE_OMEGA_IN,EDGE_TAU_IN

            INTEGER :: IE, IDX
            INTEGER :: I_SGN,ICD,ICDO,ICD_SGN,ICDO_SGN
            REAL(EB) :: DUIDXJ_USE(2),MU_DUIDXJ_USE(2)
            REAL(EB) :: OMEGA_STAGE(-2:2),TAU_STAGE(-2:2)
            LOGICAL  :: ALTERED_GRADIENT(-2:2)

            IE = EDGE_ACTIVE_LB + (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X - 1
            IF (IE < EDGE_ACTIVE_LB .OR. IE > EDGE_ACTIVE_UB) RETURN

            OMEGA_STAGE(:) = EDGE_OMEGA_IN(IE,:)
            TAU_STAGE(:)   = EDGE_TAU_IN  (IE,:)
            ALTERED_GRADIENT(:) = (OMEGA_STAGE(:) /= -1.E6_EB) .OR. (TAU_STAGE(:) /= -1.E6_EB)

            DO IDX = -2, 2
                IF (ALTERED_GRADIENT(IDX)) EXIT
            ENDDO
            IF (IDX > 2) RETURN

            DO I_SGN = -1, 1, 2
                DO ICD = 1, 2

                    ICDO = 3 - ICD

                    ICD_SGN = I_SGN*ICD

                    IF (ALTERED_GRADIENT(ICD_SGN)) THEN
                        DUIDXJ_USE(ICD)    = OMEGA_STAGE(ICD_SGN)
                        MU_DUIDXJ_USE(ICD) = TAU_STAGE(ICD_SGN)
                    ELSEIF (ALTERED_GRADIENT(-ICD_SGN)) THEN
                        DUIDXJ_USE(ICD)    = OMEGA_STAGE(-ICD_SGN)
                        MU_DUIDXJ_USE(ICD) = TAU_STAGE(-ICD_SGN)
                    ELSE
                        CYCLE
                    ENDIF

                    ICDO_SGN = I_SGN*ICDO

                    IF (ALTERED_GRADIENT(ICDO_SGN)) THEN
                        DUIDXJ_USE(ICDO)    = OMEGA_STAGE(ICDO_SGN)
                        MU_DUIDXJ_USE(ICDO) = TAU_STAGE(ICDO_SGN)
                    ELSEIF (ALTERED_GRADIENT(-ICDO_SGN)) THEN
                        DUIDXJ_USE(ICDO)    = OMEGA_STAGE(-ICDO_SGN)
                        MU_DUIDXJ_USE(ICDO) = TAU_STAGE(-ICDO_SGN)
                    ELSE
                        DUIDXJ_USE(ICDO)    = 0._EB
                        MU_DUIDXJ_USE(ICDO) = 0._EB
                    ENDIF

                    EDGE_OMEGA_IN(IE,ICD_SGN) = DUIDXJ_USE(1) - DUIDXJ_USE(2)
                    EDGE_TAU_IN  (IE,ICD_SGN) = MU_DUIDXJ_USE(1) + MU_DUIDXJ_USE(2)

                ENDDO
            ENDDO
        END SUBROUTINE VELOCITY_BC_GPU_TAU_OMEGA_KERNEL

        !--------------------------------------------------------------------------
        ! Interpolated-boundary viscosity exchange
        !--------------------------------------------------------------------------
        SUBROUTINE VISCOSITY_BC_GPU(NM,APPLY_TO_ESTIMATED_VARIABLES)
            IMPLICIT NONE

            INTEGER, INTENT(IN) :: NM
            LOGICAL, INTENT(IN) :: APPLY_TO_ESTIMATED_VARIABLES
            TYPE(OMESH_TYPE), POINTER :: OM
            TYPE(GPU_OMESH_TYPE), POINTER :: GOM_LOCAL
            INTEGER :: BLKD_1D,GRDD_1D,NOM
            REAL(EB), DEVICE, POINTER, CONTIGUOUS :: OM_MU(:,:,:),OM_KRES(:,:,:),OM_DP(:,:,:)
            REAL(EB) :: TNOW

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)

            IF (N_EXTERNAL_WALL_CELLS <= 0) RETURN

            TNOW = CURRENT_TIME()

            IF (GPU_DEBUG_H2D_D2H) THEN
                GM%MUD = MU
                GM%KRESD = KRES
                IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
                    GM%DS_D = DS
                ELSE
                    GM%D_D = D
                ENDIF
            ENDIF

            BLKD_1D = GPU_BLKD_1D
            GRDD_1D = GM%GRDD_EXTERNAL_WALL_1D

            DO NOM=1,NMESHES
                IF (.NOT.ALLOCATED(GM%OMESH)) EXIT
                GOM_LOCAL => GM%OMESH(NOM)
                IF (.NOT.GOM_LOCAL%ACTIVE) CYCLE
                OM => MESHES(NM)%OMESH(NOM)

                IF (GPU_DEBUG_H2D_D2H) THEN
                    GOM_LOCAL%MU_D = OM%MU
                    GOM_LOCAL%KRES_D = OM%KRES
                    IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
                        GOM_LOCAL%DS_D = OM%DS
                    ELSE
                        GOM_LOCAL%D_D = OM%D
                    ENDIF
                ENDIF

                OM_MU => GOM_LOCAL%MU_D
                OM_KRES => GOM_LOCAL%KRES_D
                IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
                    OM_DP => GOM_LOCAL%DS_D
                    CALL VISCOSITY_BC_EXTERNAL_WALL_KERNEL<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>(GM%MUD,GM%KRESD,GM%DS_D,OM_MU,OM_KRES,OM_DP, &
                        GM%WALL_TOPO_AOS,GM%EXTERNAL_WALL_NOMD,GM%WALL_IIO,GM%WALL_JJO,GM%WALL_KKO,GM%WALL_IIO_MAX,GM%WALL_JJO_MAX,GM%WALL_KKO_MAX, &
                        N_EXTERNAL_WALL_CELLS,NOM,LBOUND(OM_MU,1),LBOUND(OM_MU,2),LBOUND(OM_MU,3), &
                        IBP1,JBP1,SIZE(OM_MU,1),SIZE(OM_MU,2))
                ELSE
                    OM_DP => GOM_LOCAL%D_D
                    CALL VISCOSITY_BC_EXTERNAL_WALL_KERNEL<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>(GM%MUD,GM%KRESD,GM%D_D,OM_MU,OM_KRES,OM_DP, &
                        GM%WALL_TOPO_AOS,GM%EXTERNAL_WALL_NOMD,GM%WALL_IIO,GM%WALL_JJO,GM%WALL_KKO,GM%WALL_IIO_MAX,GM%WALL_JJO_MAX,GM%WALL_KKO_MAX, &
                        N_EXTERNAL_WALL_CELLS,NOM,LBOUND(OM_MU,1),LBOUND(OM_MU,2),LBOUND(OM_MU,3), &
                        IBP1,JBP1,SIZE(OM_MU,1),SIZE(OM_MU,2))
                ENDIF
            ENDDO

            IF (GPU_DEBUG_H2D_D2H) THEN
                CALL GPU_WAIT(GM%STRM)
                MU = GM%MUD
                KRES = GM%KRESD
                IF (APPLY_TO_ESTIMATED_VARIABLES) THEN
                    DS = GM%DS_D
                ELSE
                    D = GM%D_D
                ENDIF
            ENDIF

            IF (ALLOCATED(T_USED)) T_USED(4) = T_USED(4) + CURRENT_TIME() - TNOW

        END SUBROUTINE VISCOSITY_BC_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE VISCOSITY_BC_EXTERNAL_WALL_KERNEL( &
            MU_OUT,KRES_OUT,DP_OUT,OM_MU,OM_KRES,OM_DP,WALL_TOPO, &
            EXT_NOMD,EXT_IIOD,EXT_JJOD,EXT_KKOD,EXT_IIOD_MAX,EXT_JJOD_MAX,EXT_KKOD_MAX, &
            N_EXTERNAL_WALL_CELLS,NOM_FILTER,OM_LB_I,OM_LB_J,OM_LB_K,IBP1,JBP1,OM_NI,OM_NJ)
            IMPLICIT NONE

            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,NOM_FILTER,OM_LB_I,OM_LB_J,OM_LB_K,IBP1,JBP1,OM_NI,OM_NJ
            REAL(EB), DEVICE, DIMENSION(0:IBP1,0:JBP1,0:*) :: MU_OUT,KRES_OUT,DP_OUT
            REAL(EB), DEVICE, DIMENSION(0:OM_NI-1,0:OM_NJ-1,0:*) :: OM_MU,OM_KRES,OM_DP
            INTEGER, DEVICE, DIMENSION(1:*) :: EXT_NOMD
            INTEGER, DEVICE, DIMENSION(1:*) :: EXT_IIOD,EXT_JJOD,EXT_KKOD
            INTEGER, DEVICE, DIMENSION(1:*) :: EXT_IIOD_MAX,EXT_JJOD_MAX,EXT_KKOD_MAX
            INTEGER :: IW,II,JJ,KK,IIO,JJO,KKO,IIO1,IIO2,JJO1,JJO2,KKO1,KKO2,N_INT_CELLS
            REAL(EB) :: MU_OTHER,KRES_OTHER,DP_OTHER

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS) RETURN
            IF (EXT_NOMD(IW) /= NOM_FILTER) RETURN

            II = WALL_TOPO(IW)%II
            JJ = WALL_TOPO(IW)%JJ
            KK = WALL_TOPO(IW)%KK
            IIO1 = EXT_IIOD(IW)
            JJO1 = EXT_JJOD(IW)
            KKO1 = EXT_KKOD(IW)
            IIO2 = EXT_IIOD_MAX(IW)
            JJO2 = EXT_JJOD_MAX(IW)
            KKO2 = EXT_KKOD_MAX(IW)

            MU_OTHER = 0._EB
            KRES_OTHER = 0._EB
            DP_OTHER = 0._EB
            DO KKO=KKO1,KKO2
                DO JJO=JJO1,JJO2
                    DO IIO=IIO1,IIO2
                        MU_OTHER = MU_OTHER + OM_MU(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)
                        KRES_OTHER = KRES_OTHER + OM_KRES(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)
                        DP_OTHER = DP_OTHER + OM_DP(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)
                    ENDDO
                ENDDO
            ENDDO

            N_INT_CELLS = (IIO2-IIO1+1)*(JJO2-JJO1+1)*(KKO2-KKO1+1)
            IF (N_INT_CELLS <= 0) RETURN
            MU_OUT(II,JJ,KK) = MU_OTHER/REAL(N_INT_CELLS,EB)
            KRES_OUT(II,JJ,KK) = KRES_OTHER/REAL(N_INT_CELLS,EB)
            DP_OUT(II,JJ,KK) = DP_OTHER/REAL(N_INT_CELLS,EB)

        END SUBROUTINE VISCOSITY_BC_EXTERNAL_WALL_KERNEL
        !--------------------------------------------------------------------------
        ! Interpolated-boundary momentum-flux matching
        !--------------------------------------------------------------------------
        SUBROUTINE MATCH_VELOCITY_FLUX_GPU(NM)
            IMPLICIT NONE

            INTEGER, INTENT(IN) :: NM
            INTEGER :: NOM
            TYPE(OMESH_TYPE), POINTER :: OM
            TYPE(GPU_OMESH_TYPE), POINTER :: GOM_LOCAL
            REAL(EB), DEVICE, POINTER, CONTIGUOUS :: OM_FVX(:,:,:),OM_FVY(:,:,:),OM_FVZ(:,:,:)
            REAL(EB), DEVICE, POINTER, CONTIGUOUS :: DX_OTHER(:),DY_OTHER(:),DZ_OTHER(:)
            REAL(EB) :: TNOW

            IF (NMESHES==1) RETURN
            IF (SOLID_PHASE_ONLY) RETURN
            IF (CC_IBM) RETURN

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)

            IF (N_EXTERNAL_WALL_CELLS <= 0) RETURN

            TNOW = CURRENT_TIME()

            IF (GPU_DEBUG_H2D_D2H) THEN
                GM%FVXD = FVX
                GM%FVYD = FVY
                GM%FVZD = FVZ
            ENDIF

            DO NOM=1,NMESHES
                IF (.NOT.ALLOCATED(GM%OMESH)) EXIT
                GOM_LOCAL => GM%OMESH(NOM)
                IF (.NOT.GOM_LOCAL%ACTIVE) CYCLE
                OM => MESHES(NM)%OMESH(NOM)

                IF (GPU_DEBUG_H2D_D2H) THEN
                    GOM_LOCAL%FVX_D = OM%FVX
                    GOM_LOCAL%FVY_D = OM%FVY
                    GOM_LOCAL%FVZ_D = OM%FVZ
                ENDIF

                OM_FVX => GOM_LOCAL%FVX_D
                OM_FVY => GOM_LOCAL%FVY_D
                OM_FVZ => GOM_LOCAL%FVZ_D
                DX_OTHER => GPU_MESHES(NOM)%DXD
                DY_OTHER => GPU_MESHES(NOM)%DYD
                DZ_OTHER => GPU_MESHES(NOM)%DZD

                CALL MATCH_VELOCITY_FLUX_KERNEL<<<GM%GRDD_EXTERNAL_WALL_1D,GPU_BLKD_1D,0,GM%STRM>>>( &
                    GM%FVXD,GM%FVYD,GM%FVZD,OM_FVX,OM_FVY,OM_FVZ, &
                    GM%WALL_TOPO_AOS,GM%EXTERNAL_WALL_NOMD,GM%WALL_IIO,GM%WALL_JJO,GM%WALL_KKO,GM%WALL_IIO_MAX,GM%WALL_JJO_MAX,GM%WALL_KKO_MAX, &
                    DX_OTHER,DY_OTHER,DZ_OTHER,N_EXTERNAL_WALL_CELLS,IBAR,JBAR,KBAR,INTERPOLATED_BOUNDARY,NOM, &
                    LBOUND(OM_FVX,1),LBOUND(OM_FVX,2),LBOUND(OM_FVX,3),SIZE(OM_FVX,1),SIZE(OM_FVX,2))
            ENDDO

            IF (GPU_DEBUG_H2D_D2H) THEN
                CALL GPU_WAIT(GM%STRM)
                FVX = GM%FVXD
                FVY = GM%FVYD
                FVZ = GM%FVZD
            ENDIF

            IF (ALLOCATED(T_USED)) T_USED(4) = T_USED(4) + CURRENT_TIME() - TNOW

        END SUBROUTINE MATCH_VELOCITY_FLUX_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE MATCH_VELOCITY_FLUX_KERNEL( &
            FVX,FVY,FVZ,OM_FVX,OM_FVY,OM_FVZ,WALL_TOPO, &
            EXT_NOMD,WALL_IIOD,WALL_JJOD,WALL_KKOD,WALL_IIOD_MAX,WALL_JJOD_MAX,WALL_KKOD_MAX, &
            DX_OTHER,DY_OTHER,DZ_OTHER,N_EXTERNAL_WALL_CELLS,IBAR_LOCAL,JBAR_LOCAL,KBAR_LOCAL, &
            INTERPOLATED_BOUNDARY_VALUE,NOM_FILTER,OM_LB_I,OM_LB_J,OM_LB_K,OM_NI,OM_NJ)
            IMPLICIT NONE

            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,IBAR_LOCAL,JBAR_LOCAL,KBAR_LOCAL,INTERPOLATED_BOUNDARY_VALUE
            INTEGER, VALUE :: NOM_FILTER,OM_LB_I,OM_LB_J,OM_LB_K,OM_NI,OM_NJ
            REAL(EB), DEVICE, DIMENSION(0:IBAR_LOCAL+1,0:JBAR_LOCAL+1,0:*) :: FVX,FVY,FVZ
            REAL(EB), DEVICE, DIMENSION(0:OM_NI-1,0:OM_NJ-1,0:*) :: OM_FVX,OM_FVY,OM_FVZ
            REAL(EB), DEVICE, DIMENSION(0:*) :: DX_OTHER,DY_OTHER,DZ_OTHER
            INTEGER, DEVICE, DIMENSION(1:*) :: EXT_NOMD
            INTEGER, DEVICE, DIMENSION(1:*) :: WALL_IIOD,WALL_JJOD,WALL_KKOD
            INTEGER, DEVICE, DIMENSION(1:*) :: WALL_IIOD_MAX,WALL_JJOD_MAX,WALL_KKOD_MAX
            INTEGER :: IW,II,JJ,KK,IOR,IIO,JJO,KKO,BT,IIO1,IIO2,JJO1,JJO2,KKO1,KKO2
            REAL(EB) :: DA_OTHER, OTHER_VALUE, WEIGHT

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS) RETURN
            IF (EXT_NOMD(IW) /= NOM_FILTER) RETURN

            BT = WALL_TOPO(IW)%BOUNDARY_TYPE
            IF (BT /= INTERPOLATED_BOUNDARY_VALUE) RETURN

            II  = WALL_TOPO(IW)%II
            JJ  = WALL_TOPO(IW)%JJ
            KK  = WALL_TOPO(IW)%KK
            IOR = WALL_TOPO(IW)%IOR
            IIO1 = WALL_IIOD(IW)
            JJO1 = WALL_JJOD(IW)
            KKO1 = WALL_KKOD(IW)
            IIO2 = WALL_IIOD_MAX(IW)
            JJO2 = WALL_JJOD_MAX(IW)
            KKO2 = WALL_KKOD_MAX(IW)

            SELECT CASE(IOR)
                CASE( 1)
                    DA_OTHER = 0._EB
                    OTHER_VALUE = 0._EB
                    DO KKO=KKO1,KKO2
                        DO JJO=JJO1,JJO2
                            DO IIO=IIO1,IIO2
                                WEIGHT = DY_OTHER(JJO)*DZ_OTHER(KKO)
                                DA_OTHER = DA_OTHER + WEIGHT
                                OTHER_VALUE = OTHER_VALUE + OM_FVX(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)*WEIGHT
                            ENDDO
                        ENDDO
                    ENDDO
                    IF (DA_OTHER > 0._EB) OTHER_VALUE = OTHER_VALUE/DA_OTHER
                    FVX(0,JJ,KK) = 0.5_EB*(FVX(0,JJ,KK) + OTHER_VALUE)
                CASE(-1)
                    DA_OTHER = 0._EB
                    OTHER_VALUE = 0._EB
                    DO KKO=KKO1,KKO2
                        DO JJO=JJO1,JJO2
                            DO IIO=IIO1,IIO2
                                WEIGHT = DY_OTHER(JJO)*DZ_OTHER(KKO)
                                DA_OTHER = DA_OTHER + WEIGHT
                                OTHER_VALUE = OTHER_VALUE + OM_FVX(IIO-1-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)*WEIGHT
                            ENDDO
                        ENDDO
                    ENDDO
                    IF (DA_OTHER > 0._EB) OTHER_VALUE = OTHER_VALUE/DA_OTHER
                    FVX(IBAR_LOCAL,JJ,KK) = 0.5_EB*(FVX(IBAR_LOCAL,JJ,KK) + OTHER_VALUE)
                CASE( 2)
                    DA_OTHER = 0._EB
                    OTHER_VALUE = 0._EB
                    DO KKO=KKO1,KKO2
                        DO JJO=JJO1,JJO2
                            DO IIO=IIO1,IIO2
                                WEIGHT = DX_OTHER(IIO)*DZ_OTHER(KKO)
                                DA_OTHER = DA_OTHER + WEIGHT
                                OTHER_VALUE = OTHER_VALUE + OM_FVY(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)*WEIGHT
                            ENDDO
                        ENDDO
                    ENDDO
                    IF (DA_OTHER > 0._EB) OTHER_VALUE = OTHER_VALUE/DA_OTHER
                    FVY(II,0,KK) = 0.5_EB*(FVY(II,0,KK) + OTHER_VALUE)
                CASE(-2)
                    DA_OTHER = 0._EB
                    OTHER_VALUE = 0._EB
                    DO KKO=KKO1,KKO2
                        DO JJO=JJO1,JJO2
                            DO IIO=IIO1,IIO2
                                WEIGHT = DX_OTHER(IIO)*DZ_OTHER(KKO)
                                DA_OTHER = DA_OTHER + WEIGHT
                                OTHER_VALUE = OTHER_VALUE + OM_FVY(IIO-OM_LB_I,JJO-1-OM_LB_J,KKO-OM_LB_K)*WEIGHT
                            ENDDO
                        ENDDO
                    ENDDO
                    IF (DA_OTHER > 0._EB) OTHER_VALUE = OTHER_VALUE/DA_OTHER
                    FVY(II,JBAR_LOCAL,KK) = 0.5_EB*(FVY(II,JBAR_LOCAL,KK) + OTHER_VALUE)
                CASE( 3)
                    DA_OTHER = 0._EB
                    OTHER_VALUE = 0._EB
                    DO KKO=KKO1,KKO2
                        DO JJO=JJO1,JJO2
                            DO IIO=IIO1,IIO2
                                WEIGHT = DX_OTHER(IIO)*DY_OTHER(JJO)
                                DA_OTHER = DA_OTHER + WEIGHT
                                OTHER_VALUE = OTHER_VALUE + OM_FVZ(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)*WEIGHT
                            ENDDO
                        ENDDO
                    ENDDO
                    IF (DA_OTHER > 0._EB) OTHER_VALUE = OTHER_VALUE/DA_OTHER
                    FVZ(II,JJ,0) = 0.5_EB*(FVZ(II,JJ,0) + OTHER_VALUE)
                CASE(-3)
                    DA_OTHER = 0._EB
                    OTHER_VALUE = 0._EB
                    DO KKO=KKO1,KKO2
                        DO JJO=JJO1,JJO2
                            DO IIO=IIO1,IIO2
                                WEIGHT = DX_OTHER(IIO)*DY_OTHER(JJO)
                                DA_OTHER = DA_OTHER + WEIGHT
                                OTHER_VALUE = OTHER_VALUE + OM_FVZ(IIO-OM_LB_I,JJO-OM_LB_J,KKO-1-OM_LB_K)*WEIGHT
                            ENDDO
                        ENDDO
                    ENDDO
                    IF (DA_OTHER > 0._EB) OTHER_VALUE = OTHER_VALUE/DA_OTHER
                    FVZ(II,JJ,KBAR_LOCAL) = 0.5_EB*(FVZ(II,JJ,KBAR_LOCAL) + OTHER_VALUE)
                CASE DEFAULT
                    RETURN
            END SELECT

        END SUBROUTINE MATCH_VELOCITY_FLUX_KERNEL
        !--------------------------------------------------------------------------
        ! Interpolated-boundary normal-velocity matching
        !--------------------------------------------------------------------------
        SUBROUTINE MATCH_VELOCITY_GPU(NM)
            IMPLICIT NONE

            INTEGER, INTENT(IN) :: NM
            TYPE(OMESH_TYPE), POINTER :: OM
            TYPE(GPU_OMESH_TYPE), POINTER :: GOM_LOCAL
            INTEGER :: BLKD_1D, GRDD_1D, ISTAT, IW, NOM
            INTEGER, ALLOCATABLE :: BOUNDARY_TYPE_PREVIOUS_H(:)
            REAL(EB), DEVICE, POINTER, CONTIGUOUS :: OM_UU(:,:,:),OM_VV(:,:,:),OM_WW(:,:,:)
            REAL(EB), DEVICE, POINTER, CONTIGUOUS :: DX_OTHER(:),DY_OTHER(:),DZ_OTHER(:)
            REAL(EB) :: TNOW

            IF (SOLID_PHASE_ONLY) RETURN
            IF (CC_IBM) RETURN

            CALL POINT_TO_MESH(NM)
            GM => GPU_MESHES(NM)

            IF (N_EXTERNAL_WALL_CELLS <= 0) RETURN

            TNOW = CURRENT_TIME()

            BLKD_1D = GPU_BLKD_1D
            GRDD_1D = GM%GRDD_EXTERNAL_WALL_1D

            IF (PREDICTOR) THEN
                IF (GPU_DEBUG_H2D_D2H) THEN
                    GM%USD = US
                    GM%VSD = VS
                    GM%WSD = WS
                    GM%U_GHOSTD = U_GHOST
                    GM%V_GHOSTD = V_GHOST
                    GM%W_GHOSTD = W_GHOST
                    GM%UVW_SAVED = UVW_SAVE
                ENDIF

                DO NOM=1,NMESHES
                    IF (.NOT.ALLOCATED(GM%OMESH)) EXIT
                    GOM_LOCAL => GM%OMESH(NOM)
                    IF (.NOT.GOM_LOCAL%ACTIVE) CYCLE
                    OM => MESHES(NM)%OMESH(NOM)
                    IF (GPU_DEBUG_H2D_D2H) THEN
                        GOM_LOCAL%US_D = OM%US
                        GOM_LOCAL%VS_D = OM%VS
                        GOM_LOCAL%WS_D = OM%WS
                    ENDIF
                    OM_UU => GOM_LOCAL%US_D
                    OM_VV => GOM_LOCAL%VS_D
                    OM_WW => GOM_LOCAL%WS_D
                    DX_OTHER => GPU_MESHES(NOM)%DXD
                    DY_OTHER => GPU_MESHES(NOM)%DYD
                    DZ_OTHER => GPU_MESHES(NOM)%DZD
                    CALL MATCH_VELOCITY_KERNEL<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>(GM%USD,GM%VSD,GM%WSD,OM_UU,OM_VV,OM_WW,GM%U_GHOSTD,GM%V_GHOSTD,GM%W_GHOSTD, &
                        GM%UVW_SAVED,GM%WALL_BOUNDARY_TYPE_PREVIOUSD,GM%WALL_TOPO_AOS,GM%EXTERNAL_WALL_NOMD,GM%WALL_IIO,GM%WALL_JJO,GM%WALL_KKO, &
                        GM%WALL_IIO_MAX,GM%WALL_JJO_MAX,GM%WALL_KKO_MAX,GM%EXTERNAL_WALL_AREA_RATIOD,DX_OTHER,DY_OTHER,DZ_OTHER, &
                        N_EXTERNAL_WALL_CELLS,IBAR,JBAR,KBAR,INTERPOLATED_BOUNDARY,NOM, &
                        LBOUND(OM_UU,1),LBOUND(OM_UU,2),LBOUND(OM_UU,3),SIZE(OM_UU,1),SIZE(OM_UU,2))
                ENDDO

                IF (GPU_DEBUG_H2D_D2H) THEN
                    CALL GPU_WAIT(GM%STRM)
                    US = GM%USD
                    VS = GM%VSD
                    WS = GM%WSD
                    DO NOM=1,NMESHES
                        IF (.NOT.ALLOCATED(GM%OMESH)) EXIT
                        GOM_LOCAL => GM%OMESH(NOM)
                        IF (.NOT.GOM_LOCAL%ACTIVE) CYCLE
                        OM => MESHES(NM)%OMESH(NOM)
                        OM%US = GOM_LOCAL%US_D
                        OM%VS = GOM_LOCAL%VS_D
                        OM%WS = GOM_LOCAL%WS_D
                    ENDDO
                ENDIF
            ELSE
                IF (GPU_DEBUG_H2D_D2H) THEN
                    GM%UD = U
                    GM%VD = V
                    GM%WD = W
                    GM%U_GHOSTD = U_GHOST
                    GM%V_GHOSTD = V_GHOST
                    GM%W_GHOSTD = W_GHOST
                    GM%UVW_SAVED = UVW_SAVE
                ENDIF

                DO NOM=1,NMESHES
                    IF (.NOT.ALLOCATED(GM%OMESH)) EXIT
                    GOM_LOCAL => GM%OMESH(NOM)
                    IF (.NOT.GOM_LOCAL%ACTIVE) CYCLE
                    OM => MESHES(NM)%OMESH(NOM)
                    IF (GPU_DEBUG_H2D_D2H) THEN
                        GOM_LOCAL%U_D = OM%U
                        GOM_LOCAL%V_D = OM%V
                        GOM_LOCAL%W_D = OM%W
                    ENDIF
                    OM_UU => GOM_LOCAL%U_D
                    OM_VV => GOM_LOCAL%V_D
                    OM_WW => GOM_LOCAL%W_D
                    DX_OTHER => GPU_MESHES(NOM)%DXD
                    DY_OTHER => GPU_MESHES(NOM)%DYD
                    DZ_OTHER => GPU_MESHES(NOM)%DZD
                    CALL MATCH_VELOCITY_KERNEL<<<GRDD_1D,BLKD_1D,0,GM%STRM>>>(GM%UD,GM%VD,GM%WD,OM_UU,OM_VV,OM_WW,GM%U_GHOSTD,GM%V_GHOSTD,GM%W_GHOSTD, &
                        GM%UVW_SAVED,GM%WALL_BOUNDARY_TYPE_PREVIOUSD,GM%WALL_TOPO_AOS,GM%EXTERNAL_WALL_NOMD,GM%WALL_IIO,GM%WALL_JJO,GM%WALL_KKO, &
                        GM%WALL_IIO_MAX,GM%WALL_JJO_MAX,GM%WALL_KKO_MAX,GM%EXTERNAL_WALL_AREA_RATIOD,DX_OTHER,DY_OTHER,DZ_OTHER, &
                        N_EXTERNAL_WALL_CELLS,IBAR,JBAR,KBAR,INTERPOLATED_BOUNDARY,NOM, &
                        LBOUND(OM_UU,1),LBOUND(OM_UU,2),LBOUND(OM_UU,3),SIZE(OM_UU,1),SIZE(OM_UU,2))
                ENDDO

                IF (GPU_DEBUG_H2D_D2H) THEN
                    CALL GPU_WAIT(GM%STRM)
                    U = GM%UD
                    V = GM%VD
                    W = GM%WD
                    DO NOM=1,NMESHES
                        IF (.NOT.ALLOCATED(GM%OMESH)) EXIT
                        GOM_LOCAL => GM%OMESH(NOM)
                        IF (.NOT.GOM_LOCAL%ACTIVE) CYCLE
                        OM => MESHES(NM)%OMESH(NOM)
                        OM%U = GOM_LOCAL%U_D
                        OM%V = GOM_LOCAL%V_D
                        OM%W = GOM_LOCAL%W_D
                    ENDDO
                ENDIF
            ENDIF

            IF (GPU_DEBUG_H2D_D2H) THEN
                U_GHOST = GM%U_GHOSTD
                V_GHOST = GM%V_GHOSTD
                W_GHOST = GM%W_GHOSTD
                UVW_SAVE = GM%UVW_SAVED
                ALLOCATE(BOUNDARY_TYPE_PREVIOUS_H(1:N_EXTERNAL_WALL_CELLS),STAT=ISTAT)
                BOUNDARY_TYPE_PREVIOUS_H = GM%WALL_BOUNDARY_TYPE_PREVIOUSD(1:N_EXTERNAL_WALL_CELLS)
                DO IW=1,N_EXTERNAL_WALL_CELLS
                    EXTERNAL_WALL(IW)%BOUNDARY_TYPE_PREVIOUS = BOUNDARY_TYPE_PREVIOUS_H(IW)
                ENDDO
                DEALLOCATE(BOUNDARY_TYPE_PREVIOUS_H,STAT=ISTAT)
            ENDIF

            IF (ALLOCATED(T_USED)) T_USED(4) = T_USED(4) + CURRENT_TIME() - TNOW

        END SUBROUTINE MATCH_VELOCITY_GPU

        ATTRIBUTES(GLOBAL) SUBROUTINE MATCH_VELOCITY_KERNEL( &
            UU,VV,WW,OM_UU,OM_VV,OM_WW,U_GHOST,V_GHOST,W_GHOST,UVW_SAVE_D,BOUNDARY_TYPE_PREVIOUS, &
            WALL_TOPO,EXT_NOMD,WALL_IIOD,WALL_JJOD,WALL_KKOD, &
            WALL_IIOD_MAX,WALL_JJOD_MAX,WALL_KKOD_MAX,EXTERNAL_WALL_AREA_RATIO,DX_OTHER,DY_OTHER,DZ_OTHER, &
            N_EXTERNAL_WALL_CELLS,IBAR_LOCAL,JBAR_LOCAL,KBAR_LOCAL,INTERPOLATED_BOUNDARY_VALUE, &
            NOM_FILTER,OM_LB_I,OM_LB_J,OM_LB_K,OM_NI,OM_NJ)
            IMPLICIT NONE

            TYPE(WALL_TOPO_AOS_TYPE), DEVICE :: WALL_TOPO(0:*)
            INTEGER, VALUE :: N_EXTERNAL_WALL_CELLS,IBAR_LOCAL,JBAR_LOCAL,KBAR_LOCAL,INTERPOLATED_BOUNDARY_VALUE
            REAL(EB), DEVICE, DIMENSION(0:IBAR_LOCAL+1,0:JBAR_LOCAL+1,0:*) :: UU,VV,WW
            REAL(EB), DEVICE, DIMENSION(0:OM_NI-1,0:OM_NJ-1,0:*) :: OM_UU,OM_VV,OM_WW
            REAL(EB), DEVICE, DIMENSION(1:*) :: U_GHOST,V_GHOST,W_GHOST,UVW_SAVE_D
            REAL(EB), DEVICE, DIMENSION(1:*) :: EXTERNAL_WALL_AREA_RATIO
            REAL(EB), DEVICE, DIMENSION(0:*) :: DX_OTHER,DY_OTHER,DZ_OTHER
            INTEGER, DEVICE, DIMENSION(0:*) :: BOUNDARY_TYPE_PREVIOUS
            INTEGER, DEVICE, DIMENSION(1:*) :: EXT_NOMD
            INTEGER, DEVICE, DIMENSION(1:*) :: WALL_IIOD,WALL_JJOD,WALL_KKOD
            INTEGER, DEVICE, DIMENSION(1:*) :: WALL_IIOD_MAX,WALL_JJOD_MAX,WALL_KKOD_MAX
            INTEGER, VALUE :: NOM_FILTER,OM_LB_I,OM_LB_J,OM_LB_K,OM_NI,OM_NJ
            INTEGER :: IW,II,JJ,KK,IOR,IIO,JJO,KKO,BT,IIO1,IIO2,JJO1,JJO2,KKO1,KKO2
            REAL(EB) :: OLD_VALUE,OTHER_VALUE,AVG_VALUE,DA_OTHER,WEIGHT,NOM_CELLS

            IW = (BLOCKIDX%X - 1) * BLOCKDIM%X + THREADIDX%X
            IF (IW < 1 .OR. IW > N_EXTERNAL_WALL_CELLS) RETURN

            BT = WALL_TOPO(IW)%BOUNDARY_TYPE
            BOUNDARY_TYPE_PREVIOUS(IW) = BT
            IF (EXT_NOMD(IW) /= NOM_FILTER) RETURN
            IF (BT /= INTERPOLATED_BOUNDARY_VALUE) RETURN

            II  = WALL_TOPO(IW)%II
            JJ  = WALL_TOPO(IW)%JJ
            KK  = WALL_TOPO(IW)%KK
            IOR = WALL_TOPO(IW)%IOR
            IIO1 = WALL_IIOD(IW)
            JJO1 = WALL_JJOD(IW)
            KKO1 = WALL_KKOD(IW)
            IIO2 = WALL_IIOD_MAX(IW)
            JJO2 = WALL_JJOD_MAX(IW)
            KKO2 = WALL_KKOD_MAX(IW)
            NOM_CELLS = REAL((IIO2-IIO1+1)*(JJO2-JJO1+1)*(KKO2-KKO1+1),EB)

            SELECT CASE(IOR)
                CASE( 1)
                    OLD_VALUE = UU(0,JJ,KK)
                    DA_OTHER = 0._EB
                    OTHER_VALUE = 0._EB
                    DO KKO=KKO1,KKO2
                        DO JJO=JJO1,JJO2
                            DO IIO=IIO1,IIO2
                                WEIGHT = DY_OTHER(JJO)*DZ_OTHER(KKO)
                                DA_OTHER = DA_OTHER + WEIGHT
                                OTHER_VALUE = OTHER_VALUE + OM_UU(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)*WEIGHT
                            ENDDO
                        ENDDO
                    ENDDO
                    IF (DA_OTHER > 0._EB) OTHER_VALUE = OTHER_VALUE/DA_OTHER
                    AVG_VALUE = 0.5_EB*(OLD_VALUE + OTHER_VALUE)
                    UVW_SAVE_D(IW) = OLD_VALUE
                    UU(0,JJ,KK) = AVG_VALUE
                    IF (EXTERNAL_WALL_AREA_RATIO(IW) > 0.9_EB) THEN
                        DO KKO=KKO1,KKO2
                            DO JJO=JJO1,JJO2
                                DO IIO=IIO1,IIO2
                                    OM_UU(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K) = &
                                        0.5_EB*(OM_UU(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)+OLD_VALUE)
                                ENDDO
                            ENDDO
                        ENDDO
                    ENDIF
                CASE(-1)
                    OLD_VALUE = UU(IBAR_LOCAL,JJ,KK)
                    DA_OTHER = 0._EB
                    OTHER_VALUE = 0._EB
                    DO KKO=KKO1,KKO2
                        DO JJO=JJO1,JJO2
                            DO IIO=IIO1,IIO2
                                WEIGHT = DY_OTHER(JJO)*DZ_OTHER(KKO)
                                DA_OTHER = DA_OTHER + WEIGHT
                                OTHER_VALUE = OTHER_VALUE + OM_UU(IIO-1-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)*WEIGHT
                            ENDDO
                        ENDDO
                    ENDDO
                    IF (DA_OTHER > 0._EB) OTHER_VALUE = OTHER_VALUE/DA_OTHER
                    AVG_VALUE = 0.5_EB*(OLD_VALUE + OTHER_VALUE)
                    UVW_SAVE_D(IW) = OLD_VALUE
                    UU(IBAR_LOCAL,JJ,KK) = AVG_VALUE
                    IF (EXTERNAL_WALL_AREA_RATIO(IW) > 0.9_EB) THEN
                        DO KKO=KKO1,KKO2
                            DO JJO=JJO1,JJO2
                                DO IIO=IIO1,IIO2
                                    OM_UU(IIO-1-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K) = &
                                        0.5_EB*(OM_UU(IIO-1-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)+OLD_VALUE)
                                ENDDO
                            ENDDO
                        ENDDO
                    ENDIF
                CASE( 2)
                    OLD_VALUE = VV(II,0,KK)
                    DA_OTHER = 0._EB
                    OTHER_VALUE = 0._EB
                    DO KKO=KKO1,KKO2
                        DO JJO=JJO1,JJO2
                            DO IIO=IIO1,IIO2
                                WEIGHT = DX_OTHER(IIO)*DZ_OTHER(KKO)
                                DA_OTHER = DA_OTHER + WEIGHT
                                OTHER_VALUE = OTHER_VALUE + OM_VV(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)*WEIGHT
                            ENDDO
                        ENDDO
                    ENDDO
                    IF (DA_OTHER > 0._EB) OTHER_VALUE = OTHER_VALUE/DA_OTHER
                    AVG_VALUE = 0.5_EB*(OLD_VALUE + OTHER_VALUE)
                    UVW_SAVE_D(IW) = OLD_VALUE
                    VV(II,0,KK) = AVG_VALUE
                    IF (EXTERNAL_WALL_AREA_RATIO(IW) > 0.9_EB) THEN
                        DO KKO=KKO1,KKO2
                            DO JJO=JJO1,JJO2
                                DO IIO=IIO1,IIO2
                                    OM_VV(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K) = &
                                        0.5_EB*(OM_VV(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)+OLD_VALUE)
                                ENDDO
                            ENDDO
                        ENDDO
                    ENDIF
                CASE(-2)
                    OLD_VALUE = VV(II,JBAR_LOCAL,KK)
                    DA_OTHER = 0._EB
                    OTHER_VALUE = 0._EB
                    DO KKO=KKO1,KKO2
                        DO JJO=JJO1,JJO2
                            DO IIO=IIO1,IIO2
                                WEIGHT = DX_OTHER(IIO)*DZ_OTHER(KKO)
                                DA_OTHER = DA_OTHER + WEIGHT
                                OTHER_VALUE = OTHER_VALUE + OM_VV(IIO-OM_LB_I,JJO-1-OM_LB_J,KKO-OM_LB_K)*WEIGHT
                            ENDDO
                        ENDDO
                    ENDDO
                    IF (DA_OTHER > 0._EB) OTHER_VALUE = OTHER_VALUE/DA_OTHER
                    AVG_VALUE = 0.5_EB*(OLD_VALUE + OTHER_VALUE)
                    UVW_SAVE_D(IW) = OLD_VALUE
                    VV(II,JBAR_LOCAL,KK) = AVG_VALUE
                    IF (EXTERNAL_WALL_AREA_RATIO(IW) > 0.9_EB) THEN
                        DO KKO=KKO1,KKO2
                            DO JJO=JJO1,JJO2
                                DO IIO=IIO1,IIO2
                                    OM_VV(IIO-OM_LB_I,JJO-1-OM_LB_J,KKO-OM_LB_K) = &
                                        0.5_EB*(OM_VV(IIO-OM_LB_I,JJO-1-OM_LB_J,KKO-OM_LB_K)+OLD_VALUE)
                                ENDDO
                            ENDDO
                        ENDDO
                    ENDIF
                CASE( 3)
                    OLD_VALUE = WW(II,JJ,0)
                    DA_OTHER = 0._EB
                    OTHER_VALUE = 0._EB
                    DO KKO=KKO1,KKO2
                        DO JJO=JJO1,JJO2
                            DO IIO=IIO1,IIO2
                                WEIGHT = DX_OTHER(IIO)*DY_OTHER(JJO)
                                DA_OTHER = DA_OTHER + WEIGHT
                                OTHER_VALUE = OTHER_VALUE + OM_WW(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)*WEIGHT
                            ENDDO
                        ENDDO
                    ENDDO
                    IF (DA_OTHER > 0._EB) OTHER_VALUE = OTHER_VALUE/DA_OTHER
                    AVG_VALUE = 0.5_EB*(OLD_VALUE + OTHER_VALUE)
                    UVW_SAVE_D(IW) = OLD_VALUE
                    WW(II,JJ,0) = AVG_VALUE
                    IF (EXTERNAL_WALL_AREA_RATIO(IW) > 0.9_EB) THEN
                        DO KKO=KKO1,KKO2
                            DO JJO=JJO1,JJO2
                                DO IIO=IIO1,IIO2
                                    OM_WW(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K) = &
                                        0.5_EB*(OM_WW(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K)+OLD_VALUE)
                                ENDDO
                            ENDDO
                        ENDDO
                    ENDIF
                CASE(-3)
                    OLD_VALUE = WW(II,JJ,KBAR_LOCAL)
                    DA_OTHER = 0._EB
                    OTHER_VALUE = 0._EB
                    DO KKO=KKO1,KKO2
                        DO JJO=JJO1,JJO2
                            DO IIO=IIO1,IIO2
                                WEIGHT = DX_OTHER(IIO)*DY_OTHER(JJO)
                                DA_OTHER = DA_OTHER + WEIGHT
                                OTHER_VALUE = OTHER_VALUE + OM_WW(IIO-OM_LB_I,JJO-OM_LB_J,KKO-1-OM_LB_K)*WEIGHT
                            ENDDO
                        ENDDO
                    ENDDO
                    IF (DA_OTHER > 0._EB) OTHER_VALUE = OTHER_VALUE/DA_OTHER
                    AVG_VALUE = 0.5_EB*(OLD_VALUE + OTHER_VALUE)
                    UVW_SAVE_D(IW) = OLD_VALUE
                    WW(II,JJ,KBAR_LOCAL) = AVG_VALUE
                    IF (EXTERNAL_WALL_AREA_RATIO(IW) > 0.9_EB) THEN
                        DO KKO=KKO1,KKO2
                            DO JJO=JJO1,JJO2
                                DO IIO=IIO1,IIO2
                                    OM_WW(IIO-OM_LB_I,JJO-OM_LB_J,KKO-1-OM_LB_K) = &
                                        0.5_EB*(OM_WW(IIO-OM_LB_I,JJO-OM_LB_J,KKO-1-OM_LB_K)+OLD_VALUE)
                                ENDDO
                            ENDDO
                        ENDDO
                    ENDIF
                CASE DEFAULT
                    RETURN
            END SELECT

            U_GHOST(IW) = 0._EB
            V_GHOST(IW) = 0._EB
            W_GHOST(IW) = 0._EB
            DO KKO=KKO1,KKO2
                DO JJO=JJO1,JJO2
                    DO IIO=IIO1,IIO2
                        U_GHOST(IW) = U_GHOST(IW) + 0.5_EB*(OM_UU(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K) + &
                                                           OM_UU(IIO-1-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K))
                        V_GHOST(IW) = V_GHOST(IW) + 0.5_EB*(OM_VV(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K) + &
                                                           OM_VV(IIO-OM_LB_I,JJO-1-OM_LB_J,KKO-OM_LB_K))
                        W_GHOST(IW) = W_GHOST(IW) + 0.5_EB*(OM_WW(IIO-OM_LB_I,JJO-OM_LB_J,KKO-OM_LB_K) + &
                                                           OM_WW(IIO-OM_LB_I,JJO-OM_LB_J,KKO-1-OM_LB_K))
                    ENDDO
                ENDDO
            ENDDO
            IF (NOM_CELLS > 0._EB) THEN
                U_GHOST(IW) = U_GHOST(IW)/NOM_CELLS
                V_GHOST(IW) = V_GHOST(IW)/NOM_CELLS
                W_GHOST(IW) = W_GHOST(IW)/NOM_CELLS
            ENDIF

        END SUBROUTINE MATCH_VELOCITY_KERNEL
END MODULE G_VELO
