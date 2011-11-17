      subroutine brci_zeeman(jj,n_cfg,n_eigvec,weight,g_J,
     :                g_LSJdom,cfile);
      IMPLICIT NONE
      INTEGER I, J, N_CFG, N_JVALUE, N_EIGVEC
      DOUBLE PRECISION :: G_LSJ(n_cfg,n_eigvec),L(n_cfg),S(n_cfg)
      DOUBLE PRECISION :: JQ,g_J(n_eigvec);
      DOUBLE PRECISION weight(n_cfg*n_eigvec);
      DOUBLE PRECISION :: g_LSJdom(n_eigvec);
      integer dom_comp(n_eigvec);
      integer cfile,jj,ii,mx(1);
      integer IA

      rewind (cfile);
      JQ = DBLE(jj)/2.D0;
      CALL READLS(cfile,N_CFG,L,S);
      CALL GLSJ(N_CFG,N_EIGVEC,L,S,JQ,G_LSJ);
      CALL GJ(N_CFG,N_EIGVEC,JQ,WEIGHT,G_LSJ,g_j);
      IF (JQ < 0.25D0) THEN
         g_LSJdom(1:n_eigvec) = 0.d0
      ELSE
        do ii = 1,n_eigvec;
          mx = MaxLoc(abs(weight(((ii-1)*n_cfg+1):(ii*n_cfg))));
          dom_comp(ii) = mx(1)
          g_LSJdom(ii) = G_LSJ(dom_comp(ii),ii);
        end do
      END IF

      contains;

      SUBROUTINE READLS(IU,N_CFG,L,S)
      IMPLICIT NONE
      INTEGER :: IU, N_CFG, L_STRING, I
      DOUBLE PRECISION :: MAP(68:83)
      DOUBLE PRECISION :: L(n_cfg), S(n_cfg)
      CHARACTER*132 LINE
      !
      ! Define a map from the ASCII code of the spectroscopic notation 
      ! to the L quantum number
      !
      MAP(83) = 0.D0          ! S -> L = 0
      MAP(80) = 1.D0          ! P -> L = 1
      MAP(68) = 2.D0          ! D -> L = 2
      MAP(70) = 3.D0          ! F -> L = 3
      MAP(71) = 4.D0          ! G -> L = 4
      MAP(72) = 5.D0          ! H -> L = 5
      MAP(73) = 6.D0          ! I -> L = 6
      MAP(75) = 7.D0          ! K -> L = 7
      MAP(76) = 8.D0          ! L -> L = 8
      MAP(77) = 9.D0          ! M -> L = 9
      MAP(78) = 10.D0         ! N -> L = 10
      
      READ(IU,'(A)')
      READ(IU,'(A)')
      DO I = 1,N_CFG
         READ(IU,'(/A)') LINE
         L_STRING = LEN_TRIM(LINE)
!cff_02 IF (IACHAR(LINE(L_STRING-1:L_STRING-1))>57) L_STRING = L_STRING - 1   
         IA = IACHAR(LINE(L_STRING-1:L_STRING-1))
         IF (IA>57) L_STRING = L_STRING - 1
      ! To deal with configurations with only one group of equiv elec.
         S(I) = DBLE(IACHAR(LINE(L_STRING-1:L_STRING-1))-49)/2.D0
         L(I) = MAP(IACHAR(LINE(L_STRING:L_STRING)))
      END DO
      END SUBROUTINE READLS
      
      SUBROUTINE GLSJ(N_CFG,N_EIGVEC,L,S,JQ,G_LSJ)
      IMPLICIT NONE
      INTEGER :: N_CFG, N_EIGVEC, I, K
      DOUBLE PRECISION :: G_S = 2.002319304386D0
      DOUBLE PRECISION :: G_LSJ(N_CFG,N_EIGVEC)
      DOUBLE PRECISION :: L(N_CFG), S(N_CFG)
      DOUBLE PRECISION :: JQ
      DO I = 1,N_EIGVEC
         IF (JQ>0.25D0) THEN       
      ! Compute the g_j factor for all eigenvektors with J > 0
            DO K = 1,N_CFG
               G_LSJ(K,I) = 1.D0 + (G_S-1.D0)*(JQ*(JQ+1.D0)+
     :         S(K)*(S(K)+1.D0)-L(K)*(L(K)+1.D0))/(2.D0*JQ*
     :         (JQ+1.D0));
            END DO
         END IF
      END DO
      END SUBROUTINE GLSJ
      
      SUBROUTINE GJ(N_CFG,N_EIGVEC,JQ,WEIGHT,G_LSJ,G_J)
      IMPLICIT NONE
      INTEGER :: N_CFG, N_EIGVEC, I, K
      DOUBLE PRECISION :: G_J(N_EIGVEC)
      DOUBLE PRECISION ::  G_LSJ(N_CFG,N_EIGVEC)
      DOUBLE PRECISION :: JQ
      DOUBLE PRECISION :: weight(N_CFG*N_EIGVEC);
 
      DO I = 1,N_EIGVEC
         G_J(i) = 0.D0 
         IF (JQ>0.25D0) THEN
            DO K = 1,N_CFG
               G_J(i) = G_J(i) + WEIGHT(k+(i-1)*n_CFG)*
     :                     WEIGHT(k+(i-1)*n_CFG)*G_LSJ(K,I)
            END DO
         END IF

      END DO
      END SUBROUTINE GJ
      
      end subroutine brci_zeeman 
      
            
