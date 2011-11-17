!
!     ----------------------------------------------------------------
!        T R E D 2
!     ----------------------------------------------------------------
!
!
      SUBROUTINE TRED2(NM, N, A, D, E, Z) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NM 
      INTEGER , INTENT(IN) :: N 
      REAL(DOUBLE) , INTENT(IN) :: A(NM,N) 
      REAL(DOUBLE) , INTENT(INOUT) :: D(N) 
      REAL(DOUBLE) , INTENT(INOUT) :: E(N) 
      REAL(DOUBLE) , INTENT(INOUT) :: Z(NM,N) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, II, L, K, JP1 
      REAL(DOUBLE) :: H, SCALE, F, G, HH 
!-----------------------------------------------
      DO I = 1, N 
         Z(I,:I) = A(I,:I) 
      END DO 
      IF (N /= 1) THEN 
         DO II = 2, N 
            I = N + 2 - II 
            L = I - 1 
            H = 0.0 
            SCALE = 0.0 
            IF (L < 2) GO TO 130 
            DO K = 1, L 
               SCALE = SCALE + ABS(Z(I,K)) 
            END DO 
            IF (SCALE /= 0.0) GO TO 140 
  130       CONTINUE 
            E(I) = Z(I,L) 
            GO TO 290 
  140       CONTINUE 
            Z(I,:L) = Z(I,:L)/SCALE 
            H = SUM(Z(I,:L)*Z(I,:L)) 
            F = Z(I,L) 
            G = -SIGN(SQRT(H),F) 
            E(I) = SCALE*G 
            H = H - F*G 
            Z(I,L) = F - G 
            F = 0.0 
            DO J = 1, L 
               Z(J,I) = Z(I,J)/(SCALE*H) 
               G = 0.0 
               G = SUM(Z(J,:J)*Z(I,:J)) 
               JP1 = J + 1 
               IF (L >= JP1) THEN 
                  G = G + SUM(Z(JP1:L,J)*Z(I,JP1:L)) 
               ENDIF 
               E(J) = G/H 
               F = F + E(J)*Z(I,J) 
            END DO 
            HH = F/(H + H) 
            DO J = 1, L 
               F = Z(I,J) 
               G = E(J) - HH*F 
               E(J) = G 
               Z(J,:J) = Z(J,:J) - F*E(:J) - G*Z(I,:J) 
            END DO 
            Z(I,:L) = SCALE*Z(I,:L) 
  290       CONTINUE 
            D(I) = H 
         END DO 
      ENDIF 
      D(1) = 0.0 
      E(1) = 0.0 
      DO I = 1, N 
         L = I - 1 
         IF (D(I) /= 0.0) THEN 
            DO J = 1, L 
               G = SUM(Z(I,:L)*Z(:L,J)) 
               Z(:L,J) = Z(:L,J) - G*Z(:L,I) 
            END DO 
         ENDIF 
         D(I) = Z(I,I) 
         Z(I,I) = 1.0 
         IF (L < 1) CYCLE  
         Z(I,:L) = 0.0 
         Z(:L,I) = 0.0 
      END DO 
      RETURN  
      END SUBROUTINE TRED2 
