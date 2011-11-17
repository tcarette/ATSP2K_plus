!
!     ----------------------------------------------------------------
!        T Q L 2
!     ----------------------------------------------------------------
!
!
      SUBROUTINE TQL2(NM, N, D, E, Z, IERR) 
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
      INTEGER , INTENT(OUT) :: IERR 
      REAL(DOUBLE) , INTENT(INOUT) :: D(N) 
      REAL(DOUBLE) , INTENT(INOUT) :: E(N) 
      REAL(DOUBLE) , INTENT(INOUT) :: Z(NM,N) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, L, J, M, MML, II, K 
      REAL(DOUBLE) :: MACHEP, F, B, H, P, R, C, S, G 
!-----------------------------------------------
      MACHEP = 1.387878878078144568D-17 
      IERR = 0 
      IF (N /= 1) THEN 
         E(:N-1) = E(2:N) 
         F = 0.0 
         B = 0.0 
         E(N) = 0.0 
         DO L = 1, N 
            J = 0 
            H = MACHEP*(ABS(D(L))+ABS(E(L))) 
            B = DMAX1(H,B) 
            DO M = L, N 
               IF (ABS(E(M)) > B) CYCLE  
               EXIT  
            END DO 
            IF (M /= L) THEN 
  130          CONTINUE 
               IF (J == 30) GO TO 1000 
               J = J + 1 
               P = (D(L+1)-D(L))/(2.0*E(L)) 
               R = SQRT(P*P + 1.0) 
               H = D(L) - E(L)/(P + SIGN(R,P)) 
               D(L:N) = D(L:N) - H 
               F = F + H 
               P = D(M) 
               C = 1.0 
               S = 0.0 
               MML = M - L 
               DO II = 1, MML 
                  I = M - II 
                  G = C*E(I) 
                  H = C*P 
                  IF (ABS(P) >= ABS(E(I))) THEN 
                     C = E(I)/P 
                     R = SQRT(C*C + 1.0) 
                     E(I+1) = S*P*R 
                     S = C/R 
                     C = 1.0/R 
                  ELSE 
                     C = P/E(I) 
                     R = SQRT(C*C + 1.0) 
                     E(I+1) = S*E(I)*R 
                     S = 1.0/R 
                     C = C*S 
                  ENDIF 
                  P = C*D(I) - S*G 
                  D(I+1) = H + S*(C*G + S*D(I)) 
                  DO K = 1, N 
                     H = Z(K,I+1) 
                     Z(K,I+1) = S*Z(K,I) + C*H 
                     Z(K,I) = C*Z(K,I) - S*H 
                  END DO 
               END DO 
               E(L) = S*P 
               D(L) = C*P 
               IF (ABS(E(L)) > B) GO TO 130 
            ENDIF 
            D(L) = D(L) + F 
         END DO 
         GO TO 1001 
 1000    CONTINUE 
         IERR = L 
      ENDIF 
 1001 CONTINUE 
      RETURN  
      END SUBROUTINE TQL2 
