!
!     -------------------------------------------------------------------
! *** GRAD2
!     ------------------------------------------------------------------
!
! *** THE GRAD2 FUNCTION SUBPROGRAM COMPUTES THE FOLLOWING DIRECTLY
! ***        <P(I)[R.D + F(DELTA)[ P(J)>
!
      REAL(KIND(0.0D0)) FUNCTION GRAD2 (I, J) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      use param_C
      use radial_C
      use NEL_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:51:26  11/20/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE quadr_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I 
      INTEGER , INTENT(IN) :: J 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: JJ, II, K, LI, LJ, IL, MM 
      REAL(DOUBLE), DIMENSION(NOD) :: Q 
      REAL(DOUBLE) :: A1, A2, A, FACT, G, F1, F2, G0, G1, G2, F3, &
         F4, G3, G4, U, DELTA 
!-----------------------------------------------
      JJ = I 
      II = J 
      Q(:NO) = P(:NO,JJ)*R(:NO) 
      LI = L(I) 
      LJ = L(J) 
      IL = IABS(LI - LJ) 
      IF (IL==0 .OR. IL==2) THEN 
         A1 = (LI + LJ + D2)/((LI + D1)*(LJ + D1)) 
         A2 = ((LJ + D5)*(LJ + D1) + (LJ + D1 + D5)*(LI + D1))/((LI + D1)*(LJ&
             + D1)) 
         A = A1 - A2*(LI + LJ + D3)/((LI + LJ + D4)*(LJ + D5)) 
         FACT = (LJ + D5)/(LI + LJ + D3) 
         G = R(1)**2*P(1,I)*P(1,J)*FACT*(1. + A*Z*R(1)) 
         MM = MIN0(MAX_(I) + 1,MAX_(J) + 1,ND) 
         K = 2 
         F1 = D5*(P(K+1,II)-P(K-1,II)) 
         F2 = P(K+1,II) - D2*P(K,II) + P(K-1,II) 
         G0 = Q(K)*R(K) 
         G1 = D5*(Q(K+1)*R(K+1)-Q(K-1)*R(K-1)) 
         G2 = Q(K+1)*R(K+1) - D2*Q(K)*R(K) + Q(K-1)*R(K-1) 
         G = G + D2*F1*G0 + (D2*F2*G1 + F1*G2)/D3 
         DO K = 4, MM, 2 
            F1 = D5*(P(K+1,II)-P(K-1,II)) 
            F2 = P(K+1,II) - D2*P(K,II) + P(K-1,II) 
            F3 = D5*(P(K+1,II)-P(K-2,II)) - D2*F1 
            F4 = P(K+2,II) + P(K-2,II) - D4*(P(K+1,II)+P(K-1,II)) + D6*P(K,II) 
            G0 = Q(K)*R(K) 
            G1 = D5*(Q(K+1)*R(K+1)-Q(K-1)*R(K-1)) 
            G2 = Q(K+1)*R(K+1) - D2*Q(K)*R(K) + Q(K-1)*R(K-1) 
            G3 = D5*(Q(K+2)*R(K+2)-Q(K-2)*R(K-2)) - D2*G1 
            G4 = Q(K+2)*R(K+2) + Q(K-2)*R(K-2) - D4*(Q(K+1)*R(K+1)+Q(K-1)*R(K-1&
               )) + D6*Q(K)*R(K) 
            G = G + D2*F1*G0 + (D2*F2*G1 + F1*G2)/D3 - (F1*G4 - F4*G1 + D4*(F2*&
               G3 - F3*G2))/90.E0 
         END DO 
         U = QUADR(JJ,II,0) 
         G = G - D5*U 
         DELTA = LJ - LI 
         IF (DELTA <= 0.D0) THEN 
            IF (DELTA /= 0.D0) THEN 
               GRAD2 = G - (LI - D2)*U 
               RETURN  
            ENDIF 
            GRAD2 = G + (D1 + D5)*U 
            RETURN  
         ENDIF 
         GRAD2 = G + (LI + D3)*U 
         RETURN  
      ENDIF 
      WRITE (6, 101) I, J 
  101 FORMAT(5X,'L(I)-L(J) NOT=0,2 FOR I = ',I2,' AND J = ',I2) 
      STOP  
      END FUNCTION GRAD2 
