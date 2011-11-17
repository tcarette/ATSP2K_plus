!
!     ------------------------------------------------------------------
!              G R A D
!     ------------------------------------------------------------------
!
!  *****  THE GRAD FUNCTION SUBPROGRAM COMPUTES THE FOLLOWING DIRECTLY
!  *****         <P(J)^D + L(I)/R ^P(I)> WITH L(I) > L(J)
!
      DOUBLE PRECISION FUNCTION GRAD (I, J) 
      USE param_C
      USE radial_C
      USE nel_C
!
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!      PARAMETER(NOD=220)
      PARAMETER (IWRITE = 6) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  00:30:51  11/16/01  
!...Switches:                     
!      COMMON/PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
!     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
!      COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
!      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
!     :       (IQMAX,MAX_(1))
!      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
      DIMENSION :: F(NOD,4), G(NOD,4), GZ(NOD) 
      IF (IABS(L(I)-L(J)) .EQ. 1) THEN 
         LL = MAX0(L(I),L(J)) 
         II = I 
         JJ = J 
         IF (L(I) .LE. L(J)) THEN 
            II = J 
            JJ = I 
         ENDIF 
         A1 = (LL + D5)/(LL*(LL + 1)*(2*LL + 1)) 
         GRAD = R(1)*P(1,I)*P(1,J)*(D1 + A1*Z*R(1)) 
         DL = D5*P(1,I)*P(1,J)*R(1) 
         MM = MIN0(MAX_(I) + 1,MAX_(J) + 1,ND) 
         K = 2 
         F1 = D5*(P(K+1,II)-P(K-1,II)) 
         F2 = P(K+1,II) - D2*P(K,II) + P(K-1,II) 
         G0 = P(K,JJ)*R(K) 
         G1 = D5*(P(K+1,JJ)*R(K+1)-P(K-1,JJ)*R(K-1)) 
         G2 = P(K+1,JJ)*R(K+1) - D2*P(K,JJ)*R(K) + P(K-1,JJ)*R(K-1) 
         GRAD = GRAD + D2*F1*G0 + (D2*F2*G1 + F1*G2)/D3 
         DL = DL + D2*P(K,II)*P(K,JJ)*R(K) + P(K+1,II)*P(K+1,JJ)*R(K+1) 
         DO K = 4, MM, 2 
            F(K,1) = D5*(P(K+1,II)-P(K-1,II)) 
            F(K,2) = P(K+1,II) - D2*P(K,II) + P(K-1,II) 
            F(K,4) = P(K+2,II) + P(K-2,II) - D4*(P(K+1,II)+P(K-1,II)) + D6*P(K,&
               II) 
            GZ(K) = P(K,JJ)*R(K) 
            G(K,1) = D5*(P(K+1,JJ)*R(K+1)-P(K-1,JJ)*R(K-1)) 
            G(K,2) = P(K+1,JJ)*R(K+1) - D2*P(K,JJ)*R(K) + P(K-1,JJ)*R(K-1) 
            G(K,4) = P(K+2,JJ)*R(K+2) + P(K-2,JJ)*R(K-2) - D4*(P(K+1,JJ)*R(K+1)&
               +P(K-1,JJ)*R(K-1)) + D6*P(K,JJ)*R(K) 
         END DO 
         DO K = 4, MM, 2 
            F(K,3) = D5*(P(K+2,II)-P(K-2,II)) - D2*F(K,1) 
            G(K,3) = D5*(P(K+2,JJ)*R(K+2)-P(K-2,JJ)*R(K-2)) - D2*G(K,1) 
         END DO 
         DO K = 4, MM, 2 
            GZ(K) = D2*F(K,1)*GZ(K) + (D2*F(K,2)*G(K,1)+F(K,1)*G(K,2))/D3 - (F(&
               K,1)*G(K,4)-F(K,4)*G(K,1)+D4*(F(K,2)*G(K,3)-F(K,3)*G(K,2)))/&
               90.D0 
         END DO 
         DO K = 4, MM, 2 
            GRAD = GRAD + GZ(K) 
         END DO 
         DO K = 4, MM, 2 
            GZ(K) = D2*P(K,II)*P(K,JJ)*R(K) + P(K+1,II)*P(K+1,JJ)*R(K+1) 
         END DO 
         DO K = 4, MM, 2 
            DL = DL + GZ(K) 
         END DO 
         GRAD = GRAD + (LL + D5)*DL*H1 
         IF (II .EQ. I) GRAD = -GRAD 
         RETURN  
      ENDIF 
      WRITE (IWRITE, 101) I, J 
  101 FORMAT(5X,'L(I)-L(J) NOT =1 FOR I = ',I2,' AND J = ',I2) 
      STOP  
      END FUNCTION GRAD 
