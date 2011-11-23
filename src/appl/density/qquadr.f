*
*-----------------------------------------------------------------------
*     Q U A D R 
*-----------------------------------------------------------------------
*
*     QUADR INTEGRATES P(I)*P(J)*R**KK BY SIMPON'S RULE
*
      DOUBLE PRECISION FUNCTION QQUADR(I,J,KK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NWD=60)
      COMMON /PARAT/D0,D1,D2,D3,D4,D5,D6,D10,H,H1,NO,ND
      COMMON/ADATA/P(NOD,NWD),R(NOD),RR(NOD),R2(NOD),
     :       ZED,AZ(NWD),L(NWD),MAX(NWD)
      
      Z=ZED
      K = KK + 2
      LI = L(I)
      LJ = L(J)
      DEN = LI + LJ + 1 + K
      ZR = Z*R(4)
      BI = (P(4,I)/(AZ(I)*R2(4)*R(4)**LI) - D1+ZR/(LI+1))/ZR**2
      BJ = (P(4,J)/(AZ(J)*R2(4)*R(4)**LJ) - D1+ZR/(LJ+1))/ZR**2
      ALPHA= (D1/(LI + 1) + D1/(LJ + 1))/(DEN + D1)
      ZR = Z*R(1)
      BETA = (DEN+D1)*ALPHA**2 - D2*(BI+BJ+D1/((LI+1)*(LJ+1)))/(DEN+D2)
      D = P(1,I)*P(1,J)*R(1)**K*(((BETA*ZR+ALPHA)*ZR+D1)/(DEN*H1)+D5)
      M = MIN0(MAX(I),MAX(J)) - 1
      DO 1 JJ = 2,M,2
      JP = JJ + 1
    1 D=D+D2*P(JJ,I)*P(JJ,J)*R(JJ)**K+P(JP,I)*P(JP,J)*R(JP)**K
      QQUADR = D*H1
      RETURN
      END
