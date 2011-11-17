* 
*     ------------------------------------------------------------------ 
*              Z K 
*     ------------------------------------------------------------------ 
* 
*               k 
*       Stores Z (i, j; r) in the array YK. 
* 
* 
      SUBROUTINE ZK(I,J,K) 
      PARAMETER (NOD=220,NWD=30) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
      DEN = L(I) + L(J) + 3+ K 
      FACT = (D1/(L(I)+1) + D1/(L(J)+1))/(DEN + D1) 
      A = EH**K 
      A2 = A*A 
      H90 = H/90.D0 
      A3 = A2*A*H90 
      AI = H90/A 
      AN = 114.D0*A*H90 
      A34 = 34.D0*H90 
      F1 = RR(1)*P(1,I)*P(1,J) 
      F2 = RR(2)*P(2,I)*P(2,J) 
      F3 = RR(3)*P(3,I)*P(3,J) 
      F4 = RR(4)*P(4,I)*P(4,J) 
      YK(1) = F1*(D1 + Z*R(1)*FACT)/DEN 
      YK(2) = F2*(D1 + Z*R(2)*FACT)/DEN 
      YK(3) = YK(1)*A2 + H3*(F3 + D4*A*F2 + A2*F1) 
      DO 8 M = 5,NO 
      F5 = (RR(M)*P(M,I))*P(M,J) 
      YK(M-1) = YK(M-3)*A2 +    ( AN*F3 + A34*(F4+A2*F2)-F5*AI-F1*A3) 
      F1 = F2 
      F2 = F3 
      F3 = F4 
8     F4 = F5 
      YK(NO) = A*YK(NO-1) 
      IF (IABS(I-J)  +  IABS(K) .NE. 0) GO TO 2 
* 
*  *****  FOR Y0(I,I) SET THE LIMIT TO 1 AND REMOVE OSCILLATIONS 
*  *****  INTRODUCED BY THE USE OF SIMPSON'S RULE 
* 
      M1 = (NO/2)*2 - 1 
      M2 = M1 - 1 
      C1 = D1 - YK(M1) 
      C2 = D1 - YK(M2) 
      DO 3 M = 1,M1,2 
      YK(M) = YK(M) + C1 
3     YK(M+1) = YK(M+1) + C2 
      YK(NO) = D1 
      YK(NO-1) = D1 
2     RETURN 
      END 
