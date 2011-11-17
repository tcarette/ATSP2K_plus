* 
*     ------------------------------------------------------------------ 
*              Y K F 
*     ------------------------------------------------------------------ 
* 
*               k 
*       Stores Y (i, j; r) in the array YK 
* 
* 
      SUBROUTINE YKF(I,J,K,REL) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,NWD=30) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
      LOGICAL REL
      CALL ZK(I,J,K) 
      A =    EH**(K+1) 
      C = 2*K+1 
      A2 = A*A 
      H90 = C*H3/D30 
      A3 = A2*A*H90 
      AI = H90/A 
      AN = 114.D0*A*H90 
      A34 = 34.D0*H90 
      F1 = YK(NO)*EH**K 
      F2 = YK(NO) 
      F3 = YK(NO-1) 
      F4 = YK(ND) 
      DO 9 MM = 2,ND 
      M = NO -MM 
      F5 = YK(M-1) 
      YK(M) = YK(M+2)*A2 +     ( AN*F3 + A34*(F4+A2*F2)-F5*AI-F1*A3) 
      F1 = F2 
      F2 = F3 
      F3 = F4 
9     F4 = F5 
      YK(1) = YK(3)*A2+C*H3*(F4 + D4*A*F3 + A2*F2) 
      IF (.NOT.REL) RETURN
      MM = MAX0( MAX(I), MAX(J) ) 
      C = C*FINE 
      DO 10 M = 1,MM 
      YK(M) = YK(M) + C*P(M,I)*P(M,J) 
10    CONTINUE 
      RETURN 
      END 
