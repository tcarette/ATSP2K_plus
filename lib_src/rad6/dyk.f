*     Routines for MCHF_LIB_RAD6
*     (same as RAD3 except for NOD parameter)
*     Computer Physics Communication, Vol. 64, 399-405 (1991)
*
*     C O P Y R I G H T -- 1994
*
*     ------------------------------------------------------------------ 
*              D Y K 
*     ------------------------------------------------------------------ 
* 
*       Stores in YK the values of the integral of 
*              k 
*       P (s/r) (dP /ds - P /s) integrated over the interval (0,r) 
*        i         j       j 
* 
*   which enter into the spin-orbit calculation. 
* 
* 
      SUBROUTINE DYK(I,J,K) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER(NOD=220,NWD=30,NWD2=2*NWD) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID 
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS 
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD2),YK(NOD), 
     :   YR(NOD),X(NOD),AZ(NWD2),L(NWD2),MAX(NWD2),N(NWD2) 
      DEN = L(I)+ L(J)+ 2 + K 
      FACT = L(J) 
      DO 1 JJ =1,2 
1     YK(JJ) = FACT*P(JJ,I)*P(JJ,J)*R(JJ)/DEN 
      A = EH**K 
      AA = A*A 
      A = D4*A 
      F1 = FACT*P(1,I)*P(1,J)*R(1) 
      F2 = FACT*P(2,I)*P(2,J)*R(2) 
      DO 8 M =3,ND 
      F3 = (-P(M+2,J) + D8*(P(M+1,J)-P(M-1,J)) + P(M-2,J))/(D6*H) 
      F3 = D5*P(M,I)*(F3 - P(M,J))*R(M) 
      YK(M) = YK(M-2)*AA + H3*(F3+ A*F2 + AA*F1) 
      F1 = F2 
8     F2 = F3 
      A = A*(EH)**3 
      AA = A*A/D16 
      C = 2*K+3 
      HH = C*H3 
      YK(NO)= YK(ND) 
      F1 = YK(NO) 
      F2 = F1 
      DO 9 MM = 3,NO 
      M = NO -MM+1 
      F3 =YK(M) 
      YK(M) = YK(M+2)*AA + HH*(F3 +A*F2 + AA*F1) 
      F1 = F2 
9     F2 = F3 
      RETURN 
      END 
