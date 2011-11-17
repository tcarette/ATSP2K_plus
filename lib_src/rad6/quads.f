* 
*     ------------------------------------------------------------------ 
*              Q U A D S 
*     ------------------------------------------------------------------ 
* 
* 
*                                       kk 
*       Evaluates the integral of  (1/r)   YK(r) P (r) P (r)  with 
*                                                 i     j 
*   respect to r. 
* 
* 
      DOUBLE PRECISION FUNCTION  QUADS(I,J,KK) 
      PARAMETER (NOD=220,NWD=30,NWD2=2*NWD) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD2),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD2),L(NWD2),MAX(NWD2),N(NWD2)
      DEN = L(I) + L(J) + 3 
      K = 2 - KK 
      CD = D1 + Z*R(1)*(DEN-D1)/((DEN+D1)*((L(I)+1)*(L(J)+1))) 
      D = YK(1)*P(1,I)*P(1,J)*R(1)**K*( CD/(DEN*H1)+ D5) 
      DD = D0
      MX = MIN0(MAX(I),MAX(J)) - 1 
      DO 1 M = 2,MX,2 
      DD = DD + YK(M)*P(M,I)*P(M,J)*R(M)**K  
      D= D+  YK(M+1)*P(M+1,I)*P(M+1,J)*R(M+1)**K 
1     CONTINUE
      QUADS = H1*(D + D2*DD)
      RETURN 
      END 
