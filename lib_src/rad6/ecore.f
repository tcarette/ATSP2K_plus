*
*    ------------------------------------------------------------------
*      E C O R E
*    ------------------------------------------------------------------
*
*
      SUBROUTINE ECORE(EL,EC,REL) 
* 
* *** COMPUTES THE ENERGY OF THE COMMON CLOSED SHELLS 
* 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER (NOD=220,NWD=30,NWD2=2*NWD)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID 
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS 
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD2),YK(NOD), 
     :   YR(NOD),X(NOD),AZ(NWD2),L(NWD2),MAX(NWD2),N(NWD2) 
      CHARACTER EL(*)*3
      LOGICAL REL
* 
      EC = D0 
      DO 10 I = 1,NCLOSD 
         SUMI = 4*L(I)+2 
         TI   = FK(I,I,0,REL) 
         DO 20 K = 2,2*L(I),2 
            TI = TI - CA(L(I),K)*FK(I,I,K,REL) 
   20    CONTINUE 
         EC = EC + SUMI*((SUMI-1)*TI - HL(EL,I,I,REL))/D2 
         DO 30 J = 1,I-1 
            SUMJ = 4*L(J)+2 
            TIJ = FK(I,J,0,REL) 
            DO 40 K=IABS(L(I)-L(J)),L(I)+L(J),2 
               TIJ = TIJ -CB(L(I),L(J),K)*GK(I,J,K,REL) 
   40       CONTINUE 
            EC = EC + SUMI*SUMJ*TIJ 
   30    CONTINUE 
   10 CONTINUE 
      END 
