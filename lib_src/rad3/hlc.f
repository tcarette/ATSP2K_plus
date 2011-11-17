*
*     -----------------------------------------------------------------
*       H L C
*     -----------------------------------------------------------------
*
*
      DOUBLE PRECISION FUNCTION HLC(EL,I,J,REL) 
* 
* *** COMPUTES HL(I,J) MODIFIED BY THE INTERACTIONS WITH THE CLOSED 
*     SHELL 
* 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER(NOD=220,NWD=30)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID 
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS 
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD), 
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD) 
* 
      CHARACTER EL(*)*3
      LOGICAL REL
      HLC = HL(EL,I,J,REL) 
      DO 10 IP = 1,NCLOSD 
         SUMIP = 4*L(IP)+2 
         T = RK(I,IP,J,IP,0,REL) 
         DO 20 K = IABS(L(I)-L(IP)),L(I)+L(IP),2 
            T = T - CB(L(I),L(IP),K)*RK(I,IP,IP,J,K,REL) 
   20    CONTINUE 
         HLC = HLC - D2*SUMIP*T 
   10 CONTINUE 
      END 
