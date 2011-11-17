* 
*     ------------------------------------------------------------------ 
*              V K
*     ------------------------------------------------------------------ 
* 
*                  k 
*       Evaluates V (i,j) as defined by Blume and Watson (1962). 
* 
      DOUBLE PRECISION FUNCTION VK(I,J,II,JJ,K) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      CALL DYK(I,II,K) 
      VK = QUADS(J,JJ,2)*FINE
      RETURN 
      END 
