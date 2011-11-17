* 
*     ------------------------------------------------------------------ 
*              F K 
*     ------------------------------------------------------------------ 
*                             k 
*       Returns the value of F (i,j) 
* 
* 
      DOUBLE PRECISION FUNCTION FK(I,J,K,REL) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      LOGICAL REL
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      CALL YKF(I,I,K,REL) 
      FK = QUADS(J,J,1) 
      IF (MASS .EQ. 2) FK = FK*(D1 + RMASS/D2)
      RETURN 
      END 
