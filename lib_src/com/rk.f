* 
*     ------------------------------------------------------------------ 
*              R K 
*     ------------------------------------------------------------------ 
* 
*                   k 
*       Evaluates  R (i, j; ii, jj) 
* 
* 
      DOUBLE PRECISION FUNCTION RK(I,J,II,JJ,K,REL) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      LOGICAL REL
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      CALL YKF(I,II,K,REL) 
      RK = QUADS(J,JJ,1) 
      IF (MASS .GT. 0) THEN
 	 IF (MASS .EQ. 1) THEN
            IF (K .EQ. 1) RK = RK - RMASS*GRAD(I,II)*GRAD(J,JJ)
	 ELSE
	    RK = RK*(D1 + RMASS/D2)
	    IF (K .EQ. 1) RK = RK + Z*RMASS/D2*(
     :        QUADR(I,II,1)*QUADR(J,JJ,-2)+QUADR(I,II,-2)*QUADR(J,JJ,1))
	 END IF
      END IF
      RETURN 
      END 
