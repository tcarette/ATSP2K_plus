*    
*     Routines for MCHF_LIB_COM
*     Computer Physics Communication, Vol. 64, 399-405 (1991)
*
*     C O P Y R I G H T -- 1994
*
*-----------------------------------------------------------------------
*       A C U R A T
*-----------------------------------------------------------------------
*
*       Coefficients of Slater integrals are square roots of
*  rational numbers. To improve the accuracy, certain commonly
*  occurring coefficients are improved to machine precision.
*
	DOUBLE PRECISION FUNCTION ACURAT(C)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	INTEGER NUM
	INTEGER DEN(11)
	DATA DEN/2,3,7,9,15,35,49,175,189,315,441/
	DATA D1/1.D0/
*
	C2 = C*C
	ACURAT = C
	DO 1 I = 1,11
	   PROD = DEN(I)*C2
	   NUM = NINT(PROD)
	   EPS = ABS(NUM-PROD)/DEN(I)
	   IF (EPS .LE. 1.E-8) THEN
	      IF (EPS .NE. 0.) THEN
	         ACURAT = SQRT((NUM*D1)/DEN(I))
	  	 IF (C .LT. 0.) ACURAT = -ACURAT
		 RETURN
	      END IF
	   END IF
  1	CONTINUE
	END
