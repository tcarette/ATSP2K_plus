*
*-----------------------------------------------------------------------
*	A C U R A T
*-----------------------------------------------------------------------
*
*	Coefficients of Slater integrals are square roots of 
*   rational numbers.  To improve the accuaracy, certain commonly
*   occuring coefficients are improved to machine accuracy.
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
           EPS = DABS(NUM-PROD)/DEN(I)
           IF (EPS .LE. 1.D-8) THEN
              IF (EPS .NE. 0.) THEN
                 ACURAT = DSQRT((NUM*D1)/DEN(I))
                 IF (C .LT. 0.) ACURAT = -ACURAT
                 RETURN
              END IF
           END IF
  1     CONTINUE
        END
