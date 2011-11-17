*
*     -----------------------------------------------------------------
*           F A C T R L
*     -----------------------------------------------------------------
*
*
      SUBROUTINE FACTRL(NFACT) 
* 
*      GAM(I) = LOG( GAMMA(I-1) ), WHERE GAMMA(I) = FACTORIAL I-1 
* 
      IMPLICIT REAL *8(A-H,O-Z) 
* 
      COMMON/FACT/GAM(100) 
      DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/ 
* 
      GAMMA=ONE 
      GAM(1) = ZERO 
      DO 1 I=1,NFACT-1 
         GAMMA=I*GAMMA 
         GAM(I+1) = DLOG(GAMMA) 
    1 CONTINUE 
      DO 20 I = NFACT+1,(100) 
         X = I-1 
         GAM(I) = GAM(I-1) + DLOG(X) 
   20 CONTINUE 
      RETURN 
      END 
