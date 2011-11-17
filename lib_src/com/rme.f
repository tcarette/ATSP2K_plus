*
*     ------------------------------------------------------------------
*         R M E
*     ------------------------------------------------------------------
*
*
      DOUBLE PRECISION FUNCTION RME(L,LP,K)
*
      IMPLICIT REAL *8(A-H,O-Z)
*
      COMMON/FACT/GAM(100)
*
*--- EVALUATES THE REDUCED MATRIX ELEMENT (L//C(K)//LP)  -  SEE FANO
*    AND RACAH, IRREDUCIBLE TENSORIAL SETS, CHAP. 14, P. 81
*
*
      IF (MIN0(L,LP) .EQ. 0) THEN
         RME = 1.D0
       ELSE IF ( K .EQ. 0) THEN
         RME = 2*L+1
         RME = DSQRT(RME)
       ELSE 
         I2G=L+LP+K
         IG=I2G/2
         IF (I2G - 2*IG .NE. 0) THEN
             RME = 0.D0
          ELSE
            I1=IG-L
            I2=IG-LP
            I3=IG-K
            QUSQRT=(2*L+1)*(2*LP+1)
            RME=DSQRT(QUSQRT)*DEXP((GAM(2*I1+1)+GAM(2*I2+1)+GAM(2*I3+1)-
     :        GAM(I2G+2))/2.D0 +GAM(IG+1)-GAM(I1+1)-GAM(I2+1)-GAM(I3+1))
         END IF
      END IF
      RETURN
      END
