!
!     ------------------------------------------------------------------
!         R M E
!     ------------------------------------------------------------------
!
!
      REAL(KIND(0.0D0)) FUNCTION RME (L, LP, K) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE FACT_C 
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: L 
      INTEGER , INTENT(IN) :: LP 
      INTEGER , INTENT(IN) :: K 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I2G, IG, I1, I2, I3 
      REAL(DOUBLE) :: QUSQRT 
!-----------------------------------------------
!
!
!--- EVALUATES THE REDUCED MATRIX ELEMENT (L//C(K)//LP)  -  SEE FANO
!    AND RACAH, IRREDUCIBLE TENSORIAL SETS, CHAP. 14, P. 81
!
!
      IF (MIN0(L,LP) == 0) THEN 
         RME = 1.D0 
      ELSE IF (K == 0) THEN 
         RME = 2*L + 1 
         RME = DSQRT(RME) 
      ELSE IF (K == 1) THEN 
         RME = MAX0(L,LP) 
         RME = DSQRT(RME) 
      ELSE 
         I2G = L + LP + K 
         IG = I2G/2 
         IF (I2G - 2*IG /= 0) THEN 
            RME = 0.D0 
         ELSE 
            I1 = IG - L 
            I2 = IG - LP 
            I3 = IG - K 
            QUSQRT = (2*L + 1)*(2*LP + 1) 
            RME = DSQRT(QUSQRT)*DEXP((GAM(2*I1+1)+GAM(2*I2+1)+GAM(2*I3+1)-GAM(&
               I2G+2))/2.D0+GAM(IG+1)-GAM(I1+1)-GAM(I2+1)-GAM(I3+1)) 
         ENDIF 
      ENDIF 
      RETURN  
      END FUNCTION RME 
