*
*     ------------------------------------------------------------------
*                       A V E R P
*     ------------------------------------------------------------------
*
*     Add the deviations to the average energy for a partially filled
*       p- shell
*
      SUBROUTINE AVERP(J,K,N,A)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER F2PP2(3),F2PP3(3)
      DIMENSION IAVS(3),IAVV(3)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DATA IAVS/1,0,-2/
      DATA IAVV/1,0,25/
*             ... p2 coefficients
      DATA F2PP2/12,-3,3/
*             ... p3 coefficients
      DATA F2PP3/-9,6,0/
      A=ZERO 
      IF (MOD(K,2) .NE. 0) RETURN
      IF (K .GT. 2) RETURN
      IDV=0 
      IVV=25 
      KK=K+1
      IF (N .EQ. 1 .OR. N .EQ. 5) THEN
        IF (J .NE. 4) RETURN
      ELSEIF (N .EQ. 2 .OR. N .EQ. 4) THEN
        IF (K .EQ. 2) THEN
          JJ=J-5
          IDV=F2PP2(JJ)
        END IF
      ELSEIF (N .EQ. 3) THEN
        IF (K .EQ. 2) THEN
          JJ=J-2
          IDV=F2PP3(JJ)
        END IF
      ELSEIF (N .EQ. 6) THEN
        IF (J .NE. 6) RETURN
      END IF
      A=DBLE(IAVS(KK)*N*(N-1))*HALF/DBLE(IAVV(KK))+DBLE(IDV)/DBLE(IVV)
      RETURN
      END
