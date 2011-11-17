*
*     ------------------------------------------------------------------
*                       A V E R D
*     ------------------------------------------------------------------
*
*     Add the deviations to the average energy for a partially filled
*       d- shell
*
      SUBROUTINE AVERD(J,K,N,A)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER F2DD2(16),F4DD2(16),F2DD3(16),F4DD3(16),
     :F2DD4(16),F4DD4(16),F2DD5(16),F4DD5(16)
      DIMENSION IAVS(5),IAVV(5)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DATA IAVS/1,0,-2,0,-2/
      DATA IAVV/1,0,63,0,63/
*             ... d2 coefficients
      DATA F2DD2/140,0,77,3*0,-13,0,-58,3*0,50,3*0/
      DATA F4DD2/140,0,-70,3*0,50,0,5,3*0,15,3*0/
*             ... d3 coefficients
      DATA F2DD3/2*0,42,-12,105,69,2*0,-93,123,0,-57,2*0,-12,0/
      DATA F4DD3/2*0,-105,30,105,-15,2*0,-30,-45,0,55,2*0,30,0/
*             ... d4 coefficients
      DATA F2DD4/210,138,21,57,-105,39,219,111,66,12,84,-24,30,
     :48,-69,-51/
      DATA F4DD4/210,-30,70,-55,-105,-45,30,-15,45,-30,0,-10,135,
     :20,15,75/
*             ... d5 coefficients
      DATA F2DD5/-175,113,-112,320,140,104,-22,86,23,-85,59,167,
     :-85,23,-58,-76/
      DATA F4DD5/-175,-55,35,-100,140,20,-85,-40,-40,125,-25,-15,
     :-50,-5,110,50/
      A=ZERO 
      IF (MOD(K,2) .NE. 0) RETURN
      IF (K .GT. 4) RETURN
      IDV=0 
      IVV=441 
      KK=K+1
      IF (N .EQ. 1 .OR. N .EQ. 9) THEN
        IF (J .NE. 13) RETURN
      ELSEIF (N .EQ. 2 .OR. N .EQ. 8) THEN
        JJ=J-24
        IF (K .EQ. 2) THEN
          IDV=F2DD2(JJ)
        ELSEIF (K .EQ. 4) THEN
          IDV=F4DD2(JJ)
        END IF
      ELSEIF (N .EQ. 3 .OR. N .EQ. 7) THEN
        JJ=J-8
        IF (K .EQ. 2) THEN
          IDV=F2DD3(JJ)
        ELSEIF (K .EQ. 4) THEN
          IDV=F4DD3(JJ)
        END IF
      ELSEIF (N .EQ. 4 .OR. N .EQ. 6) THEN
        JJ=J-24
        IF (K .EQ. 2) THEN
          IDV=F2DD4(JJ)
        ELSEIF (K .EQ. 4) THEN
          IDV=F4DD4(JJ)
        END IF
      ELSEIF (N .EQ. 5) THEN
        JJ=J-8
        IF (K .EQ. 2) THEN
          IDV=F2DD5(JJ)
        ELSEIF (K .EQ. 4) THEN
          IDV=F4DD5(JJ)
        END IF
      ELSEIF (N .EQ. 10) THEN
        IF (J .NE. 25) RETURN
      END IF
      A=DBLE(IAVS(KK)*N*(N-1))*HALF/DBLE(IAVV(KK))+DBLE(IDV)/DBLE(IVV)
      RETURN
      END
