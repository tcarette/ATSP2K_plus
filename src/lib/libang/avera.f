*
*     ------------------------------------------------------------------
*                       A V E R A
*     ------------------------------------------------------------------
*
*     Add the deviations to the average energy for a partially filled
*       s, p, d, f- shells
*
      SUBROUTINE AVERA(L,J,K,N,A)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      IF (L .EQ. 3) THEN
        CALL AVERF(J,K,N,A)
      ELSEIF (L .EQ. 0) THEN
        A=ZERO
        IF (K .NE. 0) RETURN
        IF (N .EQ. 1) THEN
          IF (J .NE. 1) RETURN
        ELSEIF (N .EQ. 2) THEN
          IF (J .NE. 2) RETURN
        END IF
        A=DBLE(N*(N-1))*HALF
      ELSEIF (L .EQ. 1) THEN
        CALL AVERP(J,K,N,A)
      ELSEIF (L .EQ. 2) THEN
        CALL AVERD(J,K,N,A)
      END IF
      RETURN
      END
