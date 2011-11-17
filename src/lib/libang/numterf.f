*
*     -----------------------------------------------------------------
*     N U M T E R F
*     -----------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                          December 1993   *
*
      FUNCTION NUMTERF(I2N,I2S,I2L,N,I2Q)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON /MT67/ M76(238)
      NUMTERF=0
      ISL=(I2S*100)+I2L
      DO 1 IA =1,119
        IF((N/2)*2.EQ.N) THEN
          J=119+IA
        ELSE
          J=IA
        ENDIF
        LASTE=M76(J)
        LS=JTHN(LASTE,1,10000)
        IF(ISL.EQ.LS)THEN
          NR=JTHN(LASTE,4,100)
          IF(I2N.EQ.NR)GO TO 2
        ENDIF
    1 CONTINUE
      STOP
    2 NUMTERF=J
      I2Q=JTHN(LASTE,3,100)
C     WRITE(*,5) I2Q,NUMTERF
C   5 FORMAT(2X,'I2Q=',I6,'      J=',I6)
      RETURN
      END
