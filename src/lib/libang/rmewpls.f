*
*     -----------------------------------------------------------------
*      R M E W P L S
*     -----------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                          December 1993   *
*
      SUBROUTINE RMEWPLS(J1,J2,K1,K2,K3,COEF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/RIBOLS/IMPTLS(40),IMGTLS(40),IMPNLS(40),IMGNLS(40)
      DIMENSION IP120N(3,3),IP120L(3,3),IP111N(3,3),IP111L(3,3),
     *IP021N(3,3),IP021L(3,3)
      DATA IP120N/4*0,60,90,0,90,0/
      DATA IP120L/2*0,40,0,-90,0,-40,0,70/
      DATA IP111N/0,-72,0,72,108,-90,0,-90,0/
      DATA IP111L/0,72,0,-72,108,-90,0,-90,0/
      DATA IP021N/2*0,-40,0,-90,0,40,0,70/
      DATA IP021L/4*0,60,90,0,90,0/
      COEF=0.00
      JI1=IMPTLS(J1)
      JI2=IMPTLS(J2)
      IF(JI1.NE.JI2) RETURN
      IF(K1.EQ.1.AND.K2.EQ.2.AND.K3.EQ.0) THEN
        IF(J1.LT.6.AND.J2.LT.6) THEN
          I1=J1-2
          I2=J2-2
          IF(IP120N(I1,I2).GE.0) THEN
            COEF=DSQRT(DBLE(IP120N(I1,I2)))
          ELSE
            COEF=-DSQRT(-DBLE(IP120N(I1,I2)))
          ENDIF
        ELSEIF(J1.GT.5.AND.J2.GT.5) THEN
          I1=J1-5
          I2=J2-5
          IF(IP120L(I1,I2).GE.0) THEN
            COEF=DSQRT(DBLE(IP120L(I1,I2)))
          ELSE
            COEF=-DSQRT(-DBLE(IP120L(I1,I2)))
          ENDIF
        ENDIF
      ELSEIF(K1.EQ.1.AND.K2.EQ.1.AND.K3.EQ.1) THEN
        IF(J1.LT.6.AND.J2.LT.6) THEN
          I1=J1-2
          I2=J2-2
          IF(IP111N(I1,I2).GE.0) THEN
            COEF=DSQRT(DBLE(IP111N(I1,I2)))
          ELSE
            COEF=-DSQRT(-DBLE(IP111N(I1,I2)))
          ENDIF
        ELSEIF(J1.GT.5.AND.J2.GT.5) THEN
          I1=J1-5
          I2=J2-5
          IF(IP111L(I1,I2).GE.0) THEN
            COEF=DSQRT(DBLE(IP111L(I1,I2)))
          ELSE
            COEF=-DSQRT(-DBLE(IP111L(I1,I2)))
          ENDIF
        ENDIF
      ELSEIF(K1.EQ.0.AND.K2.EQ.2.AND.K3.EQ.1) THEN
        IF(J1.LT.6.AND.J2.LT.6) THEN
          I1=J1-2
          I2=J2-2
          IF(IP021N(I1,I2).GE.0) THEN
            COEF=DSQRT(DBLE(IP021N(I1,I2)))
          ELSE
            COEF=-DSQRT(-DBLE(IP021N(I1,I2)))
          ENDIF
        ELSEIF(J1.GT.5.AND.J2.GT.5) THEN
          I1=J1-5
          I2=J2-5
          IF(IP021L(I1,I2).GE.0) THEN
            COEF=DSQRT(DBLE(IP021L(I1,I2)))
          ELSE
            COEF=-DSQRT(-DBLE(IP021L(I1,I2)))
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
