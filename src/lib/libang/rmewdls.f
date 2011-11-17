*
*     -----------------------------------------------------------------
*      R M E W D L S
*     -----------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                              June 1995   *
*
      SUBROUTINE RMEWDLS(J1,J2,K1,K2,K3,COEF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/RIBOLS/IMPTLS(40),IMGTLS(40),IMPNLS(40),IMGNLS(40)
      DIMENSION IPR(16)
      DATA IPR/0,1,3,6,10,15,21,28,36,45,55,66,78,91,105,120/
      COEF=0.00
      JI1=IMPTLS(J1)
      JI2=IMPTLS(J2)
      IF(JI1.NE.JI2) RETURN
      CALL RUMT(J1,2,LQ,LS,LL)
      CALL RUMT(J2,2,LQS,LSS,LLS)
      IF(J1.GT.J2) THEN
        JI1=J2
        JI2=J1
        IFAZ=LQ+LS+LL-LQS-LSS-LLS
      ELSE
        JI1=J1
        JI2=J2
        IFAZ=4
      ENDIF
      IF(J1.GT.24) THEN
        JI1=JI1-24
        JI2=JI2-24
        L=0
      ELSE
        JI1=JI1-8
        JI2=JI2-8
        L=1
      ENDIF
      J=IPR(JI2)+JI1
      IF(K1.EQ.1.AND.K2.EQ.2.AND.K3.EQ.0) THEN
        CALL RMEWD1LS(J,L,COEF)
      ELSEIF(K1.EQ.0.AND.K2.EQ.2.AND.K3.EQ.1) THEN
        L=IABS(L-1)
        CALL RMEWD1LS(J,L,COEF)
        IF(L.EQ.0) THEN
          IFAZ2=LQ-LQS
        ELSE
          IFAZ2=LS-LSS
        ENDIF
          IF(MOD(IFAZ2,4).NE.0) COEF=-COEF
      ELSEIF(K1.EQ.1.AND.K2.EQ.1.AND.K3.EQ.1) THEN
        IF(L.EQ.0) THEN
          CALL RMEWD2LS(J,0,COEF)
        ELSE
          CALL RMEWD2LS(J,0,COEF)
          IF(MOD(LQ-LQS,4).NE.0) COEF=-COEF
        ENDIF
      ELSEIF(K1.EQ.0.AND.K2.EQ.3.AND.K3.EQ.0) THEN
        CALL RMEWD2LS(J,1,COEF)
      ELSEIF(K1.EQ.1.AND.K2.EQ.3.AND.K3.EQ.1) THEN
        IF(L.EQ.0) THEN
          CALL RMEWD2LS(J,2,COEF)
        ELSE
          CALL RMEWD2LS(J,2,COEF)
          IF(MOD(LQ-LQS,4).NE.0) COEF=-COEF
        ENDIF
      ELSEIF(K1.EQ.1.AND.K2.EQ.4.AND.K3.EQ.0) THEN
        CALL RMEWD3LS(J,L,COEF)
      ELSEIF(K1.EQ.0.AND.K2.EQ.4.AND.K3.EQ.1) THEN
        L=IABS(L-1)
        CALL RMEWD3LS(J,L,COEF)
        IF(L.EQ.0) THEN
          IFAZ2=LQ-LQS
        ELSE
          IFAZ2=LS-LSS
        ENDIF
          IF(MOD(IFAZ2,4).NE.0) COEF=-COEF
      ENDIF
      IF(MOD(IFAZ,4).NE.0) COEF=-COEF
      RETURN
      END
