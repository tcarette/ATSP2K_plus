*
*     --------------------------------------------------------------
*     D L S A 6
*     --------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE DLSA6(K,K4,K3,K5,K2,K1,J12,IRE,IAT,REC)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      REC=ZERO
      IAT=1
      IF(IRE.EQ.0) THEN
        IF(IXJTIK(K4,K2,J12,K1,K5,K3).EQ.0) IAT=0
      ELSE
        CALL SIXJ(K4,K2,J12,K1,K5,K3,0,REC)
        REC=REC*DSQRT(DBLE((J12+1)*(K3+1)))
        IF(MOD(2*K5+K2+K4-J12,4).NE.0)REC=-REC
      ENDIF
      RETURN
      END
