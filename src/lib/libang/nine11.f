*
*     -------------------------------------------------------------
*      N I N E 1 1
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT         *
*                                                                  *
*     |  J1/2  J2/2  J3/2 |                                        *
*     |  L1/2  L2/2  J3/2 |                                        *
*     |  K1/2  K1/2    1  |                                        *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE NINE11(J1,J2,J3,L1,L2,K1,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      A=ZERO
      ISKAI=DBLE((J1-L1)*(J1+L1+2)-(J2-L2)*(J2+L2+2))
      IF(ISKAI.EQ.0) RETURN
      CALL SIXJ(J1,L1,K1,L2,J2,J3,0,S)
      VARD=DSQRT(DBLE(4*K1)*DBLE(K1+2)*DBLE(K1+1)*DBLE(J3)*
     :DBLE(J3+2)*DBLE(J3+1))
      A=S*DBLE(ISKAI)/VARD
      IF(MOD(L1+J2+K1+J3,4).NE.0) A=-A 
      RETURN
      END


