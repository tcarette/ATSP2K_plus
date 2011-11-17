*
*     -------------------------------------------------------------
*      N I N E 1 2
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT         *
*                                                                  *
*     |  J1/2      J2/2  J3/2 |                                    *
*     |  L1/2      L2/2  J3/2 |                                    *
*     |  K1/2 + 1  K1/2    1  |                                    *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE NINE12(J1,J2,J3,L1,L2,K1,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      K2=K1-2
      CALL SIXJ(J1,L1,K1,L2,J2,J3,0,S1)
      SKAI1=DBLE(J2-L2+K1)*DBLE(L2-J2+K1)*DBLE(J2+L2+K2+4)*
     :DBLE(J2+L2-K2)
      S1=S1*DSQRT(SKAI1)
      CALL SIXJ(J1,L1,K2,L2,J2,J3,0,S2)
      SKAI2=DBLE(J1-L1+K1)*DBLE(L1-J1+K1)*DBLE(J1+L1+K2+4)*
     :DBLE(J1+L1-K2)
      S2=S2*DSQRT(SKAI2)
      S=S1+S2
C      VARD=DSQRT(DBLE(8*(K2+3)*(K2+2)*(K2+1)*J3*(J3+2)*(J3+1)))
      VARD=DSQRT(DBLE(8*(K2+3))*DBLE(K2+2)*DBLE(K2+1)*DBLE(J3)*
     :DBLE(J3+2)*DBLE(J3+1))
      A=S/VARD
      IF(MOD(L1+J2+K2+J3,4).NE.0) A=-A 
      RETURN
      END


