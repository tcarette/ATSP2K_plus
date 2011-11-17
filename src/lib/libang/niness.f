*
*     -------------------------------------------------------------
*      N I N E S S
*     -------------------------------------------------------------
*                                                                  *
*    THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT          *
*                                                                  *
*     |  1/2  1/2  J1  |                                           *
*     |  1/2  1/2  J2  |                                           *
*     |  J3   J4   J5  |                                           *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville               October 1996 *
*
      SUBROUTINE NINESS(J1,J2,J3,J4,J5,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      DIMENSION IS(2,2,2,2,3),IV(2,2,2,2,3)
      DATA IS/1,2*0,1,8*0,1,2*0,-1,5*0,1,-1,1,0,-1,2*1,
     :0,2*1,16*0,1/
      DATA IV/4,2*1,12,8*1,12,2*1,324,5*1,2*36,54,1,2*36,54,1,
     :2*54,16*1,81/
      K1=J1+1
      K2=J2+1
      K3=J3+1
      K4=J4+1
      K5=J5+1
      A=DBLE(IS(K1,K2,K3,K4,K5))/DSQRT(DBLE(IV(K1,K2,K3,K4,K5)))
      RETURN
      END
