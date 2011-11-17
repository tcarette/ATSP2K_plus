*
*     -------------------------------------------------------------
*      N I N E 0
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT         *
*                                                                  *
*     |  J1/2  J2/2  J3/2 |                                        *
*     |  L1/2  L2/2  L3/2 |                                        *
*     |  K1/2  K2/2    0  |                                        *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                             March 1995   *
*
      SUBROUTINE NINE0(J1,J2,J3,L1,L2,L3,K1,K2,K3,AA)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      IF(J1.EQ.0) THEN 
        CALL SIXJ(L2,K2,J2,K3,L3,L1,0,A)
        B=DBLE((J2+1)*(L1+1))
        IFA=K2+J2+L3+L1
      ELSEIF(J2.EQ.0) THEN
        CALL SIXJ(L3,K3,J3,K1,L1,L2,0,A)
        B=DBLE((J3+1)*(L2+1))
        IFA=K3+J3+L1+L2
      ELSEIF(J3.EQ.0) THEN
        CALL SIXJ(L1,K1,J1,K2,L2,L3,0,A)
        B=DBLE((J1+1)*(L3+1))
        IFA=K1+J1+L2+L3
      ELSEIF(L1.EQ.0) THEN
        CALL SIXJ(K2,J2,L2,J3,K3,K1,0,A)
        B=DBLE((L2+1)*(K1+1))
        IFA=J2+L2+K3+K1
      ELSEIF(L2.EQ.0) THEN
        CALL SIXJ(K3,J3,L3,J1,K1,K2,0,A)
        B=DBLE((L3+1)*(K2+1))
        IFA=J3+L3+K1+K2
      ELSEIF(L3.EQ.0) THEN
        CALL SIXJ(K1,J1,L1,J2,K2,K3,0,A)
        B=DBLE((L1+1)*(K3+1))
        IFA=J1+L1+K2+K3
      ELSEIF(K1.EQ.0) THEN
        CALL SIXJ(J2,J3,J1,L3,L2,K2,0,A)
        B=DBLE((J1+1)*(K2+1))
        IFA=J3+J1+L2+K2
      ELSEIF(K2.EQ.0) THEN
        CALL SIXJ(J3,J1,J2,L1,L3,K3,0,A)
        B=DBLE((J2+1)*(K3+1))
        IFA=J1+J2+L3+K3
      ELSEIF(K3.EQ.0) THEN
        CALL SIXJ(J1,J2,J3,L2,L1,K1,0,A)
        B=DBLE((J3+1)*(K1+1))
        IFA=J2+J3+L1+K1
      ELSE
        A=ZERO
        B=ONE
      ENDIF
      AA=A/DSQRT(B)
      IF(MOD(IFA,4).NE.0)AA=-AA
      RETURN
      END


