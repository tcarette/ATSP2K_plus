*
*     -------------------------------------------------------------
*      N I N E
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF 9j COEFFICIENT         *
*                                                                  *
*     |  J1/2  J2/2  J3/2 |                                        *
*     |  L1/2  L2/2  L3/2 |                                        *
*     |  K1/2  K2/2  K3/2 |                                        *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                             March 1995   *
*
      SUBROUTINE NINE(J1,J2,J3,L1,L2,L3,K1,K2,K3,I,IN,AA)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      IF(I.EQ.1) THEN
        IN=0
        IF(ITTK(J1,J2,J3).EQ.0)RETURN
        IF(ITTK(L1,L2,L3).EQ.0)RETURN
        IF(ITTK(K1,K2,K3).EQ.0)RETURN
        IF(ITTK(J1,L1,K1).EQ.0)RETURN
        IF(ITTK(J2,L2,K2).EQ.0)RETURN
        IF(ITTK(J3,L3,K3).EQ.0)RETURN
        IN=1
        RETURN
      ENDIF
      IF(J1*J2*J3*L1*L2*L3*K1*K2*K3.EQ.0) THEN
        IN=1
        CALL NINE0(J1,J2,J3,L1,L2,L3,K1,K2,K3,AA)
      ELSE
        N1=IABS(J1-K3)
        N2=IABS(L3-J2)
        N3=IABS(L1-K2)
        N4=IABS(J2-L3)
        N5=IABS(K2-L1)
        N6=IABS(J1-K3)
        MAX=MAX0(N1,N2,N3,N4,N5,N6)
        N1=J1+K3
        N2=L3+J2
        N3=J2+L3
        N4=K2+L1
        N5=J1+K3
        N6=L1+K2
        MIN=MIN0(N1,N2,N3,N4,N5,N6)
        IN=1
        AA=ZERO
        DO 1 IX=MAX,MIN,2
          CALL SIXJ(J1,J2,J3,L3,K3,IX,0,S1)
          CALL SIXJ(L1,L2,L3,J2,IX,K2,0,S2)
          CALL SIXJ(K1,K2,K3,IX,J1,L1,0,S3)
          X=S1*S2*S3*DBLE(IX+1)
          IF(MOD(IX,2).NE.0) X=-X 
          AA=X+AA
    1   CONTINUE
      ENDIF
      RETURN
      END
