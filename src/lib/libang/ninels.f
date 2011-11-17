*
*     -------------------------------------------------------------
*      N I N E L S
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
*     Vilnius, LITHUANIA                              January 1997 *
*
      SUBROUTINE NINELS(J1,J2,J3,L1,L2,L3,K1,K2,K3,I,IN,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      A=ZERO
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
      IN=1
      IF(J1*J2*J3*L1*L2*L3*K1*K2*K3.EQ.0) THEN
        CALL NINE0(J1,J2,J3,L1,L2,L3,K1,K2,K3,A)
      ELSEIF(K3.NE.2) THEN
        CALL NINE(J1,J2,J3,L1,L2,L3,K1,K2,K3,I,IN,A)
      ELSEIF(J3.EQ.L3) THEN
        IF(K1.EQ.K2) THEN
*
*    Case  J3 = L3           and   K1 = K2
*
          CALL NINE11(J1,J2,J3,L1,L2,K1,A)
        ELSEIF(K1-2.EQ.K2) THEN
*
*    Case  J3 = L3           and   K1/2 + 1 = K2/2
*
          CALL NINE12(J1,J2,J3,L1,L2,K1,A)
        ELSEIF(K1+2.EQ.K2) THEN
*
*    Case  J3 = L3           and   K1/2 = K2/2 + 1
*
          CALL NINE12(J2,J1,J3,L2,L1,K2,A)
          IF(MOD(J1+J2+J3+L1+L2+L3+K1+K2+K3,4).NE.0) A=-A 
        ELSE
          CALL NINE(J1,J2,J3,L1,L2,L3,K1,K2,K3,I,IN,A)
        ENDIF
      ELSEIF(J3-2.EQ.L3) THEN
        IF(K1.EQ.K2) THEN
*
*    Case  J3/2 + 1 = L3/2   and   K1 = K2
*
          CALL NINE12(J1,L1,K1,J2,L2,J3,A)
        ELSEIF(K1-2.EQ.K2) THEN
*
*    Case  J3/2 + 1 = L3/2   and   K1/2 + 1 = K2/2
*
          CALL NINE13(J1,J2,J3,L1,L2,K1,A)
        ELSEIF(K1+2.EQ.K2) THEN
*
*    Case  J3/2 + 1 = L3/2   and   K1/2 - 1 = K2/2
*
          CALL NINE13(J2,J1,J3,L2,L1,K2,A)
          IF(MOD(J1+J2+J3+L1+L2+L3+K1+K2+K3,4).NE.0) A=-A 
        ELSE
          CALL NINE(J1,J2,J3,L1,L2,L3,K1,K2,K3,I,IN,A)
        ENDIF
      ELSEIF(J3+2.EQ.L3) THEN
        IF(K1.EQ.K2) THEN
*
*    Case  J3/2 = L3/2 + 1   and   K1 = K2
*
          CALL NINE12(L1,J1,K1,L2,J2,L3,A)
          IF(MOD(J1+J2+J3+L1+L2+L3+K1+K2+K3,4).NE.0) A=-A 
        ELSEIF(K1-2.EQ.K2) THEN
*
*    Case  J3/2 = L3/2 + 1   and   K1/2 - 1 = K2/2
*
          CALL NINE13(L1,L2,L3,J1,J2,K1,A)
          IF(MOD(J1+J2+J3+L1+L2+L3+K1+K2+K3,4).NE.0) A=-A 
        ELSEIF(K1+2.EQ.K2) THEN
*
*    Case  J3/2 = L3/2 + 1   and   K1/2 - 1 = K2/2
*
          CALL NINE13(L2,L1,L3,J2,J1,K2,A)
        ELSE
          CALL NINE(J1,J2,J3,L1,L2,L3,K1,K2,K3,I,IN,A)
        ENDIF
      ELSE
        CALL NINE(J1,J2,J3,L1,L2,L3,K1,K2,K3,I,IN,A)
      ENDIF
      RETURN
      END
