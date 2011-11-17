*     ..........................................................   :
*                                                                  : 
*          Block                                                   : 
*                       R E C O U P L S                            : 
*                                                                  : 
*     For Calculation Recoupling Coeficients                       : 
*                                                                  : 
*     Written by G. Gaigalas,                                      : 
*                  Department  of  Computer Science,               : 
*                  Vanderbilt University,  Nashville               : 
*                                                  February 1994   : 
*                                                                  : 
*     ..........................................................   :
*
*     --------------------------------------------------------------
*     D L S A 1
*     --------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE DLSA1(K,JA1,KA,IRE,IAT,REC)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      REC=ZERO
      IA1=J1QN1(JA1,K)-1
      IB1=J1QN2(JA1,K)-1
      IF(JA1.GT.2) THEN
         K1=IHSH+JA1-2
         J1=J1QN1(K1,K)-1
         K1=K1+1
         IT1=J1QN1(K1,K)-1
         IT1S=J1QN2(K1,K)-1
      ELSE
         JJ=1
         IF(JA1.EQ.1)JJ=2
         J1=J1QN1(JJ,K)-1
         K1=IHSH+1
         IT1=J1QN1(K1,K)-1
         IT1S=J1QN2(K1,K)-1
      ENDIF
      IF(IRE.EQ.0) THEN
         IF(IXJTIK(KA,IB1,IA1,J1,IT1,IT1S).EQ.0)RETURN
         IAT=1
      ELSE
         CALL SIXJ(KA,IB1,IA1,J1,IT1,IT1S,0,A1)
         A1=A1*DSQRT(DBLE((IA1+1)*(IT1S+1)))
         IFAZ=J1+IT1+IB1+KA
         IF((IFAZ/4)*4.NE.IFAZ)A1=-A1
         REC=A1
         IAT=1
         IF(JA1.EQ.1) THEN
            IFAZ=IA1+IB1+2*J1-IT1-IT1S
            IF((IFAZ/4)*4.NE.IFAZ)REC=-REC
         ENDIF
      ENDIF
      RETURN
      END
