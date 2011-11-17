*
*     -------------------------------------------------------------
*      S I X J 4
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
*                                                                  *
*     | JC/2  JE/2  JD/2 |                                         *
*     | JB/2  JF/2    4  |                                         *
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
*
      SUBROUTINE SIXJ4(JC,JE,JD,JB,JF,ITIK,SI)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      SI=ZERO
      IF(ITIK.NE.0) THEN
C
C     CHESKED TRIANGULAR CONDITIONS
C
        IF(IXJTIK(JC,JE,JD,JB,JF,8).EQ.0)RETURN
      ENDIF
      IF(IXJTIK(JC,JE,JD,JB,JF,6).EQ.0) THEN
	CALL GRACAH1(JC,JE,JF,JB,JD,8,SI)
	IF(MOD(JC+JE+JF+JB,4).NE.0) SI=-SI
      ELSE
        A=THREE
        C=DBLE(JC)*HALF
        E=DBLE(JE)*HALF
        D=DBLE(JD)*HALF
        B=DBLE(JB)*HALF
        F=DBLE(JF)*HALF
        X1=A*DSQRT((A+B+E+TWO)*(A-B+E+ONE)*(A+B-E+ONE)*(-A+B+E)*
     :  (A+C+F+TWO)*(A-C+F+ONE)*(A+C-F+ONE)*(-A+C+F))
        X2=(A+ONE)*DSQRT((A+B+E+ONE)*(A-B+E)*(A+B-E)*(-A+B+E+ONE)*
     :  (A+C+F+ONE)*(A-C+F)*(A+C-F)*(-A+C+F+ONE))
        X3=(TWO*A+ONE)*(TWO*(A*(A+ONE)*D*(D+ONE)-B*(B+ONE)*C*(C+ONE)-
     :  E*(E+ONE)*F*(F+ONE))+
     : (A*(A+ONE)-B*(B+ONE)-E*(E+ONE))*(A*(A+ONE)-C*(C+ONE)-F*(F+ONE)))
        IF(DABS(X2).LT.EPS) THEN
	  S2=ZERO
        ELSE
          CALL SIXJ2(JC,JE,JD,JB,JF,0,S2)
        ENDIF
        IF(DABS(X3).LT.EPS) THEN
	  S3=ZERO
        ELSE
          CALL SIXJ3(JC,JE,JD,JB,JF,0,S3)
        ENDIF
        SI=(X3*S3-X2*S2)/X1
      ENDIF
      RETURN
      END
