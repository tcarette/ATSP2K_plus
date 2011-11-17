*
*     -------------------------------------------------------------
*      L M A T R I X 1 
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
*                                                  N'2 = N2        *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville           September 1997   * 
*                                                                  * 
*
      SUBROUTINE LMATRIX1(IA)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL RECOUPLS0
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      IF(JA.NE.JB)THEN
        DO 1 IN=1,3
          IF(.NOT.RECOUPLS0(IN,IA,IA,IA,IA,0))RETURN
    1   CONTINUE
      END IF
      CALL HIBFF(IA,IA,IA,IA,1)
      IF(IK1(1).NE.ID1(1)) RETURN
      A=-DBLE(ID1(4))/SQRT(DBLE(4*IK1(3)+2))
      B=A*HALF*DSQRT(DBLE(4*LJ(IA)+2))
      IF (DABS(B).GT.EPS) THEN
        if(ja.ne.jb) then
          write(*,'(A,I5,A,I5)')
     :            'Configurations',ja,'and',jb,'are identical'
          STOP
        ELSE
          LA=IJFUL(IA)
          CALL SAVENON(4,B,LJ(IA),0,LA,0,LA,JA,JB,0)
        END IF
      END IF
      RETURN
      END
