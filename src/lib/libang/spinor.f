*
*     ..........................................................   :
*                                                                  : 
*          Block                                                   : 
*                 S U B M A T R I X    E L E M E N T S             :
*                                                                  : 
*     For Calculation Angular Momentum Coefficients for            :
*     Relativistic operators                                       :
*                                                                  : 
*     Written by G. Gaigalas,                                      : 
*                  Department  of  Computer Science,               : 
*                  Vanderbilt University,  Nashville               : 
*                                                   October 1996   :  
*                                                                  : 
*     ..........................................................   :
*        i)  Spin - orbit                                          : 
*     ..........................................................   :
*
*     -------------------------------------------------------------
*      S P I N O R 
*     -------------------------------------------------------------
*
*     CALCULATE ONE-ELECTRON NUCLEAR SPIN-ORBIT OPERATOR           *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Universite Libre de Bruxelles, Belgium         October 1995  *
*
      SUBROUTINE SPINOR(L1,L2,I,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      A=ZERO
      I=5
      IF(L1.NE.L2)RETURN
      IF(L1.EQ.0)RETURN
      A=DSQRT(DBLE(6*L1*(L1+1)*(2*L1+1)))
      KK=IHSH+IHSH-1
      A=-A*DSQRT(DBLE(J1QN1(KK,2)*J1QN1(KK,3)))
      RETURN
      END
