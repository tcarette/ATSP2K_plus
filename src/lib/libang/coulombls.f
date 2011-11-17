*     ..........................................................   :
*                                                                  : 
*          Block                                                   : 
*                   N O N - R E L A T I V I S T I C                :
*                                                                  : 
*     For Calculation Angular Momentum Coefficients for Non-       :
*     Relativistic operators                                       :
*                                                                  : 
*     Written by G. Gaigalas,                                      : 
*                  Department  of  Computer Science,               : 
*                  Vanderbilt University,  Nashville               : 
*                                                 February  1994   : 
*                  Universite Libre de Bruxelles                   : 
*                                                 December  1995   : 
*                  Department  of  Computer Science,               : 
*                  Vanderbilt University,  Nashville               : 
*                                                 September 1997   : 
*                                                                  : 
*     ..........................................................   :
*
*
*     -------------------------------------------------------------
*      C O U L O M B L S
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF COULOMB INTERACTIONS BETWEEN THE ELECTRONS                *
*                                                                  *
*                          k   k+1  (k) (k)                        *
*     (n l L S  n l L S ::r  / r  ( C   C )::n l L S  n l L S )    *
*       1 1 1 1  2 2 2 2   <    >             3 3 3 3  4 4 4 4     *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      SUBROUTINE COULOMBLS(L1,L2,L3,L4,KL,AA)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      AA=ZERO
      IF(ITTK(L1,L3,KL).EQ.0)RETURN
      IF(ITTK(L2,L4,KL).EQ.0)RETURN
      AA=RME (L1,L3,KL)
      IF(DABS(AA).LT.EPS)RETURN
      AA=AA*RME (L2,L4,KL)
      RETURN
      END
