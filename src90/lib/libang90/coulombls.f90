!     ..........................................................   :
!                                                                  :
!          Block                                                   :
!                   N O N - R E L A T I V I S T I C                :
!                                                                  :
!     For Calculation Angular Momentum Coefficients for Non-       :
!     Relativistic operators                                       :
!                                                                  :
!     Written by G. Gaigalas,                                      :
!                  Department  of  Computer Science,               :
!                  Vanderbilt University,  Nashville               :
!                                                 February  1994   :
!                  Universite Libre de Bruxelles                   :
!                                                 December  1995   :
!                  Department  of  Computer Science,               :
!                  Vanderbilt University,  Nashville               :
!                                                 September 1997   :
!                                                                  :
!     ..........................................................   :
!
!
!     -------------------------------------------------------------
!      C O U L O M B L S
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF COULOMB INTERACTIONS BETWEEN THE ELECTRONS                *
!                                                                  *
!                          k   k+1  (k) (k)                        *
!     (n l L S  n l L S ::r  / r  ( C   C )::n l L S  n l L S )    *
!       1 1 1 1  2 2 2 2   <    >             3 3 3 3  4 4 4 4     *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE COULOMBLS(L1, L2, L3, L4, KL, AA) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:48:27  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ittk_I 
      USE rme_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: L1 
      INTEGER  :: L2 
      INTEGER  :: L3 
      INTEGER  :: L4 
      INTEGER  :: KL 
      REAL(DOUBLE) , INTENT(OUT) :: AA 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      AA = ZERO 
      IF (ITTK(L1,L3,KL) == 0) RETURN  
      IF (ITTK(L2,L4,KL) == 0) RETURN  
      AA = RME(L1,L3,KL) 
      IF (DABS(AA) < EPS) RETURN  
      AA = AA*RME(L2,L4,KL) 
      RETURN  
      END SUBROUTINE COULOMBLS 
