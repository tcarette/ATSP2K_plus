!
!     ..........................................................   :
!                                                                  :
!          Block                                                   :
!                 S U B M A T R I X    E L E M E N T S             :
!                                                                  :
!     For Calculation Angular Momentum Coefficients for            :
!     Relativistic operators                                       :
!                                                                  :
!     Written by G. Gaigalas,                                      :
!                  Department  of  Computer Science,               :
!                  Vanderbilt University,  Nashville               :
!                                                   October 1996   :
!                                                                  :
!     ..........................................................   :
!        i)  Spin - orbit                                          :
!     ..........................................................   :
!
!     -------------------------------------------------------------
!      S P I N O R
!     -------------------------------------------------------------
!
!     CALCULATE ONE-ELECTRON NUCLEAR SPIN-ORBIT OPERATOR           *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Universite Libre de Bruxelles, Belgium         October 1995  *
!
      SUBROUTINE SPINOR(L1, L2, I, A) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:44:41  11/14/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: L1 
      INTEGER , INTENT(IN) :: L2 
      INTEGER , INTENT(OUT) :: I 
      REAL(DOUBLE) , INTENT(OUT) :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KK 
!-----------------------------------------------
      A = ZERO 
      I = 5 
      IF (L1 /= L2) RETURN  
      IF (L1 == 0) RETURN  
      A = DSQRT(DBLE(6*L1*(L1 + 1)*(2*L1 + 1))) 
      KK = IHSH + IHSH - 1 
      A = -A*DSQRT(DBLE(J1QN1(KK,2)*J1QN1(KK,3))) 
      RETURN  
      END SUBROUTINE SPINOR 
