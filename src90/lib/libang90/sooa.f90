!
!     -------------------------------------------------------------
!      S O O A
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF   SPIN OTHER ORBIT  INTERACTIONS BETWEEN THE ELECTRONS    *
!                                                                  *
!     (n l L S  n l L S ::                 ::n l L S  n l L S )    *
!       1 1 1 1  2 2 2 2                      3 3 3 3  4 4 4 4     *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
!
      SUBROUTINE SOOA(L1, L2, L3, L4, KL, AA) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:15:20  11/16/01  
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
      INTEGER :: KK 
!-----------------------------------------------
      AA = ZERO 
      IF (KL < 0) RETURN  
      IF (ITTK(L1,L3,KL) == 0) RETURN  
      IF (ITTK(L2,L4,KL) == 0) RETURN  
      AA = RME(L1,L3,KL) 
      AA = AA*RME(L2,L4,KL) 
      KK = IHSH + IHSH - 1 
      AA = -AA*TWO*DSQRT(DBLE(J1QN1(KK,2)*J1QN1(KK,3))) 
      RETURN  
      END SUBROUTINE SOOA 
