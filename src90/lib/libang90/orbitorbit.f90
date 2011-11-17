!
!     -------------------------------------------------------------
!      O R B I T O R B I T
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF  ORBIT - ORBIT  INTERACTION                               *
!                                                                  *
!                              2                                   *
!     (n l L S  n l L S ::   -a   Sqrt(2k+1)  (2k-1)/(k+1)         *
!       1 1 1 1  2 2 2 2                                           *
!                                                                  *
!           (k-1) (1)(k)   (k-1) (k)(1)(k)(0)       k-1  k+2       *
!        [[[C   * L  ] * [[C   * L  ] ]   ]        r  / r          *
!              1    1        2    2                 <    >         *
!                                                                  *
!                                     ::n l L S  n l L S )         *
!                                        3 3 3 3  4 4 4 4          *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1997   *
!
      SUBROUTINE ORBITORBIT(L1, L2, L3, L4, KL, AA) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:22:07  11/16/01  
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
      INTEGER :: K 
      REAL(DOUBLE) :: S 
!-----------------------------------------------
      AA = ZERO 
      IF (ITTK(L1,L3,KL) == 0) RETURN  
      IF (ITTK(L2,L4,KL) == 0) RETURN  
      AA = RME(L1,L3,KL) 
      IF (DABS(AA) < EPS) RETURN  
      AA = AA*RME(L2,L4,KL) 
      IF (DABS(AA) < EPS) RETURN  
      K = KL - 1 
      S = DBLE(2*K + 1)*DBLE(L1 + L3 + K + 2)*DBLE(L1 + L3 - K)*DBLE(L1 - L3 + &
         K + 1)*DBLE(L3 - L1 + K + 1)*DBLE(L2 + L4 + K + 2)*DBLE(L2 + L4 - K)*&
         DBLE(L2 - L4 + K + 1)*DBLE(L4 - L2 + K + 1) 
      AA = AA*TWO*DSQRT(S)/DBLE(K*KL) 
      RETURN  
      END SUBROUTINE ORBITORBIT 
