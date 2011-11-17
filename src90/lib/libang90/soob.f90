!
!     -------------------------------------------------------------
!      S O O B
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
      SUBROUTINE SOOB(I1, L1, L2, KL1, KL2, AA) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:15:20  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ittk_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I1 
      INTEGER  :: L1 
      INTEGER  :: L2 
      INTEGER , INTENT(IN) :: KL1 
      INTEGER  :: KL2 
      REAL(DOUBLE) , INTENT(OUT) :: AA 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IS1, IS2, IS3 
!-----------------------------------------------
      AA = ZERO 
      IF (KL1 < 0) RETURN  
      IF (ITTK(L1,L2,KL2) == 0) RETURN  
      IF (KL1 < KL2) THEN 
         AA = ONE 
         IS1 = 2*KL1 + 1 
         IS1 = IS1*(2*KL1 + 3)*(L1 + L2 - KL1)*(KL1 + 1 - L2 + L1)*(KL1 + 1 + &
            L2 - L1)*(KL1 + L1 + L2 + 2) 
         IS2 = KL1 + 1 
      ELSE IF (KL1 == KL2) THEN 
         IF (I1 == 1) THEN 
            IF (L1 == L2) THEN 
               AA = ONE 
               IS1 = KL1*(KL1 + 1)*(2*KL1 + 1) 
               IS2 = 1 
            ELSE 
               IS1 = 2*KL1 + 1 
               IS2 = KL1*(KL1 + 1) 
               IS3 = L2*(L2 + 1) + KL1*(KL1 + 1) - L1*(L1 + 1) 
               AA = DBLE(IS3) 
            ENDIF 
         ELSE 
            AA = TWO 
            IS1 = (2*KL1 + 1)*KL1*(KL1 + 1) 
            IS2 = 1 
         ENDIF 
      ELSE 
         IF (KL1 < 1) RETURN  
         AA = ONE 
         IS1 = 2*KL1 + 1 
         IS1 = IS1*(2*KL1 - 1)*(L1 + L2 - KL1 + 1)*(KL1 - L2 + L1)*(KL1 + L2 - &
            L1)*(KL1 + L1 + L2 + 1) 
         IS2 = KL1 
      ENDIF 
      AA = AA*DSQRT(DBLE(IS1)/DBLE(IS2)) 
      RETURN  
      END SUBROUTINE SOOB 
