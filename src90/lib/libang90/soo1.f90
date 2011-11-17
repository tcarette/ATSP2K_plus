!
!     ..........................................................   :
!        iii)  Spin - other - orbit                                :
!     ..........................................................   :
!
!     -------------------------------------------------------------
!      S O O 1
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
      SUBROUTINE SOO1(L1, KL1, KS1, KL2, KS2, AA) 
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
      INTEGER , INTENT(IN) :: KL1 
      INTEGER , INTENT(IN) :: KS1 
      INTEGER  :: KL2 
      INTEGER , INTENT(IN) :: KS2 
      REAL(DOUBLE) , INTENT(OUT) :: AA 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IS, KK 
!-----------------------------------------------
      AA = ZERO 
      IF (KL2 < 0) RETURN  
      IF (ITTK(L1,L1,KL2) == 0) RETURN  
      AA = RME(L1,L1,KL2) 
      AA = AA*AA 
!   s    space
      IF (KS1 < KS2) AA = TWO*AA 
!   l    space
      IF (KL1 < KL2) THEN 
         IS = (2*KL1 + 1)*(2*KL1 + 3)*(2*L1 - KL1)*(KL1 + 1)*(KL1 + 2*L1 + 2) 
      ELSE IF (KL1 == KL2) THEN 
         IS = KL1*(KL1 + 1)*(2*KL1 + 1) 
      ELSE 
         IS = (2*KL1 + 1)*(2*KL1 - 1)*(2*L1 - KL1 + 1)*KL1*(KL1 + 2*L1 + 1) 
      ENDIF 
      AA = AA*DSQRT(DBLE(IS))*TWO 
      KK = IHSH + IHSH - 1 
      AA = -AA*DSQRT(DBLE(J1QN1(KK,2)*J1QN1(KK,3))) 
      RETURN  
      END SUBROUTINE SOO1 
