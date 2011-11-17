!
!     -------------------------------------------------------------
!      R L S P 0
!     -------------------------------------------------------------
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE RLSP0(K, JA1, JA2, KA, IAT) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE MEDEFN_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:24:23  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ittk_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: K 
      INTEGER , INTENT(IN) :: JA1 
      INTEGER , INTENT(IN) :: JA2 
      INTEGER  :: KA 
      INTEGER , INTENT(OUT) :: IAT 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KK, LK, LD, I, J, JJ, ISKR, JI 
!-----------------------------------------------
      IF (K /= 1) THEN 
         IAT = 0 
         KK = IHSH + IHSH - 1 
         LK = J1QN1(KK,K) - 1 
         LD = J1QN2(KK,K) - 1 
         IF (ITTK(LK,LD,KA) == 0) RETURN  
      ENDIF 
      IAT = 1 
      IF (IHSH == 1) RETURN  
      DO I = 1, IHSH 
         IF (JA1 == I) CYCLE  
         IF (JA2 == I) CYCLE  
         IF (J1QN1(I,K) == J1QN2(I,K)) CYCLE  
         IAT = 0 
      END DO 
      IF (IAT == 0) RETURN  
      IF (K == 1) RETURN  
      IF (IHSH <= 2) RETURN  
      IF (JA1 <= 2) RETURN  
      DO J = 3, JA1 
         JJ = IHSH - 2 + J 
         IF (J1QN1(JJ,K) /= J1QN2(JJ,K)) IAT = 0 
         IF (IAT /= 0) CYCLE  
         RETURN  
      END DO 
! cia dabar
      ISKR = IHSH - JA2 
      IF (ISKR > 0) THEN 
         DO JI = 1, ISKR 
            KK = IHSH + JA2 - 2 + JI 
            LK = J1QN1(KK,K) - 1 
            LD = J1QN2(KK,K) - 1 
            IF (ITTK(LK,LD,KA) == 0) IAT = 0 
            IF (IAT /= 0) CYCLE  
            RETURN  
         END DO 
      ENDIF 
! cia dabar
      RETURN  
      END SUBROUTINE RLSP0 
