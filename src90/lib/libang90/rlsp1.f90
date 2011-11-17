!
!     -------------------------------------------------------------
!      R L S P 1
!     -------------------------------------------------------------
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE RLSP1(K, JA1, KA, IRE, IAT, REC) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:24:23  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dlsa5_I 
      USE dlsa1_I 
      USE dlsa3_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: K 
      INTEGER  :: JA1 
      INTEGER  :: KA 
      INTEGER  :: IRE 
      INTEGER  :: IAT 
      REAL(DOUBLE) , INTENT(OUT) :: REC 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ISKR 
      REAL(DOUBLE) :: RE 
!-----------------------------------------------
      IAT = 1 
      REC = ONE/SQRT(DBLE(J1QN1(JA1,K))) 
      IF (IHSH == 1) RETURN  
      IAT = 0 
      IF (IRE /= 0) THEN 
         IF (KA == 0) THEN 
            IAT = 1 
            RETURN  
         ENDIF 
      ENDIF 
      IF (IHSH /= 2) THEN 
         CALL DLSA5 (K, JA1, KA, IRE, IAT, RE) 
         REC = RE*REC 
         IF (IAT == 0) RETURN  
         IF (JA1 == IHSH) RETURN  
         IAT = 0 
      ENDIF 
      CALL DLSA1 (K, JA1, KA, IRE, IAT, RE) 
      IF (IAT == 0) RETURN  
      REC = RE*REC 
      IF (IHSH == 2) RETURN  
      ISKR = IHSH - JA1 
      IF (JA1 == 1) ISKR = IHSH - 1 - JA1 
      IF (ISKR <= 1) RETURN  
      IAT = 0 
      CALL DLSA3 (K, JA1, IHSH, KA, IRE, IAT, RE) 
      REC = RE*REC 
      RETURN  
      END SUBROUTINE RLSP1 
