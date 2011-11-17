!
!     --------------------------------------------------------------
!     R L S P 3 2
!     --------------------------------------------------------------
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE RLSP32(K, JA1, JA2, JA3, K1, K2, K3, K4, KA, IRE, IAT, REC) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:24:23  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dlsa3_I 
      USE dlsa4_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: K 
      INTEGER  :: JA1 
      INTEGER  :: JA2 
      INTEGER  :: JA3 
      INTEGER  :: K1 
      INTEGER  :: K2 
      INTEGER  :: K3 
      INTEGER  :: K4 
      INTEGER  :: KA 
      INTEGER  :: IRE 
      INTEGER  :: IAT 
      REAL(DOUBLE) , INTENT(OUT) :: REC 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: RE 
!-----------------------------------------------
!
      REC = ONE 
      IF (JA3 - JA2 > 1) THEN 
         IAT = 0 
         CALL DLSA3 (K, JA2, JA3, K3, IRE, IAT, RE) 
         IF (IAT == 0) RETURN  
         REC = RE*REC 
      ENDIF 
!
      IAT = 0 
      CALL DLSA4 (K, JA1, JA2, K1, K2, K3, IRE, IAT, RE) 
      IF (IAT == 0) RETURN  
      REC = RE*REC 
!
      IAT = 0 
      CALL DLSA4 (K, JA2, JA3, K3, K4, KA, IRE, IAT, RE) 
      IF (IAT == 0) RETURN  
      REC = RE*REC 
      RETURN  
      END SUBROUTINE RLSP32 
