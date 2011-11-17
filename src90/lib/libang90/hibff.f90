!
!     ------------------------------------------------------------------
!     H I B F F
!     ------------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                          December 1993   *
!
      SUBROUTINE HIBFF(JA1, JA2, JA3, JA4, NSLUO) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE TRK_C, BDS1=>BD1, BDS2=>BD2, BKS1=>BK1, BKS2=>BK2, IBDS1=>ID1, IBDS2&
         =>ID2, IBKS1=>IK1, IBKS2=>IK2 
      USE TRK2_C, BDS3=>BD3, BDS4=>BD4, BKS3=>BK3, BKS4=>BK4, IBDS3=>ID3, IBDS4&
         =>ID4, IBKS3=>IK3, IBKS4=>IK4 
      USE consts_C
      USE ribof_C
      USE ribols3_C
      USE ribols_C
      USE ribolsf_C

!...Translated by Pacific-Sierra Research 77to90  4.3E  10:28:22  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE hibp31_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: JA1 
      INTEGER  :: JA2 
      INTEGER  :: JA3 
      INTEGER  :: JA4 
      INTEGER , INTENT(IN) :: NSLUO 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      CALL HIBP31 (JA1, BKS1, IBKS1, BDS1, IBDS1) 
      IF (NSLUO == 1) RETURN  
      CALL HIBP31 (JA2, BKS2, IBKS2, BDS2, IBDS2) 
      IF (NSLUO == 2) RETURN  
      CALL HIBP31 (JA3, BKS3, IBKS3, BDS3, IBDS3) 
      IF (NSLUO == 3) RETURN  
      CALL HIBP31 (JA4, BKS4, IBKS4, BDS4, IBDS4) 
      RETURN  
      END SUBROUTINE HIBFF 
