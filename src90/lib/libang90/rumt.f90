!
!     -------------------------------------------------------------
!      R U M T
!     -------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE RUMT(KNT, LL, LQ, LS, L) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE MT_C 
      USE SKMT2_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:57:38  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rumt67_I 
      USE jthn_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: KNT 
      INTEGER , INTENT(IN) :: LL 
      INTEGER  :: LQ 
      INTEGER  :: LS 
      INTEGER  :: L 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KT, KNTMIN, NR 
!-----------------------------------------------
      IF (LL < 3) THEN 
         KT = MT(KNT) 
      ELSE IF (LL == 3) THEN 
         IF (KNT > 300) THEN 
            KNTMIN = KNT - 300 
            KT = MTF(KNTMIN) 
         ELSE 
            CALL RUMT67 (KNT, NR, LQ, LS, L) 
            RETURN  
         ENDIF 
      ELSE 
         KT = MT3(KNT) 
      ENDIF 
      LQ = JTHN(KT,3,100) 
      LS = JTHN(KT,2,100) 
      L = JTHN(KT,1,100) 
      RETURN  
      END SUBROUTINE RUMT 
