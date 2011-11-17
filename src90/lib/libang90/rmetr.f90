!
!     ------------------------------------------------------------------
!       R M E T R
!     ------------------------------------------------------------------
!
! --- evaluates the angular part of the one-electron transition reduced
!     matrix element. See equations (4) and (5) of paper II.
!
      REAL(KIND(0.0D0)) FUNCTION RMETR (L1, L2) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE EMS_C 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:13:23  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rme_I 
      USE sixj_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: L1 
      INTEGER  :: L2 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LL2, L3 
      REAL(DOUBLE) :: W 
!-----------------------------------------------
      IF (IFL == 1) THEN 
!       electric multipole transitions
         RMETR = RME(L1,L2,LAM) 
         IF (MOD(L2 - L1 + LAM,4) /= 0) RMETR = -RMETR 
      ELSE 
!       magnetic multipole transitions
         RMETR = RME(L1,L2,LAM - 1) 
         IF (DABS(RMETR) > EPS) THEN 
            IF (MOD(L2 - L1 + LAM - 1,4) /= 0) RMETR = -RMETR 
            IF (IFL == 3) THEN 
               LL2 = L2 + L2 
               L3 = LAM + LAM 
               CALL SIXJ (2, LL2, LL2, L1 + L1, L3, L3 - 2, 1, W) 
               RMETR = -RMETR*W*SQRT(DBLE(L2*(L2 + 1)*(LL2 + 1)*(L3 + 1))) 
            ELSE 
               RMETR = RMETR*SQRT(HALF*THREE) 
            ENDIF 
         ENDIF 
      ENDIF 
      RETURN  
      END FUNCTION RMETR 
