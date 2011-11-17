!
!     --------------------------------------------------------------
!       F L I N E
!     --------------------------------------------------------------
!
!     evaluates the line factor for the (J,J') pair.
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville           September 1997   *
!
      REAL(KIND(0.0D0)) FUNCTION FLINE (K) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE LSJ_C 
      use ems_C
      use medefn_C
      use mult_C
      use consts_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:52:56  11/20/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE sixj_I 
      USE ninels_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: K 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IVL, IVR, I2HSH, LAM2, L4, IN 
      REAL(DOUBLE) :: F 
!-----------------------------------------------
 
      IVL = JVL(K) 
      IVR = JVR(K) 
      I2HSH = 2*IHSH - 1 
      LL1 = J1QN1(I2HSH,2) - 1 
      LL2 = J1QN2(I2HSH,2) - 1 
      IS1 = J1QN1(I2HSH,3) - 1 
      IS2 = J1QN2(I2HSH,3) - 1 
      LAM2 = LAM + LAM 
      IF (IFL /= 4) THEN 
         CALL SIXJ (LL1, IS1, IVL, IVR, LAM2, LL2, 1, F) 
         IF (MOD(LL1 + IS1 + IVR + LAM2,4) /= 0) F = -F 
         IF (IFL == 3) F = F/DBLE(LAM + 1) 
         FLINE = F 
      ELSE 
         L4 = LAM2 - 2 
         CALL NINELS (IVL, LL1, IS1, IVR, LL2, IS2, LAM2, L4, 2, 1, IN, F) 
         IF (IN == 1) THEN 
            CALL NINELS (IVL, LL1, IS1, IVR, LL2, IS2, LAM2, L4, 2, 0, IN, F) 
            FLINE = F*SQRT(DBLE(LAM2 + 1)) 
         ELSE 
            FLINE = ZERO 
         ENDIF 
      ENDIF 
      IF (MOD(LL1 + IS1 + LL2 + IS2 - IVL - IVR,4) /= 0) FLINE = -FLINE 
      RETURN  
      END FUNCTION FLINE 
