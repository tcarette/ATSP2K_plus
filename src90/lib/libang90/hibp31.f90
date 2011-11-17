!
!     ------------------------------------------------------------------
!     H I B P 3 1
!     ------------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE HIBP31(I, BK, IBK, BD, IBD) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:28:22  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE numterf_I 
      USE numter_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I 
      INTEGER  :: IBK(7) 
      INTEGER  :: IBD(7) 
      REAL(DOUBLE) , INTENT(OUT) :: BK(3) 
      REAL(DOUBLE) , INTENT(OUT) :: BD(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IREZ, I2NK, I2ND 
!-----------------------------------------------
      IBK(2) = NJ(I) 
      IBD(2) = IBK(2) 
      IBK(3) = LJ(I) 
      IBD(3) = IBK(3) 
      IBK(4) = NOSH1(I) 
      IBD(4) = NOSH2(I) 
      IBK(5) = J1QN1(I,2) - 1 
      IBD(5) = J1QN2(I,2) - 1 
      IBK(6) = J1QN1(I,3) - 1 
      IBD(6) = J1QN2(I,3) - 1 
      IREZ = 0 
      IF (IBK(3) < 3) THEN 
         IBK(7) = 2*LJ(I) + 1 - J1QN1(I,1) 
      ELSE IF (IBK(4) < 3) THEN 
         IBK(7) = 2*LJ(I) + 1 - J1QN1(I,1) 
      ELSE IF (IBK(4) < 14) THEN 
         IREZ = 1 
         I2NK = J1QN1(I,1) 
         IBK(1) = NUMTERF(I2NK,IBK(6),IBK(5),IBK(4),IBK(7)) 
      ELSE 
         IBK(7) = 2*LJ(I) + 1 - J1QN1(I,1) 
      ENDIF 
      BK(1) = HALF*DBLE(IBK(7)) 
      BK(2) = HALF*DBLE(IBK(6)) 
      BK(3) = -HALF*DBLE(2*LJ(I)+1-NOSH1(I)) 
      IF (IREZ == 0) IBK(1) = NUMTER(IBK(7),IBK(6),IBK(5),IBK(3),IBD(4),IBK(4)) 
      IREZ = 0 
      IF (IBD(3) < 3) THEN 
         IBD(7) = 2*LJ(I) + 1 - J1QN2(I,1) 
      ELSE IF (IBD(4) < 3) THEN 
         IBD(7) = 2*LJ(I) + 1 - J1QN2(I,1) 
      ELSE IF (IBD(4) < 14) THEN 
         IREZ = 1 
         I2ND = J1QN2(I,1) 
         IBD(1) = NUMTERF(I2ND,IBD(6),IBD(5),IBD(4),IBD(7)) 
      ELSE 
         IBD(7) = 2*LJ(I) + 1 - J1QN2(I,1) 
      ENDIF 
      BD(1) = HALF*DBLE(IBD(7)) 
      BD(2) = HALF*DBLE(IBD(6)) 
      BD(3) = -HALF*DBLE(2*LJ(I)+1-NOSH2(I)) 
      IF (IREZ == 0) IBD(1) = NUMTER(IBD(7),IBD(6),IBD(5),IBD(3),IBK(4),IBD(4)) 
      RETURN  
      END SUBROUTINE HIBP31 
