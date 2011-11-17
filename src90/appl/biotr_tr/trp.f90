!
!     ---------------------------------------------------------------
!        T R P
!     ---------------------------------------------------------------
!
!
      REAL(KIND(0.0D0)) FUNCTION TRP (END, LAM, DD) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      use param_C
      use inout_C
      use debug_C
      use dbg_C
!
!
!  *****  THIS FUNCTION CALCULATES THE FACTOR CONNECTING THE TRANSITION
!  *****  PROBABILITY (SEC-1) AND THE LINE STRENGTH (A.U.)
!
!  *****  WE ARE TAKING THE NUMERICAL FACTORS OF PAGES 437,439 OF SHORE
!  *****  AND MENZEL FOR M AND E TRANSITIONS, AND NOT THOSE OF PAGE 445
!  *****  DOING THAT, OUR DEFINITION OF LINE STRENGTHS DOES NOT CORRES-
!  *****  POND TO EQ.(8-15),CHAP.10 BUT RATHER TO THE SQUARE OF THE RME
!  *****  WITHOUT EXTRA FACTORS (DEF.OF SOBEL'MAN)
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:26:07  11/20/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: LAM 
      REAL(DOUBLE) , INTENT(IN) :: DD 
      CHARACTER , INTENT(IN) :: END 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L3, L4 
      REAL(DOUBLE), DIMENSION(10) :: DBLFAC 
      REAL(DOUBLE) :: C, ALPHA, RYD, PI, F, F2, R 
!-----------------------------------------------
!
!      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
!     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
!      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2),iut(2)
!      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
!      COMMON /DBG  /IBUGM
      DATA C, ALPHA, RYD, PI/ 2.99792458D10, 7.29735308D-03, 1.0973731534D05, &
         3.14159265359D00/  
!
!  *****  DBLFAC(I) = (2I+1)!!
!
      DATA DBLFAC/ 3.0000000000D00, 1.5000000000D01, 1.0500000000D02, &
         9.4500000000D02, 1.0395000000D04, 1.3513500000D05, 2.0270250000D06, &
         3.4459425000D07, 6.5472907500D08, 1.3749310575D10/  
      IF (IBUGM /= 0) WRITE (IWRITE, 501) END, LAM, DD 
  501 FORMAT(/,' ',' IN TRP ',A1,I2,'  DD =',D16.8) 
      IF (LAM > 10) THEN 
         WRITE (IWRITE, 100) 
  100    FORMAT(1X,' THE VECTOR DBLFAC IN TRP HAS TO BE EXTENDED') 
         STOP  
      ENDIF 
      L3 = LAM + LAM 
      L4 = L3 + 1 
      F = DBLE(L4*(LAM + 1))/LAM 
      F2 = DBLFAC(LAM) 
      IF (IBUGM /= 0) WRITE (IWRITE, 500) F, LAM, F2 
  500 FORMAT(8X,' F =',D16.8,' LAM =',I3,' F2 =',D16.8) 
      F2 = F2**2 
      F = F/F2 
      R = RYD**L3 
      F = F*PI*C/R 
      IF (END == 'E') THEN 
         F2 = D4/D2**L3 
         F2 = F2*ALPHA**L4 
      ELSE 
         F2 = D2**(-L3) 
         F2 = F2*ALPHA**(L3 + 3) 
      ENDIF 
      F = F*F2 
!
!  *****  MULTIPLY BY THE WAVENUMBER (CM-1) FUNCTION
!
      TRP = F*DD**L4 
      RETURN  
      END FUNCTION TRP 
