!
!     ------------------------------------------------------------------
!       N T A B 1
!     ------------------------------------------------------------------
!
      INTEGER FUNCTION NTAB1 (NELCTS, K) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  09:53:21  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NELCTS 
      INTEGER , INTENT(IN) :: K 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(0:9) :: IROW 
      INTEGER :: NPAR, L, LHALF 
!-----------------------------------------------
      DATA IROW/ 0, 2, 5, 10, 12, 14, 16, 18, 20, 22/  
!
!     THIS SUBROUTINE CALCULATES THE ROW OF NTAB CORRESPONDING TO THE
!     PARENTS WHICH MAY GIVE RISE TO THE TERM ASSOCIATED WITH SHELL
!     LAMBDA .  E.G. IF WE SEEK THE ROW OF NTAB CONTAINING THE PARENTS
!     OF ONE OF THE P**3 TERMS, THE ROW = VALUE OF NTAB1 IS THAT
!     CONTAINING THE P**2 TERMS
!
!     USE IS MADE OF THE FACT THAT THE LIST OF POSSIBLE PARENTS (SEE
!     WHITE - ATOMIC SPECTRA - APPENDIX)  IS SYMMETRICAL ABOUT THE
!     CONFIGURATION L**(2L+1)
!
!
! --- FOR ONE ELECTRON IN A TERM, THE PARENT IS ALWAYS A SINGLET S TERM
!
      IF (NELCTS == 1) THEN 
         NTAB1 = 2 
      ELSE 
         NPAR = NELCTS - 1 
         L = K - 1 
         LHALF = 2*L + 1 
         IF (NPAR > LHALF) NPAR = 2*LHALF - NPAR 
         NTAB1 = IROW(L) + NPAR 
      ENDIF 
      RETURN  
      END FUNCTION NTAB1 
