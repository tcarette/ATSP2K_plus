!
!     ------------------------------------------------------------------
!      L M A T R I X
!     ------------------------------------------------------------------
!
!     THE ROUTINE EVALUATES THE ONE-ELECTRON NON-RELATIVISTIC
!     HAMILTONIAN WITH ORTHOGONAL ORBITALS
!
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville           September 1997   *
!
      SUBROUTINE LMATRIX 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE MEDEFN_C 
      USE DIAGNL_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:57:57  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE lmatrix2_I 
      USE lmatrix1_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IX, J, N, IRHO, IRHOP, N1, K1 
!-----------------------------------------------
      IX = 0 
      DO J = 1, IHSH 
         N = NOSH1(J) - NOSH2(J) 
         IF (IABS(N) > 1) RETURN  
         IF (N == 1) THEN 
            IRHO = J 
            IX = IX + 1 
         ELSE IF (N + 1 == 0) THEN 
            IRHOP = J 
            IX = IX + 1 
         ENDIF 
      END DO 
      IF (IX > 2) RETURN  
      N1 = 2*IHSH - 1 
      IF (J1QN1(N1,2) - J1QN2(N1,2) /= 0) RETURN  
      IF (J1QN1(N1,3) - J1QN2(N1,3) /= 0) RETURN  
      IF (IX == 2) THEN 
         IF (LJ(IRHOP) /= LJ(IRHO)) RETURN  
         CALL LMATRIX2 (IRHOP, IRHO) 
      ELSE IF (IX == 0) THEN 
         DO K1 = 1, IHSH 
            IF (NOSH1(K1) == 0) CYCLE  
            CALL LMATRIX1 (K1) 
         END DO 
      ENDIF 
      RETURN  
      END SUBROUTINE LMATRIX 
