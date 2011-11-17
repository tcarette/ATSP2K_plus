!
!     ------------------------------------------------------------------
!      N O N R E L A T
!     ------------------------------------------------------------------
!
!     THE ROUTINE EVALUATES THE NON-RELATIVISTIC HAMILTONIAN
!     WITH ORTHOGONAL ORBITALS
!
!
!     THE MATRIX ELEMENT OF THE TWO-ELECTRON POTENTIAL BETWEEN TWO
!     STATES (LABELLED 1 AND 2) MAY BE EXPRESSED AS A SUM OF WEIGHTED
!     RK (SLATER) INTEGRALS.  THIS SUBROUTINE, TOGETHER WITH THOSE
!     CALLED BY IT, DETERMINES THESE WEIGHTS, WHICH ARISE FROM AN
!     INTEGRATION OVER THE ANGULAR AND SPIN CO-ORDINATES
!     FOR DETAILS, SEE   U. FANO, PHYS. REV.,140,A67,(1965)
!
!     THE =INTERACTING= SHELLS ARE DESIGNATED  IRHO,ISIG,IRHOP,ISIGP.
!     THE FIRST TWO REFER TO THE L.H.S. OF     (PSI/V/PSIP)     , WHILE
!     THE SECOND TWO REFER TO THE R.H.S.  FOR DIAGONAL AND CERTAIN OFF-
!     DIAGONAL MATRIX ELEMENTS, THESE MAY NOT BE UNIQUE, AND EACH
!     POSSIBILITY MUST BE CONSIDERED IN TURN
!     THE CONDITION =IRHO .LE. ISIG ,  IRHOP .LE. ISIGP=  IS TO BE
!     SATISFIED
!
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE NONRELAT 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE INFORM_C 
      USE DEBUG_C 
      USE MEDEFN_C 
      USE DIAGNL_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:17:07  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE coupling_I 
      USE nonrelat2_I 
      USE nonrelat4_I 
      USE nonrelat5_I 
      USE nonrelat3_I 
      USE nonrelat1_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IZERO, IX, IRHO, ISIG, IRHOP, ISIGP, IGGG, J, N, N1, IRSTO, &
         IRPSTO, K1, IINE, IIRE, ISTO, K2 
!-----------------------------------------------
    5 FORMAT(/,/,10X,' IRHO =',I3,4X,' ISIG =',I3,4X,' IRHOP =',I3,3X,&
         ' ISIGP =',I3) 
      IZERO = 0 
      IX = 0 
      IRHO = 0 
      ISIG = 0 
      IRHOP = 0 
      ISIGP = 0 
      IGGG = 0 
      DO J = 1, IHSH 
         N = NOSH1(J) - NOSH2(J) 
         IF (IABS(N) > 2) RETURN  
         IF (N > 0) THEN 
            IF (N == 1) THEN 
               ISIG = J 
               IF (IRHO == 0) IRHO = J 
               IX = IX + 1 
            ELSE 
               IRHO = J 
               IF (IABS(N) == 2) IGGG = 20 + IGGG 
               IX = IX + 2 
            ENDIF 
         ELSE IF (N < 0) THEN 
            IF (N + 1 == 0) THEN 
               ISIGP = J 
               IF (IRHOP == 0) IRHOP = J 
               IX = IX + 1 
            ELSE 
               IRHOP = J 
               IF (IABS(N) == 2) IGGG = 20 + IGGG 
               IX = IX + 2 
            ENDIF 
         ENDIF 
      END DO 
!
!     IX MEASURES THE TOTAL NUMBER OF ELECTRONS IN EITHER CONFIGURATION
!     WHICH DO NOT OCCUR IN THE OTHER.  THEN IF  IX  IS GREATER THAN 4,
!     ORTHOGONALITY OF THE ORBITALS PREVENTS A NON-ZERO MATRIX ELEMENT.
!     IF  IX  IS LESS THAN 4, THEN WE DIVIDE IX BY 2 AND NOW IX MEASURES
!     THE NUMBER OF ELECTRONS WHICH HAVE BEEN CHANGED IN GOING FROM PSI
!     TO PSIP.  IF NOW IX=0, WE HAVE A DIAGONAL MATRIX ELEMENT.  RHO AND
!     SIG MAY TAKE ON ANY VALUES LESS THAN IHSH.  IF IX=1, ONE INTER-
!     ACTING SHELL ON EACH SIDE IS FIXED, WHILE THE OTHER MAY VARY.  IF
!     IX=2, ALL INTERACTING SHELLS ARE DETERMINED
!
      IF (IX > 4) RETURN  
      CALL COUPLING (JA, JB) 
      N1 = 2*IHSH - 1 
      IF (J1QN1(N1,2) - J1QN2(N1,2) /= 0) RETURN  
      IF (J1QN1(N1,3) - J1QN2(N1,3) /= 0) RETURN  
      IX = IX/2 
      IF (IX - 1 > 0) THEN 
!
! === UNIQUE SPECIFICATION OF INTERACTING SHELLS
!
         IF (ISIG == 0) ISIG = IRHO 
         IF (ISIGP == 0) ISIGP = IRHOP 
         IF (IBUG2 > 0) WRITE (IWRITE, 5) IRHO, ISIG, IRHOP, ISIGP 
         IF (IGGG == 40) THEN 
            CALL NONRELAT2 (IRHO, IRHOP) 
         ELSE IF (IGGG == 20) THEN 
            CALL NONRELAT4 (IRHO, ISIG, IRHOP, ISIGP) 
         ELSE 
            CALL NONRELAT5 (IRHO, ISIG, IRHOP, ISIGP) 
         ENDIF 
      ELSE IF (IX - 1 == 0) THEN 
!
! === ONE INTERACTING SHELL SPECIFIED ON EACH SIDE. SUMMATION OVER OTHER
!
         IRSTO = IRHO 
         IRPSTO = IRHOP 
         DO K1 = 1, IHSH 
            IINE = 1 
            IIRE = 1 
            IF (NOSH1(K1) == 0) IIRE = 0 
            ISIG = K1 
            IF (NOSH2(K1) == 0) IIRE = 0 
            ISIGP = K1 
            IRHO = IRSTO 
            IRHOP = IRPSTO 
!
!     ORTHOGONALITY OF THE ORBITALS REQUIRES THAT THE VARYING INTER-
!     ACTING SHELL BE THE SAME ON BOTH SIDES OF THE MATRIX ELEMENT
!
! --- IRHO.LE.ISIG,   IRHOP.LE.ISIGP
!
            IF (IRHO > ISIG) THEN 
               ISTO = IRHO 
               IRHO = ISIG 
               ISIG = ISTO 
            ELSE IF (IRHO == ISIG) THEN 
               IF (NOSH1(ISIG) == 1) IIRE = 0 
               IF (NOSH1(ISIG) == 0) IINE = 0 
            ENDIF 
            IF (IRHOP > ISIGP) THEN 
               ISTO = IRHOP 
               IRHOP = ISIGP 
               ISIGP = ISTO 
            ELSE IF (IRHOP == ISIGP) THEN 
               IF (NOSH2(ISIGP) == 1) IIRE = 0 
               IF (NOSH2(ISIGP) == 0) IINE = 0 
            ENDIF 
            IF (IINE /= 1) CYCLE  
            IF (IBUG2 > 0) WRITE (IWRITE, 5) IRHO, ISIG, IRHOP, ISIGP 
            CALL NONRELAT3 (IRHO, ISIG, IRHOP, ISIGP, IIRE) 
         END DO 
      ELSE IF (IX - 1 < 0) THEN 
!
! === NO INTERACTING SHELLS SPECIFIED
!     SUMMATION OVER ALL POSSIBLE COMBINATIONS
!     IN THIS CASE, ORTHOGONALITY OF ORBITALS PRECLUDES ALL CASES
!     EXCEPT  IRHO=IRHOP    AND    ISIG=ISIGP
!
         DO K1 = 1, IHSH 
            IF (NOSH1(K1) == 0) CYCLE  
            ISIG = K1 
            DO K2 = 1, K1 
               IIRE = 1 
               IF (NOSH1(K2) == 0) CYCLE  
               IRHO = K2 
               IF (IRHO == ISIG) THEN 
                  IF (NOSH1(ISIG) == 1) IIRE = 0 
               ENDIF 
               IRHOP = IRHO 
               ISIGP = ISIG 
               IF (IBUG2 > 0) WRITE (IWRITE, 5) IRHO, ISIG, IRHOP, ISIGP 
               CALL NONRELAT1 (IRHO, ISIG, IIRE) 
            END DO 
         END DO 
      ENDIF 
      RETURN  
      END SUBROUTINE NONRELAT 
