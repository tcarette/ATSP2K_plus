!
!     --------------------------------------------------------------
!      N O N T R A N S
!     --------------------------------------------------------------
!
!     THE ROUTINE EVALUATES THE TRANSITION OPERATORS WITH          *
!     ORTHOGONAL ORBITALS                                          *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville           September 1997   *
!
      SUBROUTINE NONTRANS(KA, KB, CL, CV) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE SAVECOM_C 
      use ems_C
      use consts_C
      use medefn_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:27:41  11/20/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE transition_I 
      USE ittk_I 
      USE oneparticle2_I 
      USE oneparticle1_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: KA 
      INTEGER  :: KB 
      REAL(DOUBLE) , INTENT(OUT) :: CL 
      REAL(DOUBLE) , INTENT(OUT) :: CV 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IX, IRHO, IRHOP, J, N, K1 
!-----------------------------------------------
!      LOGICAL REL,VOK
!      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
!      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
!      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
!     :     J1QN2(31,3),IJFUL(16)
      CL2 = ZERO 
      CV2 = ZERO 
      CL = CL2 
      CV = CV2 
      IX = 0 
      IRHO = 0 
      IRHOP = 0 
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
      IF (IX == 2) THEN 
         LRHO = LJ(IRHO) 
         LSIG = LJ(IRHOP) 
         IF (IFL == 1) THEN 
            IF (ITTK(LRHO,LSIG,LAM) == 0) RETURN  
         ELSE 
            IF (ITTK(LRHO,LSIG,LAM - 1) == 0) RETURN  
         ENDIF 
         CALL ONEPARTICLE2 (KA, KB, IRHO, IRHOP, TRANSITION) 
      ELSE IF (IX == 0) THEN 
         DO K1 = 1, IHSH 
            IF (NOSH1(K1) == 0) CYCLE  
            LRHO = LJ(K1) 
            LSIG = LRHO 
            IF (IFL == 1) THEN 
               IF (ITTK(LRHO,LSIG,LAM) /= 0) CALL ONEPARTICLE1 (KA, KB, K1, &
                  TRANSITION) 
            ELSE 
               IF (ITTK(LRHO,LSIG,LAM - 1) /= 0) CALL ONEPARTICLE1 (KA, KB, K1&
                  , TRANSITION) 
            ENDIF 
         END DO 
      ENDIF 
      CL = CL2 
      CV = CV2 
      RETURN  
      END SUBROUTINE NONTRANS 
