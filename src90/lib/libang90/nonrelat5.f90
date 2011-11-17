!
!     -------------------------------------------------------------
!      N O N R E L A T 5
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :    N'1 = N1 (+-) 1        *
!                                           N'2 = N2 (+-) 1        *
!                                           N'3 = N3 (+-) 1        *
!                                           N'4 = N4 (+-) 1        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE NONRELAT5(IA, IB, IC, ID) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE MEDEFN_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:17:07  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE nonrelat51_I 
      USE nonrelat52_I 
      USE nonrelat53_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IA 
      INTEGER  :: IB 
      INTEGER  :: IC 
      INTEGER  :: ID 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      IF (IHSH <= 3) RETURN  
      IF (IB < IC) THEN 
         CALL NONRELAT51 (IA, IB, IC, ID, 1) 
      ELSE IF (IA>ID .AND. IB>ID) THEN 
         CALL NONRELAT51 (IC, ID, IA, IB, 2) 
      ELSE IF (IB>IC .AND. IB<ID .AND. IA<IC) THEN 
         CALL NONRELAT52 (IA, IC, IB, ID, 1) 
      ELSE IF (IB>IC .AND. IB>ID .AND. IA>IC) THEN 
         CALL NONRELAT52 (IC, IA, ID, IB, 2) 
      ELSE IF (IB>IC .AND. IB>ID .AND. IA<IC) THEN 
         CALL NONRELAT53 (IA, IC, ID, IB, 1) 
      ELSE IF (IB>IC .AND. IB<ID .AND. IA>IC) THEN 
         CALL NONRELAT53 (IC, IA, IB, ID, 2) 
      ELSE 
         WRITE (6, '(A)') ' KLAIDA NONRELAT5  ' 
         STOP  
      ENDIF 
      RETURN  
      END SUBROUTINE NONRELAT5 
