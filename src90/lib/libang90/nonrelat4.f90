!
!     -------------------------------------------------------------
!      N O N R E L A T 4
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :      N'1 = N1 +- 1        *
!                                             N'2 = N2 +- 1        *
!                                             N'3 = N3 -+ 2        *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!                                                                  *
      SUBROUTINE NONRELAT4(IA, IB, IC, ID) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE OPERAT_C 
      USE MEDEFN_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:17:07  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ittk_I 
      USE nonrelat41_I 
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
      IF (IHSH <= 2) RETURN  
      IF (IA == IB) THEN 
         IF (ISOTOP == 1) THEN 
            IF (ITTK(LJ(IC),LJ(IA),1) == 0) RETURN  
            IF (ITTK(LJ(ID),LJ(IA),1) == 0) RETURN  
         ENDIF 
         CALL NONRELAT41 (IC, ID, IA, 1, IA, IB, IC, ID) 
      ELSE IF (IC == ID) THEN 
         IF (ISOTOP == 1) THEN 
            IF (ITTK(LJ(IA),LJ(IC),1) == 0) RETURN  
            IF (ITTK(LJ(IB),LJ(IC),1) == 0) RETURN  
         ENDIF 
         CALL NONRELAT41 (IA, IB, IC, 2, IA, IB, IC, ID) 
      ELSE 
         WRITE (6, '(A)') ' KLAIDA NONRELAT4  ' 
         STOP  
      ENDIF 
      RETURN  
      END SUBROUTINE NONRELAT4 
