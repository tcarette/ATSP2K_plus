!
!     -------------------------------------------------------------
!      L M A T R I X 1
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                  N'2 = N2        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville           September 1997   *
!                                                                  *
!
      SUBROUTINE LMATRIX1(IA) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
      USE DIAGNL_C 
      USE TRK_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:58:20  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE recoupls0_I 
      USE hibff_I 
      USE savenon_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IA 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IN, LA 
      REAL(DOUBLE) :: A, B 
!-----------------------------------------------
      IF (JA /= JB) THEN 
         DO IN = 1, 3 
            IF (RECOUPLS0(IN,IA,IA,IA,IA,0)) CYCLE  
            RETURN  
         END DO 
      ENDIF 
      CALL HIBFF (IA, IA, IA, IA, 1) 
      IF (IK1(1) /= ID1(1)) RETURN  
      A = -DBLE(ID1(4))/SQRT(DBLE(4*IK1(3)+2)) 
      B = A*HALF*DSQRT(DBLE(4*LJ(IA)+2)) 
      IF (DABS(B) > EPS) THEN 
         IF (JA /= JB) THEN 
            WRITE (*, '(A,I5,A,I5)') 'Configurations', JA, 'and', JB, &
               'are identical' 
            STOP  
         ELSE 
            LA = IJFUL(IA) 
            CALL SAVENON (4, B, LJ(IA), 0, LA, 0, LA, JA, JB, 0) 
         ENDIF 
      ENDIF 
      RETURN  
      END SUBROUTINE LMATRIX1 
