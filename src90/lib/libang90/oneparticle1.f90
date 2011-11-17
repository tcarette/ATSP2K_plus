!     ..........................................................   :
!                                                                  :
!          Block                                                   :
!                    O N E   -   T W O                             :
!                                                                  :
!     For Calculation Angular Momentum Coefficients for            :
!     Relativistic operators                                       :
!                                                                  :
!     Written by G. Gaigalas,                                      :
!                  Department  of  Computer Science,               :
!                  Vanderbilt University,  Nashville               :
!                                                   October 1996   :
!                                                                  :
!     ..........................................................   :
!        i)   one - particle operator                              :
!     ..........................................................   :
!
!     -------------------------------------------------------------
!      O N E P A R T I C L E 1
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                  N'2 = N2        *
!     Written by G. Gaigalas,                                      *
!     Universite Libre de Bruxelles, Belgium         October 1995  *
!
      SUBROUTINE ONEPARTICLE1(K1, K2, IA, XXX) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE MEDEFN_C 
      USE DIAGNL_C 
      USE CONSTS_C 
      USE TRK_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:22:07  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
       
      USE rlsp0_I 
      USE rlsp1_I 
      USE hibff_I 
      USE w1_I 
      USE savenon_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: K1 
      INTEGER , INTENT(IN) :: K2 
      INTEGER  :: IA 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(2) :: K 
      INTEGER :: I, IAT, J, INUM, LA 
      REAL(DOUBLE) :: C1, REC, A1, W, RECLS 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
      C1 = ZERO 
      K(1) = K1 
      K(2) = K2 
      DO I = 1, 3 
         IF (I == 1) THEN 
            CALL RLSP0 (I, IA, IA, 0, IAT) 
         ELSE 
            J = I - 1 
            CALL RLSP0 (I, IA, IA, 2*K(J), IAT) 
         ENDIF 
         IF (IAT /= 0) CYCLE  
         RETURN  
      END DO 
      DO I = 2, 3 
         J = I - 1 
         CALL RLSP1 (I, IA, 2*K(J), 0, IAT, REC) 
         IF (IAT /= 0) CYCLE  
         RETURN  
      END DO 
      CALL HIBFF (IA, IA, IA, IA, 1) 
      CALL XXX (IK1(3), ID1(3), INUM, A1) 
      IF (DABS(A1) < EPS) RETURN  
      LA = IJFUL(IA) 
      CALL W1 (IK1, BK1, ID1, BD1, K(1), K(2), HALF, (-HALF), W) 
      RECLS = 1 
      DO I = 2, 3 
         J = I - 1 
         CALL RLSP1 (I, IA, 2*K(J), 1, IAT, REC) 
         RECLS = RECLS*REC 
      END DO 
      C1 = A1*W*RECLS/DSQRT(DBLE((2*K(1)+1)*(2*K(2)+1))) 
      IF (DABS(C1) > EPS) CALL SAVENON (INUM, C1, 0, 0, LA, 0, LA, JA, JB, 0) 
      RETURN  
      END SUBROUTINE ONEPARTICLE1 
