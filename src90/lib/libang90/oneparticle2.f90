!                                                                  *
!     -------------------------------------------------------------
!      O N E P A R T I C L E 2
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1 -+ 1   *
!                                                  N'2 = N2 +- 1   *
!     Written by G. Gaigalas,                                      *
!     Universite Libre de Bruxelles, Belgium         October 1995  *
!
      SUBROUTINE ONEPARTICLE2(K1, K2, IIA, IIB, XXX) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE TRK_C 
      USE MEDEFN_C 
      USE DIAGNL_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:22:07  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rlsp0_I 
      USE rlsp2_I 
      USE hibff_I 
      USE c0t5s_I 
      USE sls_I 
      USE savenon_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: K1 
      INTEGER , INTENT(IN) :: K2 
      INTEGER  :: IIA 
      INTEGER  :: IIB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(2) :: K 
      INTEGER :: IA, IB, I, IAT, J, INUM, IQMM1, IQMM2, LA, LB, JIB, IFAZ 
      REAL(DOUBLE) :: C1, REC, A1, QM1, QM2, A2, A3, S1, S2, RECLS 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
      C1 = ZERO 
      K(1) = K1 
      K(2) = K2 
      IA = MIN0(IIA,IIB) 
      IB = MAX0(IIA,IIB) 
      DO I = 1, 3 
         IF (I == 1) THEN 
            CALL RLSP0 (I, IA, IB, 0, IAT) 
         ELSE 
            J = I - 1 
            CALL RLSP0 (I, IA, IB, 2*K(J), IAT) 
         ENDIF 
         IF (IAT /= 0) CYCLE  
         RETURN  
      END DO 
      DO I = 2, 3 
         J = I - 1 
         IF (I == 3) THEN 
            CALL RLSP2 (I, IA, IB, 1, 1, 2*K(J), 0, IAT, REC) 
         ELSE 
            CALL RLSP2 (I, IA, IB, 2*LJ(IA), 2*LJ(IB), 2*K(J), 0, IAT, REC) 
         ENDIF 
         IF (IAT /= 0) CYCLE  
         RETURN  
      END DO 
      CALL HIBFF (IIA, IIB, IIA, IIA, 2) 
      CALL XXX (IK1(3), IK2(3), INUM, A1) 
      IF (DABS(A1) < EPS) RETURN  
      QM1 = HALF 
      QM2 = -HALF 
      IQMM1 = QM1 + QM1 + TENTH*QM1 
      IF (IK1(4) /= ID1(4) + IQMM1) RETURN  
      IQMM2 = QM2 + QM2 + TENTH*QM2 
      IF (IK2(4) /= ID2(4) + IQMM2) RETURN  
      CALL C0T5S (BD1(1), BD1(3), QM1, BK1(1), BK1(3), A2) 
      IF (DABS(A2) < EPS) RETURN  
      CALL C0T5S (BD2(1), BD2(3), QM2, BK2(1), BK2(3), A3) 
      IF (DABS(A3) < EPS) RETURN  
      CALL SLS (IK1(3), IK1(1), IK1(7), IK1(5), IK1(6), ID1(1), ID1(7), ID1(5)&
         , ID1(6), S1) 
      IF (DABS(S1) < EPS) RETURN  
      CALL SLS (IK2(3), IK2(1), IK2(7), IK2(5), IK2(6), ID2(1), ID2(7), ID2(5)&
         , ID2(6), S2) 
      IF (DABS(S2) < EPS) RETURN  
      LA = IJFUL(IIA) 
      LB = IJFUL(IIB) 
      A1 = S1*S2*A1*A2*A3 
      RECLS = 1 
      DO I = 2, 3 
         J = I - 1 
         IF (I == 3) THEN 
            CALL RLSP2 (I, IA, IB, 1, 1, 2*K(J), 1, IAT, REC) 
         ELSE 
            CALL RLSP2 (I, IA, IB, 2*LJ(IA), 2*LJ(IB), 2*K(J), 1, IAT, REC) 
         ENDIF 
         RECLS = RECLS*REC 
      END DO 
      C1 = A1*RECLS/DSQRT(DBLE((2*K(1)+1)*(2*K(2)+1)*(IK1(7)+1)*(IK2(7)+1))) 
      JIB = IB - 1 
      IFAZ = SUM(NOSH1(IA:JIB)) 
      IFAZ = IFAZ + 1 
      IF (MOD(IFAZ,2) /= 0) C1 = -C1 
      IF (IA /= IIA) THEN 
         IFAZ = IK1(3) + IK2(3) - K(1) - K(2) 
         IF (MOD(IFAZ,2) /= 0) C1 = -C1 
      ENDIF 
      IF (DABS(C1) > EPS) CALL SAVENON (INUM, C1, 0, 0, LA, 0, LB, JA, JB, 0) 
      RETURN  
      END SUBROUTINE ONEPARTICLE2 
