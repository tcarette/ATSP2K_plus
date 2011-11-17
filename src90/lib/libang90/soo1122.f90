!
!     -------------------------------------------------------------
!      S O O 1 1 2 2
!     -------------------------------------------------------------
!                                                                  *
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Universite Libre de Bruxelles, Belgium         October 1995  *
!
      SUBROUTINE SOO1122(IG, KL, IA, IB, XXX) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
      USE DIAGNL_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:15:20  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
       
      USE sooa_I 
      USE soob_I 
      USE savenon_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IG 
      INTEGER  :: KL 
      INTEGER  :: IA 
      INTEGER  :: IB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L1, L2, LA, LB, KLP, KLM 
      REAL(DOUBLE) :: AP1, G1, G2, A1, C1, CC1, C2, CC2, C3, CC3, C4, CC4, C55&
         , CC55, C66, CC66, C7, C8, G4, C5, C6, G3, A2 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
      L1 = LJ(IA) 
      L2 = LJ(IB) 
      CALL SOOA (L1, L1, L2, L2, KL, AP1) 
      IF (DABS(AP1) < EPS) RETURN  
      LA = IJFUL(IA) 
      LB = IJFUL(IB) 
      KLP = KL + 1 
      KLM = KL - 1 
      G1 = ZERO 
      G2 = ZERO 
      IF (KLM >= 0) THEN 
         CALL SOOB (1, L1, L2, KLM, KL, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL XXX (KLM, 1, KL, 0, 1, IA, IB, IB, IB, IB, C1, CC1) 
            CALL XXX (KLM, 0, KL, 1, 1, IA, IB, IB, IB, IB, C2, CC2) 
            C1 = C1*A1 
            C2 = TWO*C2*A1 
            G1 = C1 + C2 
         ENDIF 
      ENDIF 
      IF (KLP <= IG) THEN 
         CALL SOOB (1, L1, L2, KLP, KL, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL XXX (KLP, 1, KL, 0, 1, IA, IB, IB, IB, IB, C3, CC3) 
            CALL XXX (KLP, 0, KL, 1, 1, IA, IB, IB, IB, IB, C4, CC4) 
            C3 = C3*A1 
            C4 = TWO*C4*A1 
            G2 = C3 + C4 
         ENDIF 
      ENDIF 
      CALL XXX (KL, 1, KL, 0, 1, IA, IB, IB, IB, IB, C55, CC55) 
      CALL XXX (KL, 0, KL, 1, 1, IA, IB, IB, IB, IB, C66, CC66) 
!
!    ( k k )  tensorine struktura        V   integralas
!
      CALL SOOB (2, L1, L2, KL, KL, A1) 
      C7 = C55*A1 
      C8 = TWO*C66*A1 
      G4 = AP1*(C7 + C8) 
      IF (DABS(G4) > EPS) CALL SAVENON (7, G4, KL - 1, LA, LA, LB, LB, JA, JB, &
         0) 
!
!    ( k  k+1)  + ( k  k )
!
      CALL SOOB (1, L1, L2, KL, KL, A1) 
      C5 = C55*A1 
      C6 = TWO*C66*A1 
      G3 = C5 + C6 
      IF (KL - 1 >= 0) THEN 
         A1 = AP1*(G1 + DBLE(KL + 1)*G3) 
         IF (DABS(A1) > EPS) CALL SAVENON (6, A1, KL - 2, LA, LA, LB, LB, JA, &
            JB, 0) 
      ENDIF 
!
!    ( k  k-1)  + ( k  k )
!
      A2 = AP1*(G2 - DBLE(KL)*G3) 
      IF (DABS(A2) > EPS) CALL SAVENON (6, A2, KL, LA, LA, LB, LB, JA, JB, 0) 
      RETURN  
      END SUBROUTINE SOO1122 
