!
!     -------------------------------------------------------------
!      S O O 1 2 1 2 P
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
      SUBROUTINE SOO1212P(IG, KL, IA, IB) 
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
      USE two12_I 
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
      INTEGER :: L1, L2, LB, LA, KLM, KLP 
      REAL(DOUBLE) :: AP1, A1, C1, C2, G1, GG1, C3, C4, G2, GG2, C55, C66, CC55&
         , CC66, G4, GG4, A2, G3, GG3 
!-----------------------------------------------
      L1 = LJ(IA) 
      L2 = LJ(IB) 
      CALL SOOA (L1, L2, L1, L2, KL, AP1) 
      IF (DABS(AP1) < EPS) RETURN  
      LB = IJFUL(IB) 
      LA = IJFUL(IA) 
      KLM = KL - 1 
      IF (KLM >= 0) THEN 
         CALL SOOB (1, L1, L1, KLM, KL, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO12 (KLM, 1, KL, 0, 1, IA, IB, C1) 
            CALL TWO12 (KLM, 0, KL, 1, 1, IA, IB, C2) 
            G1 = C1 + TWO*C2 
            G1 = A1*G1 
         ELSE 
            G1 = ZERO 
         ENDIF 
!
         CALL SOOB (1, L2, L2, KLM, KL, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO12 (KL, 0, KLM, 1, 1, IA, IB, C1) 
            CALL TWO12 (KL, 1, KLM, 0, 1, IA, IB, C2) 
            GG1 = C1 + TWO*C2 
            GG1 = A1*GG1 
         ELSE 
            GG1 = ZERO 
         ENDIF 
      ELSE 
         G1 = ZERO 
         GG1 = ZERO 
      ENDIF 
      KLP = KL + 1 
      IF (KLP <= IG) THEN 
         CALL SOOB (1, L1, L1, KLP, KL, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO12 (KLP, 1, KL, 0, 1, IA, IB, C3) 
            CALL TWO12 (KLP, 0, KL, 1, 1, IA, IB, C4) 
            G2 = C3 + TWO*C4 
            G2 = A1*G2 
         ELSE 
            G2 = ZERO 
         ENDIF 
!
         CALL SOOB (1, L2, L2, KLP, KL, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO12 (KL, 0, KLP, 1, 1, IA, IB, C3) 
            CALL TWO12 (KL, 1, KLP, 0, 1, IA, IB, C4) 
            GG2 = C3 + TWO*C4 
            GG2 = A1*GG2 
         ELSE 
            GG2 = ZERO 
         ENDIF 
      ELSE 
         G2 = ZERO 
         GG2 = ZERO 
      ENDIF 
      CALL TWO12 (KL, 1, KL, 0, 1, IA, IB, C55) 
      CALL TWO12 (KL, 0, KL, 1, 1, IA, IB, C66) 
!      CALL TWO12(KL,0,KL,1,1,IA,IB,CC55)
!      CC55=-CC55
!      CALL TWO12(KL,1,KL,0,1,IA,IB,CC66)
!      CC66=-CC66
      CC55 = -C66 
      CC66 = -C55 
!
!     ( k k )  tenzorine struktura    V  integralas
!
      CALL SOOB (2, L1, L1, KL, KL, A1) 
      G4 = C55 + TWO*C66 
      G4 = AP1*A1*G4 
      IF (DABS(G4) > EPS) CALL SAVENON (7, G4, KL - 1, LA, LB, LA, LB, JA, JB, &
         0) 
      CALL SOOB (2, L2, L2, KL, KL, A1) 
      GG4 = CC55 + TWO*CC66 
      GG4 = AP1*A1*GG4 
      IF (DABS(GG4) > EPS) CALL SAVENON (7, GG4, KL - 1, LB, LA, LB, LA, JA, JB&
         , 0) 
!
!     ( k  k+1)  ( k k )
!
      CALL SOOB (1, L1, L1, KL, KL, A1) 
      CALL SOOB (1, L2, L2, KL, KL, A2) 
      G3 = C55 + TWO*C66 
      G3 = A1*G3 
      GG3 = CC55 + TWO*CC66 
      GG3 = A2*GG3 
      IF (KL - 1 >= 0) THEN 
         A1 = AP1*(G1 + DBLE(KL + 1)*G3) 
         IF (DABS(A1) > EPS) CALL SAVENON (6, A1, KL - 2, LB, LA, LB, LA, JA, &
            JB, 0) 
         A1 = AP1*(GG1 + DBLE(KL + 1)*GG3) 
         IF (DABS(A1) > EPS) CALL SAVENON (6, A1, KL - 2, LA, LB, LA, LB, JA, &
            JB, 0) 
      ENDIF 
!
!     ( k  k-1)  ( k k )
!
      IF (KLP <= IG) THEN 
         A2 = AP1*(G2 - DBLE(KL)*G3) 
         IF (DABS(A2) > EPS) CALL SAVENON (6, A2, KL, LA, LB, LA, LB, JA, JB, 0&
            ) 
         A2 = AP1*(GG2 - DBLE(KL)*GG3) 
         IF (DABS(A2) > EPS) CALL SAVENON (6, A2, KL, LB, LA, LB, LA, JA, JB, 0&
            ) 
      ENDIF 
      RETURN  
      END SUBROUTINE SOO1212P 
