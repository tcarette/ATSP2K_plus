!
!     -------------------------------------------------------------
!      S O O 1 1 2 2 P
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :      N'1 = N1 +- 1        *
!                                             N'2 = N2 +- 1        *
!                                             N'3 = N3 -+ 2        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE SOO1122P(IG, KL, IA, IB, IC, ID, IIA, IIB, IIC, IID, IREZ, XXX&
         ) 
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
      INTEGER  :: IC 
      INTEGER  :: ID 
      INTEGER , INTENT(IN) :: IIA 
      INTEGER , INTENT(IN) :: IIB 
      INTEGER , INTENT(IN) :: IIC 
      INTEGER , INTENT(IN) :: IID 
      INTEGER  :: IREZ 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L1, L2, L3, L4, LA, LB, LC, LD, KLM, KLP 
      REAL(DOUBLE) :: AP1, A1, AA1, G1, GG1, C1, CC1, C2, CC2, G2, GG2, C3, CC3&
         , C4, CC4, C55, CC55, C66, CC66, C7, C8, G4, CC7, CC8, GG4, A2, C5, &
         CC5, C6, CC6, G3, GG3 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
      L1 = LJ(IIA) 
      L2 = LJ(IIB) 
      L3 = LJ(IIC) 
      L4 = LJ(IID) 
      CALL SOOA (L1, L2, L3, L4, KL, AP1) 
      IF (DABS(AP1) < EPS) RETURN  
      LA = IJFUL(IIA) 
      LB = IJFUL(IIB) 
      LC = IJFUL(IIC) 
      LD = IJFUL(IID) 
      KLM = KL - 1 
      IF (KLM >= 0) THEN 
         CALL SOOB (1, L1, L3, KLM, KL, A1) 
         CALL SOOB (1, L2, L4, KLM, KL, AA1) 
         IF (DABS(A1)<EPS .AND. DABS(AA1)<EPS) THEN 
            G1 = ZERO 
            GG1 = ZERO 
         ELSE 
            CALL XXX (KLM, 1, KL, 0, 1, IA, IB, IC, ID, IREZ, C1, CC1) 
            CALL XXX (KLM, 0, KL, 1, 1, IA, IB, IC, ID, IREZ, C2, CC2) 
            C1 = C1*A1 
            C2 = TWO*C2*A1 
            G1 = C1 + C2 
            CC1 = CC1*AA1 
            CC2 = TWO*CC2*AA1 
            GG1 = CC1 + CC2 
         ENDIF 
      ENDIF 
      KLP = KL + 1 
      IF (KLP <= IG) THEN 
         CALL SOOB (1, L1, L3, KLP, KL, A1) 
         CALL SOOB (1, L2, L4, KLP, KL, AA1) 
         IF (DABS(A1)<EPS .AND. DABS(AA1)<EPS) THEN 
            G2 = ZERO 
            GG2 = ZERO 
         ELSE 
            CALL XXX (KLP, 1, KL, 0, 1, IA, IB, IC, ID, IREZ, C3, CC3) 
            CALL XXX (KLP, 0, KL, 1, 1, IA, IB, IC, ID, IREZ, C4, CC4) 
            C3 = C3*A1 
            C4 = TWO*C4*A1 
            G2 = C3 + C4 
            CC3 = CC3*AA1 
            CC4 = TWO*CC4*AA1 
            GG2 = CC3 + CC4 
         ENDIF 
      ENDIF 
      IF (KL > 0) THEN 
         CALL XXX (KL, 1, KL, 0, 1, IA, IB, IC, ID, IREZ, C55, CC55) 
         CALL XXX (KL, 0, KL, 1, 1, IA, IB, IC, ID, IREZ, C66, CC66) 
      ELSE 
         C55 = ZERO 
         CC55 = ZERO 
         C66 = ZERO 
         CC66 = ZERO 
      ENDIF 
!
!    ( k k )  tensorine struktura        V   integralas
!
      CALL SOOB (2, L1, L3, KL, KL, A1) 
      C7 = C55*A1 
      C8 = TWO*C66*A1 
      G4 = AP1*(C7 + C8) 
      IF (DABS(G4) > EPS) CALL SAVENON (7, G4, KL - 1, LA, LB, LC, LD, JA, JB, &
         0) 
      CALL SOOB (2, L2, L4, KL, KL, A1) 
      CC7 = CC55*A1 
      CC8 = TWO*CC66*A1 
      GG4 = AP1*(CC7 + CC8) 
      IF (DABS(GG4) > EPS) CALL SAVENON (7, GG4, KL - 1, LB, LA, LD, LC, JA, JB&
         , 0) 
!
!    ( k  k+1)  + ( k  k )
!
      CALL SOOB (1, L1, L3, KL, KL, A1) 
      CALL SOOB (1, L2, L4, KL, KL, A2) 
      C5 = C55*A1 
      CC5 = CC55*A2 
      C6 = TWO*C66*A1 
      CC6 = TWO*CC66*A2 
      G3 = C5 + C6 
      GG3 = CC5 + CC6 
      IF (KL - 1 >= 0) THEN 
         A1 = AP1*(G1 + DBLE(KL + 1)*G3) 
         IF (DABS(A1) > EPS) CALL SAVENON (6, A1, KL - 2, LB, LA, LD, LC, JA, &
            JB, 0) 
         A1 = AP1*(GG1 + DBLE(KL + 1)*GG3) 
         IF (DABS(A1) > EPS) CALL SAVENON (6, A1, KL - 2, LA, LB, LC, LD, JA, &
            JB, 0) 
      ENDIF 
!
!    ( k  k-1)  + ( k  k )
!
      IF (KLP <= IG) THEN 
         A2 = AP1*(G2 - DBLE(KL)*G3) 
         IF (DABS(A2) > EPS) CALL SAVENON (6, A2, KL, LA, LB, LC, LD, JA, JB, 0&
            ) 
         A2 = AP1*(GG2 - DBLE(KL)*GG3) 
         IF (DABS(A2) > EPS) CALL SAVENON (6, A2, KL, LB, LA, LD, LC, JA, JB, 0&
            ) 
      ENDIF 
      RETURN  
      END SUBROUTINE SOO1122P 
