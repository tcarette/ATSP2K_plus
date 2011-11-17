!
!     -------------------------------------------------------------
!      S O O 1 1 1 1 P
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
      SUBROUTINE SOO1111P(IG, KL, IA) 
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
      USE soo1_I 
      USE two11_I 
      USE savenon_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IG 
      INTEGER  :: KL 
      INTEGER  :: IA 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KLM, LA, KLP 
      REAL(DOUBLE) :: A1, C1, C2, G1, C3, C4, G2, C5, C6, G3, G4 
!-----------------------------------------------
      KLM = KL - 1 
      LA = IJFUL(IA) 
      IF (KLM >= 0) THEN 
         CALL SOO1 (LJ(IA), KLM, 1, KL, 0, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO11 (KLM, 1, KL, 0, 1, IA, C1) 
            C1 = C1*A1 
         ELSE 
            C1 = ZERO 
         ENDIF 
         CALL SOO1 (LJ(IA), KLM, 0, KL, 1, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO11 (KLM, 0, KL, 1, 1, IA, C2) 
            C2 = C2*A1 
         ELSE 
            C2 = ZERO 
         ENDIF 
         G1 = C1 + C2 
      ELSE 
         G1 = ZERO 
      ENDIF 
      KLP = KL + 1 
      IF (KLP <= IG) THEN 
         CALL SOO1 (LJ(IA), KLP, 1, KL, 0, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO11 (KLP, 1, KL, 0, 1, IA, C3) 
            C3 = C3*A1 
         ELSE 
            C3 = ZERO 
         ENDIF 
         CALL SOO1 (LJ(IA), KLP, 0, KL, 1, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO11 (KLP, 0, KL, 1, 1, IA, C4) 
            C4 = C4*A1 
         ELSE 
            C4 = ZERO 
         ENDIF 
         G2 = C3 + C4 
      ELSE 
         G2 = ZERO 
      ENDIF 
      IF (KL > 0) THEN 
         CALL SOO1 (LJ(IA), KL, 1, KL, 0, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO11 (KL, 1, KL, 0, 1, IA, C5) 
            C5 = C5*A1 
         ELSE 
            C5 = ZERO 
         ENDIF 
         CALL SOO1 (LJ(IA), KL, 0, KL, 1, A1) 
         IF (DABS(A1) > EPS) THEN 
            CALL TWO11 (KL, 0, KL, 1, 1, IA, C6) 
            C6 = C6*A1 
         ELSE 
            C6 = ZERO 
         ENDIF 
      ELSE 
         C5 = ZERO 
         C6 = ZERO 
      ENDIF 
!
!     ( k k )   tenzorine struktura    V integralas
!
      G3 = C5 + C6 
      G4 = TWO*G3 
      IF (DABS(G4) > EPS) CALL SAVENON (7, G4, KL - 1, LA, LA, LA, LA, JA, JB, &
         0) 
!
!     ( k k+1 )  +  ( k k )
!
      IF (KL - 1 >= 0) THEN 
         A1 = G1 + DBLE(KL + 1)*G3 
         IF (DABS(A1) > EPS) CALL SAVENON (6, A1, KL - 2, LA, LA, LA, LA, JA, &
            JB, 0) 
      ENDIF 
!
!     ( k k-1 )  +  ( k k )
!
      IF (KLP <= IG) THEN 
         A1 = G2 - DBLE(KL)*G3 
         IF (DABS(A1) > EPS) CALL SAVENON (6, A1, KL, LA, LA, LA, LA, JA, JB, 0&
            ) 
      ENDIF 
      RETURN  
      END SUBROUTINE SOO1111P 
