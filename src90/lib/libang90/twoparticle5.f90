!
!     -------------------------------------------------------------
!      T W O P A R T I C L E 5
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :    N'1 = N1 (+-) 1        *
!                                           N'2 = N2 (+-) 1        *
!                                           N'3 = N3 (+-) 1        *
!                                           N'4 = N4 (+-) 1        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE TWOPARTICLE5(IA, IB, IC, ID) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
      USE TRK_C 
      USE TRK2_C 
      USE PERMAT_C 
      USE KAMPAS_C 
      USE CASEOP_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:23:13  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE two51_I 
      USE two52_I 
      USE two53_I 
      USE two54_I 
      USE two55_I 
      USE two56_I 
      USE rlsp00_I 
      USE a1a2a3a4lsp_I 
      USE itrexg2_I 
      USE itrexg_I 
      USE ssc_I 
      USE sooc_I 
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
      INTEGER :: KGAL, I1, I2, IAT, IP, IKK, IG, I, KL 
      REAL(DOUBLE) :: W 
!-----------------------------------------------
      IF (IHSH <= 3) RETURN  
!
      IF (IOCASE == 1) THEN 
! Spin - spin      operator
         KGAL = 2 
      ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
         KGAL = 1 
      ENDIF 
      IRS = 0 
      IW1(1,:) = 0 
      IW1(2,:) = 0 
      IW2(1,:) = 0 
      IW2(2,:) = 0 
      IRL = 0 
      IWAA(1,1,:,:) = 0 
      IWAA(2,1,:,:) = 0 
!   i)
      IF (IB < IC) THEN 
         CALL RLSP00 (1, IA, IB, IC, ID, 0, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP00 (2, IA, IB, IC, ID, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP00 (3, IA, IB, IC, ID, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL A1A2A3A4LSP (IA, IB, IC, ID, HALF, HALF, (-HALF), (-HALF), W) 
         IF (DABS(W) < EPS) RETURN  
         IF (IOCASE == 1) THEN 
! Spin-Spin operator
            IP = ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK) + 1 
         ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
            IP = ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK) + 1 
         ENDIF 
         IF (IKK > 0) THEN 
            IG = IP + IKK - 1 
            DO I = IP, IG 
               KL = I - 1 
!
               IF (IOCASE == 1) THEN 
! Spin-spin operator
                  CALL SSC (IG, KL, IA, IB, IC, ID, IA, IB, IC, ID, 1, TWO51) 
               ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
                  CALL SOOC (IG, KL, IA, IB, IC, ID, IA, IB, IC, ID, 1, TWO51) 
               ENDIF 
!
            END DO 
         ENDIF 
         IF (IOCASE == 1) THEN 
! Spin-Spin operator
            IP = ITREXG2(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK) + 1 
         ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
            IP = ITREXG(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK) + 1 
         ENDIF 
         IF (IKK <= 0) RETURN  
         IG = IP + IKK - 1 
         DO I = IP, IG 
            KL = I - 1 
!
            IF (IOCASE == 1) THEN 
! Spin-spin operator
               CALL SSC (IG, KL, IA, IB, IC, ID, IA, IB, ID, IC, 1, TWO52) 
            ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
               CALL SOOC (IG, KL, IA, IB, IC, ID, IA, IB, ID, IC, 1, TWO52) 
            ENDIF 
!
         END DO 
!
      ELSE IF (IA>ID .AND. IB>ID) THEN 
         CALL RLSP00 (1, IC, ID, IA, IB, 0, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP00 (2, IC, ID, IA, IB, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP00 (3, IC, ID, IA, IB, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL A1A2A3A4LSP (IC, ID, IA, IB, (-HALF), (-HALF), HALF, HALF, W) 
         IF (DABS(W) < EPS) RETURN  
         IF (IOCASE == 1) THEN 
! Spin-Spin operator
            IP = ITREXG2(LJ(IC),LJ(IA),LJ(ID),LJ(IB),IKK) + 1 
         ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
            IP = ITREXG(LJ(IC),LJ(IA),LJ(ID),LJ(IB),IKK) + 1 
         ENDIF 
         IF (IKK > 0) THEN 
            IG = IP + IKK - 1 
            DO I = IP, IG 
               KL = I - 1 
!
               IF (IOCASE == 1) THEN 
! Spin-spin operator
                  CALL SSC (IG, KL, IC, ID, IA, IB, IA, IB, IC, ID, 2, TWO51) 
               ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
                  CALL SOOC (IG, KL, IC, ID, IA, IB, IA, IB, IC, ID, 2, TWO51) 
               ENDIF 
!
            END DO 
         ENDIF 
         IF (IOCASE == 1) THEN 
! Spin-Spin operator
            IP = ITREXG2(LJ(IC),LJ(IB),LJ(ID),LJ(IA),IKK) + 1 
         ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
            IP = ITREXG(LJ(IC),LJ(IB),LJ(ID),LJ(IA),IKK) + 1 
         ENDIF 
         IF (IKK <= 0) RETURN  
         IG = IP + IKK - 1 
         DO I = IP, IG 
            KL = I - 1 
!
            IF (IOCASE == 1) THEN 
! Spin-spin operator
               CALL SSC (IG, KL, IC, ID, IA, IB, IB, IA, IC, ID, 2, TWO52) 
            ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
               CALL SOOC (IG, KL, IC, ID, IA, IB, IB, IA, IC, ID, 2, TWO52) 
            ENDIF 
!
         END DO 
!  ii)
      ELSE IF (IB>IC .AND. IB<ID .AND. IA<IC) THEN 
         CALL RLSP00 (1, IA, IC, IB, ID, 0, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP00 (2, IA, IC, IB, ID, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP00 (3, IA, IC, IB, ID, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL A1A2A3A4LSP (IA, IC, IB, ID, HALF, (-HALF), HALF, (-HALF), W) 
         IF (DABS(W) < EPS) RETURN  
         IF (IOCASE == 1) THEN 
! Spin-Spin operator
            IP = ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK) + 1 
         ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
            IP = ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK) + 1 
         ENDIF 
         IF (IKK > 0) THEN 
            IG = IP + IKK - 1 
            DO I = IP, IG 
               KL = I - 1 
!
               IF (IOCASE == 1) THEN 
! Spin-spin operator
                  CALL SSC (IG, KL, IA, IC, IB, ID, IA, IB, IC, ID, 1, TWO53) 
               ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
                  CALL SOOC (IG, KL, IA, IC, IB, ID, IA, IB, IC, ID, 1, TWO53) 
               ENDIF 
!
            END DO 
         ENDIF 
         IF (IOCASE == 1) THEN 
! Spin-Spin operator
            IP = ITREXG2(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK) + 1 
         ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
            IP = ITREXG(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK) + 1 
         ENDIF 
         IF (IKK <= 0) RETURN  
         IG = IP + IKK - 1 
         DO I = IP, IG 
            KL = I - 1 
!
            IF (IOCASE == 1) THEN 
! Spin-spin operator
               CALL SSC (IG, KL, IA, IC, IB, ID, IA, IB, ID, IC, 1, TWO54) 
            ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
               CALL SOOC (IG, KL, IA, IC, IB, ID, IA, IB, ID, IC, 1, TWO54) 
            ENDIF 
!
         END DO 
!
      ELSE IF (IB>IC .AND. IB>ID .AND. IA>IC) THEN 
         CALL RLSP00 (1, IC, IA, ID, IB, 0, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP00 (2, IC, IA, ID, IB, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP00 (3, IC, IA, ID, IB, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL A1A2A3A4LSP (IC, IA, ID, IB, (-HALF), HALF, (-HALF), HALF, W) 
         IF (DABS(W) < EPS) RETURN  
         IF (IOCASE == 1) THEN 
! Spin-Spin operator
            IP = ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK) + 1 
         ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
            IP = ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK) + 1 
         ENDIF 
         IF (IKK > 0) THEN 
            IG = IP + IKK - 1 
            DO I = IP, IG 
               KL = I - 1 
!
               IF (IOCASE == 1) THEN 
! Spin-spin operator
                  CALL SSC (IG, KL, IC, IA, ID, IB, IA, IB, IC, ID, 2, TWO53) 
               ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
                  CALL SOOC (IG, KL, IC, IA, ID, IB, IA, IB, IC, ID, 2, TWO53) 
               ENDIF 
!
            END DO 
         ENDIF 
         IF (IOCASE == 1) THEN 
! Spin-Spin operator
            IP = ITREXG2(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK) + 1 
         ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
            IP = ITREXG(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK) + 1 
         ENDIF 
         IF (IKK <= 0) RETURN  
         IG = IP + IKK - 1 
         DO I = IP, IG 
            KL = I - 1 
!
            IF (IOCASE == 1) THEN 
! Spin-spin operator
               CALL SSC (IG, KL, IC, IA, ID, IB, IB, IA, IC, ID, 2, TWO54) 
            ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
               CALL SOOC (IG, KL, IC, IA, ID, IB, IB, IA, IC, ID, 2, TWO54) 
            ENDIF 
!
         END DO 
! iii)
      ELSE IF (IB>IC .AND. IB>ID .AND. IA<IC) THEN 
         CALL RLSP00 (1, IA, IC, ID, IB, 0, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP00 (2, IA, IC, ID, IB, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP00 (3, IA, IC, ID, IB, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL A1A2A3A4LSP (IA, IC, ID, IB, HALF, (-HALF), (-HALF), HALF, W) 
         IF (DABS(W) < EPS) RETURN  
         IF (IOCASE == 1) THEN 
! Spin-Spin operator
            IP = ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK) + 1 
         ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
            IP = ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK) + 1 
         ENDIF 
         IF (IKK > 0) THEN 
            IG = IP + IKK - 1 
            DO I = IP, IG 
               KL = I - 1 
!
               IF (IOCASE == 1) THEN 
! Spin-spin operator
                  CALL SSC (IG, KL, IA, IC, ID, IB, IA, IB, IC, ID, 1, TWO55) 
               ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
                  CALL SOOC (IG, KL, IA, IC, ID, IB, IA, IB, IC, ID, 1, TWO55) 
               ENDIF 
!
            END DO 
         ENDIF 
         IF (IOCASE == 1) THEN 
! Spin-Spin operator
            IP = ITREXG2(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK) + 1 
         ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
            IP = ITREXG(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK) + 1 
         ENDIF 
         IF (IKK <= 0) RETURN  
         IG = IP + IKK - 1 
         DO I = IP, IG 
            KL = I - 1 
!
            IF (IOCASE == 1) THEN 
! Spin-spin operator
               CALL SSC (IG, KL, IA, IC, ID, IB, IA, IB, ID, IC, 1, TWO56) 
            ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
               CALL SOOC (IG, KL, IA, IC, ID, IB, IA, IB, ID, IC, 1, TWO56) 
            ENDIF 
!
         END DO 
!
      ELSE IF (IB>IC .AND. IB<ID .AND. IA>IC) THEN 
         CALL RLSP00 (1, IC, IA, IB, ID, 0, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP00 (2, IC, IA, IB, ID, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL RLSP00 (3, IC, IA, IB, ID, 2*KGAL, IAT) 
         IF (IAT == 0) RETURN  
         CALL A1A2A3A4LSP (IC, IA, IB, ID, (-HALF), HALF, HALF, (-HALF), W) 
         IF (DABS(W) < EPS) RETURN  
         IF (IOCASE == 1) THEN 
! Spin-Spin operator
            IP = ITREXG2(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK) + 1 
         ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
            IP = ITREXG(LJ(IA),LJ(IC),LJ(IB),LJ(ID),IKK) + 1 
         ENDIF 
         IF (IKK > 0) THEN 
            IG = IP + IKK - 1 
            DO I = IP, IG 
               KL = I - 1 
!
               IF (IOCASE == 1) THEN 
! Spin-spin operator
                  CALL SSC (IG, KL, IC, IA, IB, ID, IA, IB, IC, ID, 2, TWO55) 
               ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
                  CALL SOOC (IG, KL, IC, IA, IB, ID, IA, IB, IC, ID, 2, TWO55) 
               ENDIF 
!
            END DO 
         ENDIF 
         IF (IOCASE == 1) THEN 
! Spin-Spin operator
            IP = ITREXG2(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK) + 1 
         ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
            IP = ITREXG(LJ(IA),LJ(ID),LJ(IB),LJ(IC),IKK) + 1 
         ENDIF 
         IF (IKK <= 0) RETURN  
         IG = IP + IKK - 1 
         DO I = IP, IG 
            KL = I - 1 
!
            IF (IOCASE == 1) THEN 
! Spin-spin operator
               CALL SSC (IG, KL, IC, IA, IB, ID, IB, IA, IC, ID, 2, TWO56) 
            ELSE IF (IOCASE == 2) THEN 
! Spin-other-orbit operator
               CALL SOOC (IG, KL, IC, IA, IB, ID, IB, IA, IC, ID, 2, TWO56) 
            ENDIF 
!
         END DO 
!
      ELSE 
         WRITE (6, '(A)') ' ERROR IN SUBROUTINE  TWOPARTICLE5  ' 
         STOP  
      ENDIF 
      RETURN  
      END SUBROUTINE TWOPARTICLE5 
