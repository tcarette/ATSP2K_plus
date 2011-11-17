!
!     -------------------------------------------------------------
!      R W L S
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                              (k1 k2 k3)          *
!     REDUCED MATRIX ELEMENT       (nl QLS::: W        :::nl QLS)  *
!                                                                  *
!     Written and extended by G. Gaigalas,                         *
!     Vilnius,  Lithuania                          December 1993   *
!     Bruxelles,  Belgium                          December 1995   *
!
      SUBROUTINE RWLS(K1, K2, K3, L, J1, J2, W) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE RIBOF_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:16:44  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rumt_I 
      USE ittk_I 
      USE rmewpls_I 
      USE rmewdls_I 
      USE sls_I 
      USE sixj_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: K1 
      INTEGER  :: K2 
      INTEGER  :: K3 
      INTEGER  :: L 
      INTEGER  :: J1 
      INTEGER  :: J2 
      REAL(DOUBLE)  :: W 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: L2QJ1, L2SJ1, L2LJ1, L2QJ2, L2SJ2, L2LJ2, KK1, KK2, KK3, LJ1, &
         IP, IG, I, L2QI, L2SI, L2LI, LAS 
      REAL(DOUBLE) :: QJ1, SJ1, D, SAKNIS, SAKNI2, S, S1, S2, SS 
!-----------------------------------------------
!GGf                                with    f-shell   ***  beginning
!GGf                                with    f-shell   ***        end
      CALL RUMT (J1, L, L2QJ1, L2SJ1, L2LJ1) 
      CALL RUMT (J2, L, L2QJ2, L2SJ2, L2LJ2) 
      W = ZERO 
      KK1 = K1*2 
      IF (ITTK(L2QJ1,L2QJ2,KK1) == 0) RETURN  
      KK2 = K2*2 
      IF (ITTK(L2LJ1,L2LJ2,KK2) == 0) RETURN  
      KK3 = K3*2 
      IF (ITTK(L2SJ1,L2SJ2,KK3) == 0) RETURN  
      QJ1 = DBLE(L2QJ1)/TWO 
      LJ1 = L2LJ1/2 
      SJ1 = DBLE(L2SJ1)/TWO 
      IF (K2 - 1 <= 0) THEN 
         IF (K2 - 1 /= 0) THEN 
            IF (K1 /= 1) THEN 
               IF (K3 /= 1) THEN 
                  D = -(TWO*DBLE(L) + ONE) 
                  GO TO 15 
               ENDIF 
               D = -TWO*SQRT(SJ1*(SJ1 + ONE)) 
               GO TO 15 
            ENDIF 
            D = -TWO*SQRT(QJ1*(QJ1 + ONE)) 
            GO TO 15 
         ENDIF 
         IF (K1/=0 .OR. K3/=0) GO TO 2 
         SAKNIS = DBLE(3*LJ1*(LJ1 + 1)) 
         SAKNI2 = DBLE(L*(L + 1)) 
         D = -SQRT(SAKNIS/SAKNI2) 
   15    CONTINUE 
         IF (J1 /= J2) RETURN  
         SAKNIS = DBLE((L2QJ1 + 1)*(L2LJ1 + 1)*(L2SJ1 + 1)) 
         SAKNI2 = DBLE(2*L + 1) 
         W = D*SQRT(SAKNIS/SAKNI2) 
         RETURN  
      ENDIF 
    2 CONTINUE 
      SELECT CASE (L)  
      CASE (1)  
         CALL RMEWPLS (J1, J2, K1, K2, K3, W) 
      CASE (2)  
         CALL RMEWDLS (J1, J2, K1, K2, K3, W) 
!GGf                                with    f-shell   ***  beginning
      CASE (3)  
         IP = IMPNF(J2) 
         IG = IMGNF(J2) 
         S = ZERO 
         DO I = IP, IG 
            CALL RUMT (I, L, L2QI, L2SI, L2LI) 
            CALL SLS (L, J1, L2QJ1, L2LJ1, L2SJ1, I, L2QI, L2LI, L2SI, S1) 
            IF (ABS(S1) > EPS) THEN 
               CALL SLS (L, I, L2QI, L2LI, L2SI, J2, L2QJ2, L2LJ2, L2SJ2, S2) 
               S1 = S1*S2 
               IF (ABS(S1) > EPS) THEN 
                  CALL SIXJ (2*L, 2*L, KK2, L2LJ2, L2LJ1, L2LI, 1, S2) 
                  S1 = S1*S2 
                  CALL SIXJ (1, 1, KK1, L2QJ2, L2QJ1, L2QI, 1, S2) 
                  S1 = S1*S2 
                  CALL SIXJ (1, 1, KK3, L2SJ2, L2SJ1, L2SI, 1, S2) 
                  S1 = S1*S2 
               ENDIF 
            ENDIF 
            S = S + S1 
         END DO 
         SS = DBLE((KK1 + 1)*(KK2 + 1)*(KK3 + 1)) 
         W = S*SQRT(SS) 
         LAS = L2QJ1 + L2LJ1 + L2SJ1 + L2QJ2 + L2SJ2 + L2LJ2 + KK1 + KK2 + KK3 
         IF (MOD(LAS,4) /= 0) W = -W 
!GGf                                with    f-shell   ***        end
      END SELECT 
      RETURN  
      END SUBROUTINE RWLS 
