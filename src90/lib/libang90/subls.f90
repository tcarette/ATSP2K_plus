!
!     ---------------------------------------------------------------
!     S U B L S
!     ---------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                          December 1993   *
!
      SUBROUTINE SUBLS(J1, J2, LL, S) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:56:58  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rumt_I 
      USE c0t5s_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: J1 
      INTEGER  :: J2 
      INTEGER  :: LL 
      REAL(DOUBLE) , INTENT(OUT) :: S 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LQ, LS, L, LQS, LSS, L1S, N, IN 
      REAL(DOUBLE) :: QQ, Q, QS, QM, QMS, A4, A1 
!-----------------------------------------------
      S = ZERO 
      CALL RUMT (J1, LL, LQ, LS, L) 
      CALL RUMT (J2, LL, LQS, LSS, L1S) 
      QQ = HALF 
      N = 2 
      Q = HALF*DBLE(LQ) 
      QS = HALF*DBLE(LQS) 
      QM = -HALF*DBLE(2*LL + 1 - N) 
      QMS = -HALF*DBLE(2*LL + 1 - (N - 1)) 
      CALL C0T5S (QS, QMS, QQ, Q, QM, A4) 
      IF (DABS(A4) < EPS) RETURN  
      A1 = DBLE(N*(LQ + 1)*(L + 1)*(LS + 1)) 
      S = DSQRT(A1)/A4 
      IN = (-N) - 1 
      IF ((IN/2)*2 /= IN) S = -S 
      RETURN  
      END SUBROUTINE SUBLS 
