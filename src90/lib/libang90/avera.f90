!
!     ------------------------------------------------------------------
!                       A V E R A
!     ------------------------------------------------------------------
!
!     Add the deviations to the average energy for a partially filled
!       s, p, d, f- shells
!
      SUBROUTINE AVERA(L, J, K, N, A) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:33:24  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE averf_I 
      USE averp_I 
      USE averd_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: L 
      INTEGER  :: J 
      INTEGER  :: K 
      INTEGER  :: N 
      REAL(DOUBLE)  :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      SELECT CASE (L)  
      CASE (3)  
         CALL AVERF (J, K, N, A) 
      CASE (0)  
         A = ZERO 
         IF (K /= 0) RETURN  
         IF (N == 1) THEN 
            IF (J /= 1) RETURN  
         ELSE IF (N == 2) THEN 
            IF (J /= 2) RETURN  
         ENDIF 
         A = DBLE(N*(N - 1))*HALF 
      CASE (1)  
         CALL AVERP (J, K, N, A) 
      CASE (2)  
         CALL AVERD (J, K, N, A) 
      END SELECT 
      RETURN  
      END SUBROUTINE AVERA 
