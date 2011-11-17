!
!     ------------------------------------------------------------------
!                       A V E R P
!     ------------------------------------------------------------------
!
!     Add the deviations to the average energy for a partially filled
!       p- shell
!
      SUBROUTINE AVERP(J, K, N, A) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:33:24  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: J 
      INTEGER , INTENT(IN) :: K 
      INTEGER , INTENT(IN) :: N 
      REAL(DOUBLE) , INTENT(OUT) :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(3) :: F2PP2, F2PP3, IAVS, IAVV 
      INTEGER :: IDV, IVV, KK, JJ 
!-----------------------------------------------
      DATA IAVS/ 1, 0, -2/  
      DATA IAVV/ 1, 0, 25/  
!             ... p2 coefficients
      DATA F2PP2/ 12, -3, 3/  
!             ... p3 coefficients
      DATA F2PP3/ -9, 6, 0/  
      A = ZERO 
      IF (MOD(K,2) /= 0) RETURN  
      IF (K > 2) RETURN  
      IDV = 0 
      IVV = 25 
      KK = K + 1 
      IF (N==1 .OR. N==5) THEN 
         IF (J /= 4) RETURN  
      ELSE IF (N==2 .OR. N==4) THEN 
         IF (K == 2) THEN 
            JJ = J - 5 
            IDV = F2PP2(JJ) 
         ENDIF 
      ELSE IF (N == 3) THEN 
         IF (K == 2) THEN 
            JJ = J - 2 
            IDV = F2PP3(JJ) 
         ENDIF 
      ELSE IF (N == 6) THEN 
         IF (J /= 6) RETURN  
      ENDIF 
      A = DBLE(IAVS(KK)*N*(N-1))*HALF/DBLE(IAVV(KK)) + DBLE(IDV)/DBLE(IVV) 
      RETURN  
      END SUBROUTINE AVERP 
