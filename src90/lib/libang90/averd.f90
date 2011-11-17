!
!     ------------------------------------------------------------------
!                       A V E R D
!     ------------------------------------------------------------------
!
!     Add the deviations to the average energy for a partially filled
!       d- shell
!
      SUBROUTINE AVERD(J, K, N, A) 
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
      INTEGER , DIMENSION(16) :: F2DD2, F4DD2, F2DD3, F4DD3, F2DD4, F4DD4, &
         F2DD5, F4DD5 
      INTEGER , DIMENSION(5) :: IAVS, IAVV 
      INTEGER :: IDV, IVV, KK, JJ 
!-----------------------------------------------
      DATA IAVS/ 1, 0, -2, 0, -2/  
      DATA IAVV/ 1, 0, 63, 0, 63/  
!             ... d2 coefficients
      DATA F2DD2/ 140, 0, 77, 3*0, -13, 0, -58, 3*0, 50, 3*0/  
      DATA F4DD2/ 140, 0, -70, 3*0, 50, 0, 5, 3*0, 15, 3*0/  
!             ... d3 coefficients
      DATA F2DD3/ 2*0, 42, -12, 105, 69, 2*0, -93, 123, 0, -57, 2*0, -12, 0/  
      DATA F4DD3/ 2*0, -105, 30, 105, -15, 2*0, -30, -45, 0, 55, 2*0, 30, 0/  
!             ... d4 coefficients
      DATA F2DD4/ 210, 138, 21, 57, -105, 39, 219, 111, 66, 12, 84, -24, 30, 48&
         , -69, -51/  
      DATA F4DD4/ 210, -30, 70, -55, -105, -45, 30, -15, 45, -30, 0, -10, 135, &
         20, 15, 75/  
!             ... d5 coefficients
      DATA F2DD5/ -175, 113, -112, 320, 140, 104, -22, 86, 23, -85, 59, 167, &
         -85, 23, -58, -76/  
      DATA F4DD5/ -175, -55, 35, -100, 140, 20, -85, -40, -40, 125, -25, -15, &
         -50, -5, 110, 50/  
      A = ZERO 
      IF (MOD(K,2) /= 0) RETURN  
      IF (K > 4) RETURN  
      IDV = 0 
      IVV = 441 
      KK = K + 1 
      IF (N==1 .OR. N==9) THEN 
         IF (J /= 13) RETURN  
      ELSE IF (N==2 .OR. N==8) THEN 
         JJ = J - 24 
         IF (K == 2) THEN 
            IDV = F2DD2(JJ) 
         ELSE IF (K == 4) THEN 
            IDV = F4DD2(JJ) 
         ENDIF 
      ELSE IF (N==3 .OR. N==7) THEN 
         JJ = J - 8 
         IF (K == 2) THEN 
            IDV = F2DD3(JJ) 
         ELSE IF (K == 4) THEN 
            IDV = F4DD3(JJ) 
         ENDIF 
      ELSE IF (N==4 .OR. N==6) THEN 
         JJ = J - 24 
         IF (K == 2) THEN 
            IDV = F2DD4(JJ) 
         ELSE IF (K == 4) THEN 
            IDV = F4DD4(JJ) 
         ENDIF 
      ELSE IF (N == 5) THEN 
         JJ = J - 8 
         IF (K == 2) THEN 
            IDV = F2DD5(JJ) 
         ELSE IF (K == 4) THEN 
            IDV = F4DD5(JJ) 
         ENDIF 
      ELSE IF (N == 10) THEN 
         IF (J /= 25) RETURN  
      ENDIF 
      A = DBLE(IAVS(KK)*N*(N-1))*HALF/DBLE(IAVV(KK)) + DBLE(IDV)/DBLE(IVV) 
      RETURN  
      END SUBROUTINE AVERD 
