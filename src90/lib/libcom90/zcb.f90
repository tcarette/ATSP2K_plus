!======================================================================
!         Z C B
!======================================================================
 
      REAL(KIND(0.0D0)) FUNCTION ZCB (K1, K2, K3) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: K1 
      INTEGER , INTENT(IN) :: K2 
      INTEGER , INTENT(IN) :: K3 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: M, N, N1, N2, N3, M1, M2, M3, I, K, IK, J 
      REAL(DOUBLE) :: A, B 
!-----------------------------------------------
!
!     CB =  3j(k1,0,k2,0,k3,0)**2
!
!
      ZCB = 0.0 
!
      M = K1 + K2 + K3 
      N = M/2 
      IF (N + N /= M) RETURN  
      M = M + 1 
!
      N1 = N - K1 
      N2 = N - K2 
      N3 = N - K3 
      IF (N1<0 .OR. N2<0 .OR. N3<0) RETURN  
      M1 = N1 + N1 
      M2 = N2 + N2 
      M3 = N3 + N3 
!
      A = 1.0 
!
      DO I = 2, M 
         K = -1                                  ! the extent of integer I in units of 1/2 
         IF (M1 >= I) K = K + 1 
         IF (M2 >= I) K = K + 1 
         IF (M3 >= I) K = K + 1 
         IF (N >= I) K = K + 2 
         IF (N1 >= I) K = K - 2 
         IF (N2 >= I) K = K - 2 
         IF (N3 >= I) K = K - 2 
!
         IK = IABS(K) 
         B = DBLE(I) 
         IF (K > 0) THEN 
            DO J = 1, IK 
               A = A*B 
            END DO 
         ENDIF 
         IF (K >= 0) CYCLE  
         DO J = 1, IK 
            A = A/B 
         END DO 
!
      END DO 
!
      ZCB = A 
      RETURN  
      END FUNCTION ZCB 
