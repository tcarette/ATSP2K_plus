!
!     ------------------------------------------------------------------
!               L I N E Q N
!     ------------------------------------------------------------------
!
!        This routine is a modification of the one in "Computer Methods
!     for Mathematical Computation" by Forsythe, Malcolm, and Moler
!     (Prentice Hall, 1975)
!
      SUBROUTINE LINEQN(NDIM, N, A, B) 
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
      INTEGER , INTENT(IN) :: NDIM 
      INTEGER , INTENT(IN) :: N 
      REAL(DOUBLE) , INTENT(INOUT) :: A(NDIM,*) 
      REAL(DOUBLE) , INTENT(INOUT) :: B(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NM1, K, KP1, M, I, J 
      REAL(DOUBLE) :: T 
!-----------------------------------------------
!
!     SOLVE A SYSTEM OF LINEAR EQUATIONS A*X = B
!
!     INPUT..
!
!        NDIM = DECLARED ROW DIMENSION OF THE ARRAY CONTAINING  A.
!        N = ORDER OF THE MATRIX
!        A =  COEFFICIENT MATRIX
!        B = CONSTANT VECTOR
!
!     OUTPUT..
!
!        B = SOLUTION VECTOR
!
      NM1 = N - 1 
!
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
      DO K = 1, NM1 
         KP1 = K + 1 
!
!        FIND PIVOT
!
         M = K 
         DO I = KP1, N 
            IF (DABS(A(I,K)) <= DABS(A(M,K))) CYCLE  
            M = I 
         END DO 
         T = A(M,K) 
         A(M,K) = A(K,K) 
         A(K,K) = T 
!
!        SKIP STEP IF PIVOT IS ZERO
!
         IF (T == 0.0D0) CYCLE  
!
!        COMPUTE MULTIPLIERS
!
         A(KP1:N,K) = -A(KP1:N,K)/T 
!
!        INTERCHANGE AND ELIMINATE BY COLUMNS
!
         DO J = KP1, N 
            T = A(M,J) 
            A(M,J) = A(K,J) 
            A(K,J) = T 
            IF (T == 0.0D0) CYCLE  
            A(KP1:N,J) = A(KP1:N,J) + A(KP1:N,K)*T 
         END DO 
         T = B(M) 
         B(M) = B(K) 
         B(K) = T 
         IF (T == 0.0D0) CYCLE  
         B(KP1:N) = B(KP1:N) + A(KP1:N,K)*T 
      END DO 
!
!     BACK SUBSTITUTION
!
      DO K = N, 2, -1 
         IF (A(K,K) == 0.D0) THEN 
            B(K) = 0.D0 
         ELSE 
            B(K) = B(K)/A(K,K) 
         ENDIF 
         T = -B(K) 
         B(:K-1) = B(:K-1) + A(:K-1,K)*T 
      END DO 
      B(1) = B(1)/A(1,1) 
      RETURN  
      END SUBROUTINE LINEQN 
