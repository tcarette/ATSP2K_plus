!
!     ------------------------------------------------------------------
!        S Y M M E Q   A N D  S Y M M S L
!     ------------------------------------------------------------------
!
!        This routine is a modification of the one in "Computer Methods
!     for Mathematical Computation" by Forsythe, Malcolm, and Moler
!     (Prentice Hall, 1975) to solve a singular system of equations.
!
      SUBROUTINE SYMMEQ(NDIM, N, A, X) 
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
      REAL(DOUBLE) , INTENT(INOUT) :: X(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NM1, K, KP1, I, J 
      REAL(DOUBLE) :: T 
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
!
!
!     SOLVE A SYSTEM OF LINEAR EQUATIONS A*X = 0
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
!        X = SOLUTION VECTOR with X(1) unchanged.
!
      NM1 = N - 1 
!
!
!      L U FACTORIZATION WITHOUT PIVOTING
!
      DO K = N, 3, -1 
         KP1 = K - 1 
         T = A(K,K) 
         IF (T == 0.0D0) CYCLE  
!
!        COMPUTE MULTIPLIERS
!
         A(KP1:2:(-1),K) = -A(KP1:2:(-1),K)/T 
!
!        INTERCHANGE AND ELIMINATE BY COLUMNS
!
         DO J = KP1, 2, -1 
            T = A(K,J) 
            IF (T == 0.0D0) CYCLE  
            A(KP1:2:(-1),J) = A(KP1:2:(-1),J) + A(KP1:2:(-1),K)*T 
         END DO 
      END DO 
!
!     At this point it is assumed that the LU factorization
!     has already been performed.
!
      ENTRY SYMMSL (NDIM, N, A, X) 
!
      X(2:N) = -A(2:N,1) 
      DO K = N, 3, -1 
         T = X(K) 
         IF (T == 0.0D0) CYCLE  
         X(K-1:2:(-1)) = X(K-1:2:(-1)) + A(K-1:2:(-1),K)*T 
      END DO 
!
!     BACK SUBSTITUTION
!
      DO K = 2, N - 1 
         IF (A(K,K) == 0.D0) THEN 
            X(K) = 0.D0 
         ELSE 
            X(K) = X(K)/A(K,K) 
         ENDIF 
         T = -X(K) 
         X(K+1:N) = X(K+1:N) + A(K+1:N,K)*T 
      END DO 
      X(N) = X(N)/A(N,N) 
      RETURN  
      END SUBROUTINE SYMMEQ 
