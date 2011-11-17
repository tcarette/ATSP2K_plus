*
*     ------------------------------------------------------------------
*               L I N E Q N
*     ------------------------------------------------------------------
*
*        This routine is a modification of the one in "Computer Methods
*     for Mathematical Computation" by Forsythe, Malcolm, and Moler
*     (Prentice Hall, 1975)
*
      SUBROUTINE LINEQN(NDIM,N,A,B)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(NDIM,*),B(*)
*
*     SOLVE A SYSTEM OF LINEAR EQUATIONS A*X = B
*
*     INPUT..
*
*        NDIM = DECLARED ROW DIMENSION OF THE ARRAY CONTAINING  A.
*        N = ORDER OF THE MATRIX
*        A =  COEFFICIENT MATRIX
*        B = CONSTANT VECTOR
*
*     OUTPUT..
*
*        B = SOLUTION VECTOR
*
      NM1 = N-1
*
*
*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
*
      DO 35 K = 1,NM1
         KP1= K+1
*
*        FIND PIVOT
*
         M = K
         DO 15 I = KP1,N
            IF (DABS(A(I,K)) .GT. DABS(A(M,K))) M = I
   15    CONTINUE
         T = A(M,K)
         A(M,K) = A(K,K)
         A(K,K) = T
*
*        SKIP STEP IF PIVOT IS ZERO
*
         IF (T .EQ. 0.0D0) GO TO 35
*
*        COMPUTE MULTIPLIERS
*
         DO 20 I = KP1,N
             A(I,K) = -A(I,K)/T
   20    CONTINUE
*
*        INTERCHANGE AND ELIMINATE BY COLUMNS
*
         DO 30 J = KP1,N
             T = A(M,J)
             A(M,J) = A(K,J)
             A(K,J) = T
             IF (T .EQ. 0.0D0) GO TO 30
             DO 25 I = KP1,N
                A(I,J) = A(I,J) + A(I,K)*T
   25        CONTINUE
   30    CONTINUE
          T = B(M)
          B(M) = B(K)
          B(K) = T
          IF (T .EQ. 0.0D0) GO TO 35
          DO 32 I = KP1, N
              B(I) = B(I) + A(I,K)*T
   32     CONTINUE
   35 CONTINUE
*
*     BACK SUBSTITUTION
*
      DO 40 K = N,2,-1
        IF (A(K,K) .EQ. 0.D0) THEN
             B(K) = 0.D0
          ELSE
             B(K) = B(K)/A(K,K)
        END IF
         T = -B(K)
         DO 41 I = 1, K-1
             B(I) = B(I) + A(I,K)*T
   41    CONTINUE
   40 CONTINUE
   50 B(1) = B(1)/A(1,1)
      RETURN
      END
