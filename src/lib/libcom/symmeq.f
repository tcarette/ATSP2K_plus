*
*     ------------------------------------------------------------------
*        S Y M M E Q   A N D  S Y M M S L
*     ------------------------------------------------------------------
*
*        This routine is a modification of the one in "Computer Methods
*     for Mathematical Computation" by Forsythe, Malcolm, and Moler
*     (Prentice Hall, 1975) to solve a singular system of equations.
*
      SUBROUTINE SYMMEQ(NDIM,N,A,X)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(NDIM,*),X(*)
*
*
*     SOLVE A SYSTEM OF LINEAR EQUATIONS A*X = 0
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
*        X = SOLUTION VECTOR with X(1) unchanged.
*
      NM1 = N-1
*
*
*      L U FACTORIZATION WITHOUT PIVOTING
*
      DO 35 K = N,3,-1
         KP1= K-1
         T = A(K,K)
         IF (T .EQ. 0.0D0) GO TO 35
*
*        COMPUTE MULTIPLIERS
*
         DO 20 I = KP1,2,-1
             A(I,K) = -A(I,K)/T
   20    CONTINUE
*
*        INTERCHANGE AND ELIMINATE BY COLUMNS
*
         DO 30 J = KP1,2,-1
             T = A(K,J)
             IF (T .EQ. 0.0D0) GO TO 30
             DO 25 I = KP1,2,-1
                A(I,J) = A(I,J) + A(I,K)*T
   25        CONTINUE
   30    CONTINUE
   35 CONTINUE
*
*     At this point it is assumed that the LU factorization
*     has already been performed.
*
      ENTRY SYMMSL(NDIM,N,A,X)
*
      DO 36 I = 2,N
        X(I) = - A(I,1)
   36 CONTINUE
      DO 50 K=N,3,-1
        T=X(K)
        IF (T.EQ.0.0D0) GO TO 50
        DO 40 I=K-1,2,-1
          X(I)=X(I)+A(I,K)*T
   40   CONTINUE
   50   CONTINUE
*
*     BACK SUBSTITUTION
*
      DO 80 K =2,N-1
        IF (A(K,K) .EQ. 0.D0) THEN
             X(K) = 0.D0
          ELSE
             X(K) = X(K)/A(K,K)
        END IF
         T = -X(K)
         DO 70 I =K+1,N
             X(I) = X(I) + A(I,K)*T
   70    CONTINUE
   80 CONTINUE
      X(N) = X(N)/A(N,N)
      RETURN
      END
