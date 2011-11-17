      SUBROUTINE DGECO(A, LDA, N, IPVT, RCOND, Z) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ddot_I 
      USE dasum_I 
      USE dgefa_I 
      USE dscal_I 
      USE daxpy_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: LDA 
      INTEGER  :: N 
      REAL(DOUBLE) , INTENT(OUT) :: RCOND 
      INTEGER  :: IPVT(1) 
      REAL(DOUBLE)  :: A(LDA,1) 
      REAL(DOUBLE)  :: Z(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: INFO, J, K, KB, KP1, L 
      REAL(DOUBLE) :: EK, T, WK, WKM, ANORM, S, SM, YNORM 
!-----------------------------------------------
!
!     DGECO FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION
!     AND ESTIMATES THE CONDITION OF THE MATRIX.
!
!     IF  RCOND  IS NOT NEEDED, DGEFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW DGECO BY DGESL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DGECO BY DGESL.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DGECO BY DGEDI.
!     TO COMPUTE  INVERSE(A) , FOLLOW DGECO BY DGEDI.
!
!     ON ENTRY
!
!        A       DOUBLE PRECISION(LDA, N)
!                THE MATRIX TO BE FACTORED.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!     ON RETURN
!
!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
!                WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
!
!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.
!
!        RCOND   DOUBLE PRECISION
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.
!
!        Z       DOUBLE PRECISION(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     LINPACK DGEFA
!     BLAS DAXPY,DDOT,DSCAL,DASUM
!     FORTRAN DABS,DMAX1,DSIGN
!
!     INTERNAL VARIABLES
!
!
!
!     COMPUTE 1-NORM OF A
!
      ANORM = 0.0D0 
      DO J = 1, N 
         ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1)) 
      END DO 
!
!     FACTOR
!
      CALL DGEFA (A, LDA, N, IPVT, INFO) 
!
!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.
!
!     SOLVE TRANS(U)*W = E
!
      EK = 1.0D0 
      Z(:N) = 0.0D0 
      DO K = 1, N 
         IF (Z(K) /= 0.0D0) EK = DSIGN(EK,(-Z(K))) 
         IF (DABS(EK - Z(K)) > DABS(A(K,K))) THEN 
            S = DABS(A(K,K))/DABS(EK - Z(K)) 
            CALL DSCAL (N, S, Z, 1) 
            EK = S*EK 
         ENDIF 
         WK = EK - Z(K) 
         WKM = (-EK) - Z(K) 
         S = DABS(WK) 
         SM = DABS(WKM) 
         IF (A(K,K) /= 0.0D0) THEN 
            WK = WK/A(K,K) 
            WKM = WKM/A(K,K) 
         ELSE 
            WK = 1.0D0 
            WKM = 1.0D0 
         ENDIF 
         KP1 = K + 1 
         IF (KP1 <= N) THEN 
            DO J = KP1, N 
               SM = SM + DABS(Z(J)+WKM*A(K,J)) 
               Z(J) = Z(J) + WK*A(K,J) 
               S = S + DABS(Z(J)) 
            END DO 
            IF (S < SM) THEN 
               T = WKM - WK 
               WK = WKM 
               Z(KP1:N) = Z(KP1:N) + T*A(K,KP1:N) 
            ENDIF 
         ENDIF 
         Z(K) = WK 
      END DO 
      S = 1.0D0/DASUM(N,Z,1) 
      CALL DSCAL (N, S, Z, 1) 
!
!     SOLVE TRANS(L)*Y = W
!
      DO KB = 1, N 
         K = N + 1 - KB 
         IF (K < N) Z(K) = Z(K) + DDOT(N - K,A(K+1,K),1,Z(K+1),1) 
         IF (DABS(Z(K)) > 1.0D0) THEN 
            S = 1.0D0/DABS(Z(K)) 
            CALL DSCAL (N, S, Z, 1) 
         ENDIF 
         L = IPVT(K) 
         T = Z(L) 
         Z(L) = Z(K) 
         Z(K) = T 
      END DO 
      S = 1.0D0/DASUM(N,Z,1) 
      CALL DSCAL (N, S, Z, 1) 
!
      YNORM = 1.0D0 
!
!     SOLVE L*V = Y
!
      DO K = 1, N 
         L = IPVT(K) 
         T = Z(L) 
         Z(L) = Z(K) 
         Z(K) = T 
         IF (K < N) CALL DAXPY (N - K, T, A(K+1,K), 1, Z(K+1), 1) 
         IF (DABS(Z(K)) <= 1.0D0) CYCLE  
         S = 1.0D0/DABS(Z(K)) 
         CALL DSCAL (N, S, Z, 1) 
         YNORM = S*YNORM 
      END DO 
      S = 1.0D0/DASUM(N,Z,1) 
      CALL DSCAL (N, S, Z, 1) 
      YNORM = S*YNORM 
!
!     SOLVE  U*Z = V
!
      DO KB = 1, N 
         K = N + 1 - KB 
         IF (DABS(Z(K)) > DABS(A(K,K))) THEN 
            S = DABS(A(K,K))/DABS(Z(K)) 
            CALL DSCAL (N, S, Z, 1) 
            YNORM = S*YNORM 
         ENDIF 
         IF (A(K,K) /= 0.0D0) THEN 
            Z(K) = Z(K)/A(K,K) 
         ELSE 
            Z(K) = 1.0D0 
         ENDIF 
         T = -Z(K) 
         CALL DAXPY (K - 1, T, A(1,K), 1, Z(1), 1) 
      END DO 
!     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1) 
      CALL DSCAL (N, S, Z, 1) 
      YNORM = S*YNORM 
      IF (ANORM /= 0.0D0) THEN 
!
         RCOND = YNORM/ANORM 
      ELSE 
         RCOND = 0.0D0 
      ENDIF 
      RETURN  
      END SUBROUTINE DGECO 
