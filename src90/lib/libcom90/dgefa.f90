      SUBROUTINE DGEFA(A, LDA, N, IPVT, INFO) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE idamax_I 
      USE dscal_I 
      USE daxpy_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: LDA 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(OUT) :: INFO 
      INTEGER , INTENT(OUT) :: IPVT(1) 
      REAL(DOUBLE)  :: A(LDA,1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, K, KP1, L, NM1 
      REAL(DOUBLE) :: T 
!-----------------------------------------------
!
!     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.
!
!     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
!     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .
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
!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
!                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO
!                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE
!                     INDICATION OF SINGULARITY.
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DSCAL,IDAMAX
!
!     INTERNAL VARIABLES
!
!
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
      INFO = 0 
      NM1 = N - 1 
      IF (NM1 >= 1) THEN 
         DO K = 1, NM1 
            KP1 = K + 1 
!
!        FIND L = PIVOT INDEX
!
            L = IDAMAX(N - K + 1,A(K,K),1) + K - 1 
            IPVT(K) = L 
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
            IF (A(L,K) /= 0.0D0) THEN 
!
!           INTERCHANGE IF NECESSARY
!
               IF (L /= K) THEN 
                  T = A(L,K) 
                  A(L,K) = A(K,K) 
                  A(K,K) = T 
               ENDIF 
!
!           COMPUTE MULTIPLIERS
!
               T = -1.0D0/A(K,K) 
               CALL DSCAL (N - K, T, A(K+1,K), 1) 
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
               DO J = KP1, N 
                  T = A(L,J) 
                  IF (L /= K) THEN 
                     A(L,J) = A(K,J) 
                     A(K,J) = T 
                  ENDIF 
                  CALL DAXPY (N - K, T, A(K+1,K), 1, A(K+1,J), 1) 
               END DO 
            ELSE 
               INFO = K 
            ENDIF 
         END DO 
      ENDIF 
      IPVT(N) = N 
      IF (A(N,N) == 0.0D0) INFO = N 
      RETURN  
      END SUBROUTINE DGEFA 
