      SUBROUTINE DPBSL(ABD, LDA, N, M, B) 
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
      USE daxpy_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: LDA 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: M 
      REAL(DOUBLE)  :: ABD(LDA,1) 
      REAL(DOUBLE)  :: B(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, KB, LA, LB, LM 
      REAL(DOUBLE) :: T 
!-----------------------------------------------
!
!     DPBSL SOLVES THE DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE
!     BAND SYSTEM  A*X = B
!     USING THE FACTORS COMPUTED BY DPBCO OR DPBFA.
!
!     ON ENTRY
!
!        ABD     DOUBLE PRECISION(LDA, N)
!                THE OUTPUT FROM DPBCO OR DPBFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        M       INTEGER
!                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!
!        B       DOUBLE PRECISION(N)
!                THE RIGHT HAND SIDE VECTOR.
!
!     ON RETURN
!
!        B       THE SOLUTION VECTOR  X .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
!        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES
!        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE
!        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED
!        CORRECTLY AND  INFO .EQ. 0 .
!
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL DPBCO(ABD,LDA,N,RCOND,Z,INFO)
!           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...
!           DO 10 J = 1, P
!              CALL DPBSL(ABD,LDA,N,C(1,J))
!        10 CONTINUE
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DDOT
!     FORTRAN MIN0
!
!     INTERNAL VARIABLES
!
!
!     SOLVE TRANS(R)*Y = B
!
      DO K = 1, N 
         LM = MIN0(K - 1,M) 
         LA = M + 1 - LM 
         LB = K - LM 
         T = DDOT(LM,ABD(LA,K),1,B(LB),1) 
         B(K) = (B(K)-T)/ABD(M+1,K) 
      END DO 
!
!     SOLVE R*X = Y
!
      DO KB = 1, N 
         K = N + 1 - KB 
         LM = MIN0(K - 1,M) 
         LA = M + 1 - LM 
         LB = K - LM 
         B(K) = B(K)/ABD(M+1,K) 
         T = -B(K) 
         CALL DAXPY (LM, T, ABD(LA,K), 1, B(LB), 1) 
      END DO 
!
!     SOLVE R*X = Y
!
      RETURN  
      END SUBROUTINE DPBSL 
