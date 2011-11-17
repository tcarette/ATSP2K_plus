      SUBROUTINE DGBSL(ABD, LDA, N, ML, MU, IPVT, B, JOB) 
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
      INTEGER , INTENT(IN) :: ML 
      INTEGER , INTENT(IN) :: MU 
      INTEGER , INTENT(IN) :: JOB 
      INTEGER , INTENT(IN) :: IPVT(1) 
      REAL(DOUBLE)  :: ABD(LDA,1) 
      REAL(DOUBLE)  :: B(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, KB, L, LA, LB, LM, M, NM1 
      REAL(DOUBLE) :: T 
!-----------------------------------------------
!
!     DGBSL SOLVES THE DOUBLE PRECISION BAND SYSTEM
!     A * X = B  OR  TRANS(A) * X = B
!     USING THE FACTORS COMPUTED BY DGBCO OR DGBFA.
!
!     ON ENTRY
!
!        ABD     DOUBLE PRECISION(LDA, N)
!                THE OUTPUT FROM DGBCO OR DGBFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!
!        N       INTEGER
!                THE ORDER OF THE ORIGINAL MATRIX.
!
!        ML      INTEGER
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
!
!        MU      INTEGER
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!
!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM DGBCO OR DGBFA.
!
!        B       DOUBLE PRECISION(N)
!                THE RIGHT HAND SIDE VECTOR.
!
!        JOB     INTEGER
!                = 0         TO SOLVE  A*X = B ,
!                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE
!                            TRANS(A)  IS THE TRANSPOSE.
!
!     ON RETURN
!
!        B       THE SOLUTION VECTOR  X .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
!        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
!        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
!        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
!        CALLED CORRECTLY AND IF DGBCO HAS SET RCOND .GT. 0.0
!        OR DGBFA HAS SET INFO .EQ. 0 .
!
!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
!           IF (RCOND IS TOO SMALL) GO TO ...
!           DO 10 J = 1, P
!              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
!        10 CONTINUE
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
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
      M = MU + ML + 1 
      NM1 = N - 1 
      IF (JOB == 0) THEN 
!
!        JOB = 0 , SOLVE  A * X = B
!        FIRST SOLVE L*Y = B
!
         IF (ML /= 0) THEN 
            IF (NM1 >= 1) THEN 
               DO K = 1, NM1 
                  LM = MIN0(ML,N - K) 
                  L = IPVT(K) 
                  T = B(L) 
                  IF (L /= K) THEN 
                     B(L) = B(K) 
                     B(K) = T 
                  ENDIF 
                  CALL DAXPY (LM, T, ABD(M+1,K), 1, B(K+1), 1) 
               END DO 
            ENDIF 
         ENDIF 
         DO KB = 1, N 
            K = N + 1 - KB 
            B(K) = B(K)/ABD(M,K) 
            LM = MIN0(K,M) - 1 
            LA = M - LM 
            LB = K - LM 
            T = -B(K) 
            CALL DAXPY (LM, T, ABD(LA,K), 1, B(LB), 1) 
         END DO 
      ELSE 
         DO K = 1, N 
            LM = MIN0(K,M) - 1 
            LA = M - LM 
            LB = K - LM 
            T = DDOT(LM,ABD(LA,K),1,B(LB),1) 
            B(K) = (B(K)-T)/ABD(M,K) 
         END DO 
!
!        NOW SOLVE TRANS(L)*X = Y
!
         IF (ML /= 0) THEN 
            IF (NM1 >= 1) THEN 
               DO KB = 1, NM1 
                  K = N - KB 
                  LM = MIN0(ML,N - K) 
                  B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1) 
                  L = IPVT(K) 
                  IF (L == K) CYCLE  
                  T = B(L) 
                  B(L) = B(K) 
                  B(K) = T 
               END DO 
            ENDIF 
         ENDIF 
      ENDIF 
      RETURN  
      END SUBROUTINE DGBSL 
