      SUBROUTINE DGBCO(ABD, LDA, N, ML, MU, IPVT, RCOND, Z) 
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
      USE dgbfa_I 
      USE dscal_I 
      USE daxpy_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: LDA 
      INTEGER  :: N 
      INTEGER  :: ML 
      INTEGER  :: MU 
      REAL(DOUBLE) , INTENT(OUT) :: RCOND 
      INTEGER  :: IPVT(1) 
      REAL(DOUBLE)  :: ABD(LDA,1) 
      REAL(DOUBLE)  :: Z(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IS, INFO, J, JU, K, KB, KP1, L, LA, LM, LZ, M, MM 
      REAL(DOUBLE) :: EK, T, WK, WKM, ANORM, S, SM, YNORM 
!-----------------------------------------------
!
!     DGBCO FACTORS A DOUBLE PRECISION BAND MATRIX BY GAUSSIAN
!     ELIMINATION AND ESTIMATES THE CONDITION OF THE MATRIX.
!
!     IF  RCOND  IS NOT NEEDED, DGBFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW DGBCO BY DGBSL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DGBCO BY DGBSL.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DGBCO BY DGBDI.
!
!     ON ENTRY
!
!        ABD     DOUBLE PRECISION(LDA, N)
!                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
!                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
!                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
!                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
!                SEE THE COMMENTS BELOW FOR DETAILS.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!                LDA MUST BE .GE. 2*ML + MU + 1 .
!
!        N       INTEGER
!                THE ORDER OF THE ORIGINAL MATRIX.
!
!        ML      INTEGER
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
!                0 .LE. ML .LT. N .
!
!        MU      INTEGER
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!                0 .LE. MU .LT. N .
!                MORE EFFICIENT IF  ML .LE. MU .
!
!     ON RETURN
!
!        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
!                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
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
!     BAND STORAGE
!
!           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT
!           WILL SET UP THE INPUT.
!
!                   ML = (BAND WIDTH BELOW THE DIAGONAL)
!                   MU = (BAND WIDTH ABOVE THE DIAGONAL)
!                   M = ML + MU + 1
!                   DO 20 J = 1, N
!                      I1 = MAX0(1, J-MU)
!                      I2 = MIN0(N, J+ML)
!                      DO 10 I = I1, I2
!                         K = I - J + M
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
!           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
!           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
!           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
!           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
!           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.
!
!     EXAMPLE..  IF THE ORIGINAL MATRIX IS
!
!           11 12 13  0  0  0
!           21 22 23 24  0  0
!            0 32 33 34 35  0
!            0  0 43 44 45 46
!            0  0  0 54 55 56
!            0  0  0  0 65 66
!
!      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN
!
!            *  *  *  +  +  +  , * = NOT USED
!            *  * 13 24 35 46  , + = USED FOR PIVOTING
!            * 12 23 34 45 56
!           11 22 33 44 55 66
!           21 32 43 54 65  *
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     LINPACK DGBFA
!     BLAS DAXPY,DDOT,DSCAL,DASUM
!     FORTRAN DABS,DMAX1,MAX0,MIN0,DSIGN
!
!     INTERNAL VARIABLES
!
!
!
!     COMPUTE 1-NORM OF A
!
      ANORM = 0.0D0 
      L = ML + 1 
      IS = L + MU 
      DO J = 1, N 
         ANORM = DMAX1(ANORM,DASUM(L,ABD(IS,J),1)) 
         IF (IS > ML + 1) IS = IS - 1 
         IF (J <= MU) L = L + 1 
         IF (J < N - ML) CYCLE  
         L = L - 1 
      END DO 
!
!     FACTOR
!
      CALL DGBFA (ABD, LDA, N, ML, MU, IPVT, INFO) 
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
      M = ML + MU + 1 
      JU = 0 
      DO K = 1, N 
         IF (Z(K) /= 0.0D0) EK = DSIGN(EK,(-Z(K))) 
         IF (DABS(EK - Z(K)) > DABS(ABD(M,K))) THEN 
            S = DABS(ABD(M,K))/DABS(EK - Z(K)) 
            CALL DSCAL (N, S, Z, 1) 
            EK = S*EK 
         ENDIF 
         WK = EK - Z(K) 
         WKM = (-EK) - Z(K) 
         S = DABS(WK) 
         SM = DABS(WKM) 
         IF (ABD(M,K) /= 0.0D0) THEN 
            WK = WK/ABD(M,K) 
            WKM = WKM/ABD(M,K) 
         ELSE 
            WK = 1.0D0 
            WKM = 1.0D0 
         ENDIF 
         KP1 = K + 1 
         JU = MIN0(MAX0(JU,MU + IPVT(K)),N) 
         MM = M 
         IF (KP1 <= JU) THEN 
            DO J = KP1, JU 
               MM = MM - 1 
               SM = SM + DABS(Z(J)+WKM*ABD(MM,J)) 
               Z(J) = Z(J) + WK*ABD(MM,J) 
               S = S + DABS(Z(J)) 
            END DO 
            IF (S < SM) THEN 
               T = WKM - WK 
               WK = WKM 
               MM = M 
               DO J = KP1, JU 
                  MM = MM - 1 
                  Z(J) = Z(J) + T*ABD(MM,J) 
               END DO 
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
         LM = MIN0(ML,N - K) 
         IF (K < N) Z(K) = Z(K) + DDOT(LM,ABD(M+1,K),1,Z(K+1),1) 
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
         LM = MIN0(ML,N - K) 
         IF (K < N) CALL DAXPY (LM, T, ABD(M+1,K), 1, Z(K+1), 1) 
         IF (DABS(Z(K)) <= 1.0D0) CYCLE  
         S = 1.0D0/DABS(Z(K)) 
         CALL DSCAL (N, S, Z, 1) 
         YNORM = S*YNORM 
      END DO 
      S = 1.0D0/DASUM(N,Z,1) 
      CALL DSCAL (N, S, Z, 1) 
      YNORM = S*YNORM 
!
!     SOLVE  U*Z = W
!
      DO KB = 1, N 
         K = N + 1 - KB 
         IF (DABS(Z(K)) > DABS(ABD(M,K))) THEN 
            S = DABS(ABD(M,K))/DABS(Z(K)) 
            CALL DSCAL (N, S, Z, 1) 
            YNORM = S*YNORM 
         ENDIF 
         IF (ABD(M,K) /= 0.0D0) THEN 
            Z(K) = Z(K)/ABD(M,K) 
         ELSE 
            Z(K) = 1.0D0 
         ENDIF 
         LM = MIN0(K,M) - 1 
         LA = M - LM 
         LZ = K - LM 
         T = -Z(K) 
         CALL DAXPY (LM, T, ABD(LA,K), 1, Z(LZ), 1) 
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
      END SUBROUTINE DGBCO 
