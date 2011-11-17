      SUBROUTINE DPBCO(ABD, LDA, N, M, RCOND, Z, INFO) 
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
      USE dpbfa_I 
      USE dscal_I 
      USE daxpy_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: LDA 
      INTEGER  :: N 
      INTEGER  :: M 
      INTEGER  :: INFO 
      REAL(DOUBLE) , INTENT(OUT) :: RCOND 
      REAL(DOUBLE)  :: ABD(LDA,1) 
      REAL(DOUBLE)  :: Z(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, J2, K, KB, KP1, L, LA, LB, LM, MU 
      REAL(DOUBLE) :: EK, T, WK, WKM, ANORM, S, SM, YNORM 
!-----------------------------------------------
!
!     DPBCO FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE
!     MATRIX STORED IN BAND FORM AND ESTIMATES THE CONDITION OF THE
!     MATRIX.
!
!     IF  RCOND  IS NOT NEEDED, DPBFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW DPBCO BY DPBSL.
!     TO COMPUTE  INVERSE(A)*C , FOLLOW DPBCO BY DPBSL.
!     TO COMPUTE  DETERMINANT(A) , FOLLOW DPBCO BY DPBDI.
!
!     ON ENTRY
!
!        ABD     DOUBLE PRECISION(LDA, N)
!                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER
!                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE
!                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE
!                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!                LDA MUST BE .GE. M + 1 .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        M       INTEGER
!                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!                0 .LE. M .LT. N .
!
!     ON RETURN
!
!        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND
!                FORM, SO THAT  A = TRANS(R)*R .
!                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
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
!                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED.
!
!        Z       DOUBLE PRECISION(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS SINGULAR TO WORKING PRECISION, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!                IF  INFO .NE. 0 , Z  IS UNCHANGED.
!
!        INFO    INTEGER
!                = 0  FOR NORMAL RETURN.
!                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
!                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
!
!     BAND STORAGE
!
!           IF  A  IS A SYMMETRIC POSITIVE DEFINITE BAND MATRIX,
!           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT.
!
!                   M = (BAND WIDTH ABOVE DIAGONAL)
!                   DO 20 J = 1, N
!                      I1 = MAX0(1, J-M)
!                      DO 10 I = I1, J
!                         K = I-J+M+1
!                         ABD(K,J) = A(I,J)
!                10    CONTINUE
!                20 CONTINUE
!
!           THIS USES  M + 1  ROWS OF  A , EXCEPT FOR THE  M BY M
!           UPPER LEFT TRIANGLE, WHICH IS IGNORED.
!
!     EXAMPLE..  IF THE ORIGINAL MATRIX IS
!
!           11 12 13  0  0  0
!           12 22 23 24  0  0
!           13 23 33 34 35  0
!            0 24 34 44 45 46
!            0  0 35 45 55 56
!            0  0  0 46 56 66
!
!     THEN  N = 6 , M = 2  AND  ABD  SHOULD CONTAIN
!
!            *  * 13 24 35 46
!            * 12 23 34 45 56
!           11 22 33 44 55 66
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     LINPACK DPBFA
!     BLAS DAXPY,DDOT,DSCAL,DASUM
!     FORTRAN DABS,DMAX1,MAX0,MIN0,DREAL,DSIGN
!
!     INTERNAL VARIABLES
!
!
!
!     FIND NORM OF A
!
      DO J = 1, N 
         L = MIN0(J,M + 1) 
         MU = MAX0(M + 2 - J,1) 
         Z(J) = DASUM(L,ABD(MU,J),1) 
         K = J - L 
         IF (M < MU) CYCLE  
         DO I = MU, M 
            K = K + 1 
            Z(K) = Z(K) + DABS(ABD(I,J)) 
         END DO 
      END DO 
      ANORM = 0.0D0 
      DO J = 1, N 
         ANORM = DMAX1(ANORM,Z(J)) 
      END DO 
!
!     FACTOR
!
      CALL DPBFA (ABD, LDA, N, M, INFO) 
      IF (INFO == 0) THEN 
!
!        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
!        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .
!        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!        SOLVE TRANS(R)*W = E
!
         EK = 1.0D0 
         Z(:N) = 0.0D0 
         DO K = 1, N 
            IF (Z(K) /= 0.0D0) EK = DSIGN(EK,(-Z(K))) 
            IF (DABS(EK - Z(K)) > ABD(M+1,K)) THEN 
               S = ABD(M+1,K)/DABS(EK - Z(K)) 
               CALL DSCAL (N, S, Z, 1) 
               EK = S*EK 
            ENDIF 
            WK = EK - Z(K) 
            WKM = (-EK) - Z(K) 
            S = DABS(WK) 
            SM = DABS(WKM) 
            WK = WK/ABD(M+1,K) 
            WKM = WKM/ABD(M+1,K) 
            KP1 = K + 1 
            J2 = MIN0(K + M,N) 
            I = M + 1 
            IF (KP1 <= J2) THEN 
               DO J = KP1, J2 
                  I = I - 1 
                  SM = SM + DABS(Z(J)+WKM*ABD(I,J)) 
                  Z(J) = Z(J) + WK*ABD(I,J) 
                  S = S + DABS(Z(J)) 
               END DO 
               IF (S < SM) THEN 
                  T = WKM - WK 
                  WK = WKM 
                  I = M + 1 
                  DO J = KP1, J2 
                     I = I - 1 
                     Z(J) = Z(J) + T*ABD(I,J) 
                  END DO 
               ENDIF 
            ENDIF 
            Z(K) = WK 
         END DO 
         S = 1.0D0/DASUM(N,Z,1) 
         CALL DSCAL (N, S, Z, 1) 
!
!        SOLVE  R*Y = W
!
         DO KB = 1, N 
            K = N + 1 - KB 
            IF (DABS(Z(K)) > ABD(M+1,K)) THEN 
               S = ABD(M+1,K)/DABS(Z(K)) 
               CALL DSCAL (N, S, Z, 1) 
            ENDIF 
            Z(K) = Z(K)/ABD(M+1,K) 
            LM = MIN0(K - 1,M) 
            LA = M + 1 - LM 
            LB = K - LM 
            T = -Z(K) 
            CALL DAXPY (LM, T, ABD(LA,K), 1, Z(LB), 1) 
         END DO 
         S = 1.0D0/DASUM(N,Z,1) 
         CALL DSCAL (N, S, Z, 1) 
!
         YNORM = 1.0D0 
!
!        SOLVE TRANS(R)*V = Y
!
         DO K = 1, N 
            LM = MIN0(K - 1,M) 
            LA = M + 1 - LM 
            LB = K - LM 
            Z(K) = Z(K) - DDOT(LM,ABD(LA,K),1,Z(LB),1) 
            IF (DABS(Z(K)) > ABD(M+1,K)) THEN 
               S = ABD(M+1,K)/DABS(Z(K)) 
               CALL DSCAL (N, S, Z, 1) 
               YNORM = S*YNORM 
            ENDIF 
            Z(K) = Z(K)/ABD(M+1,K) 
         END DO 
         S = 1.0D0/DASUM(N,Z,1) 
         CALL DSCAL (N, S, Z, 1) 
         YNORM = S*YNORM 
!
!        SOLVE  R*Z = W
!
         DO KB = 1, N 
            K = N + 1 - KB 
            IF (DABS(Z(K)) > ABD(M+1,K)) THEN 
               S = ABD(M+1,K)/DABS(Z(K)) 
               CALL DSCAL (N, S, Z, 1) 
               YNORM = S*YNORM 
            ENDIF 
            Z(K) = Z(K)/ABD(M+1,K) 
            LM = MIN0(K - 1,M) 
            LA = M + 1 - LM 
            LB = K - LM 
            T = -Z(K) 
            CALL DAXPY (LM, T, ABD(LA,K), 1, Z(LB), 1) 
         END DO 
!        MAKE ZNORM = 1.0
         S = 1.0D0/DASUM(N,Z,1) 
         CALL DSCAL (N, S, Z, 1) 
         YNORM = S*YNORM 
         IF (ANORM /= 0.0D0) THEN 
!
            RCOND = YNORM/ANORM 
         ELSE 
            RCOND = 0.0D0 
         ENDIF 
      ENDIF 
      RETURN  
      END SUBROUTINE DPBCO 
