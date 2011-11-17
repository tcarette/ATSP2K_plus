      SUBROUTINE DGBFA(ABD, LDA, N, ML, MU, IPVT, INFO) 
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
      INTEGER , INTENT(IN) :: ML 
      INTEGER , INTENT(IN) :: MU 
      INTEGER , INTENT(OUT) :: INFO 
      INTEGER , INTENT(INOUT) :: IPVT(1) 
      REAL(DOUBLE)  :: ABD(LDA,1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, I0, J, JU, JZ, J0, J1, K, KP1, L, LM, M, MM, NM1 
      REAL(DOUBLE) :: T 
!-----------------------------------------------
!
!     DGBFA FACTORS A DOUBLE PRECISION BAND MATRIX BY ELIMINATION.
!
!     DGBFA IS USUALLY CALLED BY DGBCO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
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
!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
!                     INDICATE THAT DGBSL WILL DIVIDE BY ZERO IF
!                     CALLED.  USE  RCOND  IN DGBCO FOR A RELIABLE
!                     INDICATION OF SINGULARITY.
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
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DSCAL,IDAMAX
!     FORTRAN MAX0,MIN0
!
!     INTERNAL VARIABLES
!
!
!
      M = ML + MU + 1 
      INFO = 0 
!
!     ZERO INITIAL FILL-IN COLUMNS
!
      J0 = MU + 2 
      J1 = MIN0(N,M) - 1 
      IF (J1 >= J0) THEN 
         DO JZ = J0, J1 
            I0 = M + 1 - JZ 
            ABD(I0:ML,JZ) = 0.0D0 
         END DO 
      ENDIF 
      JZ = J1 
      JU = 0 
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
      NM1 = N - 1 
      IF (NM1 >= 1) THEN 
         DO K = 1, NM1 
            KP1 = K + 1 
!
!        ZERO NEXT FILL-IN COLUMN
!
            JZ = JZ + 1 
            IF (JZ <= N) THEN 
               IF (ML >= 1) THEN 
                  ABD(:ML,JZ) = 0.0D0 
               ENDIF 
            ENDIF 
!
!        FIND L = PIVOT INDEX
!
            LM = MIN0(ML,N - K) 
            L = IDAMAX(LM + 1,ABD(M,K),1) + M - 1 
            IPVT(K) = L + K - M 
!
!        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
!
            IF (ABD(L,K) /= 0.0D0) THEN 
!
!           INTERCHANGE IF NECESSARY
!
               IF (L /= M) THEN 
                  T = ABD(L,K) 
                  ABD(L,K) = ABD(M,K) 
                  ABD(M,K) = T 
               ENDIF 
!
!           COMPUTE MULTIPLIERS
!
               T = -1.0D0/ABD(M,K) 
               CALL DSCAL (LM, T, ABD(M+1,K), 1) 
!
!           ROW ELIMINATION WITH COLUMN INDEXING
!
               JU = MIN0(MAX0(JU,MU + IPVT(K)),N) 
               MM = M 
               IF (JU >= KP1) THEN 
                  DO J = KP1, JU 
                     L = L - 1 
                     MM = MM - 1 
                     T = ABD(L,J) 
                     IF (L /= MM) THEN 
                        ABD(L,J) = ABD(MM,J) 
                        ABD(MM,J) = T 
                     ENDIF 
                     CALL DAXPY (LM, T, ABD(M+1,K), 1, ABD(MM+1,J), 1) 
                  END DO 
               ENDIF 
            ELSE 
               INFO = K 
            ENDIF 
         END DO 
      ENDIF 
      IPVT(N) = N 
      IF (ABD(M,N) == 0.0D0) INFO = N 
      RETURN  
      END SUBROUTINE DGBFA 
