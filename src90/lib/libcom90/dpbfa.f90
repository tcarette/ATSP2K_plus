      SUBROUTINE DPBFA(ABD, LDA, N, M, INFO) 
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
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: LDA 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: M 
      INTEGER , INTENT(OUT) :: INFO 
      REAL(DOUBLE)  :: ABD(LDA,1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IK, J, JK, K, MU 
      REAL(DOUBLE) :: T, S 
!-----------------------------------------------
!
!     DPBFA FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE
!     MATRIX STORED IN BAND FORM.
!
!     DPBFA IS USUALLY CALLED BY DPBCO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
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
!
!        INFO    INTEGER
!                = 0  FOR NORMAL RETURN.
!                = K  IF THE LEADING MINOR OF ORDER  K  IS NOT
!                     POSITIVE DEFINITE.
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
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DDOT
!     FORTRAN MAX0,DSQRT
!
!     INTERNAL VARIABLES
!
!     BEGIN BLOCK WITH ...EXITS TO 40
!
!
      DO J = 1, N 
         INFO = J 
         S = 0.0D0 
         IK = M + 1 
         JK = MAX0(J - M,1) 
         MU = MAX0(M + 2 - J,1) 
         IF (M >= MU) THEN 
            DO K = MU, M 
               T = ABD(K,J) - DDOT(K - MU,ABD(IK,JK),1,ABD(MU,J),1) 
               T = T/ABD(M+1,JK) 
               ABD(K,J) = T 
               S = S + T*T 
               IK = IK - 1 
               JK = JK + 1 
            END DO 
         ENDIF 
         S = ABD(M+1,J) - S 
!     ......EXIT
         IF (S <= 0.0D0) GO TO 40 
         ABD(M+1,J) = DSQRT(S) 
      END DO 
      INFO = 0 
   40 CONTINUE 
      RETURN  
      END SUBROUTINE DPBFA 
