      SUBROUTINE DGEMV(TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE lsame_I 
      USE xerbla_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: M 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: LDA 
      INTEGER , INTENT(IN) :: INCX 
      INTEGER , INTENT(IN) :: INCY 
      REAL(DOUBLE) , INTENT(IN) :: ALPHA 
      REAL(DOUBLE) , INTENT(IN) :: BETA 
      CHARACTER  :: TRANS 
      REAL(DOUBLE) , INTENT(IN) :: A(LDA,*) 
      REAL(DOUBLE) , INTENT(IN) :: X(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: Y(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ONE = 1.0D+0 
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY 
      REAL(DOUBLE) :: TEMP 
!-----------------------------------------------
!   I n t r i n s i c  F u n c t i o n s
!-----------------------------------------------
      INTRINSIC MAX 
!-----------------------------------------------
!     .. Scalar Arguments ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
!     .. Local Scalars ..
!     .. External Functions ..
!     .. External Subroutines ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0 
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND. .NOT.LSAME(&
         TRANS,'C')) THEN 
         INFO = 1 
      ELSE IF (M < 0) THEN 
         INFO = 2 
      ELSE IF (N < 0) THEN 
         INFO = 3 
      ELSE IF (LDA < MAX(1,M)) THEN 
         INFO = 6 
      ELSE IF (INCX == 0) THEN 
         INFO = 8 
      ELSE IF (INCY == 0) THEN 
         INFO = 11 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DGEMV ', INFO) 
         RETURN  
      ENDIF 
!
!     Quick return if possible.
!
      IF (M==0 .OR. N==0 .OR. ALPHA==ZERO .AND. BETA==ONE) RETURN  
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF (LSAME(TRANS,'N')) THEN 
         LENX = N 
         LENY = M 
      ELSE 
         LENX = M 
         LENY = N 
      ENDIF 
      IF (INCX > 0) THEN 
         KX = 1 
      ELSE 
         KX = 1 - (LENX - 1)*INCX 
      ENDIF 
      IF (INCY > 0) THEN 
         KY = 1 
      ELSE 
         KY = 1 - (LENY - 1)*INCY 
      ENDIF 
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF (BETA /= ONE) THEN 
         IF (INCY == 1) THEN 
            IF (BETA == ZERO) THEN 
               Y(:LENY) = ZERO 
            ELSE 
               Y(:LENY) = BETA*Y(:LENY) 
            ENDIF 
         ELSE 
            IY = KY 
            IF (BETA == ZERO) THEN 
               Y(IY:(LENY-1)*INCY+IY:INCY) = ZERO 
            ELSE 
               Y(IY:(LENY-1)*INCY+IY:INCY) = BETA*Y(IY:(LENY-1)*INCY+IY:INCY) 
            ENDIF 
         ENDIF 
      ENDIF 
      IF (ALPHA == ZERO) RETURN  
      IF (LSAME(TRANS,'N')) THEN 
!
!        Form  y := alpha*A*x + y.
!
         JX = KX 
         IF (INCY == 1) THEN 
            DO J = 1, N 
               IF (X(JX) /= ZERO) THEN 
                  TEMP = ALPHA*X(JX) 
                  Y(:M) = Y(:M) + TEMP*A(:M,J) 
               ENDIF 
               JX = JX + INCX 
            END DO 
         ELSE 
            DO J = 1, N 
               IF (X(JX) /= ZERO) THEN 
                  TEMP = ALPHA*X(JX) 
                  IY = KY 
                  Y(IY:(M-1)*INCY+IY:INCY) = Y(IY:(M-1)*INCY+IY:INCY) + TEMP*A(&
                     :M,J) 
               ENDIF 
               JX = JX + INCX 
            END DO 
         ENDIF 
      ELSE 
!
!        Form  y := alpha*A'*x + y.
!
         JY = KY 
         IF (INCX == 1) THEN 
            DO J = 1, N 
               TEMP = ZERO 
               TEMP = TEMP + SUM(A(:M,J)*X(:M)) 
               Y(JY) = Y(JY) + ALPHA*TEMP 
               JY = JY + INCY 
            END DO 
         ELSE 
            DO J = 1, N 
               TEMP = ZERO 
               IX = KX 
               TEMP = TEMP + SUM(A(:M,J)*X(IX:(M-1)*INCX+IX:INCX)) 
               Y(JY) = Y(JY) + ALPHA*TEMP 
               JY = JY + INCY 
            END DO 
         ENDIF 
      ENDIF 
!
      RETURN  
!
!     End of DGEMV .
!
      END SUBROUTINE DGEMV 
