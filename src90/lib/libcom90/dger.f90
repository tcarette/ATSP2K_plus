      SUBROUTINE DGER(M, N, ALPHA, X, INCX, Y, INCY, A, LDA) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE xerbla_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: M 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: INCX 
      INTEGER , INTENT(IN) :: INCY 
      INTEGER , INTENT(IN) :: LDA 
      REAL(DOUBLE) , INTENT(IN) :: ALPHA 
      REAL(DOUBLE) , INTENT(IN) :: X(*) 
      REAL(DOUBLE) , INTENT(IN) :: Y(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: A(LDA,*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, INFO, IX, J, JY, KX 
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
!  DGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
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
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
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
!     .. External Subroutines ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0 
      IF (M < 0) THEN 
         INFO = 1 
      ELSE IF (N < 0) THEN 
         INFO = 2 
      ELSE IF (INCX == 0) THEN 
         INFO = 5 
      ELSE IF (INCY == 0) THEN 
         INFO = 7 
      ELSE IF (LDA < MAX(1,M)) THEN 
         INFO = 9 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DGER  ', INFO) 
         RETURN  
      ENDIF 
!
!     Quick return if possible.
!
      IF (M==0 .OR. N==0 .OR. ALPHA==ZERO) RETURN  
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF (INCY > 0) THEN 
         JY = 1 
      ELSE 
         JY = 1 - (N - 1)*INCY 
      ENDIF 
      IF (INCX == 1) THEN 
         DO J = 1, N 
            IF (Y(JY) /= ZERO) THEN 
               TEMP = ALPHA*Y(JY) 
               A(:M,J) = A(:M,J) + X(:M)*TEMP 
            ENDIF 
            JY = JY + INCY 
         END DO 
      ELSE 
         IF (INCX > 0) THEN 
            KX = 1 
         ELSE 
            KX = 1 - (M - 1)*INCX 
         ENDIF 
         DO J = 1, N 
            IF (Y(JY) /= ZERO) THEN 
               TEMP = ALPHA*Y(JY) 
               IX = KX 
               A(:M,J) = A(:M,J) + X(IX:(M-1)*INCX+IX:INCX)*TEMP 
            ENDIF 
            JY = JY + INCY 
         END DO 
      ENDIF 
!
      RETURN  
!
!     End of DGER  .
!
      END SUBROUTINE DGER 
