      SUBROUTINE DSPR2(UPLO, N, ALPHA, X, INCX, Y, INCY, AP) 
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
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: INCX 
      INTEGER , INTENT(IN) :: INCY 
      REAL(DOUBLE) , INTENT(IN) :: ALPHA 
      CHARACTER  :: UPLO 
      REAL(DOUBLE) , INTENT(IN) :: X(*) 
      REAL(DOUBLE) , INTENT(IN) :: Y(*) 
      REAL(DOUBLE) , INTENT(INOUT) :: AP(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: ZERO = 0.0D+0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, INFO, IX, IY, J, JX, JY, K, KK, KX, KY 
      REAL(DOUBLE) :: TEMP1, TEMP2 
!-----------------------------------------------
!     .. Scalar Arguments ..
!     .. Array Arguments ..
!     ..
!
!  Purpose
!  =======
!
!  DSPR2  performs the symmetric rank 2 operation
!
!     A := alpha*x*y' + alpha*y*x' + A,
!
!  where alpha is a scalar, x and y are n element vectors and A is an
!  n by n symmetric matrix, supplied in packed form.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the matrix A is supplied in the packed
!           array AP as follows:
!
!              UPLO = 'U' or 'u'   The upper triangular part of A is
!                                  supplied in AP.
!
!              UPLO = 'L' or 'l'   The lower triangular part of A is
!                                  supplied in AP.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
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
!  AP     - DOUBLE PRECISION array of DIMENSION at least
!           ( ( n*( n + 1 ) )/2 ).
!           Before entry with  UPLO = 'U' or 'u', the array AP must
!           contain the upper triangular part of the symmetric matrix
!           packed sequentially, column by column, so that AP( 1 )
!           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!           and a( 2, 2 ) respectively, and so on. On exit, the array
!           AP is overwritten by the upper triangular part of the
!           updated matrix.
!           Before entry with UPLO = 'L' or 'l', the array AP must
!           contain the lower triangular part of the symmetric matrix
!           packed sequentially, column by column, so that AP( 1 )
!           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!           and a( 3, 1 ) respectively, and so on. On exit, the array
!           AP is overwritten by the lower triangular part of the
!           updated matrix.
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
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0 
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN 
         INFO = 1 
      ELSE IF (N < 0) THEN 
         INFO = 2 
      ELSE IF (INCX == 0) THEN 
         INFO = 5 
      ELSE IF (INCY == 0) THEN 
         INFO = 7 
      ENDIF 
      IF (INFO /= 0) THEN 
         CALL XERBLA ('DSPR2 ', INFO) 
         RETURN  
      ENDIF 
!
!     Quick return if possible.
!
      IF (N==0 .OR. ALPHA==ZERO) RETURN  
!
!     Set up the start points in X and Y if the increments are not both
!     unity.
!
      IF (INCX/=1 .OR. INCY/=1) THEN 
         IF (INCX > 0) THEN 
            KX = 1 
         ELSE 
            KX = 1 - (N - 1)*INCX 
         ENDIF 
         IF (INCY > 0) THEN 
            KY = 1 
         ELSE 
            KY = 1 - (N - 1)*INCY 
         ENDIF 
         JX = KX 
         JY = KY 
      ENDIF 
!
!     Start the operations. In this version the elements of the array AP
!     are accessed sequentially with one pass through AP.
!
      KK = 1 
      IF (LSAME(UPLO,'U')) THEN 
!
!        Form  A  when upper triangle is stored in AP.
!
         IF (INCX==1 .AND. INCY==1) THEN 
            DO J = 1, N 
               IF (X(J)/=ZERO .OR. Y(J)/=ZERO) THEN 
                  TEMP1 = ALPHA*Y(J) 
                  TEMP2 = ALPHA*X(J) 
                  K = KK 
                  AP(K:J-1+K) = AP(K:J-1+K) + X(:J)*TEMP1 + Y(:J)*TEMP2 
               ENDIF 
               KK = KK + J 
            END DO 
         ELSE 
            DO J = 1, N 
               IF (X(JX)/=ZERO .OR. Y(JY)/=ZERO) THEN 
                  TEMP1 = ALPHA*Y(JY) 
                  TEMP2 = ALPHA*X(JX) 
                  IX = KX 
                  IY = KY 
                  AP(KK:J-1+KK) = AP(KK:J-1+KK) + X(IX:(J-1)*INCX+IX:INCX)*&
                     TEMP1 + Y(IY:(J-1)*INCY+IY:INCY)*TEMP2 
               ENDIF 
               JX = JX + INCX 
               JY = JY + INCY 
               KK = KK + J 
            END DO 
         ENDIF 
      ELSE 
!
!        Form  A  when lower triangle is stored in AP.
!
         IF (INCX==1 .AND. INCY==1) THEN 
            DO J = 1, N 
               IF (X(J)/=ZERO .OR. Y(J)/=ZERO) THEN 
                  TEMP1 = ALPHA*Y(J) 
                  TEMP2 = ALPHA*X(J) 
                  K = KK 
                  AP(K:N-J+K) = AP(K:N-J+K) + X(J:N)*TEMP1 + Y(J:N)*TEMP2 
               ENDIF 
               KK = KK + N - J + 1 
            END DO 
         ELSE 
            DO J = 1, N 
               IF (X(JX)/=ZERO .OR. Y(JY)/=ZERO) THEN 
                  TEMP1 = ALPHA*Y(JY) 
                  TEMP2 = ALPHA*X(JX) 
                  IX = JX 
                  IY = JY 
                  AP(KK:N-J+KK) = AP(KK:N-J+KK) + X(IX:(N-J)*INCX+IX:INCX)*&
                     TEMP1 + Y(IY:(N-J)*INCY+IY:INCY)*TEMP2 
               ENDIF 
               JX = JX + INCX 
               JY = JY + INCY 
               KK = KK + N - J + 1 
            END DO 
         ENDIF 
      ENDIF 
!
      RETURN  
!
!     End of DSPR2 .
!
      END SUBROUTINE DSPR2 
