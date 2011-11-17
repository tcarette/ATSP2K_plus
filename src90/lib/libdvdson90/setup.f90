!=======================================================================
      SUBROUTINE SETUP(A, INDROW, INDCOL, IUPPER, NZER, TM, TP, N, LIM, NUME, &
         HIEND, DIAG, MINELEM, BASIS, AB, S, NIV) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!=======================================================================
!       Subroutine for setting up (i) the initial BASIS if not provided,
!       (ii) the product of the matrix A with the Basis into matrix AB,
!       and (iii) the small matrix S=B^TAB. If no initial estimates are
!       available, the BASIS =(e_i1,e_i2,...,e_iNUME), where i1,i2,...,
!       iNUME are the indices of the NUME lowest diagonal elements, and
!       e_i the i-th unit vector. (ii) and (iii) are handled by ADDABS.
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  18:49:20  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE addabs_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NZER 
      INTEGER  :: N 
      INTEGER  :: LIM 
      INTEGER , INTENT(IN) :: NUME 
      INTEGER  :: NIV 
      LOGICAL  :: IUPPER 
      LOGICAL  :: HIEND 
      INTEGER  :: INDROW(NZER) 
      INTEGER  :: INDCOL(N) 
      REAL(DOUBLE) , INTENT(INOUT) :: MINELEM(LIM) 
      REAL(DOUBLE)  :: A(NZER) 
      REAL(DOUBLE)  :: TM(N) 
      REAL(DOUBLE)  :: TP(N) 
      REAL(DOUBLE) , INTENT(IN) :: DIAG(N) 
      REAL(DOUBLE)  :: BASIS(N*LIM) 
      REAL(DOUBLE)  :: AB(N*LIM) 
      REAL(DOUBLE)  :: S(LIM*(LIM + 1)/2) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, IMIN, KPASS 
!-----------------------------------------------
!-----------------------------------------------------------------------
!   on entry
!   --------
!   OP          The block matrix-vector operation, passed to ADDABS
!   N           the order of the matrix A
!   LIM         The limit on the size of the expanding Basis
!   NUME        Largest index of the wanted eigenvalues.
!   HIEND       Logical. True only if the highest eigenpairs are needed.
!   DIAG        Array of size N with the diagonal elements of A
!   MINELEM     Array keeping the indices of the NUME lowest diagonals.
!
!   on exit
!   -------
!   BASIS       The starting basis.
!   AB, S       The starting D=AB, and small matrix S=B^TAB
!   NIV         The starting dimension of BASIS.
!-----------------------------------------------------------------------
 
      IF (NIV>LIM .OR. NIV<NUME) THEN 
!
!          ..Initial estimates are not available. Give as estimates unit
!          ..vectors corresponding to the NUME minimum diagonal elements
!          ..First find the indices of these NUME elements (in MINELEM).
!          ..Array AB is used temporarily as a scratch array.
!
         CALL DINIT (N, -1.D0, AB, 1) 
         DO I = 1, NUME 
!             ..imin= the first not gotten elem( NUME<=N )
            DO J = 1, N 
               IF (AB(J) >= 0) CYCLE  
               EXIT  
            END DO 
            IMIN = J 
            DO J = IMIN + 1, N 
               IF (AB(J)>=0 .OR. DIAG(J)>=DIAG(IMIN)) CYCLE  
               IMIN = J 
            END DO 
            MINELEM(I) = IMIN 
            AB(IMIN) = 1.D0 
         END DO 
!
!          ..Build the Basis. B_i=e_(MINELEM(i))
!
         CALL DINIT (N*LIM, 0.D0, BASIS, 1) 
         DO J = 1, NUME 
            I = (J - 1)*N + MINELEM(J) 
            BASIS(I) = 1 
         END DO 
 
         NIV = NUME 
      ENDIF 
!
! Find the matrix AB by matrix-vector multiplies, as well as the
! small matrix S = B^TAB.
!
      KPASS = 0 
      CALL ADDABS (A, INDROW, INDCOL, IUPPER, NZER, TM, TP, N, LIM, HIEND, &
         KPASS, NIV, BASIS, AB, S) 
 
      RETURN  
      END SUBROUTINE SETUP 
