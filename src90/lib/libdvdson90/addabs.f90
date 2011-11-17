!=======================================================================
      SUBROUTINE ADDABS(A, INDROW, INDCOL, IUPPER, NZER, TM, TP, N, LIM, HIEND&
         , KPASS, NNCV, BASIS, AB, S) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!=======================================================================
!       Called by: DVDRVR, SETUP
!
!       Calculates the new column in the D matrix and the new column
!       in the S matrix. The new D column is D(new)=AB(new). S has a
!       new row and column, but being symmetric only the new column is
!       stored. S(i,kpass+1)=B(i)^T D(kpass+1) for all i.
!
!       subroutines called:
!       OP, DDOT, DSCAL
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:22  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dssbmv_I 
      USE dscal_I 
      USE ddot_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NZER 
      INTEGER  :: N 
      INTEGER , INTENT(IN) :: LIM 
      INTEGER , INTENT(IN) :: KPASS 
      INTEGER  :: NNCV 
      LOGICAL  :: IUPPER 
      LOGICAL , INTENT(IN) :: HIEND 
      INTEGER  :: INDROW(NZER) 
      INTEGER  :: INDCOL(N) 
      REAL(DOUBLE)  :: A(NZER) 
      REAL(DOUBLE)  :: TM(N) 
      REAL(DOUBLE)  :: TP(N) 
      REAL(DOUBLE)  :: BASIS(N*LIM) 
      REAL(DOUBLE)  :: AB(N*LIM) 
      REAL(DOUBLE) , INTENT(OUT) :: S(LIM*(LIM + 1)/2) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IDSTART, ISSTART, IV, IBSTART, IBV 
      REAL(DOUBLE) :: SS 
!-----------------------------------------------
!        EXTERNAL OP
!-----------------------------------------------------------------------
!   on entry
!   -------
!   N           The order of the matrix A
!   kpass       The current dimension of the expanding sub-basis
!   NNCV        Number of new basis vectors.
!   Basis       the basis vectors, including the new NNCV ones.
!   on exit
!   -------
!   AB          The new matrix D=AB. (with new NNCV columns)
!   S           The small matrix with NNCV new columns at the last part
!-----------------------------------------------------------------------
!
! The user specified matrix-vector routine is called with the new
! basis vector B(*,kpass+1) and the result is assigned to AB(idstart)
!
      IDSTART = KPASS*N + 1 
      CALL DSSBMV (A, INDROW, INDCOL, IUPPER, NZER, TM, TP, N, NNCV, BASIS(&
         IDSTART), AB(IDSTART)) 
!
! If highest pairs are sought, use the negative of the matrix
!
      IF (HIEND) CALL DSCAL (N*NNCV, -1.D0, AB(IDSTART), 1) 
!
! The new S is calculated by adding the new last columns
! S(new)=B^T D(new).
!
      ISSTART = KPASS*(KPASS + 1)/2 
      DO IV = 1, NNCV 
         IBSTART = 1 
         DO IBV = 1, KPASS + IV 
            SS = DDOT(N,BASIS(IBSTART),1,AB(IDSTART),1) 
            S(ISSTART+IBV) = SS 
            IBSTART = IBSTART + N 
         END DO 
         ISSTART = ISSTART + KPASS + IV 
         IDSTART = IDSTART + N 
      END DO 
 
      RETURN  
      END SUBROUTINE ADDABS 
