!=======================================================================
      SUBROUTINE DVDRVR(A, INDROW, INDCOL, IUPPER, NZER, TM, TP, N, HIEND, LIM&
         , MBLOCK, DIAG, NUME, NIV, NEIG, ISELEC, CRITE, CRITC, CRITR, ORTHO, &
         MAXITER, EIGVAL, BASIS, AB, S, TEMPS, SVEC, SCRA1, ISCRA2, INCV, ICV, &
         OLDVAL, NLOOPS, NMV, IERR) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!=======================================================================
!       called by DVDSON
!
!       Driver routine implementing Davidson's main loop. On entry it
!       is given the Basis, the work matrix D=AB and the small symmetric
!       matrix to be solved, S=B^TAB (as found by SETUP). In each step
!       the small problem is solved by calling DSPEVX.
!       TSTSEL tests for eigenvalue convergence and selects the next
!       pairs to be considered for targeting (as a block).
!       NEWVEC computes the new vectors (block) to be added in the
!       expanding basis, and tests for residual convergence.
!       ADDABS is the critical step of matrix multiplication. The new
!       vectors of D are found Dnew=ABnew, and the new small problem S,
!       is calculated. The algorithm is repeated.
!       In case of a large expanding basis (KPASS=LIM) the Basis, AB,
!       SVEC and S are collapsed.
!       At the end the current eigenvector estimates are computed as
!       well as the residuals and eigenvalue differences.
!
!       Subroutines called:
!       DSPEVX, MULTBC, TSTSEL, OVFLOW, NEWVEC, ADDABS,
!       DCOPY, DDOT, DAXPY
!-----------------------------------------------------------------------
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:22  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE tstsel_I 
      USE dcopy_I 
!      USE dspevx_I 
      USE multbc_I 
!      USE ovflow_I 
      USE newvec_I 
      USE addabs_I 
      USE daxpy_I 
      USE ddot_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NZER 
      INTEGER  :: N 
      INTEGER  :: LIM 
      INTEGER  :: MBLOCK 
      INTEGER  :: NUME 
      INTEGER , INTENT(IN) :: NIV 
      INTEGER  :: NEIG 
      INTEGER , INTENT(IN) :: MAXITER 
      INTEGER , INTENT(INOUT) :: NLOOPS 
      INTEGER , INTENT(INOUT) :: NMV 
      INTEGER , INTENT(OUT) :: IERR 
      REAL(DOUBLE)  :: CRITE 
      REAL(DOUBLE)  :: CRITC 
      REAL(DOUBLE)  :: CRITR 
      REAL(DOUBLE)  :: ORTHO 
      LOGICAL  :: IUPPER 
      LOGICAL  :: HIEND 
      INTEGER  :: INDROW(NZER) 
      INTEGER  :: INDCOL(N) 
      INTEGER  :: ISELEC(NEIG) 
      INTEGER  :: ISCRA2(5*LIM) 
      INTEGER  :: INCV(LIM) 
      INTEGER  :: ICV(NUME) 
      REAL(DOUBLE)  :: A(NZER) 
      REAL(DOUBLE)  :: TM(N) 
      REAL(DOUBLE)  :: TP(N) 
      REAL(DOUBLE)  :: DIAG(N) 
      REAL(DOUBLE)  :: EIGVAL(LIM) 
      REAL(DOUBLE)  :: BASIS(N*LIM) 
      REAL(DOUBLE)  :: AB(N*LIM) 
      REAL(DOUBLE)  :: S(LIM*(LIM + 1)/2) 
      REAL(DOUBLE)  :: TEMPS(LIM*(LIM + 1)/2) 
      REAL(DOUBLE)  :: SVEC(LIM*NUME) 
      REAL(DOUBLE)  :: SCRA1(8*LIM) 
      REAL(DOUBLE)  :: OLDVAL(NUME) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, KPASS, NNCV, NFOUND, INFO 
      LOGICAL :: RESTART, FIRST, DONE 
!-----------------------------------------------
 
!-----------------------------------------------------------------------
!
!   on entry
!   -------
!
!   OP          The user specified block-matrix-vector routine
!   N           The order of the matrix A
!   HIEND       Logical. True only if the highest eigenpairs are needed.
!   LIM         The limit on the size of the expanding Basis
!   MBLOCK      Number of vectors to be targeted in each iteration.
!   DIAG        Array of size N with the diagonal elements of A
!   NUME        The largest index of the eigenvalues wanted.
!   NIV         Starting dimension of expanding basis.
!   NEIG        Number of eigenvalues wanted.
!   ISELEC      Array containg the indices of those NEIG eigenpairs.
!   CRITE       Convergence thresholds for eigenvalues, coefficients
!   CRITC,CRITR and residuals.
!   BASIS       Array with the basis vectors.
!   AB          Array with the vectors D=AB
!   S           Array keeping the symmetric matrix of the small problem.
!   TEMPS       scratch array
!   SVEC        Array for holding the eigenvectors of S
!   SCRA1       Srcatch array used by DSPEVX.
!   ISCRA2      Integer Srcatch array used by DSPEVX.
!   INCV        Srcatch array used in DSPEVX. Also used in TSTSEL and
!               NEWVEC where it holds the Indices of uNConVerged pairs
!   ICV         It contains "1" to the locations of ConVerged eigenpairs
!   OLDVAL      Array keeping the previous' step eigenvalue estimates.
!
!   on exit
!   -------
!
!   EIGVAL      Array containing the NUME lowest eigenvalues of the
!               the matrix A (or -A if the highest are sought).
!   Basis       On exit Basis stores the NUME corresponding eigenvectors
!   OLDVAL      On exit it stores the final differences of eigenvalues.
!   SCRA1       On exit it stores the NUME corresponding residuals.
!   NLOOPS      Number of loops taken by the algorithm
!   NMV         Number of matrix-vector products performed.
!
!-----------------------------------------------------------------------
      EIGVAL(:NUME) = 1.D30 
      ICV(:NUME) = 0 
      FIRST = .TRUE. 
      KPASS = NIV 
      NNCV = KPASS 
 
   10 CONTINUE 
      CALL DCOPY (NUME, EIGVAL, 1, OLDVAL, 1) 
      CALL DCOPY ((KPASS*(KPASS + 1))/2, S, 1, TEMPS, 1) 
      CALL DSPEVX ('Vectors also', 'In a range', 'Upper triangular', KPASS, &
         TEMPS, -1, -1., 1, NUME, 0.D0, NFOUND, EIGVAL, SVEC, KPASS, SCRA1, &
         ISCRA2, INCV, INFO) 
      IERR = -ABS(INFO) 
      IF (IERR == 0) THEN 
!
! TeST for convergence on the absolute difference of eigenvalues between
! successive steps. Also SELect the unconverged eigenpairs and sort them
! by the largest magnitude in the last added NNCV rows of Svec.
!
         DONE = TSTSEL(KPASS,NUME,NEIG,ISELEC,SVEC,EIGVAL,ICV,CRITE,CRITC,SCRA1&
            ,ISCRA2,OLDVAL,NNCV,INCV) 
         IF (.NOT.(DONE .OR. KPASS>=N)) THEN 
 
            IF (KPASS == LIM) THEN 
! Maximum size for expanding basis. Collapse basis, D, and S, Svec
! Consider the basis vectors found in TSTSEL for the newvec.
!
               CALL MULTBC (N, LIM, NUME, SVEC, SCRA1, BASIS) 
               CALL MULTBC (N, LIM, NUME, SVEC, SCRA1, AB) 
               CALL OVFLOW (NUME, LIM, S, SVEC, EIGVAL) 
               KPASS = NUME 
            ENDIF 
!
! Compute and add the new vectors. NNCV is set to the number of new
! vectors that have not converged. If none, DONE=true, exit.
!
            CALL NEWVEC (N, NUME, LIM, MBLOCK, KPASS, CRITR, ORTHO, NNCV, INCV&
               , DIAG, SVEC, EIGVAL, AB, BASIS, ICV, RESTART, DONE) 
 
!          ..An infinite loop is avoided since after a collapsing Svec=I
!          ..=> Res=Di-lBi which is just computed and it is orthogonal.
!          ..The following is to prevent an improbable infinite loop.
            IF (.NOT.RESTART) THEN 
               FIRST = .TRUE. 
            ELSE IF (FIRST) THEN 
               FIRST = .FALSE. 
               CALL MULTBC (N, KPASS + NNCV, NUME, SVEC, SCRA1, BASIS) 
               CALL MULTBC (N, KPASS + NNCV, NUME, SVEC, SCRA1, AB) 
               CALL OVFLOW (NUME, KPASS + NNCV, S, SVEC, EIGVAL) 
               KPASS = NUME 
               GO TO 10 
            ELSE 
               IERR = IERR + 1024 
               GO TO 30 
            ENDIF 
 
            IF (.NOT.DONE) THEN 
!
! Add new columns in D and S, from the NNCV new vectors.
!
               CALL ADDABS (A, INDROW, INDCOL, IUPPER, NZER, TM, TP, N, LIM, &
                  HIEND, KPASS, NNCV, BASIS, AB, S) 
 
               NMV = NMV + NNCV 
               KPASS = KPASS + NNCV 
               NLOOPS = NLOOPS + 1 
 
               IF (NLOOPS <= MAXITER) GO TO 10 
               IERR = IERR + 2048 
               NLOOPS = NLOOPS - 1 
               KPASS = KPASS - NNCV 
            ENDIF 
         ENDIF 
   30    CONTINUE 
         DO I = 1, NUME 
            OLDVAL(I) = ABS(OLDVAL(I)-EIGVAL(I)) 
         END DO 
 
         CALL MULTBC (N, KPASS, NUME, SVEC, SCRA1, BASIS) 
         CALL MULTBC (N, KPASS, NUME, SVEC, SCRA1, AB) 
!
! i=1,NUME residual(i)= DCi-liBCi= newDi-linewBi
! temporarily stored in AB(NUME*N+1)
!
         DO I = 1, NUME 
            CALL DCOPY (N, AB((I-1)*N+1), 1, AB(NUME*N+1), 1) 
            CALL DAXPY (N, (-EIGVAL(I)),BASIS((I-1)*N+1), 1, AB(NUME*N+1), 1) 
            SCRA1(I) = DDOT(N,AB(NUME*N+1),1,AB(NUME*N+1),1) 
            SCRA1(I) = SQRT(SCRA1(I)) 
         END DO 
      ENDIF 
 
      RETURN  
      END SUBROUTINE DVDRVR 
