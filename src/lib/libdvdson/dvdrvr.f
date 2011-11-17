*=======================================================================
        SUBROUTINE DVDRVR(A,INDROW,INDCOL,IUPPER,NZER,TM,TP,
     :                    N,HIEND,LIM,MBLOCK,DIAG,
     :                    NUME,NIV,NEIG,ISELEC,
     :                    CRITE,CRITC,CRITR,ORTHO,MAXITER,
     :                    EIGVAL,BASIS,AB,S,TEMPS,SVEC,
     :                    SCRA1,ISCRA2,INCV,ICV,OLDVAL,
     :                    NLOOPS,NMV,IERR)
*=======================================================================
*       called by DVDSON
*
*       Driver routine implementing Davidson's main loop. On entry it
*       is given the Basis, the work matrix D=AB and the small symmetric
*       matrix to be solved, S=B^TAB (as found by SETUP). In each step
*       the small problem is solved by calling DSPEVX.
*       TSTSEL tests for eigenvalue convergence and selects the next
*       pairs to be considered for targeting (as a block).
*       NEWVEC computes the new vectors (block) to be added in the
*       expanding basis, and tests for residual convergence.
*       ADDABS is the critical step of matrix multiplication. The new
*       vectors of D are found Dnew=ABnew, and the new small problem S,
*       is calculated. The algorithm is repeated.
*       In case of a large expanding basis (KPASS=LIM) the Basis, AB,
*       SVEC and S are collapsed.
*       At the end the current eigenvector estimates are computed as
*       well as the residuals and eigenvalue differences.
*
*       Subroutines called:
*       DSPEVX, MULTBC, TSTSEL, OVFLOW, NEWVEC, ADDABS,
*       DCOPY, DDOT, DAXPY
*-----------------------------------------------------------------------

        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        LOGICAL IUPPER
        DIMENSION A(NZER),INDROW(NZER),INDCOL(N)
        DIMENSION TM(N),TP(N)
        DIMENSION DIAG(N)
        DIMENSION S(LIM*(LIM+1)/2),TEMPS(LIM*(LIM+1)/2)
        DIMENSION SVEC(LIM*NUME),EIGVAL(LIM)
        DIMENSION ISELEC(NEIG)
        DIMENSION BASIS(N*LIM),AB(N*LIM)
        DIMENSION SCRA1(8*LIM),ISCRA2(5*LIM),INCV(LIM)
        DIMENSION ICV(NUME),OLDVAL(NUME)
        LOGICAL RESTART,FIRST,DONE,HIEND,TSTSEL

*-----------------------------------------------------------------------
*
*   on entry
*   -------
*
*   OP          The user specified block-matrix-vector routine
*   N           The order of the matrix A
*   HIEND       Logical. True only if the highest eigenpairs are needed.
*   LIM         The limit on the size of the expanding Basis
*   MBLOCK      Number of vectors to be targeted in each iteration.
*   DIAG        Array of size N with the diagonal elements of A
*   NUME        The largest index of the eigenvalues wanted.
*   NIV         Starting dimension of expanding basis.
*   NEIG        Number of eigenvalues wanted.
*   ISELEC      Array containg the indices of those NEIG eigenpairs.
*   CRITE       Convergence thresholds for eigenvalues, coefficients
*   CRITC,CRITR and residuals.
*   BASIS       Array with the basis vectors.
*   AB          Array with the vectors D=AB
*   S           Array keeping the symmetric matrix of the small problem.
*   TEMPS       scratch array
*   SVEC        Array for holding the eigenvectors of S
*   SCRA1       Srcatch array used by DSPEVX.
*   ISCRA2      Integer Srcatch array used by DSPEVX.
*   INCV        Srcatch array used in DSPEVX. Also used in TSTSEL and
*               NEWVEC where it holds the Indices of uNConVerged pairs
*   ICV         It contains "1" to the locations of ConVerged eigenpairs
*   OLDVAL      Array keeping the previous' step eigenvalue estimates.
*
*   on exit
*   -------
*
*   EIGVAL      Array containing the NUME lowest eigenvalues of the
*               the matrix A (or -A if the highest are sought).
*   Basis       On exit Basis stores the NUME corresponding eigenvectors
*   OLDVAL      On exit it stores the final differences of eigenvalues.
*   SCRA1       On exit it stores the NUME corresponding residuals.
*   NLOOPS      Number of loops taken by the algorithm
*   NMV         Number of matrix-vector products performed.
*
*-----------------------------------------------------------------------
        DO 5 I=1,NUME
           EIGVAL(I)=1.D30
  5        ICV(I)=0
        FIRST =.true.
        KPASS =NIV
        NNCV  =KPASS

10      CONTINUE
*       (iterations for kpass=NUME,LIM)
*
* Diagonalize the matrix S. Find only the NUME smallest eigenpairs
*
           CALL DCOPY(NUME,EIGVAL,1,OLDVAL,1)
           CALL DCOPY((KPASS*(KPASS+1))/2,S,1,TEMPS,1)
           CALL DSPEVX('Vectors also','In a range','Upper triangular',
     :          KPASS,TEMPS,-1.,-1.,1,NUME,0.D0,
     :          NFOUND,EIGVAL,SVEC,KPASS,SCRA1,ISCRA2,INCV,INFO)
           IERR=-ABS(INFO)
           IF (IERR.NE.0) GOTO 60
*
* TeST for convergence on the absolute difference of eigenvalues between
* successive steps. Also SELect the unconverged eigenpairs and sort them
* by the largest magnitude in the last added NNCV rows of Svec.
*
           DONE=TSTSEL(KPASS,NUME,NEIG,ISELEC,SVEC,EIGVAL,ICV,
     :            CRITE,CRITC,SCRA1,ISCRA2,OLDVAL,NNCV,INCV)
           IF ((DONE).OR.(KPASS.GE.N)) GOTO 30

           IF (KPASS.EQ.LIM) THEN
* Maximum size for expanding basis. Collapse basis, D, and S, Svec
* Consider the basis vectors found in TSTSEL for the newvec.
*
              CALL MULTBC(N,LIM,NUME,SVEC,SCRA1,BASIS)
              CALL MULTBC(N,LIM,NUME,SVEC,SCRA1,AB)
              CALL OVFLOW(NUME,LIM,S,SVEC,EIGVAL)
              KPASS=NUME
           ENDIF
*
* Compute and add the new vectors. NNCV is set to the number of new
* vectors that have not converged. If none, DONE=true, exit.
*
           CALL NEWVEC(N,NUME,LIM,MBLOCK,KPASS,CRITR,ORTHO,NNCV,INCV,
     :                 DIAG,SVEC,EIGVAL,AB,BASIS,ICV,RESTART,DONE)

*          ..An infinite loop is avoided since after a collapsing Svec=I
*          ..=> Res=Di-lBi which is just computed and it is orthogonal.
*          ..The following is to prevent an improbable infinite loop.
           IF (.NOT.RESTART) THEN
              FIRST=.true.
           ELSEIF (FIRST) THEN
              FIRST=.false.
              CALL MULTBC(N,KPASS+NNCV,NUME,SVEC,SCRA1,BASIS)
              CALL MULTBC(N,KPASS+NNCV,NUME,SVEC,SCRA1,AB)
              CALL OVFLOW(NUME,KPASS+NNCV,S,SVEC,EIGVAL)
              KPASS=NUME
              GOTO 10
           ELSE
              IERR=IERR+1024
              GOTO 30
           ENDIF

           IF (DONE) GOTO 30
*
* Add new columns in D and S, from the NNCV new vectors.
*
           CALL ADDABS(A,INDROW,INDCOL,IUPPER,NZER,TM,TP,
     :                 N,LIM,HIEND,KPASS,NNCV,BASIS,AB,S)

           NMV=NMV+NNCV
           KPASS=KPASS+NNCV
           NLOOPS=NLOOPS+1

        IF (NLOOPS.LE.MAXITER) GOTO 10
        IERR=IERR+2048
        NLOOPS=NLOOPS-1
        KPASS=KPASS-NNCV
 30     CONTINUE
*
* Calculate final results. EIGVAL contains the eigenvalues, BASIS the
* eigenvectors, OLDVAL the eigenvalue differences, and SCRA1 residuals.
*
        DO 40 I=1,NUME
 40        OLDVAL(I)=ABS(OLDVAL(I)-EIGVAL(I))

        CALL MULTBC(N,KPASS,NUME,SVEC,SCRA1,BASIS)
        CALL MULTBC(N,KPASS,NUME,SVEC,SCRA1,AB)
*
* i=1,NUME residual(i)= DCi-liBCi= newDi-linewBi
* temporarily stored in AB(NUME*N+1)
*
        DO 50 I=1,NUME
           CALL DCOPY(N,AB((I-1)*N+1),1,AB(NUME*N+1),1)
           CALL DAXPY(N,-EIGVAL(I),BASIS((I-1)*N+1),1,AB(NUME*N+1),1)
           SCRA1(I)=DDOT(N,AB(NUME*N+1),1,AB(NUME*N+1),1)
           SCRA1(I)=SQRT(SCRA1(I))
 50     CONTINUE

 60     RETURN
        END
