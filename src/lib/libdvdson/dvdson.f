*=======================================================================
        SUBROUTINE DVDSON(A,INDROW,INDCOL,IUPPER,NZER,TM,TP,
     :                   N,LIM,DIAG,
     :                   ILOW,IHIGH,ISELEC,NIV,MBLOCK,
     :                   CRITE,CRITC,CRITR,ORTHO,MAXITER,
     :                   WORK,IWRSZ,IWORK,IIWSZ,
     :                   HIEND,NLOOPS,NMV,IERR)
*=======================================================================
*
*       Author: Andreas Stathopoulos, Charlotte F. Fischer
*
*       Computer Science Department
*       Vanderbilt University
*       Nashville, TN 37212
*       andreas@vuse.vanderbilt.edu
*       cff@vuse.vanderbilt.edu                       DECEMBER 1993
*
*       Copyright (c) by Andreas Stathopoulos and Charlotte F. Fischer
*
*       DVDSON is a Fortran77 program that finds a few selected
*       eigenvalues and their eigenvectors at either end of spectrum of
*       a large, symmetric (and usually sparse) matrix, denoted as A.
*       The matrix A is only referenced indirectly through the user
*       supplied routine OP which implements a block matrix-vector
*       operation(see below). Either the range of the eigenvalues wanted
*       or an array of the indices of selected ones can be specified.
*       DVDSON is a front-end routine for setting up arrays, and initial
*       guess (calling SETUP). It also performs detailed error checking.
*       DVDRVR is the driver routine that implements a version of the
*       Davidson algorithm. The characteristics of this version are:
*        o  All arrays used by the program are stored in MEMORY.
*        o  BLOCK method (many vectors may be targeted per iteration.)
*        o  Eigenvectors are targeted in an optimum way without
*           the need to compute all unconverged residuals,
*        o  It REORTHOGONILIZES the basis in case of orthogonality loss.
*        o  Finds HIGHEST eigenpairs by using the negative of the A.
*        o  Finds SELECTED eigenpairs specified by the user.
*        o  It accepts INITIAL eigenvector ESTIMATES or it can
*           CREATE INITIAL ESTIMATES from the diagonal elements.
*        o  It uses a USER SUPPLIED block matrix-vector operation, OP.
*           Depending on the implementation, OP can operate in either
*           memory or on disc, and for either sparse or dense matrix.
*        o  The user can provide STOPPING CRITERIA for eigenvalues,
*           and residuals. The user can also CONTROL reorthogonalization
*            and block size.
*        o  On exit INFORMATION is given about the convergence status
*           of eigenpairs and the number of loops and OP operations.
*
*       The program consists of the following routines:
*       DVDSON, SETUP, DVDRVR, ADDABS, TSTSEL,
*       MULTBC, OVFLOW, NEWVEC, ORTHNRM.

*       It also calls some basic BLAS routines:
*       DCOPY, DSCAL, DDOT, DAXPY, IDAMAX, DGEMV, DINIT

*       For solving the small eigenproblem, the routine DSPEVX from
*       LAPACK is used. DSPEVX is obtainable from NETLIB, together
*       with a series of subroutines that it calls.
*
*       All the routines have IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
*-----------------------------------------------------------------------
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        LOGICAL IUPPER
        DIMENSION A(NZER),INDROW(NZER),INDCOL(N)
        DIMENSION TM(N),TP(N)
        DIMENSION DIAG(N),WORK(IWRSZ),IWORK(IIWSZ)
        DIMENSION ISELEC(LIM)
        LOGICAL HIEND
*-----------------------------------------------------------------------
*  (Important to the following is the concept of NUME, the distance of
*   the index of the eigenpair wanted which is farthest from the
*   extremes,i.e.,
*      if  lowest  eigepairs i1<i2<...<ik are wanted, NUME=ik
*      if highest eigenpairs i1<i2<...<ik are wanted, NUME=N-i1+1
*   where i1,...,ik are the indices of the wanted eigenpairs.
*   Obviously, NUME.GE.(No. of EiGenpairs wanted). )

*   on entry
*   -------
*   OP          User supplied routine with calling sequence OP(N,M,B,C).
*               B and C are N x M matrices and C stores the result AxB.
*               It should be declared external in the main program.
*   N           Order of the matrix.
*   LIM         The upper limit on the dimension of the expanding basis.
*               NUME.LT.LIM.LE.N must hold. The case LIM=NUME is allowed
*               only for LIM=NUME=N. The choice of LIM depends on the
*               available workspace (see below). If the space is
*               available it is preferable to have a large LIM, but not
*               larger than NUME$+$40.
*   DIAG        Array of size N with the diagonal elements of the
*               matrix A.
*   ILOW        The index of the lowest eigepair to be computed. If
*               (ILOW.LE.0).or.(ILOW.GT.N), the selected eigenpairs
*               to be computed should be contained in array ISELEC.
*               (Modified on exit).
*   IHIGH       The index of the highest eigenpair to be computed.
*               Considered ONLY when ILOW is in the range
*               (0.LT.ILOW.LE.N). (Modified on exit).
*   ISELEC      Array of size LIM holding the user specified indices
*               for the eigenpairs to be computed. Considered only when
*               (ILOW.LE.0).or.(ILOW.GT.N). The indices are read from
*               the first position until a non positive integer is met.
*                  Example: if N=500, ILOW=0, and ISELEC(1)=495,
*                  ISELEC(2)=497, ISELEC(3)=-1, the program will find
*                  2 of the highest eigenpairs, pairs 495 and 497.
*               Any order of indices is acceptable (Modified on exit).
*   NIV         Number of Initial Vector estimates provided by the user.
*               If NIV is in the range:  (NUME).LE.(NIV).LE.(LIM),
*               the first NIV columns of size N of WORK should contain
*               the estimates (see below). In all other cases of NIV,
*               the program generates initial estimates.
*   MBLOCK      Number of vectors to be targeted in each iteration.
*               1.LE.MBLOCK.LE.(No. EiGenpairs wanted) should hold.
*               Large block size reduces the number of iterations
*               (matrix acceses) but increases the matrix-vector
*               multiplies. It should be used when the matrix accese
*               is expensive (disc, recomputed or distributed).
*   CRITE       Convergence threshold for eigenvalues.
*               If ABS(EIGVAL-VALOLD) is less than CRITE for all wanted
*               eigenvalues, convergence is signaled.
*   CRITC       Convergence threshold for the coefficients of the last
*               added basis vector(s). If all of those corresponding to
*               unconverged eigenpairs are less than CRITC convergence
*               is signaled.
*   CRITR       Convergence threshold for residual vector norms. If
*               all the residual norms ||Ax_i-l_ix_i|| of the targeted
*               x_i are less than CRITR convergence is signaled.
*               If ANY of the criteria are satisfied the algorithm stops
*   ORTHO       The threshold over which loss of orthogonality is
*               assumed. Usually ORTHO.LE.CRITR*10 but the process can
*               be skipped by setting ORTHO to a large number(eg,1.D+3).
*   MAXITER     Upper bound on the number of iterations of the
*               algorithm. When MAXITER is exceeded the algorithm stops.
*               A typical MAXITER can be MAX(200,NUME*40), but it can
*               be increased as needed.
*   WORK        Real array of size IWRSZ. Used for both input and output
*               If NIV is in ((NUME).LE.(NIV).LE.(LIM)), on input, WORK
*               must have the NIV initial estimates. These NIV N-element
*               vectors start from WORK(1) and continue one after the
*               other. They must form an orthonormal basis.
*   IWRSZ       The size of the real workspace. It must be at least as
*               large as:
*
*                       2*N*LIM + LIM*LIM + (NUME+10)*LIM + NUME
*
*   IWORK       Integer work array of size IIWSZ. Used as scrath array
*               for indices and for use in the LAPACK routines.
*   IIWSZ       The size of the integer workspace. It must be at least
*               as large as:
*                                    6*LIM + NUME
*
*               If LIM or NUME needs to be increased, the space should
*               also be increased accordingly. For given IWRSZ and
*               IIWSZ one can calculate how big a problem one can
*               solve (LIM,NUME).
*
*   on exit
*   -------
*   WORK(1)     The first NUME*N locations contain the approximations to
*               the NUME extreme eigenvectors. If the lowest eigenpairs
*               are required, (HIEND=false), eigenvectors appear in
*               ascending order, otherwise (HIEND=false), they appear in
*               descending order. If only some are requested, the order
*               is the above one for all the NUME extreme eigenvectors,
*               but convergence has been reached only for the selected
*               ones. The rest are the current approximations to the
*               non-selected eigenvectors.
*   WORK(NUME*N+1)
*               The next NUME locations contain the approximations to
*               the NUME extreme eigenvalues, corresponding to the above
*               NUME eigenvectors. The same ordering and convergence
*               status applies here as well.
*   WORK(NUME*N+NUME+1)
*               The next NUME locations contain the corresponding values
*               of ABS(EIGVAL-VALOLD) of the NUME above eigenvalues, of
*               the last step of the algorithm.
*   WORK(NUME*N+NUME+NUME+1)
*               The next NUME locations contain the corresponding
*               residual norms of the NUME above eigenvectors, of the
*               last step.
*   HIEND       Logical. If .true. on exit the highest eigenpairs are
*               found in descending order. Otherwise, the lowest
*               eigenpairs are arranged in ascending order.
*   NLOOPS      The number of iterations it took to reach convergence.
*               This is also the number of matrix references.
*   NMV         The number of Matrix-vector(M-V) multiplies. Each matrix
*               reference can have up to size(block) M-V multiplies.
*   IERR        An integer denoting the completions status:
*               IERR = 0        denotes normal completion.
*               IERR = -k       denotes error in DSPEVX (k eigenpairs
*                               not converged)
*               0<IERR<=2048    denotes some inconsistency as follows:
*        If (INT( MOD(IERR,  2)/1  ) N < LIM
*        If (INT( MOD(IERR,  4)/2  ) LIM < 1
*        If (INT( MOD(IERR,  8)/4  ) ISELEC(1)<1, and no range specified
*        If (INT( MOD(IERR, 16)/8  ) IHIGH > N (in range or ISELEC)
*        If (INT( MOD(IERR, 32)/16 ) IHIGH < ILOW (Invalid range)
*        If (INT( MOD(IERR, 64)/32 ) NEIG >= LIM (Too many wanted)
*        If (INT( MOD(IERR,128)/64 ) Probable duplication in ISELEC
*        If (INT( MOD(IERR,256)/128) NUME >= LIM (max eigen very far)
*        If (INT( MOD(IERR,512)/256) MBLOCK is out of bounds
*        If (INT( MOD(IERR,1024)/512) IWRSZ or IIWSZ is not enough
*        If (INT( MOD(IERR,2048)/1024) Orthogonalization Failed
*        If (INT( MOD(IERR,4096)/2048) NLOOPS > MAXITER
*
*               The program will also print an informative message to
*               the standard output when NIV is not proper but it will
*               continue by picking initial estimates internally.
*-----------------------------------------------------------------------
*
* Checking user input errors, and setting up the problem to solve.
*
        IERR=0
        IF (LIM.GT.N) IERR=IERR+1
        IF (LIM.LE.0) IERR=IERR+2

        HIEND=.false.

        IF ((ILOW.LE.0).OR.(ILOW.GT.N)) THEN
*          ..Look for user choice of eigenpairs in ISELEC
           IF (ISELEC(1).LE.0) THEN
*             ..Nothing is given in ISELEC
              IERR=IERR+4
           ELSE
*             ..Find number of eigenpairs wanted, and their
*             ..min/max indices
              NEIG=1
              ILOW=ISELEC(1)
              IHIGH=ISELEC(1)
              DO 10 I=2,LIM
                 IF (ISELEC(I).LE.0) GOTO 20
                 ILOW=MIN(ILOW,ISELEC(I))
                 IHIGH=MAX(IHIGH,ISELEC(I))
                 NEIG=NEIG+1
 10           CONTINUE
*             ..Check if a very large index is asked for
 20           IF (IHIGH.GT.N) IERR=IERR+8
           ENDIF
        ELSE
*          ..Look for a range between ILOW and IHIGH
*          ..Invalid range. IHIGH>N
           IF (IHIGH.GT.N) IERR=IERR+8
           NEIG=IHIGH-ILOW+1
*          ..Invalid range. IHIGH<ILOW
           IF (NEIG.LE.0) IERR=IERR+16
           IF (NEIG.GT.LIM) THEN
*             ..Not enough Basis space. Increase LIM or decrease NEIG
              IERR=IERR+32
           ELSE
*             ..Fill in the ISELEC with the required indices
              DO 40 I=1,NEIG
 40              ISELEC(I)=ILOW+I-1
           ENDIF
        ENDIF

        IF (IERR.NE.0) RETURN

        NUME=IHIGH
*       ..Identify if few of the highest eigenpairs are wanted.
        IF ((ILOW+IHIGH-1).GT.N) THEN
           HIEND=.true.
           NUME=N-ILOW+1
*          ..Change the problem to a minimum eipenpairs one
*          ..by picking the corresponding eigenpairs on the
*          ..opposite side of the spectrum.
           DO 50 I=1,NEIG
 50           ISELEC(I)=N-ISELEC(I)+1
        ENDIF
*       ..duplications in ISELEC
        IF (NEIG.GT.NUME) IERR=IERR+64
*       ..Not enough Basis space. Increase LIM or decrease NUME
        IF ((NUME.GT.LIM).OR.((NUME.EQ.LIM).AND.(NUME.NE.N)))
     :     IERR=IERR+128
*       ..Size of Block out of bounds
        IF ( (MBLOCK.LT.1).OR.(MBLOCK.GT.NEIG) ) IERR=IERR+256

*       ..Check for enough workspace for Dvdson
        IF ((IWRSZ.LT.(LIM*(2*N+LIM+(NUME+10))+NUME)).OR.
     :      (IIWSZ.LT.(6*LIM+NUME))) IERR=IERR+512

        IF (IERR.NE.0) RETURN

        IF (NIV.GT.LIM) THEN
*          ..Check number of initial estimates NIV is lower than LIM.
           PRINT*,'WARNING: Too many initial estimates.?'
           PRINT*,'The routine will pick the appropriate number'
        ELSEIF ((NIV.LT.NUME).AND.(NIV.GT.0)) THEN
*          ..check if enough initial estimates.
*          ..(NIV<1 => program chooses)
           PRINT*,'WARNING: Not enough initial estimates'
           PRINT*,'The routine will pick the appropriate number'
        ENDIF
ctc
*        print *, 'l295 dvdson'
ctc


*
* Assigning space for the real work arrays
*
        iBasis    =1
        ieigval   =iBasis  +N*LIM
        iAB       =ieigval +LIM
        iS        =iAB     +N*LIM
        itempS    =iS      +LIM*(LIM+1)/2
        iSvec     =itempS  +LIM*(LIM+1)/2
        iscra1    =iSvec   +LIM*NUME
        ioldval   =iscra1  +8*LIM
*
* Assigning space for the integer work arrays
*
        iscra2    =1
        iscra3    =iscra2  +5*LIM
        iIcv      =iscra3  +LIM
ctc
*        print *, 'l317 dvdson'
ctc
        IF (HIEND) CALL DSCAL(N,-1.D0,DIAG,1)
ctc
*        print *, 'diag :',(diag(i), i=1,n) 
ctc

        iSTART=NIV
        CALL SETUP(A,INDROW,INDCOL,IUPPER,NZER,TM,TP,
     :             N,LIM,NUME,HIEND,DIAG,WORK(iscra1),
     :             WORK(iBasis),WORK(iAB),WORK(iS),iSTART)

ctc
*       do i=1,n*lim
*        if (abs(work(i)) .gt. 1.0e-9) then
*          print *, 'Basis :',i,work(i)
*        end if
*       end do
ctc

        NLOOPS=1
        NMV=ISTART
ctc
*        print *, 'l328 dvdson'
ctc
        CALL DVDRVR(A,INDROW,INDCOL,IUPPER,NZER,TM,TP,
     :             N,HIEND,LIM,MBLOCK,DIAG,
     :             NUME,iSTART,NEIG,ISELEC,
     :             CRITE,CRITC,CRITR,ORTHO,MAXITER,
     :             WORK(ieigval),WORK(iBasis),WORK(iAB),
     :             WORK(iS),WORK(itempS),WORK(iSvec),
     :             WORK(iscra1),IWORK(iscra2),IWORK(iscra3),
     :             IWORK(iIcv),WORK(ioldval),
     :             NLOOPS,NMV,IERR)

ctc
*       do i=1,n*lim,1000
*        if (abs(work(iBasis+i)) .gt. 1.0e-9) then
*          print *, 'Eigv :',i,work(iBasis+i)
*        end if
*       end do
*       do i=1,n*lim
*        if (abs(work(i)) .gt. 1.0e-9) then
*          print *, 'Basis :',i,work(i)
*        end if
*       end do
ctc
ctc
*        print *, 'l340 dvdson'
ctc
        IF (HIEND) THEN
           CALL DSCAL(N,-1.D0,DIAG,1)
           CALL DSCAL(NUME,-1.D0,WORK(ieigval),1)
        endif
*
* -Copy the eigenvalues after the eigenvectors
* -Next, copy the difference of eigenvalues between the last two steps
* -Next, copy the residuals for the first NUME estimates
*
ctc
*        print *, 'l352 dvdson'
ctc
        CALL DCOPY(NUME,WORK(ieigval),1,WORK(iBasis+N*NUME),1)
        CALL DCOPY(NUME,WORK(ioldval),1,WORK(iBasis+(N+1)*NUME),1)
        CALL DCOPY(NUME,WORK(iscra1),1,WORK(iBasis+(N+2)*NUME),1)

 100    RETURN
        END
