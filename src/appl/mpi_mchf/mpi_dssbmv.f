************************************************************************
        SUBROUTINE DSSBMV(A,INDROW,INDCOL,iupper,nzer,tempb,tempc,
     :                    n,m,B,C)
************************************************************************
*
*       Computes the product of matrix A with a block of vectors B(N,M)
*                               C=A B
*       where A(NxN) is a Symmetric Sparse matrix. Only the nonzero
*       parts of the matrix are kept and this is managed by Indices
*
*       Subroutines called:
*       dgathr,dscatr,ddot,daxpy,dini
************************************************************************
*
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        LOGICAL iupper
        DIMENSION A(NZER),INDCOL(N),INDROW(NZER)
        DIMENSION TEMPB(N),TEMPC(N)
        DIMENSION B(N*M),C(N*M)
*
*     MPI stuff ***********************************************
*
	INCLUDE 'mpif.h'
	parameter (MAXPROC=100)
	common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)	
	common /PVM/ istart,ifinish
************************************************************************
*
*   on entry
*   --------
*   N		the dimension (rank) of the matrix A
*   B	    	the matrix (block of vectors) to multiply A with
*   A		Linear array keeping the the nonzero elements of
*		the matrix A. It stores columns one after the other
*		starting from the first. It is either the upper 
*		triangular part or the lower part depending on logical
*		iupper. (max elements that may contain= n^2/(2*nodes))
*   INDCOL	It is the index showing for each column where that 
*		column ends in array A.
*   INDROW	It is the index showing for each element of A to its
*		row number in the square matrix
*   iupper	logical. If .true. the upper part of A is used.
*   TEMPB,TEMPC Linear scratch arrays (dim=N)
*   
*   on exit
*   -------
*   
*   C 		the result of the multiplication (dim=NxM)
* 
************************************************************************

*
*     	calculates C=A*B matrix-vector multiplication
*
	ioffs=1
	if (iupper) ioffs=0
        istart = -1;
	mycol=0
	jstart=1
	call dinit(N*M,0.d0,C,1)

	if (istart.eq.-1) then
*	   ..interleaved
	   i1 = myid+1
	   i2 = N
	   i3 =  nprocs
	else
*	   ..block
	   i1 = istart
	   i2 = ifinish
	   i3 = 1
	endif
      	do 20 icol=i1,i2,i3
	   mycol=mycol+1
	   Num_elem = Indcol(mycol)-jstart+1
	   icur=1
	   do 10 iv=1,M
		 DL = 0.0
              if (IOFFS .eq. 1) then
*
* for the lower case, ioffs = 1, which means daxpy skips the first
* element the A(JSTART)*B(ICUR-1...) above is the dot of the first
* element done manually
* DL, upon return from dmerge has the dot of the rest of the elements.
*
*
                 DIAG=C(ICUR-1+ICOL) +
     :                A(JSTART) * B(ICUR-1+INDROW(JSTART))
                 CALL DMERGE (NUM_ELEM-1,B(ICUR), C(ICUR),
     :                INDROW(JSTART+1), A(JSTART+1),
     :                B(ICUR-1+ICOL),
     :                DL )
                 C(ICUR-1+ICOL)= DIAG + DL
              else
*
*
*     for upper case, ioffs = 0, which means daxpy skips the last element.
*     so we manually dot the last element before calling dmerge.
*
*
                 DIAG = C(ICUR+1+ICOL) +
     :                  A(JSTART-1+NUM_ELEM) *
     :                  B(ICUR-1+INDROW(JSTART-1+NUM_ELEM))

                 CALL DMERGE (NUM_ELEM-1,B(ICUR), C(ICUR),
     :                INDROW(JSTART), A(JSTART),
     :                B(ICUR-1+ICOL),
     :                DL )
                 C(ICUR-1+ICOL) = DIAG + DL
              end if
              ICUR=ICUR+N
   10      continue

	    jstart=Indcol(mycol)+1
   20   continue

*
*	Sums the contibutions on C from all nodes
*     
	icur=1
	do 30 iv=1,M
	   !call gdsummpi(C(icur),N,tempb)
           call mpi_allr_dp(C(icur),N)
	   icur=icur+N
   30	continue


	return
	end

      subroutine dmerge ( n, db, dc, idy, da, dconst, dl )
C
C  this merge version has the advantage of loading da(i)
C  and idy(i) only once.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DA(N), DB(*), DC(*), IDY(N)

      dsum = 0.0
      do 300 i = 1, n
         dsum = dsum + da(i) * db(idy(i))
         dc(idy(i)) = dc(idy(i)) + dconst * da(i)
 300  continue
      dl = dsum
      return
      end

