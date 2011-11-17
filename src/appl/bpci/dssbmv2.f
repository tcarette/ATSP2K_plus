************************************************************************
        SUBROUTINE DSSBMV(A,INDROW,INDCOL,IUPPER,NZER,TEMPB,TEMPC,
     :			  N,M,B,C)
************************************************************************
*
*	Computes the product of matrix A with a block of vectors B(N,M)
*			        C=A B
*	where A(NxN) is a Symmetric Sparse matrix. Only the nonzero 
*	parts of the matrix are kept and this is managed by Indices.
*       The matrix may be either in memory or on disk. If
*       idisk =0, the matrix must be in memory; otherwise on disk
*
*	Subroutines called:
*	dinit, dmerge 
************************************************************************
* 
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	LOGICAL IUPPER 
	DIMENSION A(NZER),INDCOL(N),INDROW(NZER) 
	DIMENSION TEMPB(N),TEMPC(N)
	DIMENSION B(N*M),C(N*M)
        COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk
*
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
	IOFFS=1
	IF (IUPPER) IOFFS=0

	ISTART=1
	CALL DINIT(N*M,0.D0,C,1)
	IF (idisk .eq. 1) then
	  rewind(iouhm)
	END IF

      	DO 20 ICOL=1,N
    	   ICUR=1
	   IF (idisk .eq. 1) then
	     Read(iouhm) jb,mr,(a(i),i=1,mr),(indrow(i),i=1,mr)
	     numelem = mr
	     if (iLS .eq. 1) then
*              .. LS matrix on disk not shifted
	       if (icol .eq. 1) shift = a(1)
	       a(1) = a(1)- shift
	     end if
	   else
	     numelem = indcol(icol) -istart + 1
	   end if

    	   DO 10 IV=1,M
              DL = 0.0
              if (IOFFS .eq. 1) then
*
* for the lower case, ioffs = 1, which means daxpy skips the first 
* element the A(ISTART)*B(ICUR-1...) above is the dot of the first 
* element done manually 
* DL, upon return from dmerge has the dot of the rest of the elements.
*
*
                 DIAG=C(ICUR-1+ICOL) + 
     :                A(ISTART) * B(ICUR-1+INDROW(ISTART)) 
                 CALL DMERGE (NUMELEM-1,B(ICUR), C(ICUR), 
     :                INDROW(ISTART+1), A(ISTART+1), 
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
     :                  A(ISTART-1+NUMELEM) *
     :                  B(ICUR-1+INDROW(ISTART-1+NUMELEM))

                 CALL DMERGE (NUMELEM-1,B(ICUR), C(ICUR), 
     :                INDROW(ISTART), A(ISTART), 
     :                B(ICUR-1+ICOL),
     :                DL ) 
                 C(ICUR-1+ICOL) = DIAG + DL
              end if 
    	      ICUR=ICUR+N
   10      CONTINUE

	   IF (idisk .eq. 0) istart = indcol(icol) + 1

   20   CONTINUE
	RETURN
	END

      subroutine dmerge ( n, db, dc, idy, da, dconst, dl )
C 
C  this merge version has the advantage of loading da(i)
C  and idy(i) only once.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DA(N), DB(*), DC(*), IDY(N)

      dsum = 0.0
      do 30 i = 1, n
         dsum = dsum + da(i) * db(idy(i))
         dc(idy(i)) = dc(idy(i)) + dconst * da(i)
 30   continue
      dl = dsum
      return
      end
      
