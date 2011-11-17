*
*     ------------------------------------------------------------------
*       O U T L S
*     ------------------------------------------------------------------
*
      SUBROUTINE OUTLS(intgrl,nint,intptr,cnn,ncoef,jann,jbnn)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*     The list of unique integrals L(j,i) is formed, in the order of 
*     increasing symmetry, with j .le. i.  With each integral (INTGRL)
*     there is a pointer (INTPTR) that points to the last element in the
*     array of coefficients for that integral.
*
      INTEGER intgrl(*), intptr(*), jann(*), jbnn(*)
      DOUBLE PRECISION cnn(*)
*
* On Exit
* -------
*     INTGRL -- array of integrals (packed form of (l, j, i))
*     NINT   -- number of integrals
*     INTPTR -- array of pointers to last element in list of coefficients
*     CNN    -- array of coefficients
*     NCOEF  -- number of coefficients
*     JANN    -- row of matrix
*     JBNN    -- column of matrix
*
      parameter (LSDIM=30000)
*
      COMMON /CLOSED/B1ELC(4),NCLOSD,IBK
      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON /fout/lcount,nrec,iflag,lij,nij
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(4),IALL,JSC(3),ISCW
      PARAMETER (LSTACK=128)
      CHARACTER*1 INT(8),ff
*
*     ..local arrays
      POINTER (qcn,cn(1)),(qpack,ipackn(1)),(qjan,jan(1)),
     :        (qjbn,jbn(1)),(qipt,ipt(1))
      INTEGER ISTACK(LSTACK) , ncc
*
      DATA INT/'F','G','R','R','R','R','R','L'/
*
*===  Begin processing data
*
      ncc = lcount + nrec*lsdim
      iou = isc(4)

      write(iscw,'(A/(4I15))') 'Number of L-integrals',ncc

      call alloc(qcn,ncc,8)
      call alloc(qpack,ncc,4)
      call alloc(qjan,ncc,4)
      call alloc(qjbn,ncc,4)
      call alloc(qipt,ncc,4)
*
      ii     = 0
*
      nint = 0
*        .. read data
         ist = 0
         do 130 jj = 1,nrec
           read(iou) (cn(j+ist),j=1,lsdim),(ipackn(j+ist),j=1,lsdim),
     :               (jan(j+ist),j=1,lsdim),(jbn(j+ist),j=1,lsdim)
           ist = ist + lsdim
  130    continue
         m = lcount
         if (m .ne. 0) read(iou) (cn(j+ist),j=1,m),
     :      (ipackn(j+ist),j=1,m),(jan(j+ist),j=1,m),(jbn(j+ist),j=1,m)
         close(iou)
*
         CALL QSORT(NCC,IPACKN,IPT,ISTACK,LSTACK,IERR)
         IF (IERR .EQ. 1) THEN
            WRITE(iscw,*) ' Stack dimension not large enough for sort',
     :           'CASE = ',icase, 'N = ',nc
            stop 
         END IF
*
*        Form the list of integrals with pointers to the data
*
         LAST = 0
  110    J = LAST +1
         LAST = J
         IF (J .LE. NCC) THEN
*
*           Find  last item in the list with this integral
*           
  120      LAST = LAST + 1
           IF (LAST .LE. NCC) THEN
             IF (IPACKN(IPT(J)) 
     :           .EQ. IPACKN(IPT(LAST)))  GO TO 120
           END IF
           LAST = LAST -1
           NINT = NINT + 1
           intgrl(nint) = ipackn(ipt(last) )
           intptr(nint) = last
           GO TO 110
        END IF
	ncoef = ncc
*
*       .. transfer the sorted data
*
	do 150 j =1, ncoef
	  k = ipt(j)
	  cnn(j) = cn(k)
	  jann(j) = jan(k)
	  jbnn(j) = jbn(k)
  150   continue
*
*      .. deallocate the local arrays.
      call dalloc(qcn,ncc)
      call dalloc(qpack,ncc)
      call dalloc(qjan,ncc)
      call dalloc(qjbn,ncc)
      call dalloc(qipt,ncc)

      END
