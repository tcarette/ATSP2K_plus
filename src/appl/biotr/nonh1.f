*
*     ------------------------------------------------------------------
*	N O N H 1  
*     ------------------------------------------------------------------
*
      SUBROUTINE nonh1(IK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (LSDIM=30000, NWD=128)
*
      POINTER (qcn,cn(lsdim)), (qpackn,ipackn(lsdim)),
     :        (qjan,jan(lsdim)),(qjbn,jbn(lsdim)) 
      COMMON /buffer/qcn, qpackn, qjan, qjbn
      POINTER (qintgrl1,intgrl1(1)),(qintptr1,intptr1(1)),
     :        (qcnn1,cnn1(1)),(qjann1,jann1(1)),(qjbnn1,jbnn1(1)),
     :        (qintgrl2,intgrl2(1)),(qintptr2,intptr2(1)),
     :        (qcnn2,cnn2(1)),(qjann2,jann2(1)),(qjbnn2,jbnn2(1))
      COMMON /CSF/qintgrl1,nint1,qintptr1,qcnn1,ncoef1,qjann1,qjbnn1,
     :            qintgrl2,nint2,qintptr2,qcnn2,ncoef2,qjann2,qjbnn2
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/DIMEN/KFL1,KFL2,KFL3,KFL4,KFL5,KFL6,KFL7,MXIHSH
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,
     : IALL,JSC(3),ISCW
      COMMON /DIAGNL/IDIAG,JA,JB
      COMMON /fout/lcount,nrec,iflag,lij,nij
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      COMMON /DBG  /IBUGM
      character*1 ff
      common/nlorb/nq1(nwd),lq1(nwd),nq2(nwd),lq2(nwd)
*
    1 FORMAT(//' IOUT =  FGR.LST (OUTPUT FILE)'/
     :         ' IBUG1  =',I3,' (DEBUG IN WEIGHTS OF 1-EL PART)'/
     :         ' IBUG2  =',I3,' (DEBUG IN WEIGHTS OF 2-EL PART)'/
     :         ' IBUG3  =',I3,' (DEBUG IN RECOUPLING PACKAGE)'//)
*
*     ...  THE FOLLOWING SECTION CONCERNS INPUT/OUTPUT AND MAY BE
*          SYSTEM DEPENDENT.  CHECK ALLOWED UNIT NUMBERS AND
*          FILE NAME CONVENTIONS - MODIFY, IF NECESSARY.
*
      IREAD = IK
      IOUT = 14+IK
      ISC3=17
*
      OPEN(UNIT=ISC3,STATUS='SCRATCH',FORM='UNFORMATTED')
*
*     ... END OF MACHINE DEPENDENT SECTION
*
      IBUG1 = 0
      IBUG2 = 0
      IBUG3 = 0
*     WRITE(ISCW, '(A)') ' IBUG1,IBUG2,IBUG3 ?  FORMAT(3(I1,1X)) '
*     READ(5,'(I1,1X,I1,1X,I1)') IBUG1,IBUG2,IBUG3
*     WRITE(IWRITE,1) IBUG1,IBUG2,IBUG3
*     WRITE(ISCW, '(A)') ' FULL PRINT-OUT ? (Y/N) '
      IFULL = 0
*     READ(5,'(A2)') ANS
*     IF ( ANS .EQ. 'Y' .OR. ANS .EQ. 'y') IFULL = 1
*
      NEW = 0
      NZER0 = 0
      NONE = 0
*
*  ---  Determine input data; non-orthogonal case
*
      CALL CFGN1(NCLOSD)
      CALL SETSUPRAS(IK,NCLOSD)
*
      IF ( NEW .EQ. 0 ) NEW = NCFG
      IF (NZERO .EQ. 0) THEN
         NZERO = NCFG
         IFIRST = 0
      ELSE
         IFIRST = 1
      END IF
*
      IALL = 1
*
*     Initialize parameters for output
*
      NIJ = 0
      lcount = 0
      nrec = 0
*     .. allocate memory for buffered i/o
      call alloc(qcn,lsdim,8)
      call alloc(qpackn,lsdim,4)
      call alloc(qjan,lsdim,4)
      call alloc(qjbn,lsdim,4)
*
      CALL ANGMOM(NEW,NZERO,IFIRST)
*
*     .. clear arrays which have not been written to disk
      iou = 9+ik
      n = lcount
      if (n .ne. 0) then
            write(isc3) (cn(j),j=1,n), (ipackn(j),j=1,n),
     :      (jan(j),j=1,n),(jbn(j),j=1,n)
      end if
      rewind(isc3)
*
*     deallocate the arrays that will not be used any longer
      call dalloc(qcn,lsdim)
      call dalloc(qpackn,lsdim)
      call dalloc(qjan,lsdim)
      call dalloc(qjbn,lsdim)

      nl = nrec*lsdim + lcount

*     .. we are not taking advantage of the fact that the number
*     of integrals could be up to an order of magnitude less than
*     the number of coefficients.  This could be done with a
*     "realloc" after the call to outls if deemed important.

      if (ik.eq. 1) then
	call alloc(qintgrl1,nl,4)
	call alloc(qintptr1,nl,4)
        call alloc(qcnn1,nl,8)
        call alloc(qjann1,nl,4)
        call alloc(qjbnn1,nl,4)
        call outls(intgrl1,nint1,intptr1,cnn1,ncoef1,jann1,jbnn1)
	nint = nint1
	ncoef = ncoef1
      else
	call alloc(qintgrl2,nl,4)
	call alloc(qintptr2,nl,4)
        call alloc(qcnn2,nl,8)
        call alloc(qjann2,nl,4)
        call alloc(qjbnn2,nl,4)
        call outls(intgrl2,nint2,intptr2,cnn2,ncoef2,jann2,jbnn2)
	nint = nint2
	ncoef = ncoef2
      end if

      write(iscw,'(I10, A)') nij, ' non-zero matrix elements'
      write(iscw,'(I10, A)') ncoef, ' number of L(i,j) coefficients'
      write(iscw,'(I10, A)') nint, ' number of L(i,j) integrals'

      CLOSE (ISC3)

* Obtain n and l values of the orbital for the initial (or the final state)
* and check if the order found by nonh satisfies the built-in RAS order.
* (We might be in trouble with the (j.leq.i) selection of savels otherwise.)

      print*,' ik = ',ik,' maxorb = ',maxorb
      print*,'      njcomp = ',njcomp
      print*,'      ljcomp = ',ljcomp
      if (ik.eq.1) then
	lmx1 = 0
        do i = 1,maxorb
          nq1(i) = njcomp(i)
          lq1(i) = ljcomp(i)
	  if (lq1(i).gt.lmx1) lmx1 = lq1(i)
          write(*,*) ' Initial Orbital n l',nq1(i),lq1(i)
        end do
        call rscheck(nq1,lq1,lmx1,maxorb)
      else
	lmx2 = 0
        do i = 1,maxorb
          nq2(i) = njcomp(i)
          lq2(i) = ljcomp(i)
	  if (lq2(i).gt.lmx2) lmx2 = lq1(i)
          write(*,*) ' Final   Orbital n l',nq2(i),lq2(i)
        end do
      call rscheck(nq2,lq2,lmx2,maxorb)
      endif

* Here, one can now deallocate

      print*,' dalloc, qnjcomp: maxorb = ',maxorb
      call dalloc(qnjcomp,maxorb)
      print*,' dalloc, qljcomp: maxorb = ',maxorb
      call dalloc(qljcomp,maxorb)


*     Lets do some debugging on files 'debug1/2'
*     (but only if needed!)
 
      !if(ibugm.eq.0) go to 6
      if (ik .eq. 1) then
	!open(unit=21,file='debug1',status='unknown')
        !write (21,*) ' Data for l.h.s.'
        !write (21,*) nint1, ' Integrals'
        do i = 1,nint1
	lv   = intgrl1(i)
	iv   = mod(lv,64)
	lv   = lv/64
	jv   = mod(lv,64)
	lv   = lv/64
	!write(21,'(4I6)')   lv,jv,iv, intptr1(I)
        end do
        !write(22,*) ncoef2, ' Coefficients'
        !write(22,'(F12.8,2I6)')
!     :         (cnn1(j),jann1(j),jbnn1(j),j=1,ncoef1)
	!close(21)

!        open(unit=50,file='st1.lst',form='unformatted');
!        write(50) ncoef
!        write(50) nint1
!        write(50) (intptr1(I),intgrl1(i),i=1,nint1);
!        write(50) (cnn1(j),jann1(j),jbnn1(j),j=1,ncoef)
!        close(50)
      else
	!open(unit=22,file='debug2',status='unknown')
        !write(22,*)
        !write (22,*) ' Data for r.h.s.'
        !write (22,*) nint2, ' Integrals'
        do i = 1,nint2
           lv   = intgrl2(i)
           iv   = mod(lv,64)
           lv   = lv/64
           jv   = mod(lv,64)
           lv   = lv/64
!           write(22,'(4I6)')   lv,jv,iv, intptr2(I)
        end do
!        write(22,*) ncoef2, ' Coefficients'
!        write(22,'(F12.8,2I6)')
!     :     (cnn2(j),jann2(j),jbnn2(j),j=1,ncoef2)
!        close(22)
!        open(unit=50,file='st2.lst',form='unformatted');
!        write(50) ncoef
!        write(50) nint2
!        write(50) (intptr2(I),intgrl2(i),i=1,nint2);
!        write(50) (cnn2(j),jann2(j),jbnn2(j),j=1,ncoef2)
!        close(50)
      end if

       !write(ff,'(A1)') ik

6     write(iscw,'(//A/A//)') ' END OF NON',' =========='
      END
