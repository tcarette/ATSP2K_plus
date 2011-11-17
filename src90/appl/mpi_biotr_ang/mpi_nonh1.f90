!
!     ------------------------------------------------------------------
!     N O N H 1
!     ------------------------------------------------------------------
!
      subroutine nonh1(ik) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use parameters_biotr_C
!      USE dimen_C 
!      USE diagnl_C 
      USE fout_C 
      USE nlorb_C 
      use buffer_C
      use csf_C
      use debug_C
      use inform_C
      use ndims_C, ONLY: NCFG
      use non30_C
      use dbg_C
      use closed_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:50:15  11/17/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use cfgn1_I 
      use setsupras_I 
      use angmom_I 
      use outls_I 
      use rscheck_I 
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ik 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: new, nzer0, none, nzero, ifirst, iou, n, j, nl, &
         nint, ncoef, lmx1, i, lmx2, njcom, ljcom&
         , lv, iv, jv, ierr 
      character :: st 
!-----------------------------------------------
!
!      POINTER (qcn,cn(lsdim)), (qpackn,ipackn(lsdim)),
!     :        (qjan,jan(lsdim)),(qjbn,jbn(lsdim))
!      COMMON /buffer/qcn, qpackn, qjan, qjbn
!      POINTER (qintgrl1,intgrl1(1)),(qintptr1,intptr1(1)),
!     :        (qcnn1,cnn1(1)),(qjann1,jann1(1)),(qjbnn1,jbnn1(1)),
!     :        (qintgrl2,intgrl2(1)),(qintptr2,intptr2(1)),
!     :        (qcnn2,cnn2(1)),(qjann2,jann2(1)),(qjbnn2,jbnn2(1))
!      COMMON /CSF/qintgrl1,nint1,qintptr1,qcnn1,ncoef1,qjann1,qjbnn1,
!     :            qintgrl2,nint2,qintptr2,qcnn2,ncoef2,qjann2,qjbnn2
!      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
!      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,
!     : IALL,JSC(3),ISCW
!      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
!     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
!      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
!     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
!     :       (QLJCLSD,LJCLSD(1))
!      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
!      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
!      COMMON /DBG  /IBUGM

! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, istat
      character*4 idstring
      character*128 startdir, permdir, tmpdir
      integer lclst,lentmp,lenstart,lenperm
! >>>>>>     <<<<<<<


    1 format(/,/,' IOUT =  FGR.LST (OUTPUT FILE)'/,' IBUG1  =',i3,&
         ' (DEBUG IN WEIGHTS OF 1-EL PART)'/,' IBUG2  =',i3,&
         ' (DEBUG IN WEIGHTS OF 2-EL PART)'/,' IBUG3  =',i3,&
         ' (DEBUG IN RECOUPLING PACKAGE)'/,/) 
!
!     ...  THE FOLLOWING SECTION CONCERNS INPUT/OUTPUT AND MAY BE
!          SYSTEM DEPENDENT.  CHECK ALLOWED UNIT NUMBERS AND
!          FILE NAME CONVENTIONS - MODIFY, IF NECESSARY.
!
      iread = ik 
      iout = 14 + ik 
      isc3 = 17 
!
!      OPEN(UNIT=ISC3,STATUS='SCRATCH',FORM='FORMATTED')
      open(unit=isc3, status='scratch', form='UNFORMATTED', position='asis') 
 
!
!     ... END OF MACHINE DEPENDENT SECTION
!
      ibug1 = 0 
      ibug2 = 0 
      ibug3 = 0 
!     WRITE(ISCW, '(A)') ' IBUG1,IBUG2,IBUG3 ?  FORMAT(3(I1,1X)) '
!     READ(5,'(I1,1X,I1,1X,I1)') IBUG1,IBUG2,IBUG3
!     WRITE(IWRITE,1) IBUG1,IBUG2,IBUG3
!      WRITE(ISCW, '(A)') ' FULL PRINT-OUT ? (Y/N) '
!      IFULL = 1
!      READ(5,'(A2)') ANS
!      IF ( ANS .EQ. 'Y' .OR. ANS .EQ. 'y') IFULL = 1
!
      new = 0 
      nzero = 0 
      none = 0 
      nzer0 = 0
!
!  ---  Determine input data; non-orthogonal case
!
      call cfgn1 () 
      call setsupras (ik,nclosd) 
!
      new = ncfg

      if (new == 0) new = ncfg 
      if (nzero == 0) then 
         nzero = ncfg 
         ifirst = 0 
      else 
         ifirst = 1 
      endif 
!
      iall = 1 
!
!     Initialize parameters for output
!
      nij = 0 
      lcount = 0 
      nrec = 0 
!     .. allocate memory for buffered i/o
 
      allocate (cn(lsdim),stat=ierr) 
      if (ierr.ne.0) call mem_fail(6,lsdim,'nonh::cn',ierr);
      allocate (ipackn(lsdim),stat=ierr) 
      if (ierr.ne.0) call mem_fail(6,lsdim,'nonh::ipackn',ierr);
      allocate (jan(lsdim),stat=ierr) 
      if (ierr.ne.0) call mem_fail(6,lsdim,'nonh::jan',ierr);
      allocate (jbn(lsdim),stat=ierr) 
      if (ierr.ne.0) call mem_fail(6,lsdim,'nonh::jbn',ierr);
      qcn=>cn;
      qpackn=>ipackn;
      qjan=>jan;
      qjbn=>jbn;      
 
!
      call angmom (new, nzero, ifirst) 
!
!     .. clear arrays which have not been written to disk
      iou = 9 + ik 
      n = lcount 
      if (n /= 0) write (isc3) (cn(j),j=1,n), (ipackn(j),j=1,n), (jan(j),j=1,n)&
         , (jbn(j),j=1,n) 
      rewind (isc3) 
!
!     deallocate the arrays that will not be used any longer
!
      deallocate (cn) 
      deallocate (ipackn) 
      deallocate (jan) 
      deallocate (jbn) 
 
      nl = nrec*lsdim + lcount 
 
      if (ik == 1) then 
         allocate (intgrl1(nl)) 
         if (ierr.ne.0) call mem_fail(6,lsdim,'nonh::intgrl1',ierr);
         allocate (intptr1(nl)) 
         if (ierr.ne.0) call mem_fail(6,lsdim,'nonh::intptr1',ierr);
         allocate (cnn1(nl)) 
         if (ierr.ne.0) call mem_fail(6,lsdim,'nonh::cnn1',ierr);
         allocate (jann1(nl)) 
         if (ierr.ne.0) call mem_fail(6,lsdim,'nonh::jann1',ierr);
         allocate (jbnn1(nl)) 
         if (ierr.ne.0) call mem_fail(6,lsdim,'nonh::jbnn1',ierr);
         qintgrl1=>intgrl1;
         qintptr1=>intptr1;
         qcnn1=>cnn1;
         qjann1=>jann1;
         qjbnn1=>jbnn1;
         call outls (intgrl1, nint1, intptr1, cnn1, ncoef1, jann1, jbnn1) 
         nint = nint1 
         ncoef = ncoef1 
         size_int_1 = nint1;
         size_c_1 = ncoef1;
      else 
         allocate (intgrl2(nl)) 
         if (ierr.ne.0) call mem_fail(6,lsdim,'nonh::intgrl2',ierr);
         allocate (intptr2(nl)) 
         if (ierr.ne.0) call mem_fail(6,lsdim,'nonh::intptr2',ierr);
         allocate (cnn2(nl)) 
         if (ierr.ne.0) call mem_fail(6,lsdim,'nonh::cnn2',ierr);
         allocate (jann2(nl)) 
         if (ierr.ne.0) call mem_fail(6,lsdim,'nonh::jann2',ierr);
         allocate (jbnn2(nl)) 
         if (ierr.ne.0) call mem_fail(6,lsdim,'nonh::jbnn2',ierr);
         qintgrl2=>intgrl2;
         qintptr2=>intptr2;
         qcnn2=>cnn2;
         qjann2=>jann2;
         qjbnn2=>jbnn2;
         call outls (intgrl2, nint2, intptr2, cnn2, ncoef2, jann2, jbnn2) 
         nint = nint2 
         ncoef = ncoef2 
         size_int_2 = nint2;
         size_c_2 = ncoef2
      endif 
 
      if (myid == 0) then
         write (iscw, '(I10, A)') nij, ' non-zero matrix elements' 
         write (iscw, '(I10, A)') ncoef, ' number of L(i,j) coefficients' 
         write (iscw, '(I10, A)') nint, ' number of L(i,j) integrals' 
      end if
 
      close(isc3) 
 
! Obtain n and l values of the orbital for the initial (or the final state)
! and check if the order found by nonh satisfies the built-in RAS order.
! (We might be in trouble with the (j.leq.i) selection of savels otherwise.)

      if (myid == 0) then 
         write (6, *) ' ik = ', ik, ' maxorb = ', maxorb 
         write (6, *) '      njcomp = ', njcomp 
         write (6, *) '      ljcomp = ', ljcomp 
      end if 

      if (ik == 1) then 
         lmx1 = 0 
         do i = 1, maxorb 
            nq1(i) = njcomp(i) 
            lq1(i) = ljcomp(i) 
            lmx1 = max0(lq1(i),lmx1) 
            if (myid == 0) &
                write (*, *) ' Initial Orbital n l', nq1(i), lq1(i) 
         end do 
         call rscheck (nq1, lq1, lmx1, maxorb) 
      else 
         lmx2 = 0 
         do i = 1, maxorb 
            nq2(i) = njcomp(i) 
            lq2(i) = ljcomp(i) 
            if (lq2(i) > lmx2) lmx2 = lq1(i) 
            if (myid == 0) write (*, *) &
              ' Final   Orbital n l', nq2(i), lq2(i) 
         end do 
         call rscheck (nq2, lq2, lmx2, maxorb) 
      endif 
 
! Here, one can now deallocate
 
      if (myid == 0) then 
         write (6, *) ' dalloc, qnjcomp: maxorb = ', maxorb 
         write (6, *) ' dalloc, qljcomp: maxorb = ', maxorb 
      end if
         deallocate (njcomp) 
         deallocate (ljcomp) 
 
!     Lets do some debugging on files 'debug1/2'
!     (but only if needed!)
 
      if (ik == 1) st = 'I' 
      if (ik == 2) st = 'F' 
 
!      open(unit=50, file='bc.lst.'//st//'.form.s', status='unknown', form=&
!         'FORMATTED', position='asis') 
 
!      if(ibugm.eq.0) go to 6
      if (ik == 1) then 
!      !open(unit=21,file='debug1',status='unknown')
!         write (50, *) ' Data for l.h.s.' 
!         write (50, *) nint1, ' Integrals' 
         do i = 1, nint1 
            lv = intgrl1(i) 
            iv = mod(lv,64) 
            lv = lv/64 
            jv = mod(lv,64) 
            lv = lv/64 
!            write (50, '(4I6)') lv, jv, iv, intptr1(i) 
         end do 
!         write (50, *) ncoef1, ' Coefficients' 
!         write (50, '(F12.8,2I6)') (cnn1(j),jann1(j),jbnn1(j),j=1,ncoef1) 
!         close(50) 
      else 
!      !open(unit=22,file='debug2',status='unknown')
!         write (50, *) 
!         write (50, *) ' Data for r.h.s.' 
!         write (50, *) nint2, ' Integrals' 
         do i = 1, nint2 
            lv = intgrl2(i) 
            iv = mod(lv,64) 
            lv = lv/64 
            jv = mod(lv,64) 
            lv = lv/64 
!            write (50, '(4I6)') lv, jv, iv, intptr2(i) 
         end do 
!         write (50, *) ncoef2, ' Coefficients' 
!         write (50, '(F12.8,2I6)') (cnn2(j),jann2(j),jbnn2(j),j=1,ncoef2) 
!         close(50) 
      endif 
      if (myid==0) write (iscw, '(2/A/A,2/)') ' END OF NON', ' ==========' 
      return  
      end subroutine nonh1 
