* ======================================================================
*
*     GENERAL PROGRAM TO COMPUTE MATRIX ELEMENTS OF THE  NON-
*     RELATIVISTIC HAMILTONIAN UNDER THE FOLLOWING ASSUMPTIONS -
*         (1) LS-COUPLING
*         (2) ORTHO-NORMAL CONFIGURATION STATE FUNCTIONS
*         (3) ALL ORBITALS ARE ORTHOGONAL
*     VERSION ADAPTED FOR INTERFACING WITH STOCK
*
*     WRITTEN BY -
*     C. FROESE FISCHER, DEP'T OF COMPUTER SCIENCE
*        VANDERBILT UNIVERISTY
*     FEBRUARY, 1994
*
*     MODIFIED:
*     AUGUST,   1994  (A Ynnerman for unsorted lists)
*     DECEMBER, 1995                                 ( f-sell included )
*               1997  (G. Gaigalas, Vilnius for new angular codes)
*     DECEMBER, 1998  (C. Froese Fischer and G. Tachiev - block version)
*
*     C O P Y R I G H T  2000
*
*     MODIFIED:
*     NOVEMBER, 2011  (T. Carette - interface with Stock)
* ======================================================================
*
      PROGRAM NONH
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWCD=20,NBD=20,LSDIM=30000,NWD=70)

      POINTER (qcn,cn(LSDIM)),(qinptr,inptr(lsdim)),
     :        (qpackn,ipackn(lsdim)),
     :        (qnijptr,nijptr(lsdim)),(qjan,jan(lsdim)),
     :        (qjbn,jbn(lsdim)),(qintptr,idummy(lsdim)),
     :        (qlused,lused(lsdim)),
     :        (qico,ico(lsdim)) 
      COMMON /buffer/qcn,qinptr,qpackn,qlused,qintptr,lmax,qnijptr,
     :               qjan,qjbn,qico
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/DIMEN/KFL1,KFL2,KFL3,KFL4,KFL5,KFL6,KFL7,MXIHSH
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,
     : IALL,JSC(3),ISCW, state
      COMMON /DIAGNL/IDIAG,JA,JB
      POINTER  (qjptr, jptr(lsdim))
      COMMON /fout/ncoff,ntot,idum(6),nrec(8),iflag,lij,nij,qjptr
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG, nlines, endc
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      COMMON /CLOSED/B1ELC(4),NCLOSD,IBK
      COMMON /OPERAT/ ICOLOM,ISOTOP,IORBORB
      DIMENSION noint(4)
      INTEGER NBsize(NBD)
      REAL TIME(2), ETIME, ELAPSE
      LOGICAL lused, yclist, endc 
      Character*72 string
      CHARACTER LINE*72, EL(NWD)*3
      character*3 term
      character*5 su 
      integer l_clist
      character*2 ih_file, ls

      CHARACTER*16 INPUTP,INPUTR
      EXTERNAL INITT

    1 FORMAT(//' IOUT =  FGR.LST (OUTPUT FILE)'/
     :         ' IBUG1  =',I3,' (DEBUG IN WEIGHTS OF 1-EL PART)'/
     :         ' IBUG2  =',I3,' (DEBUG IN WEIGHTS OF 2-EL PART)'/
     :         ' IBUG3  =',I3,' (DEBUG IN RECOUPLING PACKAGE)'//)
*
    2 FORMAT(///20X,  '   ==============================='/
     :            20X,'       N O N H_E X P      2011    ',/
     :            20X,'   ==============================='//)
  3   FORMAT(A30,I3,I4,I6,I8,I8,I8,2x,a5)

      su = 'snonh'
      i = iargc()
      if (i .eq. 0) then
         INPUTP = 'parents'
         inputr = 'ground'
      else
         call getarg(1,inputp)
         if(i .eq. 2) then
           call getarg(2,inputr)
         else
           inputr='ground'
         endif
      end if
      inquire( FILE=inputp, exist=yclist) 
      if (yclist) then 
        print *, 'parents input file is ', trim(inputp)
      else 
        print *, trim(adjustl(inputp)),' not found: nonh is exiting!...'
        call exit(0)
      endif	
      inquire( FILE=inputr, exist=yclist) 
      if (yclist) then 
        print *, '"Ground states" input file is ', trim(inputr)
      else 
        print *, trim(adjustl(inputr)),' not found: nonh is exiting!...'
        call exit(0)
      endif	

* ...  THE FOLLOWING SECTION CONCERNS INPUT/OUTPUT AND MAY BE
*      SYSTEM DEPENDENT.  CHECK ALLOWED UNIT NUMBERS AND
*      FILE NAME CONVENTIONS - MODIFY, IF NECESSARY.

      IREAD=15
      IREADP=IREAD+1
      IREADR=IREAD+2
      ISCW = 0
      IWRITE = 6
      IOUT=8
      WRITE(IWRITE,2)
*
      OPEN(UNIT=IREADP,FORM='FORMATTED',FILE=inPUTp,STATUS='UNKNOWN')
      OPEN(UNIT=IREADR,FORM='FORMATTED',FILE=inPUTr,STATUS='UNKNOWN')
      OPEN(UNIT=IREAD,FORM='FORMATTED',FILE='cfg.inp',STATUS='UNKNOWN')
      
      OPEN(UNIT=IOUT, FILE='yint.lst',STATUS='UNKNOWN',
     :     FORM='unformatted')
      OPEN(UNIT=12, FILE='ico.lst',STATUS='UNKNOWN',
     :     FORM='unformatted')
      OPEN (UNIT=50,FILE='c.lst',STATUS='UNKNOWN',
     :     FORM='UNFORMATTED')
      OPEN(UNIT=25,FILE='cfg.h',STATUS='UNKNOWN')
*
*     ... END OF MACHINE DEPENDENT SECTION
*
*
Cww  Read the CSFs in the multireference and put them first in clist
      call merge(10,'2Pe')
      close(unit=ireadp)
      close(unit=ireadr)
      rewind(unit=iread)

      ICOLOM=1
      ISOTOP=0
      IORBORB=0
      IBUG1 = 0
      IBUG2 = 0
      IBUG3 = 0
      ist=0
      IFIRST = 0
      IALL = 1

*
*  ---  Determine input data; non-orthogonal case
*
      call inita

      call analy_blk(IREAD,IWRITE,NCLOSD,MAXORB,NB,NBsize,EL)

      call orbitals(maxorb,el,qiajcmp,qljcomp,
     :                    qnjcomp,qiajcld,qljclsd,nb)

*  .. find maximum l-value
      lmax = 0
      do i=1,maxorb
         lmax = max (lmax,ljcomp(i))
      end do

*  .. allocate memory for buffered i/o
      call alloc(qcn,lsdim,8)
      call alloc(qinptr,lsdim,4)
      call alloc(qnijptr,lsdim,4)
      call alloc(qjan,lsdim,4)
      call alloc(qjbn,lsdim,4)
      call alloc(qico,lsdim,4)

*  .. generate list of integrals
      call genint(maxorb,lmax,qljcomp,qintptr,qpackn,qlused,noint)

*  .. initialize lused
      do i=1,noint(4)
       lused(i) = .FALSE.
      end do

*  .. write global information to cfg.inp
      nint = noint(4)
      write(25, '(I3,2I8,3x,A5)') nb, nint, lsdim, su

      rewind (iread)
 
*  .. for all blocks
      do N = 1, NB
         write(ih_file,'(I2.2)') N 
         OPEN(UNIT=11, FILE='ih.'//ih_file//'.lst',STATUS='UNKNOWN',
     :     FORM='unformatted')
 
         ncfg = NBsize(N)
         NEW = NCFG
         NZERO = NCFG
         nrec = 0

*     .. allocate memory for this block
         call alloc(qnoc,ncfg,4)
         call alloc(qnelcsh,8*ncfg,4)
         call alloc(qnocorb,8*ncfg,4)
         call alloc(qj1,15*ncfg,4)
         call alloc(qjptr, ncfg,4)
         call alloc(qjptr, ncfg,4) 

*        read CFG's for this block
         CALL CFGBLK(ncfg,maxorb,QIAJCMP,QLJCOMP,QNJCOMP,QNOC,
     :                  QNELCSH,QNOCORB,QJ1,QIAJCLD,QLJCLSD,term)
 
*
*        Initialize parameters for output
*
         lij = 0
         nij = 0
         ncol = 0
         ntot = 0
         ncoff = 0
         max_nze = 0;
         do jb = 1,ncfg
           nih = 0
           CALL SHELLSJB(jb)
           CALL ANGMOMG(NEW,NZERO,IFIRST,nih)
           write(11) nih, (jan(i),i=1,nih)
           write(12) nih, (ico(i),i=1,nih)
           max_nze = max(max_nze,nih)
           ncol = ncol + 1
           jptr(ncol) = nij
         end do

*     .. finish writing the coefficient data, if non empty arrays
	if (ncoff.eq.lsdim) then
	   write(50) lsdim,(cn(j),j=1,lsdim),(inptr(j),j=1,lsdim)
           ncoff=0
           cn(ncoff)=0
           inptr(ncoff)=0
         end if
         write(50) ncoff,(cn(j),j=1,ncoff),(inptr(j),j=1,ncoff)

         write(iout) ncol, (jptr(i),i=1,ncol)

*     .. deallocate memory for buffered i/o

      nf = nrec(1)
      ng = nrec(2)
      nr = nrec(3)
      nl = nrec(4)
      ITOTAL = NF+NG+NR+NL
      write(iscw,220) nij,nf,ng,nr,nl,itotal
  220 FORMAT( I8, ' non-zero matrix elements'/
     :      I8,' NF',I8,' NG',I8,' NR',I8,' NL'/
     :       I8,' Total number of integrals')

*    .. write block information to cfg.inp
      write(25,'(2x,A3,I6,I10,I15,I10)') term,ncfg,nij,ITOTAL,max_nze
      
         call dalloc(qnoc,ncfg)
         call dalloc(qnelcsh,8*ncfg)
         call dalloc(qnocorb,8*ncfg)
         call dalloc(qj1,15*ncfg)
         close(unit=11);
      end do
*    ..end loop on all blocks
*<<<<< END BLOCKS >>>

        iscase = 1
        do icase=1,4
          write(iout) icase,noint(icase)
          write(iout) (ipackn(i),i=iscase,noint(icase)),
     :               (lused(i),i=iscase,noint(icase))
          iscase = noint(icase) + 1
        end do

*     .. end the processing
      call dalloc(qljclsd,nwcd)
      call dalloc(qpackn,noint)
      call dalloc(qintptr,2*lmax+1)
      call dalloc(qcn,lsdim)
      call dalloc(qinptr,lsdim)
      call dalloc(qnijptr,lsdim)
      call dalloc(qjan,lsdim)
      call dalloc(qjbn,lsdim)
      call dalloc(qico,lsdim)
      call dalloc(qnjcomp,nwfd)
      call dalloc(qljcomp,nwfd)

*      elapse = etime(time)
*      elapse = etime_(time)
*6     write(iscw,'(//A/A//A,F8.3,A//)') ' END OF CASE',' ===========',
*     :      ' Total CPU time was ', TIME(1)/60, ' minutes'

      print *
      print *,  'end-of-file clist!!!'
      endfile 25
      close(unit=25)
      close(unit=11)
      close(unit=12)
      CLOSE(UNIT=50)
      CLOSE(UNIT=IOUT)

      END
