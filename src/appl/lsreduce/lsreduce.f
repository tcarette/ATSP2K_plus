* ======================================================================
*
*     GENERAL PROGRAM TO DETERMINE WHICH CSFs IN CLIST DO
*     INTERACT WITH THE MULTIREFERENCE IN MRLIST
*
*     WRITTEN BY -
*
*     P. JONSSON, DEP'T OF COMPUTER SCIENCE
*        VANDERBILT UNIVERISTY, JULY 1995
*
*     THE PROGRAM IS A SMALL MODIFICATION OF THE GDUNONH
*     PROGRAM BY
*
*     G. GAIGALAS, INSTITUTE OF THEOERETICAL PHYSICS
*        AND ASTRONOMY, VILNIUS, LITHUANIA
*
*     C. FROESE FISCHER, DEP'T OF COMPUTER SCIENCE
*        VANDERBILT UNIVERISTY
*
*      ON THE INPUT.
*      mrlist;    CONTAINS THE MULTIREFERENCE CONFIGURATIONS
*                 WITH ANGULAR COUPLINGS
*      clist:     CONTAINS THE CFSs FROM LSGEN THAT SHOULD BE ANALYZED
*
*      ON THE OUTPUT
*      clist.out; CONTAINS 1) THE MULTIREFERENCE
*                          2) THE CSFs FROM clist THAT INTERACTS WITH THE
*                             MULTIREFERENCE
*
*     C O P Y R I G H T  2000
* ======================================================================
*
      PROGRAM LSREDUCE
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
ctc      PARAMETER (NWCD=20,NBD=20,LSDIM=30000,NWD=60)
      PARAMETER (NWCD=20,NBD=20,LSDIM=1000000,NWD=70)

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
      pointer(qmrsdint,mrsdint(1))
      common/mrsdcom/qmrsdint
      DIMENSION noint(4)
      INTEGER NBsize(NBD)
      REAL TIME(2), ETIME, ELAPSE
      LOGICAL lused, yclist, endc 
      Character*72 string

      CHARACTER LINE*72, EL(NWD)*3
      character*3 term
      character*5 su 
      integer l_clist
      character*2 ih_file

      CHARACTER*16 INPUT
      EXTERNAL INITT

*   1 FORMAT(//' IOUT =  FGR.LST (OUTPUT FILE)'/
*    :         ' IBUG1  =',I3,' (DEBUG IN WEIGHTS OF 1-EL PART)'/
*    :         ' IBUG2  =',I3,' (DEBUG IN WEIGHTS OF 2-EL PART)'/
*    :         ' IBUG3  =',I3,' (DEBUG IN RECOUPLING PACKAGE)'//)
*
    2 FORMAT(///20X,'   ==============================='/
     :            20X,'     L S R E D U C E    2000',/
     :            20X,'   ==============================='//)
  3   FORMAT(A30,I3,I4,I6,I8,I8,I8,2x,a5)

      i = iargc()
      if (i .eq. 0) then
         INPUT = 'cfg.inp'
         inquire( FILE=input, exist=yclist) 
         if (yclist) then 
            print *, 'input file is clist...'
         else 
            print *, 'cfg.inp not found: nonh is exiting!...'
            call exit(0)
          endif	
       end if

* ...  THE FOLLOWING SECTION CONCERNS INPUT/OUTPUT AND MAY BE
*      SYSTEM DEPENDENT.  CHECK ALLOWED UNIT NUMBERS AND
*      FILE NAME CONVENTIONS - MODIFY, IF NECESSARY.

      IREAD=4
      ISCW = 0
      IWRITE = 6
      IOUT=8
      WRITE(IWRITE,2)
*
      OPEN(UNIT=17,FILE=inPUT,STATUS='UNKNOWN')
      OPEN(UNIT=IREAD,FILE='slask',STATUS='UNKNOWN')
      OPEN(UNIT=15,FILE='mrlist',STATUS='UNKNOWN')

*
*     ... END OF MACHINE DEPENDENT SECTION
*
Cww  Read the CSFs in the multireference and put them first in clist
      call merge(ncfgmr)
      close(unit=17)
      close(unit=15)
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
* 
*     IN lsreduce there is only one block
      ncfg = nbsize(1)
      call alloc(qjptr, ncfg,4)
      call alloc(qmrsdint, ncfg,4)
      do i=1,ncfg
         mrsdint(i) = 0
      end do


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
*      do N = 1, NB
       do N = 1, 1
*        write(ih_file,'(I2.2)') N 
*        OPEN(UNIT=11, FILE='ih.'//ih_file//'.lst',STATUS='UNKNOWN',
*    :     FORM='unformatted')
 
         ncfg = NBsize(N)
         NEW = NCFG
      IF (NZERO .EQ. 0) THEN
         NZERO = NCFG
         IFIRST = 0
      ELSE
         IFIRST = 1
      END IF
      IALL = 1

         nrec = 0

*     .. allocate memory for this block
         call alloc(qnoc,ncfg,4)
         call alloc(qnelcsh,8*ncfg,4)
         call alloc(qnocorb,8*ncfg,4)
         call alloc(qj1,15*ncfg,4)
cgd      call alloc(qjptr, ncfg,4)
cgd      call alloc(qjptr, ncfg,4) 

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
*
*        Loop only over ncfgmr
         do jb = 1,ncfgmr
           nih = 0
           CALL SHELLSJB(jb)
           CALL ANGMOMG(NEW,NZERO,IFIRST,nih)
*          write(11) nih, (jan(i),i=1,nih)
*          write(12) nih, (ico(i),i=1,nih)
*          max_nze = max(max_nze,nih)
           ncol = ncol + 1
           jptr(ncol) = nij
         end do

*     .. deallocate memory for buffered i/o

      nf = nrec(1)
      ng = nrec(2)
      nr = nrec(3)
      nl = nrec(4)
      ITOTAL = NF+NG+NR+NL
*     write(iscw,220) nij,nf,ng,nr,nl,itotal
* 220 FORMAT( I8, ' non-zero matrix elements'/
*    :      I8,' NF',I8,' NG',I8,' NR',I8,' NL'/
*    :       I8,' Total number of integrals')

*    .. write block information to cfg.inp
*     write(25,'(2x,A3,I6,I10,I15,I10)') term,ncfg,nij,ITOTAL,max_nze
      
         call dalloc(qnoc,ncfg)
         call dalloc(qnelcsh,8*ncfg)
         call dalloc(qnocorb,8*ncfg)
         call dalloc(qj1,15*ncfg)
      end do
*    ..end loop on all blocks
*<<<<< END BLOCKS >>>


*     .. end the processing
      call dalloc(qljclsd,nwcd)
      call dalloc(qpackn,noint)
      call dalloc(qintptr,2*lmax+1)
      call dalloc(qcn,lsdim)
      call dalloc(qinptr,lsdim)
cgd   call dalloc(qnijptr,lsdim)
cgd   call dalloc(qjan,lsdim)
cgd   call dalloc(qjbn,lsdim)
cgd   call dalloc(qico,lsdim)
cgd   call dalloc(qnjcomp,nwfd)
cgd   call dalloc(qljcomp,nwfd)

*      elapse = etime(time)
*      elapse = etime_(time)
*6     write(iscw,'(//A/A//A,F8.3,A//)') ' END OF CASE',' ===========',
*     :      ' Total CPU time was ', TIME(1)/60, ' minutes'

      print *
      print *,  'end-of-file clist!!!'
      call intcfg
      CLOSE(UNIT=IREAD,STATUS='DELETE')



      END
