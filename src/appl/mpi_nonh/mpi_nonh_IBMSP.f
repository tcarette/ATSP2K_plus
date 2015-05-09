* ======================================================================
*
*     GENERAL PROGRAM TO COMPUTE MATRIX ELEMENTS OF THE  NON-
*     RELATIVISTIC HAMILTONIAN UNDER THE FOLLOWING ASSUMPTIONS -
*         (1) LS-COUPLING
*         (2) ORTHO-NORMAL CONFIGURATION STATE FUNCTIONS
*         (3) ALL ORBITALS ARE ORTHOGONAL
*
*     WRITTEN BY -
*     G. GAIGALAS, INSTITUTE OF THEOERETICAL PHYSICS
*        AND ASTRONOMY, VILNIUS, LITHUANIA
*
*     C. FROESE FISCHER, DEP'T OF COMPUTER SCIENCE
*        VANDERBILT UNIVERISTY
*
*
*     FEBRUARY, 1994
*     AUGUST,   1994  (A Ynnerman for unsorted lists)
*     DECEMBER, 1995                                 ( f-sell included )
*     DECEMBER, 1998  (C. Froese Fischer and G. Tachiev - block version)
*
*
* ======================================================================
*
      PROGRAM snonh
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWCD=20,NBD=20,LSDIM=30000,NWD=94)

*
*     MPI stuff ***********************************************
*
        INCLUDE 'mpif.h'
        parameter (MAXPROC=100)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
        Character*2 idstring
        Character*72 program,arch,hostn,odir,output,outc
        character*128 NAME(2),workpath
        logical :: f_out
        integer tids, nnn(6), working_procs, cf_tot
        integer group_L
	integer, allocatable, dimension(:) :: proc_ncoff
        integer ncfg_left, new_group, new_id, color,MPI_GROUP_WORLD
        integer comm_L, comm_last,itot_ng,itot_nr,itot_nf,itot_nl 
        integer, allocatable, dimension(:) :: nij_buff
        double precision :: message(200)
        integer :: max_buffer
        real*4 speed(0:MAXPROC),timarr(2),rstart,rfin,total,etime
        real :: timing(4)
        data speed / 101*1000.0 /
        character*(128) mpi_dir,cwd,sh_command,tmpdir
        character*(128) file_cl,file_c,file_y,file_cfg,p_name

        integer lmpi_dir,lcwd,size,llc,lly,llcf,llcl,lpn
        integer*2 serr
****************************************************************

      POINTER (qcn,cn(1)),(qinptr,inptr(1)),
     :        (qpackn,ipackn(1)),
     :        (qnijptr,nijptr(1)),(qjan,jan(1)),
     :        (qjbn,jbn(1)),(qintptr,idummy(1)),
     :        (qlused,lused(1)),
     :        (qico,ico(1)) 
      COMMON /buffer/qcn,qinptr,qpackn,qlused,qintptr,lmax,qnijptr,
     :               qjan,qjbn,qico
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/DIMEN/KFL1,KFL2,KFL3,KFL4,KFL5,KFL6,KFL7,MXIHSH
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,
     : IALL,JSC(3),ISCW, state
      COMMON /DIAGNL/IDIAG,JA,JB
      POINTER  (qjptr, jptr(1))
      pointer (qcn_g, cn_g(1)) ,(qinptr_g, inptr_g(1))
      COMMON /fout/ncoff,ntot,idum(6),nrec(8),iflag,lij,nij,qjptr,cf_tot
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
      POINTER (qltused,ltused(1))
      logical ltused
      INTEGER NBsize(NBD), n_sum
*      REAL TIME(2), ETIME, ELAPSE
      LOGICAL lused, yclist, endc 
      Character*72 string
      CHARACTER LINE*72, EL(NWD)*3
      character*3 term
      character*5 su 
      integer l_clist

      CHARACTER*16 INPUT
      EXTERNAL INITT

    1 FORMAT(//' IOUT =  FGR.LST (OUTPUT FILE)'/
     :         ' IBUG1  =',I3,' (DEBUG IN WEIGHTS OF 1-EL PART)'/
     :         ' IBUG2  =',I3,' (DEBUG IN WEIGHTS OF 2-EL PART)'/
     :         ' IBUG3  =',I3,' (DEBUG IN RECOUPLING PACKAGE)'//)
*
    2 FORMAT(///20X,'   ==============================='/
     :            20X,'         S N O N H_M P I   2000',/
     :            20X,'   ==============================='//)
  3   FORMAT(A30,I3,I4,I6,I8,I8,I8,2x,a5)

* ...  THE FOLLOWING SECTION CONCERNS INPUT/OUTPUT AND MAY BE
*      SYSTEM DEPENDENT.  CHECK ALLOWED UNIT NUMBERS AND
*      FILE NAME CONVENTIONS - MODIFY, IF NECESSARY.

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr) 
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      timing(1) = MPI_WTIME()
 
      write(idstring,'(I2.2)') myid
      input = 'cfg.inp'

      su = 'snonh'
      i = iargc()
      if (i .eq. 0) then
         INPUT = 'cfg.inp'
         inquire( FILE=input, exist=yclist)
         if (yclist) then
            if(myid == 0) write (iscw,*) 'input file is cfg.inp ...'
         else
            write (iscw,*) 'cfg.inp not found: nonh is exiting!...'
            call exit(0)
          endif
       end if

      IREAD=4
      IOUT=8

      if (myid.eq.0) then
         ISCW = 0
         IWRITE = 6
      else
         iscw = 40
         iwrite = 40
      end if

      WRITE(IWRITE,2)
      write(iscw,'(A,i4,A)') 
     :      '                 ...snonh_mpi running on ',
     :                  nprocs,' processors...'
      write(iscw,*)


*>>>>>>>> specify files for input output >>>>>>>>>>>>>>>>>

      cwd = " "; mpi_dir = " "; tmpdir = " ";
      call mpi_work_dir(cwd, mpi_dir, tmpdir);
      lcwd = LEN_TRIM(cwd)  - 1
      lenperm = LEN_TRIM(mpi_dir) - 1
      lentmp = LEN_TRIM(tmpdir)  - 1

      file_c = tmpdir(1:lentmp)//'/c.lst.'//idstring
      file_y = tmpdir(1:lentmp)//'/yint.lst.'//idstring
      file_cl = cwd(1:lcwd)//'/cfg.inp'
      file_cfg = cwd(1:lcwd)//'/cfg.h'
      llc = len_trim(file_c)
      lly = len_trim(file_y)
      llcl = len_trim(file_cl)
      llcf = len_trim(file_cfg)

      !print*, file_c,file_y,file_cl,file_cfg
      open(unit=39,file=file_c(1:llc),status='unknown',
     :     form='unformatted');
      !serr = chmod(file_c(1:llc),511)
      !if (serr.ne.0) print *, 'can''t chmod to 4777 of', file_c
      open(unit=38,file=file_y(1:lly),status='unknown',
     :     form='unformatted');
      !serr = chmod(file_y(1:lly),511)
      !if (serr.ne.0) print *, 'can''t chmod to 4777 of', file_y
      OPEN(UNIT=IREAD,FILE=file_cl(1:llcl),STATUS='UNKNOWN')
      if (myid == 0) then
         OPEN(UNIT=25,FILE=file_cfg(1:llcf),STATUS='UNKNOWN')
      end if
       
*>>>>>>>>>>>> end files <<<<<<<<<<<<<<<<<<<<<<<<<<
*
*     ... END OF MACHINE DEPENDENT SECTION
*

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

*  ..   find maximum l-value
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

      call genint(maxorb,lmax,qljcomp,qintptr,qpackn,qlused,noint,iscw)

*  .. initialize lused
      lused(1:noint(4)) = .FALSE.

*  .. write global information to cfg.inp
      nint = noint(4)

      if (myid==0) then
      write(25, '(I3,2I8,3x,A5)') nb, nint, lsdim, su
      end if

      rewind (iread)
*      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      timing(2) = MPI_WTIME() 
      ncfg_total = sum(NBsize)
      allocate (proc_ncoff(NB*nprocs));
      proc_ncoff = 0;
      max_buffer = 0
      cf_tot = 0
*  .. for all blocks
      do NLB = 1, NB
         ncfg = NBsize(NLB)
         NEW = NCFG
         NZERO = NCFG
         nrec = 0
         cf_tot = 0
*     .. allocate memory for this block
         call alloc(qnoc,ncfg,4)
         call alloc(qnelcsh,8*ncfg,4)
         call alloc(qnocorb,8*ncfg,4)
         call alloc(qj1,15*ncfg,4)
         call alloc(qjptr,ncfg,4)

*        read CFG's for this block
         CALL CFGBLK(ncfg,maxorb,QIAJCMP,QLJCOMP,QNJCOMP,QNOC,
     :                  QNELCSH,QNOCORB,QJ1,QIAJCLD,QLJCLSD,term)
        write(iscw,'(A,A,A,i8,A)')'processing ',term,' with ',
     :                            ncfg,' configurations'
*
*        Initialize parameters for output
*
         lij = 0; nij = 0; mycol = 0; ntot = 0; ncoff = 0;
         n_start = 1; nj_start = 1; nij_tot = 0; njptr_start = 1;
         ntot_tot = 0;
         
*        ... create communicator for the last ncfg_last 
         ncfg_last = modulo(ncfg,nprocs)
!         if (myid < ncfg_last) then
!            new_s = 1
!         else
!            new_s = 0 
!         end  if
!         call MPI_COMM_SPLIT(group_L,new_s,myid,comm_last,ierr)
!         call MPI_COMM_SIZE(comm_last,nsize,ierr )
!         call MPI_COMM_RANK(comm_last,new_id,ierr)

*        ... allocate memeory for buffered output
         if (myid==0) then
           call alloc(qcn_g,lsdim,8)
           allocate(nij_buff(ncfg))
           nij_buff = 0
         end if

         ncfg_left = ncfg
         f_out = .false.

         do jb = myid+1,ncfg,nprocs
           if (mod(jb,1000) .eq. 0) write(0,*) '   jb =',jb
!           if (mod(jb,10) .eq. 0) write(0,*) '   jb =',jb
           if (jb.eq.ncfg) write(0,*) '   jb =',jb
           if (jb > (ncfg-nprocs)) f_out = .true.
!           if (ncfg_left >=  nprocs) then
!              working_procs = MPI_COMM_WORLD 
!              npw = nprocs - 1 
!           else 
!              working_procs = comm_last !new_group !group_L
!              working_procs = MPI_COMM_WORLD
!              npw = ncfg_last - 1
!           end if 
          
           CALL SHELLSJB(jb)
           CALL ANGMOMG(NEW,NZERO,IFIRST)

           if (myid < ncfg_left) then
!           call MPI_SEND(lij,1,MPI_INTEGER,0,94,MPI_COMM_WORLD,ierr)
!           call MPI_SEND(nij,1,MPI_INTEGER,0,97,MPI_COMM_WORLD,ierr)
           end if

           mycol = mycol + 1
           jptr(mycol) = nij
           ncfg_left = ncfg_left - nprocs
        end do
*     call mpi_barrier(MPI_COMM_WORLD,ierr)
*>>>>
      call MPI_Reduce(nij,nij_tot,1,MPI_INTEGER,MPI_MAX,0,
     :                 MPI_COMM_WORLD,ierr)

*>>>>>
*     .. finish writing the coefficient data, if non empty arrays
*      ..... write jptr contained in nij_buff

       write(38) lij,(jan(i),i=1,lij), (ico(i),i=1,lij)
*       if (myid == 0) then
*          print*, lij,(jan(i),i=1,lij), (ico(i),i=1,lij)
*       endif
       write(38) mycol,(jptr(i),i=1,mycol)
*>>>>
        if (ncoff.eq.lsdim) then
           write(39) lsdim,(cn(j),j=1,lsdim),(inptr(j),j=1,lsdim)
           cf_tot = cf_tot + lsdim
           ncoff=0
           cn(ncoff)=0
           inptr(ncoff)=0
	 else 
            cf_tot = cf_tot + ncoff
         end if
         write(39) ncoff, (cn(j),j=1,ncoff),(inptr(j),j=1,ncoff)
*         if (myid == 0) then
*            print*, ncoff, (cn(j),j=1,ncoff),(inptr(j),j=1,ncoff) 
*         end if
      call MPI_GATHER(cf_tot,1,MPI_INTEGER,
     :                 proc_ncoff(((NLB-1)*nprocs)+myid+1),
     :                 1,MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)

*>>>>

      timing(3) = MPI_WTIME()
*     .. deallocate memory for buffered i/o

      nf = nrec(1)
      ng = nrec(2)
      nr = nrec(3)
      nl = nrec(4)

*      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call MPI_Reduce(nf,itot_nf,1,MPI_INTEGER,MPI_SUM,0,
     :                 MPI_COMM_WORLD,ierr)
      call MPI_Reduce(ng,itot_ng,1,MPI_INTEGER,MPI_SUM,0,
     :                 MPI_COMM_WORLD,ierr)
      call MPI_Reduce(nr,itot_nr,1,MPI_INTEGER,MPI_SUM,0,
     :                 MPI_COMM_WORLD,ierr)
      call MPI_Reduce(nl,itot_nl,1,MPI_INTEGER,MPI_SUM,0,
     :                 MPI_COMM_WORLD,ierr)

      ITOTAL = NF+NG+NR+NL
      itot_tot = itot_nf+itot_ng+itot_nr+itot_nl
      write(iscw,220) nij_tot,itot_nf,itot_ng,itot_nr,itot_nl,itot_tot
  220 FORMAT( I8, ' non-zero matrix elements'/
     :      I8,' NF',I8,' NG',I8,' NR',I8,' NL'/
     :       I8,' Total number of integrals')

*    .. write block information to cfg.inp
      if(myid == 0) then
         write(25,'(3x,A3,I8,I8,i20)') term, ncfg, nij_tot
      end if 

         call dalloc(qnoc,ncfg)
         call dalloc(qnelcsh,8*ncfg)
         call dalloc(qnocorb,8*ncfg)
         call dalloc(qj1,15*ncfg)
         if (myid == 0) deallocate(nij_buff)
         if (myid == 0) call dalloc(qcn_g,lsdim);
      end do
*    ..end loop on all blocks
     
      if (myid == 0) then
        do ix1 = 1, NB
          write(25,'(A6,i3)') "Block ", ix1
          do ix2 = 1, nprocs
            write(25,'(i20)') proc_ncoff(((ix1-1)*nprocs)+ix2)
          end do 
        end do
      end if

       call alloc(qltused,noint(4),4)
       ltused(1:noint(4)) = .false.

      call MPI_ALLREDUCE(lused,ltused,noint(4),MPI_LOGICAL,MPI_LOR,
     :                MPI_COMM_WORLD,ierr)

*      if (myid == 0) then
        iscase = 1
        do icase=1,4
          write(38) icase,noint(icase)
          write(38) (ipackn(i),i=iscase,noint(icase)),
     :               (ltused(i),i=iscase,noint(icase))
          iscase = noint(icase) + 1
        end do
*      end if

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
      deallocate (proc_ncoff)

*      elapse = etime(time)
*      elapse = etime_(time)
*6     write(iscw,'(//A/A//A,F8.3,A//)') ' END OF CASE',' ===========',
*     :      ' Total CPU time was ', TIME(1)/60, ' minutes'

      write (iscw,*)  'end-of-file cfg.inp!!!'
      endfile 25
      if (myid == 0) then
        close(50)
        close(8)
      end if
      close(39)
      close(38)
      close(unit=25)
      close(unit=39)
*      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      timing(4) = MPI_WTIME()
*      write(iscw,*) 'time 1= ', timing(2)-timing(1)
*      write(iscw,*) 'time 2= ', timing(3)-timing(2)
*      write(iscw,*) 'time 3= ', timing(4)-timing(3) 
*      write(iscw,*) 'the time per 1 cfg =',   
*     :                (timing(4)-timing(1))/ncfg_total
*      write(iscw,*) 'total =', (timing(4)-timing(1))
      call MPI_FINALIZE(ierr)

      END PROGRAM SNONH
