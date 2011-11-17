!***********************************************************************
      subroutine mpi_work_dir(startdir, permdir, tmpdir)
      implicit none
      character(len=*), intent(out):: permdir, tmpdir, startdir
      character*(60) en,ev,path,uname,mpi_dir
      integer*2  mode,serr,hostnm
      integer lpn, plstr,lmpi_dir,lstring
      character*(10) p_name


!      startdir - path where the current node started from.
!      permdir  - path where node-0 performs serial i/o.
!      tmpdir   - path where the current node goes to.
!
      include 'mpif.h'
      integer  myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      COMMON /mpi/ myid, nprocs, ierr, istat
      integer, parameter:: lendisk0 = 128, lenidstring = 3
      character(len=lendisk0) disk
      character(len=lenidstring) idstring
      integer   lendisk, i, len_cwd, iii

!      serr = hostnm(p_name);
!      p_name = trim(p_name);
!      lpn = len_trim(p_name);

!      if (serr.ne.0) print*, 'could''nt get hostname, myid = ',
!     :                       myid, ' exiting..'
!      if (serr.ne.0)  call exit(21);

!     ... get user name
      en = "USER";
      call getenv(en,ev);
      lstring = len_trim(ev);
      uname = trim(ev);

!     ... get the path where mpiruns are supposed to run:
      en = "MPI_TMP";
      call getenv(en,ev);
      lmpi_dir = len_trim(ev);
      mpi_dir = trim(ev);

      startdir = '';
      call sys_getwd (startdir,len_cwd) ! get startdir
      print*, startdir, 'inside mpi_workdir'

      permdir = startdir
      iii = len_trim(startdir) ! - 1
      tmpdir = mpi_dir
      iii = len_trim(tmpdir)
      if (myid == 0) call sys_chdir (tmpdir, iii, ierr)
      if((myid == 0).and.(ierr.ne.0)) call sys_mkdir (tmpdir,iii,ierr)

      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      if(myid > 0) call sys_chdir (tmpdir, iii, ierr)
      
      if (ierr .ne. 0) then
         print *, 'cpath failed, myid = ', myid
         stop
      endif

      return
      end

