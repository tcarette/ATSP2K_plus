!***********************************************************************
      subroutine mpi_work_dir(startdir, permdir, tmpdir)
      implicit none
      character(len=*), intent(out):: permdir, tmpdir, startdir

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

      call sys_getwd (startdir,len_cwd) ! get startdir

      permdir = startdir
      iii = len_trim(startdir) - 1
      tmpdir = startdir(1:iii)//'/tmp_mpi'//'\0'
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

