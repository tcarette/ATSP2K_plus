      subroutine mpi_allr_dp(array_dp, size)
      INCLUDE 'mpif.h'
      parameter (MAXPROC=9999)
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer size;
      POINTER(qtmp,tmp(size));
      double precision array_dp(size), tmp;
       
      call alloc(qtmp,size,8);
      tmp(1:size) = 0.0;
 
      call MPI_ALLREDUCE(array_dp,tmp,size,
     :     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
      array_dp(1:size) = tmp(1:size)
      if(qtmp.ne.0) call dalloc(qtmp,size)
     
      return
      end
    
  
