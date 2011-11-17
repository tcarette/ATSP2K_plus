      subroutine mpi_util_allr_dp(array_dp, size)
      INCLUDE 'mpif.h'
      parameter (MAXPROC=9999)
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer size;
      double precision array_dp(size);
      POINTER(qtmp,tmp(size));
       
      call alloc(qtmp,size,8);
      tmp(1:size) = 0.0;
 
      call MPI_ALLREDUCE(value,tmp_value,n,
     :     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
      array_dp(1:n) = tmp(1:size)
      if(qtmp_value.ne.0) call dalloc(qtmp_value,n)
     
      return
      end
    
  
