      subroutine mpi_allr_int(array_int, size)
      INCLUDE 'mpif.h'
      parameter (MAXPROC=9999)
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer size, tmp, array_int(size);
      POINTER(qtmp,tmp(size));
 
      call alloc(qtmp,size,4);
      tmp(1:size) = 0;
 
      call MPI_ALLREDUCE(array_int,tmp,size,
     :     MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror)
      array_int(1:size) = tmp(1:size)
      if(qtmp.ne.0) call dalloc(qtmp,size)
     
      return
      end
    
  
