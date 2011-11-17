      subroutine getwd(myid,mpi_dir,lmpi_dir,p_name,lpn)

      character*(128) en,ev,path,uname,mpi_dir 
      integer*2  mode,serr
      integer myid,plstr,lmpi_dir,lstring
      character*(10) p_name

      serr = hostnm(p_name);
      p_name = trim(p_name);
      lpn = len_trim(p_name);
     
      if (serr.ne.0) print*, 'could''nt get hostname, myid = ', 
     :                       myid, ' exiting..'
      if (serr.ne.0)  call exit(21);

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
 
      if (mpi_dir == '') then
         print*, 'Error! $MPI_TMP not set, exiting... myid: ', 
     :            p_name,':',myid;
         call exit(22);
      else 
      if (myid == 0) print*, 'mpifiles will be in directory ', 
     :            mpi_dir
         print'(A8,A5,A1,A6,I2)', 'myname: ', 
     :            p_name,':','myid: ',myid;
      end if;

      end subroutine getwd 

      
       

      
