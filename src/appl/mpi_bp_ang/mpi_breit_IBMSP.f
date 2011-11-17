*     ==================================================================
*
*      THIS PROGRAM EVALUATES THE BREIT-PAULI OPERATORS
*
*     ONE-ELECTRON OPERATOR
*     ELECTROSTATIC INTERACTION
*     SPIN-ORBIT INTERACTION
*     SPIN-OTHER-ORBIT INTERACTION
*     SPIN-SPIN INTERACTION
*
*       by C. Froese Fischer and G. Tachiev
*          Vanderbilt University
*          Nashville, TN 37235 USA
*          May, 1983
*
*       Modified for Gediminas g-angular structure
*       Modified for Unsorted version
*       Modified for restarting
*                               [Sep 1994]
*
*
*     ==================================================================
*
*       The PARAMETER values in this program define the following:
*               LSDIM - Number of coeffficients in a block
*               NTERMD- Maximum number of terms
*
*     ------------------------------------------------------------------
      PROGRAM Breit_Matrix
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NTERMD=31)
*
*
* >>>>>> mpi <<<<<<
      INCLUDE 'mpif.h'
      parameter (MAXPROC=9999)
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      character*4 idstring
      character*128 startdir, permdir, tmpdir
      character*128 file_cint,file_ih,file_ico,file_clst,file_atomc
* >>>>>>     <<<<<<<

      CHARACTER*26 STRING,ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW
      CHARACTER QREL, QMASS
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(8),ISCW

      CHARACTER*128 input
      EXTERNAL INITT
      logical REL
*
*     Set machine dependent unit numbers
*               IREAD - The unit number of standard input
*               IWRITE- The unit number of standard output
*               ISCW  - The unit number of standard error (screen)
*               iuc   - The unit number for the <name>.c file
*
      iread = 4
      IOUT = 8
      in = 5

      if (myid.eq.0) then
         ISCW = 0
         IWRITE = 6
      else
         iscw = 40
         iwrite = 40
      end if

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
     
      write(idstring,'(I4.4)') myid
      
      startdir = ' '; permdir = ' '; tmpdir = ' ';
      call mpi_work_dir(startdir, permdir, tmpdir)
      lenstart = LEN_TRIM (startdir) - 1
      lenperm = LEN_TRIM (permdir) - 1
      lentmp = LEN_TRIM (tmpdir) - 1

      if (myid == 0) then
        call inp_atom(ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW,
     :      REL,MASS,iscw,in)
      end if

      call MPI_BCAST(ATOMC,26,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(MASS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(REL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr);


      file_cint  = tmpdir(1:lentmp)//'/cint.lst.'//idstring
      file_ih    = tmpdir(1:lentmp)//'/ih.lst.'//idstring
      file_ico   = tmpdir(1:lentmp)//'/ico.lst.'//idstring
      file_clst  = tmpdir(1:lentmp)//'/c.lst.'//idstring
      file_atomc = permdir(1:lenperm)//'/'//atomc
      lcint  = len_trim(file_cint);
      lih    = len_trim(file_ih);
      lico   = len_trim(file_ico);
      lclst  = len_trim(file_clst);
      latomc = len_trim(file_atomc);

      OPEN(UNIT=IREAD,FILE=file_atomc(1:latomc),STATUS='UNKNOWN')
      OPEN(UNIT=IOUT, FILE=file_cint(1:lcint), STATUS='UNKNOWN',
     :     FORM='unformatted')
      OPEN(UNIT=11, FILE=file_ih(1:lih), STATUS='UNKNOWN',
     :     FORM='unformatted')
      OPEN(UNIT=12, FILE=file_ico(1:lico), STATUS='UNKNOWN',
     :     FORM='unformatted')
      OPEN (UNIT=50,FILE=file_clst(1:lclst), STATUS='UNKNOWN',
     :     FORM='UNFORMATTED')
*
*     Let's postpone the restart for the time being
*
      CALL INITA
      CALL BREVALA(MASS)

      !close(IREAD)
      !close(IOUT)
      !close(11)
      !close(12)
      !close(50)

      if(myid == 0) print *, 'Finished with the file'

*      call etime(time)
*      write(iscw,'(//A/A//A,F8.3,A//)')
*     :       ' END OF CASE',' ===========',
*     :       ' Total time was ', time(1)/60, ' minutes' 
      call MPI_FINALIZE(ierr)
      END program breit_matrix
