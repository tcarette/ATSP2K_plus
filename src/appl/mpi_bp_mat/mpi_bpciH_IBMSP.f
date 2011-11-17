*     ==================================================================
*
*       A CONFIGURATION INTERACTION PROGRAM EITHER NON-RELATIVISTIC
*             OR IN THE BREIT-PAULI APPROXIMATION
*
*       by C. Froese Fischer and G. Tachiev
*          Vanderbilt University
*          Nashville, TN 37235 USA
*          May, 2000
*

*       Modified for Gediminas g-angular structure
*       Modified for Unsorted version
*       Modified for restarting
*       Modified for iso-electronic processing
*                               [May, 2000]
*
*
*     ==================================================================
*
*       The PARAMETER values in this program define the following:
*               NOD   - Maximum number of points in the range of a
*                         - function
*               LIMD  - Maximum number of basis vectors
*
*       Variables for these dimensions are
*               NWF   - Number of functions (or electrons)
*               NO    - Number of points in the range of the function
*               NCFG  - Number of configuration state functions
*               NUME  - Number of eigenvalues/eigenvectors
*               NZE   - Number of non-zero matrix elements
*
*     ------------------------------------------------------------------
      PROGRAM BRCI
*	WITH RESTARTING
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220, NWD=60)
* >>>>>> mpi <<<<<<
      INCLUDE 'mpif.h'
      parameter (MAXPROC=9999)
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      character*4 idstring
      character*128 NAME(2), startdir, permdir, tmpdir
      character*128 file_cint,file_ih,file_ico,file_clst
      character*128 file_atomc,file_atoml,file_atomj,file_atomw
      character*128 file_atomnew,file_iuhnr,file_iouhz,file_iouhs
      integer lenl,lenw,lenc,lenj,lenn
* >>>>>>     <<<<<<<

*
      CHARACTER*26 STRING,ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW
      CHARACTER QREL, QMASS
      LOGICAL REL, SKIP
      logical onlydvd,jhere,lhere,mhere,zhere,shere,restart
*
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,iout,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      INTEGER NJEIG(20)
      REAL TIME(2)
*
*     Set machine dependent unit numbers
*               IREAD - The unit number of standard input
*               IWRITE- The unit number of standard output
*               ISCW  - The unit number of standard error (screen)
*               iuc   - The unit number for the <name>.c file
*               iuw   - The unit number for the <name>.w file
*               ioul  - The unit number for the <name>.l file
*               iouj  - The unit number for the <name>.j file
*               iouhn  - The unit number for storing the LS matrix
*               iouhz - The unit number for storing the ZNV data
*               iouhs - The unit number for storing the S data
*               iouhm - The unit number for writing the scratch matrix
*               iout  - The unit number for cint.lst
*
      IREAD = 5
      iuc    = 4
      iuw    = 2
      ioul   = 7
      iouj   = 8
      iounew = 13

      iouhn  = 9
      iouhz  = 10
      iouhs  = 11
      iouhm  = 12
      iout   = 14

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

      if (myid.eq.0) then
         ISCW = 0
         IWRITE = 6
      else
         iscw = 40
         iwrite = 40
      end if

      if (myid == 0 )  then
         write(iscw,'(A/A)')
         write(iscw,*) '                 ...mpi_bpci_H running on ',
     :                  nprocs,' processors...'
         write(iscw,*)
      end if

      write(idstring,'(I4.4)') myid

      startdir = ' '; permdir = ' '; tmpdir = ' ';
      call mpi_work_dir(startdir, permdir, tmpdir)
      lenstart = LEN_TRIM (startdir) - 1
      lenperm = LEN_TRIM (permdir) - 1
      lentmp = LEN_TRIM (tmpdir) - 1

      CALL INITA
      CALL INITR
      REL = .true.
      
      if(myid == 0) then
        call inp_atom(QREL,QMASS,ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW,
     :      REL,MASS,iscw,iread)
      end if

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_BCAST(QREL,1,MPI_character,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(QMASS,1,MPI_character,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(ATOMC,26,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(REL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(MASS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(restart,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(onlydvd,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(jhere,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(lhere,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(mhere,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr);
      lenc = LEN_TRIM(atomc)
      lenl = LEN_TRIM(atoml)
      lenj = LEN_TRIM(atomj)
      lenw = LEN_TRIM(atomw)
      lenn = LEN_TRIM(atomnew) 
*
*    Units iuc and iuw will be CLOSED by BREVAL
*
      file_cint    = tmpdir(1:lentmp)//'/cint.lst.'//idstring
      file_ih      = tmpdir(1:lentmp)//'/ih.lst.'//idstring
      file_ico     = tmpdir(1:lentmp)//'/ico.lst.'//idstring
      file_clst    = tmpdir(1:lentmp)//'/c.lst.'//idstring
      file_atomc   = permdir(1:lenperm)//'/'//atomc(1:lenc)
      file_atoml   = permdir(1:lenperm)//'/'//atoml(1:lenl)
      file_atomj   = permdir(1:lenperm)//'/'//atomj(1:lenj)
      file_atomw   = permdir(1:lenperm)//'/'//atomw(1:lenw)
      file_atomnew = permdir(1:lenperm)//'/'//atomnew(1:lenn)
      file_iuhnr   = tmpdir(1:lentmp)//'/hnr.lst.'//idstring
      file_iouhz   = tmpdir(1:lentmp)//'/hzeta.lst.'//idstring
      file_iouhs   = tmpdir(1:lentmp)//'/hspin.lst.'//idstring

      if (myid == 0) then
        OPEN(UNIT=iuc,FILE=file_atomc,STATUS='OLD')
        OPEN(UNIT=iuw,FILE=file_atomw,STATUS='OLD',
     :                FORM='UNFORMATTED' )
        OPEN(UNIT=ioul,FILE=file_atoml,STATUS='UNKNOWN',
     :                FORM='UNFORMATTED')
!        OPEN(UNIT=iounew,FILE=file_atomnew,STATUS='UNKNOWN',
!     :                FORM='UNFORMATTED')
      end if

      OPEN(unit=50,file=file_clst,form='UNFORMATTED', status='OLD')
      OPEN(unit=51,file=file_ih,form='UNFORMATTED', status='OLD')
      OPEN(unit=52,file=file_ico,form='UNFORMATTED', status='OLD')
      OPEN(UNIT=iouhn,FILE=file_iuhnr,STATUS='unknown',
     :        form='unformatted')
      OPEN(UNIT=iout,FILE=file_cint,status='OLD',Form='UNFORMATTED')

      if (rel) then
         OPEN(UNIT=iouhz,FILE=file_iouhz,STATUS='unknown',
     :        form='unformatted')
         OPEN(UNIT=iouhs,FILE=file_iouhs,STATUS='unknown',
     :        form='unformatted')
         if (myid == 0) then
            OPEN(UNIT=iouj,FILE=file_atomj,STATUS='UNKNOWN',
     :                 FORM='FORMATTED')
         end if
      end if

!      if (myid == 0) print *, 'Calling Brevalr'
      CALL BREVALR(rel,skip,ec,nze,restart,onlydvd,iounew)
      if (myid == 0) then
!        Print *, 'returned with rel,skip,ec,nze,restart,onlydvd,ioune', 
!     :	rel,skip,ec,nze,restart,onlydvd,iounew 
      end if

        close(iouhn)
        close(iouhz)
        close(iouhs)
      if (myid == 0) print *, 'Finished with the file'

*      call etime(time)
*      write(iscw,'(//A/A//A,F8.3,A//)')
*     :       ' END OF CASE',' ===========',
*     :       ' Total time was ', time(1)/60, ' minutes' 
      call MPI_FINALIZE(MPI_COMM_WORLD,ierr)
      END
