*     ==================================================================
*
*      Program to generate the interaction matrix files in sparse matrix
*      format given files of BREIT-PAULI angular integraion data.
*      
*
*       by C. Froese Fischer and G. Tachiev
*          Vanderbilt University
*          Nashville, TN 37235 USA
*          May, 2000
*
*                  C O P Y R I G H T   2000
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
      PROGRAM bp_mat 
*	WITH RESTARTING
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220, NWD=60)
      character*1 idstring
      character*128 NAME(2), startdir, permdir, tmpdir
      character*128 file_cint,file_ih,file_ico,file_clst
      character*128 file_atomc,file_atoml,file_atomj,file_atomw
      character*128 file_atomnew,file_iuhnr,file_iouhz,file_iouhs
      integer lenl,lenw,lenc,lenj,lenn

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
!      REAL TIME(2)
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

      myid = 0
      nprocs = 1;

         ISCW = 0
         IWRITE = 6
         write(iscw,'(A/A)')
         write(iscw,*) '                 ...bp_mat 2000 '
         write(iscw,*)

!      write(idstring,'(I4.4)') myid
      idstring = 's';

      startdir = '.'; permdir = '.'; tmpdir = '.';
      lenstart = 1; lenperm = 1; lentmp = 1;
      !call work_dir(startdir, permdir, tmpdir)
      !lenstart = LEN_TRIM (startdir) 
      !lenperm = LEN_TRIM (permdir) 
      !lentmp = LEN_TRIM (tmpdir)

      CALL INITA
      CALL INITR
      REL = .true.
      
      if(myid == 0) then
        call inp_atom(QREL,QMASS,ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW,
     :      REL,MASS,iscw,iread)
      end if

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

      CALL BPevalf(rel,skip,ec,nze,restart,onlydvd,iounew)

        close(iouhn)
        close(iouhz)
        close(iouhs)
      print *, 'Finished with the file'

*      call etime(time)
*      write(iscw,'(//A/A//A,F8.3,A//)')
*     :       ' END OF CASE',' ===========',
*     :       ' Total time was ', time(1)/60, ' minutes' 
      END
