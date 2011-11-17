*     ==================================================================
*
*       A PROGRAM TO FIND SELECTED EIGENVALUES AND EIGENVECTORS 
*       GIVEN FILES FOR THE INTERACTION MATRIX OF NON-RELATIVISTIC
*             THE BREIT-PAULI HAMILTONIANS
*
*       by C. Froese Fischer and G. Tachiev
*          Vanderbilt University
*          Nashville, TN 37235 USA
*          May, 2000
*
*                C O P Y R I G H T   2000
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
      PROGRAM bp_eiv 
*	WITH RESTARTING
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220, NWD=60, ntermd = 31)
*
      character*1 idstring
      character*128 NAME(2), startdir, permdir, tmpdir
      character*128 file_cint,file_ih,file_ico,file_clst
      character*128 file_atomc,file_atoml,file_atomj,file_atomw
      character*128 file_atomnew,file_iuhnr,file_iouhz,file_iouhs
      integer lenl,lenw,lenc,lenj,lenn

      CHARACTER*26 STRING,ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW
      CHARACTER QREL, QMASS
      LOGICAL REL, SKIP
      logical onlydvd,jhere,lhere,mhere,zhere,shere
      logical leigen(NTERMD,NTERMD)
      integer nume(ntermd),njeig(ntermd),termsh(ntermd),indx(ntermd)
      pointer(pflsj,nze_flsj(1,1,1,1));
      logical lhmx_memory, lhmx_disk
      common/hmx_nze/nze_t,nze_c
      double precision flsj;

      COMMON /NCFGS/ IQLSP,index,IJK,termsh,nterm
*
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,iout,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
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

      myid = 0; nprocs = 1

      IREAD = 5
      iuc    = 4
      iuw    = 2
      ioul   = 7
      iouj   = 8
      ionew = 13;

      iouhn  = 9
      iouhz  = 10
      iouhs  = 11
      iouhm  = 12
      iout   = 14

         ISCW = 0
         IWRITE = 6
         write(iscw,'(A/A)')
         write(iscw,*) '             ...bp_E 2000 '
         write(iscw,*)

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

      if (myid == 0) then
         call inp_case(QREL,QMASS,ATOMC,ATOML,ATOMJ,ATOMW,
     :      ATOMNEW,REL,MASS,iscw,iread)
         call inp_eig(njv,MAXJ,MINJ,leigen,nume,iscw,iread,onlydvd,
     :          termsh,indx,maxnum)
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
        OPEN(UNIT=iuw,FILE=file_atomw,STATUS='OLD',
     :                FORM='UNFORMATTED' )
        OPEN(UNIT=ioul,FILE=file_atoml,STATUS='UNKNOWN',
     :                FORM='UNFORMATTED')
!        OPEN(UNIT=iounew,FILE=file_atomnew,STATUS='UNKNOWN',
!     :                FORM='UNFORMATTED')
      end if

      OPEN(UNIT=iuc,FILE=file_atomc,STATUS='OLD')
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

      call inp_data(rel,skip,ec,nze,njv,minj,maxj,pflsj)

      if(.false.) call term_shift(iscw,nterm,ncfg,indx,lsp,termsh);

      if (skip) then
	iLS = 1
	close(unit=iouhz,status='delete')
	close(unit=iouhs,status='delete')
	if (idisk .eq. 1) iouhm = iouhn
      else
	iLS = 0
	if (idisk .eq. 1)
     :         open(iouhm,status='scratch',form='unformatted')
      end if

      ij = 1
      do j = maxj, minj,-2
        call lsjmat(j,nume(ij),skip,nclosd,nwf,ncfg,nze,ec,
     :    onlydvd,jhere,lhere,leigen(1,ij),pflsj,njv,maxj)
        ij = ij+1
      end do

      if (myid == 0 ) then
      onlydvd = .false.
!      print *, 'onlydvd', onlydvd
      if (onlydvd) then
	print *, 'iounew =',iounew
CGG        write(iounew,*) '**'
        write(iounew,*) 'END '
	close(iounew)
      else
!	print *, 'ILS =', ILS
        if (iLS .eq. 1) then
CGG         write(ioul,*) '**'
          write(ioul,*) 'END '
	 close(ioul)
        else
CGG          write(iouj,*) '**'
          write(iouj,*) 'END '
	  close(iouj)
        end if
      end if
      print *, 'Finished with the file'
      end if

*      call etime(time)
*      write(iscw,'(//A/A//A,F8.3,A//)')
*     :       ' END OF CASE',' ===========',
*     :       ' Total time was ', time(1)/60, ' minutes' 
      END
