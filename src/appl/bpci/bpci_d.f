*     ==================================================================
*
*       A CONFIGURATION INTERACTION PROGRAM EITHER NON-RELATIVISTIC
*             OR IN THE BREIT-PAULI APPROXIMATION
*
*       by C. Froese Fischer
*          Vanderbilt University
*          Nashville, TN 37235 USA
*          May, 1983
*
*       Modified August, 1992 to combine BREIT and the sparse
*       matrix Dvdson eigensolver.

*       Modified by Gediminas Gaigalas for new angular structure 
*       Modified for Unsorted version
*       Modified for restarting [Sep 1994]
*       Modified by Gediminas for new angular libraries and orbit-orbit (1997)
*
*              C O P Y R I G H T   2000
*     ==================================================================
*
*       The PARAMETER values in this program define the following:
*               NOD   - Maximum number of points in the range of a
*                         - function
*               LIMD  - Maximum number of basis vectors
*               LCDIM - Number of array elements in a block
*                       of the direct access file
*
*       Variables for these dimensions are
*               NWF   - Number of functions (or electrons)
*               NO    - Number of points in the range of the function
*               NCFG  - Number of configuration state functions
*               NUME  - Number of eigenvalues/eigenvectors
*               NZE   - Number of non-zero matrix elements
*               leigen- logical array ntermdXntermd for requested eigenvalues
*     ------------------------------------------------------------------
      PROGRAM BPCI
*	WITH RESTARTING
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220, NWD=60, NTERMD=31)
*
      CHARACTER*26 STRING,ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW
      CHARACTER QREL, QMASS,str*72
      LOGICAL REL, SKIP
      logical onlydvd,jhere,lhere,mhere,zhere,shere,restart
      logical leigen(NTERMD,NTERMD)
      integer nume(NTERMD)
*     .. set big to be 2**28
      integer big
      data big/268435456/
*
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,
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
*               iounew- The unit number for the <name>.new file
*               iouhn  - The unit number for storing the LS matrix
*               iouhz - The unit number for storing the ZNV data
*               iouhs - The unit number for storing the S data
*               iouhm - The unit number for writing the scratch matrix
*
      ISCW=0
      IREAD = 5
      IWRITE = 6
      in = 5
      iuc    = 4
      iuw    = 2
      ioul   = 7
      iouj   = 8
      iounew = 13

      iouhn  = 9
      iouhz  = 10
      iouhs  = 11
      iouhm  = 12
 
      CALL INITA
      CALL INITR
      REL = .true.
1     WRITE(ISCW,'(/A,A)') ' Enter ATOM, relativistic (Y/N)',
     :               ' with mass correction (Y/N)'
      READ(IREAD,'(A)') STRING
      I = INDEX(STRING,',')
      IF (I .eq. 0)     THEN
         WRITE(ISCW,*) ' Separate with commas'
         GO TO 1
      ELSE IF (I .GT. 22) THEN
         WRITE(ISCW,*) ' ATOM name can have at most 22 characters'
         GO TO 1
      END IF
      QREL  = STRING(I+1:I+1)
      QMASS = STRING(I+3:I+3)
      ATOMC = STRING(1:I-1)//'.c'
      ATOML = STRING(1:I-1)//'.l'
      ATOMJ = STRING(1:I-1)//'.j'
      ATOMW = STRING(1:I-1)//'.w'
      ATOMNEW = STRING(1:I-1)//'.new'
      IF (QREL.EQ.'N' .OR. QREL.EQ.'n') REL = .FALSE.
      MASS = 0
      IF (QMASS.EQ.'Y' .OR. QMASS.EQ.'y') THEN
         WRITE(ISCW,'(A)') ' Gradient or Slater form? (G/S):'
         READ(IREAD,'(A1)') QMASS
         MASS = 1
         IF (QMASS.EQ.'S' .OR. QMASS.EQ.'s') MASS = 2
      END IF

*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*       ... this is a disk version.  No counting of nze
        idisk = 1
	restart = .false.
	onlydvd = .false.
	WRITE(ISCW,*) ' Restarting (Y/y) ?'
        READ(IREAD,'(A)') STRING
        IF(STRING(1:1).EQ.'Y'.OR.STRING(1:1).EQ.'y') THEN
	  restart = .true.
        ELSE
          WRITE(ISCW,*) ' Use existing Matrix and ',
     :                  '<atom>.l/j initial guess (Y/y)?'
          READ(IREAD,'(A)') STRING
          if (STRING(1:1).eq.'Y'.or.STRING(1:1).eq.'y') then
            onlydvd = .true.
            OPEN(UNIT=iounew,FILE=ATOMNEW,STATUS='UNKNOWN')
            write(ISCW,*) ' The file <>.new will contain',
     :          ' the new <>.l or <>.j info.'
          endif
	ENDIF
*        * * * * * * Test if matrix exists * * * * * * * * * * *
        inquire(file='hnr.lst',EXIST=mhere)
        inquire(file=ATOML,EXIST=lhere)
        inquire(file=ATOMJ,EXIST=jhere)
        if (onlydvd .or. restart) then
           if (.not.mhere) then
	     WRITE(ISCW,*) 'hnr.lst file not found'
	     stop
	   end if
	   if (rel) then
             inquire(file='hzeta.lst',EXIST=zhere)
	     if (.not. zhere) then
	       WRITE(ISCW,*) 'hzeta.lst file not found'
	       stop
	     end if
             inquire(file='hspin.lst',EXIST=shere)
	     if (.not. shere) then
	       WRITE(ISCW,*) 'hzeta.lst file not found'
	       stop
	     end if
           endif
        endif

*     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

*
*    Units iuc and iuw will be CLOSED by BREVAL
*
      OPEN(UNIT=iuc,FILE=ATOMC,STATUS='OLD')
      OPEN(UNIT=iuw,FILE=ATOMW,STATUS='OLD',FORM='UNFORMATTED' )
      OPEN(UNIT=iouhn,FILE='hnr.lst',STATUS='unknown',
     :        form='unformatted')
      OPEN(UNIT=ioul,FILE=ATOML,STATUS='UNKNOWN')
      if (rel) then
         OPEN(UNIT=iouhz,FILE='hzeta.lst',STATUS='unknown',
     :        form='unformatted')
         OPEN(UNIT=iouhs,FILE='hspin.lst',STATUS='unknown',
     :        form='unformatted')
         OPEN(UNIT=iouj,FILE=ATOMJ,STATUS='UNKNOWN')
* 
      end if
      maxj = 0
      minj = 0
      If (rel) then
90       WRITE(ISCW,'(/A)')'  Enter Maximum and minimum values of 2*J'
         READ(IREAD,*) MAXJ,MINJ
         njv = (maxj-minj)/2 + 1
         if (njv .gt. 20) then
	   write(iscw, *)
     :       ' Current dimensions only allow upto 20 J-values: Re-enter'
	   go to 90
	 end if
      end if
      ij = 1
      maxnum = 1
*  determine eigenvalues
**      do 94 j = maxj,minj,-2
**        Write(iscw,'(A,I3)')
**     :      '  Enter Number of eigenvalues for 2*J =',j
**        Read(iread,*) nume
**        njeig(ij) = nume
**	ij = ij + 1
**        if (nume .gt. maxnum) maxnum = nume
**   94 continue
       
      nb = 1
      leigen = .false.
      write(iscw,*)
      write(iscw,*) 'Enter eigenvalues: ',
     :             'one line per term, eigenvalues separated by commas'
      do nj = maxj,minj,-2
        Write(iscw,'(A,I3)')
     :      '2*J =',nj
         read(in,'(A)') str
         nch = 1
         len = len_trim(str)
         do while (nch <= len)
            ipos = index(str(nch:len),',')
            if (ipos .eq. 0) ipos = len+2-nch
            read (str(nch:nch+ipos-2),*) keigv
            if (keigv .gt. ntermd) then
              write(0,*) 'Too high an eigenvalue requested:',
     :                   'Maximum for current dimensions is',ntermd
              stop
            end if
            leigen(keigv,nb) = .true.
            nch = nch + ipos
         end do
         nume(nb) = keigv
         njeig(nb) = nume(nb)
         if (nume(nb).gt. maxnum) maxnum = nume(nb)
         nb = nb + 1
      end do

*
      CALL BREVAL(rel,skip,ec,nze,restart,onlydvd,iounew,idisk)
      if (.not.onlydvd) then
        endfile(iouhn)
        endfile(iouhz)
        endfile(iouhs)
      else
        if (idisk .eq. 0) then
           nze = 0
           do i=1,NCFG
              read(iouhn) jb,m
              if (nze + m .gt. big) then
                 nze = ncfg
                 idisk = 1
              else
                nze = nze+m
              end if
           enddo
         end if
      endif
   
      rewind(iouhn)
      rewind(iouhz)
      rewind(iouhs)
*     .. allocate the memory
      call alcmat(ncfg,nze,maxnum)
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
      do 100 j = maxj, minj,-2
	print *, 'J =', j, njeig(ij)
*	nume = njeig(ij)
	call lsjmat(j,nume(ij),skip,nclosd,nwf,ncfg,nze,ec,
     :              onlydvd,jhere,lhere,iounew,leigen(1,ij))
	ij = ij+1
100   continue
      print *, 'onlydvd', onlydvd
      if (onlydvd) then
	print *, 'iounew =',iounew
        write(iounew,*) '**'
	close(iounew)
      else
	print *, 'ILS =', ILS
        if (iLS .eq. 1) then
         write(ioul,*) '**'
	 close(ioul)
        else
          write(iouj,*) '**'
	  close(iouj)
        end if
      end if
      print *, 'Finished with the file'

******************************************************
*      call etime(time)
*      write(iscw,'(//A/A//A,F8.3,A//)')
*     :       ' END OF CASE',' ===========',
*     :       ' Total time was ', time(1)/60, ' minutes' 
      END
