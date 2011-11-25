!
! Read the lists produced by nonh and show the angularly integrated 
! Hamiltonian to a human
!
! Written by Thomas Carette
!                   November, 2011
!
!

      Program show_H_ang
      IMPLICIT NONE

      INTEGER, PARAMETER :: NWCD=20,NBD=20,LSDIM=30000,NWD=70

!**********************************************************************
!
!     IO
!

      INTEGER      IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,
     :             IALL,JSC(3),ISCW, state
      CHARACTER*16 INPUT
      LOGICAL      yclist
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,
     :              IALL,JSC,ISCW, state

!**********************************************************************
!
!    CSF list variables
!
      INTEGER      MAXORB,NB,LMAX,ncfg,nclosd
      INTEGER      NBsize(NBD)
      CHARACTER    EL(NWD)*3
      character*2  ih_file

      INTEGER jan,nijptr,jbn,ico
      DOUBLE PRECISION cn
      POINTER (qnijptr,nijptr(lsdim)),(qjan,jan(lsdim)),
     :        (qjbn,jbn(lsdim)),(qico,ico(lsdim)) 

!**********************************************************************
!
!     orbital variables
!

      INTEGER IAJCMP,LJCOMP,NJCOMP,IAJCLD,LJCLSD
      POINTER(QIAJCMP,IAJCMP(nwcd)),(QLJCOMP,LJCOMP(nwcd)),
     :       (QNJCOMP,NJCOMP(nwcd)),(QIAJCLD,IAJCLD(nwcd)),
     :       (QLJCLSD,LJCLSD(nwcd))
      COMMON /buffer/qcn,qinptr,qpackn,qlused,qintptr,lmax,qnijptr,
     :               qjan,qjbn,qico
!**********************************************************************
!
!     integral variables
!

      INTEGER noint(4)

      INTEGER l,intptr,ipackn,idummy,inptr
      LOGICAL lused
      POINTER (ql,l(lsdim)),(qintptr,idummy(lsdim)),
     :        (qinptr,inptr(lsdim)),(qcn,cn(lsdim)),
     :        (qpackn,ipackn(lsdim)),(qlused,lused(lsdim))

!**********************************************************************
!
!     dummy variables
!

      INTEGER i,n

!**********************************************************************


!**********************************************************************

      INPUT = 'cfg.inp'

      i = iargc()
      if (i .eq. 0) then
       INPUT = 'cfg.inp'
       inquire( FILE=input, exist=yclist) 
       if (yclist) then 
          print *, 'input file is cfg.inp ...'
       else 
          print *, 'cfg.inp not found: nonh is exiting!...'
          call exit(0)
        endif	
      end if

      IREAD=9
      IWRITE=6
      IOUT=8

      OPEN(UNIT=IREAD,FILE=INPUT,STATUS='UNKNOWN')

      OPEN(UNIT=IOUT, FILE='yint.lst',STATUS='UNKNOWN',
     :     FORM='unformatted')
      OPEN(UNIT=12, FILE='ico.lst',STATUS='UNKNOWN',
     :     FORM='unformatted')
      OPEN(UNIT=50,FILE='c.lst',STATUS='UNKNOWN',
     :     FORM='UNFORMATTED')


!
!  ---  Determine input data; non-orthogonal case
!

      call analy_blk(NCLOSD,MAXORB,NB,NBsize,EL)

      call orbitals(maxorb,el,qiajcmp,qljcomp,
     :                    qnjcomp,qiajcld,qljclsd,nb,nclosd)

!  .. find maximum l-value
      lmax = 0
      do i=1,maxorb
         lmax = max (lmax,ljcomp(i))
      end do

!  .. allocate memory for buffered i/o
      call alloc(qcn,lsdim,8)
      call alloc(qinptr,lsdim,4)
      call alloc(qnijptr,lsdim,4)
      call alloc(qjan,lsdim,4)
      call alloc(qjbn,lsdim,4)
      call alloc(qico,lsdim,4)

!  .. generate list of integrals
      call genint(maxorb,lmax,qljcomp,qintptr,qpackn,qlused,noint)

!  .. initialize lused
      do i=1,noint(4)
       lused(i) = .FALSE.
      end do

!  .. for all blocks
      rewind (iread)
      do N=1,NB
        write(ih_file,'(I2.2)') N 
        OPEN(UNIT=11, FILE='ih.'//ih_file//'.lst',STATUS='UNKNOWN',
     :       FORM='unformatted')

        ncfg = NBsize(N)

      enddo

      end program show_H_ang
