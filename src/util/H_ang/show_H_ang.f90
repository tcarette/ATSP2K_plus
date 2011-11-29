!
! Read the lists produced by nonh and show the angularly integrated 
! Hamiltonian to a human
!
! Written by Thomas Carette
!                   November, 2011
!
!

      Program show_H_ang

      use safeio
      use mdata

      IMPLICIT NONE

!**********************************************************************
!
!     IO
!

      CHARACTER*16 INPUT
      INTEGER      fp
      LOGICAL yclist

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


      call data


      call unsparse(.true.,.false.)


      end program show_H_ang
