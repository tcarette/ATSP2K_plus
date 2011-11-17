!***********************************************************************
c      subroutine sys_mkdir (dir, lendir, ierr, machine)
      subroutine sys_mkdir (dir, lendir, ierr)
      implicit none
      character(len=*), intent(in):: dir
      integer, intent(in):: lendir
      integer, intent(out):: ierr
c      character(len=*), parameter::iam = 'T3E'  ! soly for machine
c      character(len=len_trim (iam)), optional, intent(out):: machine

!   This routine makes a sub-dir under the current working directory.
!   lendir is the length of character string dir; 
!   ierr will be zero if successful, otherwise non-zero;
!   machine is an optional parameter specifying the name of the system
!
!   Xinghong He  98-08-21
!
!***********************************************************************

      integer ishell

c      machine = iam

      ierr = ishell ('mkdir -p -m 775 ' // dir(1:lendir) )

      return
      end
