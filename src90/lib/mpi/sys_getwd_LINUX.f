!***********************************************************************
      subroutine sys_getwd (dir, dir_len)
      implicit none
      character(len=*), intent(out):: dir
      integer, intent(out):: dir_len 
      character(len=*), parameter::iam = 'SUN'  ! soly for machine
c      character(len=len_trim (iam)), optional, intent(out):: machine

!   This routine gets current working directory and assigns it to dir.
!   ierr will be zero if successful, otherwise non-zero;
!   machine is an optional parameter specifying the name of the system
!
!   Xinghong He  98-08-24
!
!***********************************************************************

      dir = ''
      call GETCWD(dir)
      dir_len = len_trim(dir) 

      return
      end
