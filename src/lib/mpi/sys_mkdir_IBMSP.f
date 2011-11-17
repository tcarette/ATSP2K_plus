!***********************************************************************
      subroutine sys_mkdir (dir, lendir, ierr)
      implicit none
      character(len=80) :: dir
      integer lendir, ierr, system

      lendir = len_trim(dir) 
      ierr = system ('mkdir '//dir(1:lendir) )
      if (ierr.ne.0) print*, ierr 

      end 
