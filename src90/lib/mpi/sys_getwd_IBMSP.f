      subroutine sys_getwd (dir,len_cwd) 
      CHARACTER (len=*), intent(out):: dir
      INTEGER GETWD, len_cwd

      call GETCWD(dir)
      len_cwd = len_trim(dir) 
      if (len_cwd.le. 0) then
         print *, 'Can''t get the current directory: exit!'
         call exit(1);
      endif

      end 
        
