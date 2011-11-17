      subroutine print_memory(icycle,iblock,clst_memory,ico_memory,
     :  hmx_memory,diag_hmx_memory,nblock,term,ncfg)

      logical clst_memory,ico_memory,hmx_memory,diag_hmx_memory
      integer iblock,ncfg,free_mem
      character*3 term
       if(iblock == 1)  then
         write(0,*)
         write(0,*) ' >>> MEMORY AND DISK USE: <<<'
           if (clst_memory) then  
             write(0,*) ' ALL IN MEMORY : c.lst ih ico hmx '
           end if
           if (hmx_memory.and.(.not.clst_memory).and.ico_memory) then
             write(0,*) ' For All Blocks:'
             write(0,*) '   -> IN MEMORY: ih ico hmx '
             write(0,*) '   -> ON DISK  : c.lst '
           end if 
           if (hmx_memory.and.(.not.ico_memory)) then
             write(0,*) ' For All Blocks:'
             write(0,*) '   -> IN MEMORY: ih hmx '
             write(0,*) '   -> ON DISK  : c.lst ico '
           end if
           if(.not.hmx_memory)  then
              write(0,*) ' For Each Block:'
              write(0,'(A6,A6,A10,A8,A8,A8)') 'Block','Term',
     :                 'ncfg','Memory','Disk','DVDSON'
           end if
        end if

      if (.not.hmx_memory) then
        if (diag_hmx_memory) then
          write(0,'(I6,A6,I10,A8,A8,A8)') iblock,term,ncfg,
     :    'hmx', ' c.lst ', ' MEMORY '
          write(0,'(A30)') 'ih'
          write(0,'(A30)') 'ico'
        else
          write(0,'(I6,A6,I10,A8,A8,A8)') iblock,term,ncfg,
     :    'ih', ' c.lst', ' DISK '
          write(0,'(A38)') 'hmx'
          write(0,'(A38)') 'ico'
        end if
      end if

      if ((icycle == 0).and.(iblock == nblock)) then
        write(0,*) ' >>> END MEMORY AND DISK USE<<< '
        write(0,*)
      end if
     
*      write (0,*), free_memory, ' kB remaining free in block ', iblock

      end subroutine print_memory


