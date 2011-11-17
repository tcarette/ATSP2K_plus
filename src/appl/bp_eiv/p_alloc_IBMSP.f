      subroutine p_alloc(pointer,size_block,size)

      pointer(pointer,array(1))
      integer size_block, iblock, ierr_mem,size;

      pointer=malloc(%val(size_block*size));

      end subroutine p_alloc 


      subroutine p_dalloc(pointer,sz)
      pointer(pointer,array(1))
      integer sz
      if (pointer == 0) then
        call free(%val(pointer));
        call disclaim(%val(pointer));
        pointer = 0;
      end if
      end subroutine p_dalloc 
