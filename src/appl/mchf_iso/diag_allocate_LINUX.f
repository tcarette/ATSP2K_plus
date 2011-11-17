      subroutine diag_allocate(pointer,size_block,size_column,
     :            logical_success,size)

      pointer(pointer,array(1))
      logical logical_success;
      integer size_block, iblock, size_column,ierr_mem,size;

      logical_success = .true.

cgd   pointer=malloc(size_block*size);
      call alloc(pointer,size_block,size);
      if (pointer == 0) then
        call free(pointer)
cgd     pointer = malloc(size_column*size);
        call alloc(pointer,size_column,size);
        if (pointer.ne.0) logical_success = .false.
        if (pointer == 0) then
            write(0,*) 'Insufficient Memory: Stop in diag_alloc_memory'
        end if
      end if

      end subroutine diag_allocate
