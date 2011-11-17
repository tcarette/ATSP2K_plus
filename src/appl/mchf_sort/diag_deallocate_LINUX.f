      subroutine diag_deallocate(pointer)
      pointer(pointer,array(1))
      integer ierr_mem;

      call free(pointer)
      pointer = 0

      end subroutine diag_deallocate

