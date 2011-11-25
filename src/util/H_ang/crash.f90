  subroutine crash()
    integer, allocatable :: i(:)
    deallocate(i)
    return
  end subroutine crash
