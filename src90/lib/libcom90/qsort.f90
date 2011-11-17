!
!-----------------------------------------------------------------------
!        Q S O R T
!-----------------------------------------------------------------------
 
!     The method use to sort the data is quick sort with a pivot value
!     that is the larger value of the first 2 different value from the
!     the sublist to be sorted.
!     This sorting method used a stack to maintain the unsorted section,
!     and sorting will be finished when the stack is empty.
      subroutine qsort(n, key, pt, stack, nstack, ierr) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      integer , intent(in) :: nstack 
      integer , intent(out) :: ierr 
      integer , intent(in) :: key(*) 
      integer , intent(inout) :: pt(*) 
      integer , intent(inout) :: stack(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: top, from, to, pivot, left, right, i, k, ismall, kp 
!-----------------------------------------------
!
!     Set the initial pointer values
!
      do i = 1, n 
         pt(i) = i 
      end do 
!
!     Initialize the stack and error indicator
!
      top = 1 
      stack(top) = 1 
      ierr = 0 
!
!     Repeat Until the Stack is empty
!
  100 continue 
      if (top /= 0) then 
         from = stack(top) 
         if (top /= 1) then 
            to = stack(top-1) - 1 
         else 
            to = n 
         endif 
         top = top - 1 
!
!        Find the position k of the pivot value that partitions
!        the current section. Return a value of k=0 when there is
!        no disinct value.
!
         if (from == to) then 
            k = 0 
         else 
            k = 0 
            ismall = key(pt(from)) 
            do i = from + 1, to 
               if (key(pt(i)) == ismall) cycle  
               k = i 
               exit  
            end do 
         endif 
         if (k /= 0) then 
            if (ismall > key(pt(k))) k = from 
         endif 
         if (k /= 0) then 
!
!         Rearange the section of the keys such that all values
!         smaller than the pivot value are stored on the left and
!         larger values on the right.
!
            pivot = key(pt(k)) 
            left = from 
            right = to 
  300       continue 
  310       continue 
            if (key(pt(left)) < pivot) then 
               left = left + 1 
               go to 310 
            endif 
  320       continue 
            if (key(pt(right)) >= pivot) then 
               right = right - 1 
               go to 320 
            endif 
            if (left < right) then 
               i = pt(left) 
               pt(left) = pt(right) 
               pt(right) = i 
               go to 300 
            endif 
            kp = left 
            if (top + 2 <= nstack) then 
               stack(top+1) = kp 
               stack(top+2) = from 
               top = top + 2 
            else 
               ierr = 1 
               return  
            endif 
         endif 
         go to 100 
      endif 
!
!     Keys are sorted
!
      return  
      end subroutine qsort 
