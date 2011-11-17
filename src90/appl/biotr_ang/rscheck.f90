!
!     ------------------------------------------------------------------
!     R S C H E C K
!     ------------------------------------------------------------------
!
!. check if the order found by nonh is consistent with the built-in
!  ras order: we want to use the (j.leq.i) restriction adopted in savels.
!
      subroutine rscheck(nq, lq, lmx, mx) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  00:09:14  11/17/01  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: nwd = 128 
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: lmx 
      integer , intent(in) :: mx 
      integer , intent(in) :: nq(nwd) 
      integer , intent(in) :: lq(nwd) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: l, nst, i 
!-----------------------------------------------
!
      do l = 0, lmx 
         nst = 0 
         do i = 1, mx 
            if (lq(i) /= l) cycle  
            if (nq(i) < nst) then 
               write (6, *) &
                  ' Problem: orbital order found by nonh inconsistent' 
               write (6, *) '          with the built-in ras order in savels!' 
               write (6, *) ' Sorry about that. ' 
               stop  
            endif 
            write (6, *) ' l = ', l, ' i = ', i, ' lq(i) = ', lq(i), &
               ' nq(i) = ', nq(i) 
            nst = nq(i) 
         end do 
      end do 
      write (6, *) 
      write (6, *) ' rscheck report: Congratulations: RAS test satisfied.' 
      write (6, *) 
      return  
      end subroutine rscheck 
