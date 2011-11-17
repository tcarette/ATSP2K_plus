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
!   MPI data
!-----------------------------------------------
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
      character*4 idstring
      character*128 startdir, permdir, tmpdir
      integer lclst,lentmp,lenstart,lenperm
! >>>>>>     <<<<<<<

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
            if (myid == 0) then
               write (6, *) &
                  ' Problem: orbital order found by nonh inconsistent' 
               write (6, *) '          with the built-in ras order in savels!' 
               write (6, *) ' Sorry about that. ' 
               stop  
            end if
            endif 
            if (myid == 0) then
            write (6, *) ' l = ', l, ' i = ', i, ' lq(i) = ', lq(i), &
               ' nq(i) = ', nq(i) 
            end if
            nst = nq(i) 
         end do 
      end do 
      if (myid == 0) then
         write (6, *) 
         write (6, *) &
           ' rscheck report: Congratulations: RAS test satisfied.' 
         write (6, *) 
      end if
      return  
      end subroutine rscheck 
