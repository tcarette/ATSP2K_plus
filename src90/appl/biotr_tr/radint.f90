!
!     ------------------------------------------------------------
!     R A D I N T
!     ------------------------------------------------------------
!
! --- Calculates the radial overlap and transition integrals
!
      subroutine radint(im) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use debug_C
      use DBG_C
      use nel_C
      use rdint_C
      use non30_C
      use ems_C
      use ras_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:23:54  11/20/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use quadr_I 
      use grad_I 
      use grad2_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character , intent(in) :: im 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j, il, ierr 
!-----------------------------------------------
!
! --- allocate the /RDINT/
!
      write (6, *) ' maxorb = ', maxorb 
      if (ibugm /= 0) then 
         write (6, *) ' iqrl    allocation: maxorb**2 = ', maxorb**2 
         write (6, *) ' iqov    allocation: maxorb**2 = ', maxorb**2 
         write (6, *) ' iqrv    allocation: maxorb**2 = ', maxorb**2 
         write (6, *) ' transition ', im, lam 
      endif 
      allocate (rlint(maxorb,maxorb),STAT=ierr) 
      if (ierr.ne.0) call mem_fail(6,maxorb*maxorb,'radint::rlint',ierr);
      rlint = 0;
      IQRL=>RLINT;
      allocate (rvint(maxorb,maxorb),STAT=ierr) 
      if (ierr.ne.0) call mem_fail(6,maxorb*maxorb,'rvint::rlint',ierr);
      rvint = 0;
      IQRV=>RVINT;
      allocate (ovrlp(maxorb,maxorb),STAT=ierr) 
      if (ierr.ne.0) call mem_fail(6,maxorb*maxorb,'ovrlap::rlint',ierr);
      ovrlp = 0;
      IQOV=>ovrlp;
 
      do i = 1, maxorb 
         do j = 1, maxorb 
            if (itab(i,1)==0 .or. itab(j,2)==0) cycle  
            il = iabs(l(itab(i,1))-l(itab(j,2))) 
            ovrlp(i,j) = quadr(itab(i,1),itab(j,2),0) 
            if (im == 'E') then 
               rlint(i,j) = quadr(itab(i,1),itab(j,2),lam) 
!            if (.not.rel) then
               if (lam==1.and.il==1) rvint(i,j) = grad(itab(i,1),itab(j,2)) 
               if (lam==2.and.(il==0.or.il==2)) rvint(i,j) &
                     = grad2(itab(i,1),itab(j,2)) 
!            end if
            else if (im == 'M') then 
               rlint(i,j) = quadr(itab(i,1),itab(j,2),lam-1) 
            endif 
         end do 
      end do 
      return  
      end subroutine radint 
