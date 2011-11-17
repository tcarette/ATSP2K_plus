 
!     ------------------------------------------------------------------
!       R E A D W 2
!     ------------------------------------------------------------------
!
      subroutine readw2 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      USE dbg_C 
      use ndims_C
      use non30_C
      use nor_C
      use inout_C
      use param_C
      use ras_C
      use nel_C
      use fo_C
      use elt_C
      use dbg_C
      use radial_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:49:02  11/18/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use eptr_I 
      use lval_I 
      use initm2_I 
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j1, j2, ncf, i, k, iw, ik, m, j, ij,ierr 
      real(double), dimension(nod) :: pt 
      real(double) :: zz, eti, eki, azd 
      logical , dimension(nwd,2) :: found 
      logical :: check 
      character :: atom*6, term*6, el1*3 
      character, dimension(2) :: label*7 
!-----------------------------------------------
      data label/ 'Initial', 'Final  '/  
!
! --- allocate the /NEL/
!
      if (ibugm /= 0) write (6, *)&
           ' iqp  allocation: nod*maxnfo = ',nod* maxnfo 
      allocate (p(nod,maxnfo),STAT=ierr) 
      if (ierr.ne.0) call mem_fail(6,nod*maxnfo,'readw2::p',ierr);
      if (ibugm /= 0) write (6, *) ' iqn  allocation: maxnfo = ', maxnfo 
      allocate (n(maxnfo),STAT=ierr) 
      if (ierr.ne.0) call mem_fail(6,maxnfo,'readw2::n',ierr);
      if (ibugm /= 0) write (6, *) ' iql  allocation: maxnfo = ', maxnfo 
      allocate (l(maxnfo),STAT=ierr) 
      if (ierr.ne.0) call mem_fail(6,maxnfo,'readw2::l',ierr);
      if (ibugm /= 0) write (6, *) ' iqaz allocation: maxnfo = ', maxnfo 
      allocate (az(maxnfo),STAT=ierr) 
      if (ierr.ne.0) call mem_fail(6,maxnfo,'readw2::az',ierr);
      if (ibugm /= 0) write (6, *) ' iqmax allocation: maxnfo = ', maxnfo 
      allocate (max_(maxnfo),STAT=ierr) 
      if (ierr.ne.0) call mem_fail(6,maxnfo,'readw2::max_',ierr);

! associated pointers, pointers not used though
      IQP=>P;
      IQN=>N;
      IQL=>L;
      IQAZ=>AZ;
      IQMAX_=>MAX_ 
 
!
! --- read the radial wavefunctions on units iuw(1) and iuw(2)
!     and sort them according the RAS order
!
!
! --- transfer NDIMS/NON30 -> PARAM
!
      nwf = maxnfo 
!      ncfg = ncf 
!      write (6, *) ' in readw2, nwf = ', nwf, ' ncfg = ', ncfg 
!
      do i = 1, nwf 
         if (i <= iwf(1)) then 
            elrasi(i) = elras(i) 
            found(i,1) = .FALSE. 
         else 
            k = i - iwf(1) 
            elrasf(k) = elras(i) 
            found(k,2) = .FALSE. 
         endif 
      end do 
      write (6, *) ' elrasi = ', (elrasi(i),i=1,iwf(1)) 
      write (6, *) ' elrasf = ', (elrasf(i),i=1,iwf(2)) 
!
!  *****  READ THE RADIAL FUNCTIONS
!
      iw = 1 
   14 continue 
      if (iw <= iwf(1)) then 
         ik = 1 
         read (iuw(1)) atom, term, el1, m, zz, eti, eki, azd, (pt(j),j=1,m) 
         call eptr (elrasi, el1, i, j1) 
         if (j1 == 1) go to 999 
      else 
         ik = 2 
         read (iuw(2)) atom, term, el1, m, zz, eti, eki, azd, (pt(j),j=1,m) 
         call eptr (elrasf, el1, i, j2) 
         if (j2 == 1) go to 999 
      endif 
      if (.not.found(i,ik)) then 
         if (ik == 1) then 
            ij = i 
         else if (ik == 2) then 
            ij = i + iwf(1) 
         endif 
         write (6, *) label(ik), ' electron ', el1, ' P index = ', ij 
         p(:m,ij) = pt(:m) 
         n(ij) = ichar(el1(2:2)) - ichar('1') + 1 
         l(ij) = lval(el1(3:3)) 
         az(ij) = azd 
         max_(ij) = m 
         p(m+1:no,ij) = d0 
         found(i,ik) = .TRUE. 
      endif 
      z = zz 
      iw = iw + 1 
  999 continue 
      if (iw <= nwf) go to 14 
      close(iuw(1)) 
      close(iuw(2)) 
!
! ---  check that all functions were found
!
      check = .TRUE. 
      check = check .and. all(found(:iwf(1),1)) 
      check = check .and. all(found(:iwf(2),2)) 
      if (.not.check) then 
         write (iscw, *) ' Not all radial functions were found' 
         write (iscw, '(1X,A3,2X,L1)') (elrasi(i),found(i,1),i=1,iwf(1)) 
         write (iscw, '(1X,A3,2X,L1)') (elrasf(i),found(i,2),i=1,iwf(2)) 
         stop  
      endif 
      do i = 1, maxorb 
         do j = 1, iwf(1) 
            if (eltens(i) /= elrasi(j)) cycle  
            itab(i,1) = j 
         end do 
         do j = 1, iwf(2) 
            if (eltens(i) /= elrasf(j)) cycle  
            itab(i,2) = j + iwf(1) 
         end do 
      end do 
      write (6, *) 
      write (6, *) 
      write (6, *) ' Position of biorthonormal shells in P vector : ' 
      write (6, *) ' ---------------------------------------------- ' 
      write (6, *) label(1), ' :  ', (itab(i,1),i=1,maxorb) 
      write (6, *) label(2), ' :  ', (itab(i,2),i=1,maxorb) 
      write (6, *) 
!     .. initialize R and Rydberg constant
      call initm2 (1) 
      return  
      end subroutine readw2 
