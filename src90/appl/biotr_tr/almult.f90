!
!     -------------------------------------------------------------
!      A L M U L T
!     -------------------------------------------------------------
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville           September 1997   *
!
      subroutine almult(npair) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use dbg_C
      use ems_C
      use state_C
      use mult_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:32:18  11/20/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use tritst_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(out) :: npair 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nmult, ierr, i, j 
      real(double) :: stat, atst 
!-----------------------------------------------
! --- allocate the /MULT/
!
      nmult = nvc(1)*nvc(2) 
      !if (ibugm /= 0) then 
         write (6, *) ' nvc(1) = ',nvc(1), ' nvc(2) = ', & 
                        nvc(2), ' nmult = ', nmult 
         write (6, *) ' qsl     allocation: nmult    = ', nmult 
         write (6, *) ' qsv     allocation: nmult    = ', nmult 
         write (6, *) ' qil     allocation: nmult    = ', nmult 
         write (6, *) ' qir     allocation: nmult    = ', nmult 
         write (6, *) ' qjvl    allocation: nmult    = ', nmult 
         write (6, *) ' qjvr    allocation: nmult    = ', nmult 
      !endif 

      allocate (sl(nmult), stat=ierr) 
      if (ierr.ne.0) call mem_fail(6,nmult,'almult::SL',ierr);
      allocate (sv(nmult), stat=ierr) 
      if (ierr.ne.0) call mem_fail(6,nmult,'almult::SV',ierr);
      allocate (il(nmult), stat=ierr) 
      if (ierr.ne.0) call mem_fail(6,nmult,'almult::IL',ierr);
      allocate (ir(nmult), stat=ierr) 
      if (ierr.ne.0) call mem_fail(6,nmult,'almult::IR',ierr);
      allocate (jvl(nmult), stat=ierr) 
      if (ierr.ne.0) call mem_fail(6,nmult,'almult::JVL',ierr);
      allocate (jvr(nmult), stat=ierr) 
      if (ierr.ne.0) call mem_fail(6,nmult,'almult::JVR',ierr);
      qsl => sl 
      qsv => sv 
      qil => il 
      qir => ir 
      qjvl => jvl 
      qjvr => jvr 
 
!
! --- determine the number of (J,J') pairs satisfying selection rules
!
      npair = 0 
      do i = 1, nvc(1) 
         do j = 1, nvc(2) 
            atst = tritst(jv1(i),jv2(j),lam+lam) 
            if (atst /= 0.0) cycle  
            npair = npair + 1 
            il(npair) = i 
            ir(npair) = j 
            jvl(npair) = jv1(i) 
            jvr(npair) = jv2(j) 
            sl(npair) = 0.0 
            sv(npair) = 0.0 
            write (6, *) ' npair = ', npair, ' jvl = ', jvl(npair), ' jvr = ', &
               jvr(npair) 
         end do 
      end do 
      return  
      end subroutine almult 
