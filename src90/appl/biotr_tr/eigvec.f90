!
!     -------------------------------------------------------------
!      E I G V E C
!     -------------------------------------------------------------
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville           September 1997   *
!
      subroutine eigvec(ici, configi, configf, lsj) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
      use dbg_C
      use inout_C
      use ndims_C
      use ems_C
      use nor_C
      use state_C
      use param_C
      use fo_C
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:08:05  11/18/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use rdegvc_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ici 
      character  :: configi*64 
      character  :: configf*64 
      integer , intent(in) :: lsj(2) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(8) :: q 
      integer :: lwf, k, j, nocc, jd, ja, jb, jc, mc, iu, k1, k2, &
                 k3, jv, nb, m , ierr
      real(double) :: max_wt1, max_wt2 
      character :: ligne*80, line*70 
      character , dimension(8) :: elc*3 
      character, dimension(15) :: couple*3 
!-----------------------------------------------
!


    1 format(15x,f14.7,1x,a) 
   14 format(/,/,8x,i4,10x,i4) 
   15 format(/,i6,f16.9,8a8,/(7f11.8)) 
      if (ici == 0) then 
         nvc(1) = 1 
         nvc(2) = 1 
         lgth(1) = ncf(1) 
         lgth(2) = ncf(2) 
      else 
         call rdegvc 
      endif 
      if (ibugm /= 0) then 
!
! --- allocate the /STATE/
!
         write (6, *) ' nvc(1) = ', nvc(1), ' nvc(2) = ', nvc(2) 
         write (6, *) ' lgth(1)= ', lgth(1), ' lgth(2)= ', lgth(2) 
         write (6, *) ' qet1    allocation: nvc(1)   = ', nvc(1) 
      endif 
      if (ibugm ==0) ibugm = 99;
      if (ibugm /= 0) write (6, *) ' qlbl1   allocation: nvc(1)   = ', nvc(1) 
      if (ibugm /= 0) write (6, *) ' qwt1    allocation: lgth(1)  = ', lgth(1) 
      if (ibugm /= 0) write (6, *) ' qjv1    allocation: nvc(1)   = ', nvc(1) 
      if (ibugm /= 0) write (6, *) ' qcfg1   allocation: 8*nvc(1) = ', 8*nvc(1) 
      if (ibugm /= 0) write (6, *) ' qet2    allocation: nvc(2)   = ', nvc(2) 
      if (ibugm /= 0) write (6, *) ' qlbl2   allocation: nvc(2)   = ', nvc(2) 
      if (ibugm /= 0) write (6, *) ' qwt2    allocation: lgth(2)  = ', lgth(2) 
      if (ibugm /= 0) write (6, *) ' qjv2    allocation: nvc(2)   = ', nvc(2) 
      if (ibugm /= 0) write (6, *) ' qcfg2   allocation: 8*nvc(2) = ', 8*nvc(2) 
      if (ibugm==99) ibugm = 0;

      allocate(et1(nvc(1)),stat=ierr);
      if (ierr.ne.0) call mem_fail(iscw,nvc(1),'eigvec::et1',ierr);
      allocate(lbl1(nvc(1)),stat=ierr);
      if (ierr.ne.0) call mem_fail(iscw,nvc(1),'eigvec::lbl1',ierr);
      allocate(wt1(lgth(1)),stat=ierr);
      if (ierr.ne.0) call mem_fail(iscw,lgth(1),'eigvec::wt1',ierr);
      allocate(jv1(nvc(1)),stat=ierr);
      if (ierr.ne.0) call mem_fail(iscw,nvc(1),'eigvec::jv1',ierr);
      allocate(cfg1(8,nvc(1)),stat=ierr);
      if (ierr.ne.0) call mem_fail(iscw,8*nvc(1),'eigvec::cfg1',ierr);
      allocate(et2(nvc(2)),stat=ierr);
      if (ierr.ne.0) call mem_fail(iscw,nvc(2),'eigvec::et2',ierr);
      allocate(lbl2(nvc(2)),stat=ierr);
      if (ierr.ne.0) call mem_fail(iscw,nvc(2),'eigvec::lbl2',ierr);
      allocate(wt2(lgth(2)),stat=ierr);
      if (ierr.ne.0) call mem_fail(iscw,lgth(2),'eigvec::wt2',ierr);
      allocate(jv2(nvc(2)),stat=ierr);
      if (ierr.ne.0) call mem_fail(iscw,nvc(2),'eigvec::jv2',ierr);
      allocate(cfg2(8,nvc(2)),stat=ierr);
      if (ierr.ne.0) call mem_fail(iscw,8*nvc(2),'eigvec::cfg2',ierr);

!associate pointers, pointers not used
      qlbl1=>lbl1;
      qwt1=>wt1;
      qjv1=>jv1;
      qcfg1=>cfg1;
      qet2=>et2;
      qlbl2=>lbl2;
      qwt2=>wt2;
      qjv2=>jv2;
      qcfg2=>cfg2;

!
! --- read in all the eigenvectors of the initial state
!
! ja = seniority; jb = 2*L+1 ; jc = 2*S+1
      if (ici == 0) then 
         read (iuc(1), 1) et1(1), ligne(1:36) 
         read (iuc(1), '(A)') configi 
!           .. closed shells
         lwf = iwf(1) - nclos(1) 
         if (ibugm /= 0) write (6, *) ' lwf = ', lwf 
    2    continue 
         read (iuc(1), '(A)') ligne 
!           .. orbitals (could have 20)
         lwf = lwf - 20 
         if (lwf > 0) go to 2 
 
         max_wt1 = 0.0D0 
         do k = 1, ncf(1) 
!          if(k.eq.1) then
            read (iuc(1), '(8(1X,A3,1X,I2,1X),F15.12)') (elc(j),q(j),j=1,8), &
               wt1(k) 
            read (iuc(1), '(15(1X,A3))') (couple(j),j=1,15) 
 
            if (abs(wt1(k)) <= max_wt1) cycle  
            max_wt1 = abs(wt1(k)) 
            nocc = 0 
    4       continue 
            if (elc(nocc+1) /= '   ') then 
               nocc = nocc + 1 
               if (nocc < 8) go to 4 
            endif 
            call pack (nocc, elc, q, couple, configi) 
 
!          else
!            read(iuc(1),'(t65,f15.12)') wt1(k)
!            read(iuc(1),'(A)')
!          end if
 
         end do 
         jd = j1qnrd(2*noccsh(1) - 1,1) 
         ja = mod(jd,64) 
         jd = jd/64 
         jb = mod(jd,64) 
         jc = jd/64 
         jv1(1) = jb + jc - 2 
         write (6, *) ' seniority = ', ja, ' (2L+1) = ', jb, ' (2S+1) = ', jc, &
            ' 2*J = ', jv1(1) 
         read (iuc(2), 1) et2(1), ligne(1:36) 
         read (iuc(2), '(A)') configf 
!           .. closed shells
         lwf = iwf(2) - nclos(2) 
         if (ibugm /= 0) write (6, *) ' lwf = ', lwf 
         write (6, *) ' lwf = ', lwf 
    5    continue 
         read (iuc(2), '(A)') ligne 
!           .. orbitals (could have 20)
         lwf = lwf - 20 
         if (lwf > 0) go to 5 
 
         max_wt2 = 0.0D0 
         do k = 1, ncf(2) 
!          if(k.eq.1) then
 
            read (iuc(2), '(8(1X,A3,1X,I2,1X),F15.12)') (elc(j),q(j),j=1,8), &
               wt2(k) 
            read (iuc(2), '(15(1X,A3))') (couple(j),j=1,15) 
            if (abs(wt2(k)) <= max_wt2) cycle  
            max_wt2 = abs(wt2(k)) 
            nocc = 0 
    7       continue 
            if (elc(nocc+1) /= '   ') then 
               nocc = nocc + 1 
               if (nocc < 8) go to 7 
            endif 
            call pack (nocc, elc, q, couple, configf) 
 
!          else
!            read(iuc(2),'(t65,f15.12)') wt2(k)
!            read(iuc(2),'(A)')
!          end if
         end do 
!GG        mc = mcfg+1
         mc = ncf(1) + 1 
         jd = j1qnrd(2*noccsh(mc) - 1,mc) 
         ja = mod(jd,64) 
         jd = jd/64 
         jb = mod(jd,64) 
         jc = jd/64 
         jv2(1) = jb + jc - 2 
         write (6, *) ' seniority = ', ja, ' (2L+1) = ', jb, ' (2S+1) = ', jc, &
            ' 2*J = ', jv2(1) 
      else 
         if (rel) then 
            iu = iuj(1) 
         else 
            iu = iul(1) 
         endif 
         k1 = 1 
         k2 = ncf(1) 
         k3 = 0 
         rewind (iu) 
         read (iu, '(A)') line 
    8    continue 
         read (iu, 14, end=10) jv, nb 
         if (ibugm /= 0) write (6, *) '     jv = ', jv, ' nb = ', nb 
         if (.not.rel) jv = lsj(1) 
         do m = 1, nb 
            k3 = k3 + 1 
            jv1(k3) = jv 
            read (iu, '(/I6,F16.9,2X,8A8)') lbl1(k3), et1(k3), (cfg1(k,k3),k=1,&
               8) 
            read (iu, '(7F11.8)') (wt1(k),k=k1,k2) 
            k1 = k1 + ncf(1) 
            k2 = k2 + ncf(1) 
         end do 
         go to 8 
!
! --- read in all the eigenvectors of the final state
!
   10    continue 
         close(iu) 
         if (ibugm /= 0) write (6, *) ' jv1 = ', (jv1(m),m=1,nvc(1)) 
         if (rel) then 
            iu = iuj(2) 
         else 
            iu = iul(2) 
         endif 
         k1 = 1 
         k2 = ncf(2) 
         k3 = 0 
         read (iu, '(A)') line 
   11    continue 
         read (iu, 14, end=13) jv, nb 
         if (ibugm /= 0) write (6, *) '     jv = ', jv, ' nb = ', nb 
         if (.not.rel) jv = lsj(2) 
         do m = 1, nb 
            k3 = k3 + 1 
            jv2(k3) = jv 
            read (iu, '(/I6,F16.9,2X,8A8)') lbl2(k3), et2(k3), (cfg2(k,k3),k=1,&
               8) 
            read (iu, '(7F11.8)') (wt2(k),k=k1,k2) 
            k1 = k1 + ncf(2) 
            k2 = k2 + ncf(2) 
         end do 
         go to 11 
   13    continue 
         close(iu) 
      endif 
      if (ibugm /= 0) write (6, *) ' jv2 = ', (jv2(m),m=1,nvc(2)) 
      return  
      end subroutine eigvec 
