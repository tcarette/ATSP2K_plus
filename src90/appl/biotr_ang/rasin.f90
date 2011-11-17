!
!     ------------------------------------------------------------------
!     R A S I N
!     ------------------------------------------------------------------
!
! --- Sets up the RAS information of shells.
!     It is basically the original rasin module, but with all the
!     the active shells in ras2, assuming that both CSF expansions
!     satisfy c.u.d. where de-excitation means n->n' with n' .leq. n,
!     for a given l-value.
!
      subroutine rasin 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: double 
      USE ras_C 
      USE fo_C 
      use inout_C
      use non30_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  00:08:29  11/17/01  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, l, is, io, j, n, ic, k 
      character :: bl 
      character , dimension(2) :: label*7 
      character :: set*22 
!-----------------------------------------------
      data label/ 'Initial', 'Final  '/  
!      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2),iut(2)
!      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
!     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
!     :       (QLJCLSD,LJCLSD(1))
!      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      data set/ 'spdfghiklmnSPDFGHIKLMN'/  
   10 format(' inactive : ',11i5) 
   11 format(' RAS1     : ',11i5) 
   14 format(' Number of shells per l for ',a7,' state') 
   15 format('          : ',11i5) 
   16 format(' Labels of shells for ',a7,' state') 
      bl = ' ' 
      write (iscw, *) 
      write (6, *) ' RAS-type calculation for initial state: ' 
      write (6, *) ' --------------------------------------: ' 
      write (6, *) &
         '               s    p    d    f    g    h    i    k    l    m    n ' 
      write (iwrite, 10) (ninac(i,1),i=1,lmax(1)) 
      write (6, '(A,11I5)') ' GAS      : ', (nras2(l,1),l=1,lmax(1)) 
      write (6, *) 
      write (iscw, *) 
      write (6, *) ' RAS-type calculation for final   state: ' 
      write (6, *) ' --------------------------------------: ' 
      write (6, *) &
         '               s    p    d    f    g    h    i    k    l    m    n ' 
      write (iwrite, 10) (ninac(i,2),i=1,lmax(2)) 
      write (6, '(A,11I5)') ' GAS      : ', (nras2(l,2),l=1,lmax(2)) 
      write (6, *) 
!
      do is = 1, 2 
         write (iwrite, 14) label(is) 
         nl(:lmax(is),is) = ninac(:lmax(is),is) + nras2(:lmax(is),is) 
         write (6, *) &
      '               s    p    d    f    g    h    i    k    l    m    n ' 
         write (iwrite, 15) (nl(i,is),i=1,lmax(is)) 
      end do 
      write (6, *) 
!
      io = 0 
      do i = 1, 2 
         write (iwrite, 16) label(i) 
         do j = 1, lmax(i) 
            n = j 
! Inactive
            ic = ninac(j,i) 
    3       continue 
            if (ic /= 0) then 
               ic = ic - 1 
               io = io + 1 
               elras(io) = bl//char(ichar('0') + n)//set(j:j) 
               write (6, *) ' inactive: ', n, set(j:j) 
               n = n + 1 
               go to 3 
            endif 
!. GAS
            ic = nras2(j,i) 
    5       continue 
            if (ic == 0) cycle  
            ic = ic - 1 
            io = io + 1 
            elras(io) = bl//char(ichar('0') + n)//set(j:j) 
            write (6, *) '   GAS     : ', n, set(j:j) 
            n = n + 1 
            go to 5 
!
         end do 
      end do 
      write (6, *) (elras(k),k=1,maxnfo) 
      return  
      end subroutine rasin 
