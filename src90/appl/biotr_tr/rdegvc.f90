!     ------------------------------------------------------------
!     R D E G V C
!     ------------------------------------------------------------
!
! --- determine the number of different J-values (njv)
!               the total number of eigenvectors (nvc)
!               the total length of the wt vector to allocate (lgth)
!
      subroutine rdegvc 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double1=>double 
      use state_C
      use inout_C
      use ems_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:09:36  11/18/01  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ii, iu, nl, jv, nb, m, k 
!      logical :: rel, vok 
!-----------------------------------------------
!      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
!      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2),iut(2)
!      COMMON /STATE/QET1,QLBL1,QWT1,QJV1,QET2,QLBL2,QWT2,QJV2,NJV(2),
!     : NVC(2),LGTH(2),NCF(2),QCFG1,QCFG2
!      POINTER(QET1,ET1(1)),(QLBL1,LBL1(1)),(QWT1,WT1(1)),(QJV1,JV1(1)),
!     :       (QET2,ET2(1)),(QLBL2,LBL2(1)),(QWT2,WT2(1)),(QJV2,JV2(1)),
!     :       (QCFG1,CFG1(8,1)),(QCFG2,CFG2(8,1))
!
   10 format(/,/,8x,i4,10x,i4) 
!
!dbg  print*,' ncf(1) = ',ncf(1),' ncf(2) = ',ncf(2)
      do ii = 1, 2 
         if (rel) then 
            iu = iuj(ii) 
         else 
            iu = iul(ii) 
         endif 
         nvc(ii) = 0 
         njv(ii) = 0 
         nl = (ncf(ii)-1)/7 + 1 
!
         read (iu, '(A)') 
    2    continue 
         read (iu, 10, end=5) jv, nb 
         njv(ii) = njv(ii) + 1 
         nvc(ii) = nvc(ii) + nb 
         do m = 1, nb 
            read (iu, '(A)') 
            read (iu, '(A)') 
            do k = 1, nl 
               read (iu, '(A)') 
            end do 
         end do 
         go to 2 
    5    continue 
         lgth(ii) = nvc(ii)*ncf(ii) 
!dbg  print*,' ii = ',ii,' nvc = ',nvc(ii),'ncf = ',ncf(ii)
!dbg  print*,'             lgth = ',lgth(ii),' njv = ',njv(ii)
         rewind (unit=iu) 
      end do 
      return  
      end subroutine rdegvc 
