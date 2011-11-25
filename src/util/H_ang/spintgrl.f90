!----------------------------------------------------------------------
!       S P I N T G R L 
!
!                from ATSP2K [CPC 176 (2007) 559]
!                modfied November 2011
!----------------------------------------------------------------------
!       This routine reads the yint.lst data storing the concatenated
!       arrays for the different blocks

      SUBROUTINE SPINTGRL

      use safeio

      IMPLICIT NONE
!
      integer n_read, n_left
      INTEGER  maxorb,ncl, mwf, nbb, lsd,iblock,itype
      integer ib1,ib2,n,nij,lasti,ibe,ncol,icount,int,i
      character*2 ih_file
! IO variables
      INTEGER fp1,fp2,fp3,ierr

      fp1=safe_open('yint.lst',ST='O',FO='U',AC='R')
      rewind(fp1)
      read(fp1) ncl, mwf, nbb, lsd
      maxorb = nwf-nclosd
      if (ncl.ne.nclosd.or.mwf.ne.maxorb.or.nbb.ne.nblock) then
        write(0,*) 'Yint.lst data not consistent with cfg.h data'
        write(0,*) 'Yint.lst :', ncl, mwf, nbb
        write(0,*) 'cfg.h    :', nclosd, nwf, nblock
        stop
      end if
      read(fp1) 
      read(fp1)

!     Allocate stuffs

      allocate(cn(ncn_tot))
      allocate(inptr(ncn_tot))

      allocate(jptr(ncfg_tot))

      allocate(jan(nze_tot))
      allocate(ico(nze_tot))    ! it's too much!(?)! TODO: clean

      allocate(lused(nint))
      allocate(ipackn(nint))

      ico = 0
      jan = 0
      jptr = 0
      ipackn = 0
      ib1 = 0
      ib2 = 0
!    read for each block
      fp3=safe_open('ico.lst',ST='O',FO='U',AC='R')
      do iblock = 1, nblock
        write(ih_file,'(I2.2)') iblock 
        fp2=safe_open('ih.'//ih_file//'.lst',ST='O',FO='U',AC='R')
        nij = 0;
        n = lsd;
        do while(nij < nze_bl(iblock));
          read(fp2,iostat=ierr) n, (jan(i+nij),i=ib1+1,ib1+n);
          read(fp3,iostat=ierr) n, (ico(i+nij),i=ib1+1,ib1+n);
          nij = nij + n;
        end do
        read(fp1,iostat=ierr) ncol,(jptr(i),i=ib2+1,(ib2+ncol))
        ib1 = ib1 + nij
        ib2 = ib2 + ncol
        if(safe_close(fp2)/=0) call crash
      end do
      if(safe_close(fp3)/=0) call crash
 
        nij = 0
!
! ***** READ  THE LIST OF INTEGRALS
!
      lasti = 0
!          ...F, G, R, or L integrals....
      DO INT = 1,4
        icount = 1
        read(fp1) itype, nint

!       .. integral pack number and logical variable indicating usage
        read(fp1) (ipackn(i),i=lasti+1,nint), &
                  (lused(i), i=lasti+1,nint)
        lasti = nint
        icount = icount + 1
        intptr(int) = lasti
      end do

      if(safe_close(fp1)/=0) call crash
      fp1=safe_open('c.lst',ST='O',FO='U',AC='R')

      ibe = 0;
      do iblock = 1, nblock
        n_left = ncn_bl(iblock)
        if (ncn_bl(iblock) > lsdim) then
          n_read = lsdim; 
        else
          n_read = ncn_bl(iblock);
        end if


        do while (n_left > 0)
          if (n_left<lsdim) n_read = n_left
          read(fp1) n,cn(ibe+1:ibe+n_read),inptr(ibe+1:ibe+n_read)
          ibe = ibe + n_read;
          n_left = n_left - n_read
        end do
      end do

      if(safe_close(fp1)/=0) call crash

      return
      end subroutine SPINTGRL

