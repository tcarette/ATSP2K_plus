*----------------------------------------------------------------------
*       S P I N T G R L
*----------------------------------------------------------------------
*       This routine reads the yint.lst data storing the concatenated
*       arrays for the different blocks

	SUBROUTINE SPINTGRL(cf_tot)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
        PARAMETER (NWD=60,NOD=220,NOFFD=800)
        PARAMETER(MEIG=20,MTERM=20)
*
        CHARACTER EL*3,ATOM*6,TERM*6
        INTEGER       IN,OUT,PRI,OUC,OUF,OUH,OUD
        COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,IUD,OUC,OUF,OUD,OUH,ISCW
        COMMON/LABEL/ EL(60),ATOM,TERM
        LOGICAL       FAIL,OMIT,EZERO,REL,ALL,TRACE
        COMMON/TEST/  FAIL,OMIT,EZERO,REL,ALL,TRACE
       COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
        COMMON/PARAM/  H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,
     :  NCFG,IB,IC,ID,
     :                D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,
     :                NSCF,NCLOSD,RMASS
        COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      POINTER (QWT,WT(1)),(QWP,WP(1)), (IQWPTR,W(1))
      COMMON /CFGS/ETOTAL,QWT,QWP,IQWPTR
*
      POINTER (pkval,kval(1)),(pvalue,value(1)),
     :        (pih,ih(1)),(pjh,jh(1)),
     :        (pcoeff,coeff(1)),(pnijptr,nijptr(1)),
     :        (pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PVALUE,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR

*     .. Logical array passed from nonh to give the used integrals
      POINTER (plused,lused(1))
      COMMON/USED/plused
      LOGICAL lused

      POINTER (pico,ico(1))
      COMMON/ICOFF/pico

      POINTER (qjptr,jptr(1))
      COMMON/COLS/qjptr

      integer cf_tot(nblock), n_read, n_left;
      logical leigen(meig,mterm), lguess
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm),
     :        nze_max(mterm)
      character term_bl(mterm)*3
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess, nze_max
*
      logical   :: hmx_memory, ico_memory, ih_memory, clst_memory
      common/memory_use/hmx_memory, ico_memory, ih_memory, clst_memory

      CHARACTER*72 string
      character*2 ih_file

       rewind(iud)
       read(iud) ncl, mwf, nbb, lsd
       maxorb = nwf-nclosd
       if (ncl.ne.nclosd.or.mwf.ne.maxorb.or.nbb.ne.nblock) then
         write(0,*) 'Yint.lst data not consistent with cfg.h data'
         write(0,*) 'Yint.lst :', ncl, mwf, nbb
         write(0,*) 'cfg.h    :', nclosd, nwf, nblock
         stop
       end if
       read(iud) string

101   read(iud) string
      maxorb = maxorb - 24
      if (maxorb .gt. 0) go to 101

*  Allocate memory for c.lst and pointer data
       call alcsts(cf_tot);

*      ico = 0
*      ih = 0;
      jptr = 0;
      kval = 0;
      ib1 = 0
      ib2 = 0
*    read for each block
      do iblock = 1, nblock
        if (hmx_memory) then
          write(ih_file,'(I2.2)') iblock 
          open(unit=11,file='ih.'//ih_file//'.lst',status='old',
     :         form='unformatted');
          nij = 0;
          n = lsd;
          do while(nij < nze_bl(iblock));
            read(11,iostat=ierr) n, (ih(i+nij),i=ib1+1,ib1+n);
            if (ico_memory) then
              read(12,iostat=ierr) n, (ico(i+nij),i=ib1+1,ib1+n);
            end if;
            nij = nij + n;
          end do;
        end if;
      read(iud,iostat=ierr) ncol,(jptr(i),i=ib2+1,(ib2+ncol))
      ib1 = ib1 + nij
      ib2 = ib2 + ncol
      close (11);
       !print '(8i10)',(ih(i),i=1,ib1)
       !print '(8i10)',(ico(i),i=1,ib1) 
       !print '(8i10)',(jptr(i),i=1,ib2)
      end do
 
        nij = 0
C
C ***** READ  THE LIST OF INTEGRALS
C
        lasti = 0
        
*          ...F, G, R, or L integrals....
      DO INT = 1,4
        icount = 1
        read(iud) itype, nint

*       .. ipackn of nonh is called kval here 
*       .. integral pack number and logical variable indicating usage
        read(iud) (kval(i),i=lasti+1,nint),
     :           (lused(i),i=lasti+1,nint)
        lasti = nint
        icount = icount + 1
        intptr(int) = lasti
      end do

*      .. adjust for the fact that overlaps have been removed

      intptr(6) = intptr(4)
      intptr(5) = intptr(3)
      intptr(4) = intptr(2)
      intptr(3) = intptr(2)
	
      if (clst_memory) then
        ibe = 0;
        do iblock = 1, nblock
          n_left = cf_tot(iblock)
          if (cf_tot(iblock) > ncodim) then
            n_read = ncodim; 
          else
            n_read = cf_tot(iblock);
          end if

          do while (n_left > 0)
            if (n_left<ncodim) n_read = n_left;
            read(39) n, (coeff(j),j = ibe+1,ibe+n_read),
     :                    (inptr(j), j = ibe+1,ibe+n_read)
            ibe = ibe + n_read;
            n_left = n_left - n_read
          end do
!       print *, 'c.lst'
!       print '(8i10)',(inptr(i),i=1,n_read)
!       print '(4F12.8)', (coeff(i),i=1,n_read) 
        end do
      end if
      return
      end

