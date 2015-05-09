      SUBROUTINE DIAG(ECONV,ACFG,CFGTOL,LAST,icycle,eigst_weight,cf_tot)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER (NWD=94,NOD=220,NOFFD=800)
      parameter(MEIG=20,MTERM=20)
      CHARACTER EL*3,ATOM*6,TERM*6,config*64, couple*80
      INTEGER       IN,OUT,PRI,OUC,OUF,OUH,OUD
      COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,iud,OUC,OUF,OUD,OUH,ISCW
      COMMON/LABEL/ EL(NWD),ATOM,TERM
      LOGICAL       FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
      COMMON/TEST/  FAIL,OMIT,EZERO,REL,ALL,TRACE
      COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
      COMMON/PARAM/  H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,
     :                NCFG,IB,IC,ID,
     :                D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,
     :                NSCF,NCLOSD,RMASS
      COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      POINTER (QWT,WT(1)),(QWP,WP(1)),  (IQWPTR,W(1))
      COMMON /CFGS/ETOTAL,QWT,QWP,IQWPTR
      POINTER (IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :   (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :   (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :             QMETH,QIEPTR
*
      POINTER (pkval,kval(1)),(pvalue,value(1)),
     :        (pih,ih(1)),(pjh,jh(1)),
     :        (pcoeff,coeff(1)),(pnijptr,nijptr(1)),
     :        (pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PVALUE,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR

      POINTER (plused,lused(1))
      COMMON/USED/plused
      LOGICAL lused

      logical   :: hmx_memory, ico_memory, ih_memory, clst_memory
      common/memory_use/hmx_memory, ico_memory, ih_memory, clst_memory
      common/dsk/iblock,ih_total,diag_hmx_memory,diag_ih_memory
      integer ih_total
      logical diag_ih_memory,diag_ico_memory,diag_hmx_memory;     
      double precision, allocatable, dimension(:,:) :: isom_shift

      POINTER (pico,ico(1))
      COMMON/ICOFF/pico

      POINTER (qjptr,jptr(1))
      COMMON/COLS/qjptr
    
      integer idispl(2), icount(2)

      POINTER (qhmx,hmx(1)),(qtm,tm(1)),(qtp,tp(1)),
     :        (qdiag,hii(1)),(qiwork,iwork(1))
      common/spd/qhmx,qtm,qtp,qdiag,qiwork 

      integer cf_tot(nblock), nze_count
 
      double precision eigst_weight(meig,mterm)
      logical   ldused
      character*64  string
      logical leigen(meig,mterm), lguess, leig_out
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm),
     :        nze_max(mterm)
      character term_bl(mterm)*3
      character*2 ih_file
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess,nze_max
      INTEGER iselec(20)
      LOGICAL ECONV,LAST,iguess,iupper,hiend,lfirst
      DATA iguess,iupper,lfirst/.false.,.false.,.true./
      DATA crite,critc,critr,trhold,ortho/
     :       1.d-14,1.d-8,1.d-8,1.d0,1.d0/


      idisk = 0
      leig_out = .false.
      nze = maxval(nze_bl(1:nblock))
      lim = maxval(iiws(1:nblock))
      iiwsz = 6*lim + maxval(nume(1:nblock))
      iworksz = maxval(iws(1:nblock))
      rewind 39
      etl = 0
      sum_energy = 0
      nz = 0
      nj = 0
      ijp = 1
      ievstart = 1
      ih_start = 1;
      ies = 1
      ibe = 0;
      ncoef = 0;
      rewind (11);
      rewind (12);

*<><><><><><><> start blocks <><><><><><><><><> 
      do iblock = 1,nblock
*       .. set parameters to current block
        ncfg = ncfg_bl(iblock)
        nze = nze_bl(iblock)
        ih_total = nze_bl(iblock)
        lim = iiws(iblock)
        iworksz = iws(iblock)
        iiwsz = 6*lim + nume(iblock)
        nijcurr = 1
        max_col = 1; 
        jjh = 1
        diag_ih_memory  = .true.
        diag_ico_memory = .true.
        diag_hmx_memory = .true.


*  if hmx_mem = .false. allocate memory for current block
        if (.not.hmx_memory) then
          call diag_allocate(pih,nze_bl(iblock),nze_max(iblock),
     :                           diag_ih_memory,4);
!          diag_ih_memory = .false.   ! DBG code
          if (.not.diag_ih_memory) ih_total = nze_max(iblock)
          call diag_allocate(qhmx,nze,nze_max(iblock),
     :                          diag_hmx_memory,8);
!          diag_hmx_memory = .false.  ! DBG code
!         if (.not.diag_ih_memory) diag_hmx_memory = .false.  ! DBG code

          if (.not.diag_hmx_memory) diag_ico_memory = .false.
          nze_count = 0;

          if (diag_ih_memory) then
          write(ih_file,'(I2.2)') iblock
            open(unit=11,file='ih.'//ih_file//'.lst',status='old',
     :          form='unformatted');
            do while (nze_count.lt.nze_bl(iblock))
              read(11) nze_ih, (ih(jj),jj=1+nze_count,nze_count+nze_ih)
              nze_count = nze_count + nze_ih
            end do
            ih_start = 1
            close (11)
          end if

          if (diag_hmx_memory.and.diag_ih_memory) then
            call diag_disk_ico(ncoef,ico,hmx,coeff,value,inptr,
     :         jptr,jjh,nijcurr,ncodim,ncfg,nze_bl(iblock),
     :         nze_max(iblock))
          else 
            open(unit=13,status='scratch',form='unformatted');
            call diag_disk_hmx(ncoef,ico,hmx,coeff,value,inptr,
     :        jptr(nj+1),jjh,nijcurr,ncodim,ncfg,nze_max(iblock),
     :        nze_bl(iblock),hii,shift)
          end if
        end if

* if hmx_mem = .true.
        if (clst_memory) then
          call diag_memory_all(cf_tot(iblock),ncoef,ico,nz,hmx,coeff,
     :         value,inptr,jptr,jjh,nijcurr,nze_bl(iblock))
        end if

        if (ico_memory.and.(.not.clst_memory)) then 
          call diag_disk_clst(ncoef,ico,hmx,coeff,value,inptr,
     :         nijcurr,ncodim,nze_bl(iblock),nz)
        end if
  
        if (hmx_memory.and.(.not.ico_memory)) then
          call diag_disk_ico(ncoef,ico,hmx,coeff,value,inptr,
     :       jptr(nj+1),jjh,nijcurr,ncodim,ncfg,nze_bl(iblock),
     :       nze_max(iblock))
        end if          

      if (icycle == 0) call print_memory(icycle,iblock,clst_memory,
     :  ico_memory,hmx_memory,diag_hmx_memory,nblock,term_bl(iblock),
     :  ncfg);

*     .. compute diagonals
      if (ico_memory.or.diag_hmx_memory) then
        shift = hmx(1)
        hmx(1) = 0.d0
        hii(1) = 0.d0
        do i = 2,ncfg
          hii(i) = hmx(jptr(i-1+nj)+1) - shift
          hmx(jptr(i-1+nj)+1) = hii(i)
        end do
*        ih_start = nze + 1
      else 
        nze = nze_max(iblock);
        ih_start = 1
      end if
*
****** COMPUTE THE EIGENVALUES AND EIGENVECTORS
*
      hiend = .false. 
      ilow = -1
      ihigh = -1
      ie = 1
      do i = 1,meig
         if (leigen(i,iblock)) then
           iselec(ie) = i
           ie = ie + 1
         end if
      end do
      iselec(ie) = -1

      if (.not. lguess) then
        niv = 0
      else
        niv = niv_bl(iblock)
        call dcopy(ncfg*niv,eigvec(ievstart),1,wt(1),1)
      end if
      mblock = ie-1

*      print *, 'hmx col 1/10000:'
*      print '(3F20.16)', (hmx(jptr(i)),i=1,ncfg,10000)
*      print *, 'ih col 1/10000:'
*      print '(4I6)', (ih(jptr(i)+nz),i=1,ncfg,10000)
*      print *, 'jptr:', (jptr(i+nj),i=1,ncfg)
*      print *, 'iupper,nze,ncfg,lim:',iupper,nze,ncfg,lim
*      write(iscw,*), 'shift =',shift
*      write(iscw,*), 'diag :'
*      write(iscw,'3F20.16'),(hii(i),i=1,ncfg)
*     print *, 'ilow,ihigh,iselec(1),iselec(2),iselec(3)',
*    :          ilow,ihigh,iselec(1),iselec(2),iselec(3)

      if (clst_memory.or.ico_memory) nze = nze_bl(iblock) 

      if(.not.diag_ih_memory) then
        write(ih_file,'(I2.2)') iblock
          open(unit=11,file='ih.'//ih_file//'.lst',status='old',
     :       form='unformatted');
      end if

      maxiter = max0(30*nume(iblock),100)

ctc      print *, 'lim =', lim !ctc

      CALL dvdson(hmx,ih(ih_start),jptr(nj+1),iupper,nze,tm,tp,
     :        Ncfg,lim,hii,ilow,ihigh,iselec,niv,mblock,
     :        crite,critc,critr,ortho,maxiter,
     :        wt,iworksz,iwork,iiwsz,hiend,nloops,nmv,ierr)

      if (ierr .ne. 0) then
            write(0,*) ' Dvdson returned ierr=',ierr
            write(0,*) ' ... continuing'
      end if
      if (.not.diag_hmx_memory) close (11);
      close (13);
 
      if (hiend) then
         print *, 'exit: dvdson returned hiend = true; enter 0 to stop'
         read *, i
         if (i.eq.0) stop
      end if
ctc      print *, 'wt :',(wt(i), i=1,ncfg) !ctc
     
      if (lguess .eqv. .false.) lguess = .true.
      do icc = 1, nume(iblock)
        if (leigen(icc,iblock)) then
           ioffw = ncfg*(icc-1)
           ioffe = icc
           j = idamax(ncfg,wt(ioffw+1),1)
           if (wt(j+ioffw) .lt. 0) then
             wt(ioffw+1:ioffw+ncfg) = -wt(ioffw+1:ioffw+ncfg)
           end if
           etl =  wt(ncfg*nume(iblock)+ioffe) + ec + shift
           wt(ncfg*nume(iblock)+ioffe) = etl
           sum_energy = sum_energy +etl*eigst_weight(icc,iblock)
         end if
      end do 
ctc
*        print *, 'l265 diag'
ctc
*        .. save eigenvectors and eigenvalues
         call dcopy(ncfg*nume(iblock),wt(1),1,eigvec(ievstart),1)
         call dcopy(nume(iblock), wt(ncfg*nume(iblock)+1),1,en(ies),1)
*        .. print Energies and some component of eigenvectors
*         ioffw = ncfg*nume(iblock)
         do ie = 1,nume(iblock)
         ioffw = ncfg*nume(iblock)
           if( leigen(ie,iblock)) then
            write (out,'(/A,F15.8,3X,A,I3,2(1PD11.3))')
     :      '         ETOTAL=',wt(ioffw+ie),
     :      'Loops,DeltaE,Res.: ', nloops,
     :      wt(ioffw+nume(iblock)+ie),
     :      wt(ioffw+2*nume(iblock)+ie)
            ioffw = ncfg*(ie-1)
            write(out,'(4(4X,I3,F11.7))') (i,wt(ioffw+i),
     : i=1,min(8,ncfg))
           end if
         end do
ctc
*      print *, 'stop l285 diag'
*      STOP
ctc
*        .. we now have nume starting vectors for the next iteration
         niv_bl(iblock) = nume(iblock)
*        .. update the terminating/starting locations
         ievstart = ievstart + ncfg*nume(iblock)
         ies = ies +nume(iblock)
         nj = nj + ncfg   !nijcurr 
         nz = nz + nze_bl(iblock) 
         if (clst_memory.or.ico_memory) ih_start = nz + 1
         if (hmx_memory.and.(.not.ico_memory)) ih_start = nz + 1
         if (.not.hmx_memory) then 
           call diag_deallocate(qhmx);
           call diag_deallocate(pih);
         end if
      end do
        
*<><><><><><><><><> end blocks <><><><><><><><><><><
      if (last) then
        if (hmx_memory) call dalloc(qhmx,nze)
        call dalloc(qtm,ncfg)
        call dalloc(qtp,ncfg)
        call dalloc(qdiag,ncfg)
        call dalloc(qiwork,iiwsz)
      end if
  
      allocate(isom_shift(MEIG,nblock)); 
      isom_shift(1:MEIG,1:nblock) = 0;
      lguess = .true.
      deltae = sum_energy - etotal
      if (abs(deltae) .lt. cfgtol) econv = .TRUE.
      etotal = sum_energy
      write(out, '(/9x,A,F15.8)') 'Sum of ETOTAL :  ',etotal
      write(iscw,*)' DeltaE =', deltae, 'Sum_Energy = ' , sum_energy
      ievstart = 1 
      !rewind(11)
      rewind(12)
      rewind(39)
      nze_count = 0
      ijp = 1
      nz  = 1
      ibe = 0;
      ncoef = 0;
*<<<<<<<<<<<<<<<<<<<updatc>>>>>>>>>>>>>>>>>>>>>>>>>>>.
      do iblock = 1,nblock
         ncfg = ncfg_bl(iblock)

        if (.not.hmx_memory) then
          call diag_allocate(pih,nze_bl(iblock),nze_max(iblock),
     :                           diag_ih_memory,4);
!         diag_ih_memory = .false.   ! DBG code

          nze_count = 0
          if (diag_ih_memory) then
            write(ih_file,'(I2.2)') iblock
            open(unit=11,file='ih.'//ih_file//'.lst',status='old',
     :          form='unformatted');
            do while (nze_count.lt.nze_bl(iblock))
              read(11) nze_ih, (ih(jj),jj=1+nze_count,nze_count+nze_ih)
              nze_count = nze_count + nze_ih
            end do
            ih_start = 1
            CALL UPDATC_disk_ico(iblock,eigvec(ievstart),ncfg,
     :          nume(iblock),jptr(ijp),ih(1),eigst_weight,
     :          max_col,nijcurr,last,isom_shift)
             close(11)
          else 
            CALL UPDATC_disk_ih(iblock,eigvec(ievstart),ncfg,
     :          nume(iblock),jptr(ijp),ih(1),eigst_weight,
     :          max_col,nijcurr,last,isom_shift)
          end if
            call diag_deallocate(pih);
        end if

         if (clst_memory) then
           CALL UPDATC_memory_all(iblock,eigvec(ievstart),ncfg,
     :          nume(iblock),jptr(ijp),ih(nz), ico(nz), eigst_weight,
     :          max_col,nijcurr,cf_tot,last,isom_shift)
         end if

         if (ico_memory.and.(.not.clst_memory)) then
           CALL UPDATC_disk_clst(iblock,eigvec(ievstart),ncfg,
     :          nume(iblock),jptr(ijp),ih(nz),ico(nz),eigst_weight,
     :          max_col,nijcurr,last,isom_shift)
         end if

         if (hmx_memory.and.(.not.ico_memory)) then
           CALL UPDATC_disk_ico(iblock,eigvec(ievstart),ncfg,
     :          nume(iblock),jptr(ijp),ih(nz),eigst_weight,
     :          max_col,nijcurr,last,isom_shift)
         end if
 
         ievstart = ievstart + ncfg*nume(iblock)
         ijp = ijp + ncfg_bl(iblock) 
         nz = nz + nze_bl(iblock)
      end do
*<<<<<<<<<<<<<<<<<<end updatc>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      DO I=1,nwf
        IF(SUM(I).EQ.0.0)write(iscw,*)'DIAG: SUM=',SUM(I),'  for ',EL(I)
      end do

      if ((mod(icycle,10)).eq.0) then
        leig_out = .true.
      end if

      if (last.or.leig_out) then
          rewind(60)
          rewind(30)
          ievstart =1
          ies = 1
          do iblock = 1,nblock
            call eig_out(iblock,eigvec(ievstart),en(ies),z,isom_shift)
            ievstart = ievstart + ncfg_bl(iblock)*nume(iblock)
            ies = ies + nume(iblock)
          end do
      endif
      deallocate(isom_shift);

      end subroutine diag;
