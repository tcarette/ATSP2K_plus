      SUBROUTINE DIAG(ECONV,ACFG,CFGTOL,LAST,icycle,eigst_weight,
     :                tmp,cf_tot)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER (NWD=94,NOD=220,NOFFD=800)
      parameter(MEIG=20,MTERM=20)
      INCLUDE 'mpif.h'
      parameter (MAXPROC=100)
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      common /PVM/ istart,ifinish

*
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
*
      POINTER (IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :  (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :  (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :             QMETH,QIEPTR
*
      POINTER (pkval,kval(1)),(pvalue,value(1)),(pih,ih(1)),(pjh,jh(1)),
     :  (pcoeff,coeff(1)),(pnijptr,nijptr(1)),(pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PVALUE,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR

      POINTER (plused,lused(1))
      COMMON/USED/plused
      LOGICAL lused

      logical                   ::clst_disk, clst_memory
      common/use_disk/clst_disk, clst_memory
      double precision, allocatable, dimension(:,:) :: isom_shift

      POINTER (pico,ico(1))
      COMMON/ICOFF/pico

      POINTER (qjptr,jptr(1))
      COMMON/COLS/qjptr
    
      integer idispl(2), icount(2)

      POINTER (qhmx,hmx(1)),(qtm,tm(1)),(qtp,tp(1)),
     :        (qdiag,hii(1)),(qiwork,iwork(1))
      common/spd/qhmx,qtm,qtp,qdiag,qiwork 

      dimension tmp(*)
      integer cf_tot(nblock)
 
      double precision eigst_weight(meig,mterm)
      logical   ldused
      character*64  string
      logical leigen(meig,mterm), lguess, leig_out
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm)
      character term_bl(mterm)*3
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess
*
      INTEGER iselec(20)
      LOGICAL ECONV,LAST,iguess,iupper,hiend,lfirst
      DATA iguess,iupper,lfirst/.false.,.false.,.true./
      DATA crite,critc,critr,trhold,ortho/
     :       1.d-14,1.d-8,1.d-8,1.d0,1.d0/

*     .. this is a matrix in memory version
      idisk = 0
      leig_out = .false.
******
      nze = maxval(nze_bl(1:nblock))
!!!      lim = maxval(iiws(1:nblock))
      lim = 1000;
      iiwsz = 6*lim + maxval(nume(1:nblock))
      iworksz = maxval(iws(1:nblock))
Ctc 19/03/2008 : crash of large calculations on several nodes of hydra@ulb
Ctc   print *, 'IWSZ, iworksz = ', iiwsz, iworksz
      WRITE (80+myid,*) 'IWSZ, iworksz = ', iiwsz, iworksz
Ctc
*  .. for block calculation, we use block information ??????

*     In generating the matrix, we need to read in the data,


      rewind 39
      etl = 0
      sum_energy = 0

*<<<<<<<<<<<<<<<<< for all blocks >>>>>>>>>>>>>>>
*
*    .. we have global arrays from yint.lst
*      ih, ico -- these increment according to the number of nze
*      jptr    -- this increments by ncfg
*    .. arrays from c.lst are buffered with dimensions ncodim (also lsdim)
*
*     nz is a pointer to the end of the previous block of ih, ico
*     nj is a pointer to the end of the previous block of jptr
*     ievstart is the beginning of the current block of eigenvectors
*     ies is start of current set of eigenvalue
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      nz = 0
      nj = 0
      ijp = 1
      ievstart = 1
      ies = 1
      ibe = 0;
      ncoef = 0;
      do iblock = 1,nblock
      istart =  -1;
*       .. set parameters to current block
        ncfg = ncfg_bl(iblock)
        nze = nze_bl(iblock)
        lim = iiws(iblock)
	iworksz = iws(iblock)
        iiwsz = 6*lim + nume(iblock)
        call dfill(nze,0.0d0,hmx,1)
        nijcurr = 1
        max_col = 1; jjh = 1

        if (clst_memory) then
           call coeff_memory(cf_tot(iblock),ncoef,ico,nz,hmx,coeff,
     :                    value,inptr,jptr,jjh,nijcurr)

        else 
           call coeff_disk(ncoef,ico,nz,hmx,coeff,
     :                  value,inptr,jptr,jjh,nijcurr,ncodim)
        end if
  
*     .. compute diagonals

        if (myid.eq.0) then
           shift = hmx(1)
        end if

        call MPI_BCAST(shift,1,MPI_DOUBLE_PRECISION,0,
     :                  MPI_COMM_WORLD,ierr)

        call DINIT(NCFG,0.d0,hii,1)

*        ..all the nodes should shift their first hmx(1)
*        ..interleaved
 
         hii(myid+1)=hmx(1)-shift
         hmx(1)=hii(myid+1)
         m_col=1
         do irow=myid+1+nprocs,NCFG,nprocs
            m_col=m_col+1
            hii(irow)=hmx(jptr(m_col-1+nj)+1) - shift
            hmx(jptr(m_col-1+nj)+1)= hii(irow)
         end do

       call gdsummpi(hii,NCFG,tm)
       ! call mpi_allr_dp(hii,ncfg)
        
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
*      if (myid == 0) then
*      print *, 'hmx:'
*      print '(3F20.16)', (hmx(i),i=1,jptr(ncfg))
*      end if
*      print *, 'ih:'
*      print '(4I6)', (ih(i+nz),i=1,jptr(ncfg))
*      print *, 'jptr:', (jptr(i+nj),i=1,ncfg)
*      print *, 'iupper,nze,ncfg,lim:',iupper,nze,ncfg,lim
*      if (myid == 0 ) then
*      print *, 'diag :'
*      print '(3F20.16)',(hii(i),i=1,ncfg)
*      end if
*     print *, 'ilow,ihigh,iselec(1),iselec(2),iselec(3)',
*    :          ilow,ihigh,iselec(1),iselec(2),iselec(3)
ctc avoid ierr=2048 in CI
ctc      maxiter = max0(30*nume(iblock),100)
      NIT=NWF-IB+1
      IF (NIT==0) THEN
        write(0,*) 'CI calculation, maxiter increased'
        maxiter = max0(30*nume(iblock),400)
      ELSE
        maxiter = max0(30*nume(iblock),100)
      END IF
ctc
*     print *, 'niv,crite,critc,critr,trhold,ortho,maxiter',
*    :          niv,crite,critc,critr,trhodl,ortho,maxiter
*     print *,  'ncfg,lim', ncfg,lim
*     write(0,*) ' Calling eigen: worksize=',iworksz,n,nze

      CALL dvdson(hmx,ih(nz+1),jptr(nj+1),iupper,nze,tm,tp,
     :        Ncfg,lim,hii,ilow,ihigh,iselec,niv,mblock,
     :        crite,critc,critr,ortho,maxiter,
     :        wt,iworksz,iwork,iiwsz,hiend,nloops,nmv,ierr)

****** deallocate when not needed

      if (ierr .ne. 0) then
         if (myid == 0) then
            write(0,*) ' Dvdson returned ierr=',ierr
            write(0,*) ' ... continuing'
         end if
      end if
   
      if (hiend) then
Ctc 19/03/2008 : crash of large calculations on several nodes of hydra@ulb
Ctc      print *, 'exit: dvdson returned hiend = true; enter 0 to stop'
         WRITE (80+myid,*)
     :             'exit: dvdson returned hiend = true; enter 0 to stop'
Ctc
         read *, i
         if (i.eq.0) stop
      end if
      
*      print *, (wt(i), i=1,ncfg) 
     
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

*        .. save eigenvectors and eigenvalues
         call dcopy(ncfg*nume(iblock),wt(1),1,eigvec(ievstart),1)
         call dcopy(nume(iblock), wt(ncfg*nume(iblock)+1),1,en(ies),1)
*        .. print Energies and some component of eigenvectors
*	 ioffw = ncfg*nume(iblock)
	 do ie = 1,nume(iblock)
         ioffw = ncfg*nume(iblock)
	   if( leigen(ie,iblock)) then
            if (myid == 0) then
            write (out,'(/A,F15.8,3X,A,I3,2(1PD11.3))')
     :      '         ETOTAL=',wt(ioffw+ie),
     :      'Loops,DeltaE,Res.: ', nloops,
     :      wt(ioffw+nume(iblock)+ie),
     :      wt(ioffw+2*nume(iblock)+ie)
            end if
	    ioffw = ncfg*(ie-1)
            if (myid == 0 ) then
	    write(out,'(4(4X,I3,F11.7))') (i,wt(ioffw+i),i=1,min(8,ncfg))
            end if
	   end if
	 end do
*        .. we now have nume starting vectors for the next iteration
         niv_bl(iblock) = nume(iblock)

*        .. update the terminating/starting locations
         ievstart = ievstart + ncfg*nume(iblock)
         ies = ies +nume(iblock)
         nj = nj + m_col   !nijcurr 
	 nz = nz + nijcurr !jptr(nj) 
      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
*<<<<<<<<<<<< end blocks >>>>>>>>>>>>

      if (last) then
        if (qhmx.ne.0) call dalloc(qhmx,nze)
        if (qtm.ne.0) call dalloc(qtm,ncfg)
        if (qtp.ne.0) call dalloc(qtp,ncfg)
        if (qdiag.ne.0) call dalloc(qdiag,ncfg)
        if (qiwork.ne.0) call dalloc(qiwork,iiwsz)
      end if

      allocate(isom_shift(MEIG,MTERM));
      isom_shift(1:MEIG,1:MTERM) = 0.0D0;
      lguess = .true.
      deltae = sum_energy - etotal
      if (abs(deltae) .lt. cfgtol) econv = .TRUE.
      etotal = sum_energy
      if (myid.eq.0) then
        write(out, '(/9x,A,F15.8)') 'Sum of ETOTAL :  ',etotal
        write(iscw,*)' DeltaE =', deltae, 'Sum_Energy = ' , sum_energy
      end if
      ievstart = 1 
      rewind(39)
      ijp = 1
      nz  = 1
      ibe = 0;
      ncoef = 0;
      do iblock = 1,nblock
         ncfg = ncfg_bl(iblock)
         if(clst_memory) then
           CALL UPDATC_memory(iblock,eigvec(ievstart),ncfg,nume(iblock),
     :     jptr(ijp),ih(nz), ico(nz), eigst_weight,tmp,max_col,nijcurr,
     :     cf_tot,last,isom_shift)
         else
           CALL UPDATC_disk(iblock,eigvec(ievstart),ncfg,nume(iblock),
     :     jptr(ijp),ih(nz), ico(nz), eigst_weight,tmp,max_col,nijcurr,
     :     last,isom_shift)
         end if
         ievstart = ievstart + ncfg*nume(iblock)
         ijp = ijp + max_col
         nz = nz + nijcurr
      end do
*
      if (myid.eq.0) then
        DO I=1,nwf
        IF(SUM(I).EQ.0.0)write(iscw,*)'DIAG: SUM=',SUM(I),'  for ',EL(I)
        end do
      end if


      if ((mod(icycle,10)).eq.0) then
        leig_out = .true.
      end if
 
      if (myid == 0) then
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
      end if
      deallocate(isom_shift); 
!       call MPI_Barrier(MPI_COMM_WORLD,ierr)
      end
