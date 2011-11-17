*-----------------------------------------------------------------------
*               U P D A T C
*-----------------------------------------------------------------------
*
*     Evaluate all coeffients of integrals.  This involves a
*     sequential pass through all the c.lst data.
*
*
      SUBROUTINE UPDATC_memory(iblock,wt,ncfg,maxev,jptr,ihh,ico,
     :    eigst_weight,tmp,max_col,nijcurr,cf_tot,last,isom_shift)
      IMPLICIT  double precision(A-H,O-Z)
      PARAMETER (NWD=70,NWC=20,NOD=220,NOFFD=800)
      parameter(MEIG=20,MTERM=20)
*
*     MPI stuff ***********************************************
*
        INCLUDE 'mpif.h'
        parameter (MAXPROC=100)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
****************************************************************

      POINTER (pkval,kval(1)),(pcoef,coef(1)),(pih,ih(1)),(pjh,jh(1)),
     :  (pcoeff,coeff(1)),(pnijptr,nijptr(1)),(pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PCOEF,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
*
      INTEGER       IN,OUT,PRI,OUC,OUF,OUH,OUD
      COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,iud,OUC,OUF,OUD,OUH,ISCW
      COMMON/LABEL/ EL(NWD),ATOM,TERM
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :  (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :  (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :           QMETH,QIEPTR

      COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
      dimension wt(ncfg*maxev), jptr(*), ico(*), ihh(*)
      character string*64

      double precision eigst_weight(meig,mterm)
      double precision isom_shift(meig,mterm)
      integer cf_tot(nblock)
      integer, allocatable, dimension(:) :: jptr_tmp,ico_tmp
      logical leigen(meig,mterm),lguess,last
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm)
      character term_bl(mterm)*3
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess
      dimension tmp(*)
      double precision, allocatable, dimension(:) :: tmp_coef
      double precision, allocatable, dimension(:) :: tmp_tc
      save ncoef
   
*     .. put zeros in coef
      if ((iblock.eq.1)) then
        nint = intptr(6)
        coef(1:nint) = 0
      end if

*  allocate temp_coef for computing isotope shift
      if (last) then
        itmpu = 0;
        do im1 = 1,maxev
          if (leigen(im1,iblock)) itmpu = itmpu + 1;
        end do
        allocate(tmp_coef(itmpu*(intptr(6))));
        allocate(tmp_tc(itmpu*intptr(6)));
        tmp_coef(1:itmpu*intptr(6)) = 0.0;
        tmp_tc(1:itmpu*intptr(6)) = 0.0D0;
      end if

      nijcurr = 1; 
      iih = 1 + myid;
      nij_tmp = 1; 
      jjh = 1 + myid 
      max_col = 1;
      n_count_tmp = 0; 
      if (iblock == 1) ncoef = 0;
*        .. add contribution from this record to coef
*        .. make use of the fact that all coefficients are read in
*           increasing order of the nze elements of the matrix
        do i = 1, cf_tot(iblock);
          n_count_tmp = ncoef+i
          if (i.gt.ico(nijcurr)) then
            nijcurr = nijcurr + 1
            if (nijcurr.gt.jptr(max_col)) then
              jjh = jjh + nprocs;
              max_col = max_col+1
            end if
          end if
          iih = ihh(nijcurr)
          im1 = 0;
          do j = 1,maxev
            ioffw = (j-1)*ncfg
            if (leigen(j,iblock)) then
               wcoef = eigst_weight(j,iblock)
               W = wcoef*wt(ioffw+iih)*wt(ioffw+jjh)
               T = W*coeff(n_count_tmp)
               IF (IIH .NE. JJH) T = T+T
               coef(inptr(n_count_tmp)) = coef(inptr(n_count_tmp)) + T
               if (last) then
                 W0 = wt(ioffw+iih)*wt(ioffw+jjh);
                 T0 = W0*coeff(n_count_tmp);
                 IF (IIH .NE. JJH) T0 = T0 + T0;
                 itmp_s = (inptr(n_count_tmp))+intptr(6)*im1;
                 tmp_coef(itmp_s) = tmp_coef(itmp_s) + T0;
                 im1 = im1 + 1;
               end if
            end if
          end do
*	     if (inptr(i) .gt. intptr(5)) then
*         print '(2I4,3F15.10,I6,F15.10,i6)',
*     :  iih,jjh,wt(iih+ioffw),wt(jjh+ioffw),coeff(n_count_tmp),
*     :  inptr(n_count_tmp), coef(inptr(n_count_tmp)),nijcurr
*	end if
      end do
 
      ncoef = ncoef + cf_tot(iblock);

      if (last) then
         call MPI_ALLREDUCE(tmp_coef,tmp_tc,itmpu*intptr(6),
     :          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
         if (myid == 0) then
           im2 = 0;
           do im1 = 1,maxev
             if (leigen(im1,iblock)) then
               istart = 1+intptr(6)*im2;
               call prprty(im1,iblock,tmp_tc(istart),intptr(6),
     :                         isom_shift(im1,iblock));
               im2 = im2 + 1;
             end if
           end do
         end if
         deallocate(tmp_coef);
         deallocate(tmp_tc);
      end if

      if (iblock == nblock) then   
       ! call gdsummpi(coef, intptr(6), tmp)
        call mpi_allr_dp(coef, intptr(6))

*
*  *****  DEFINE SUM(I)
*
         IBEGIN = INTPTR(5)+1
         IEND = INTPTR(6)

         DO I = IBEGIN,IEND
           call unpacki(6,i,kv,iel1,iel2,iel3,iel4)
           IF (IEL1 .EQ. IEL2) SUM(IEL1) = -2*COEF(I)
         end do 

*         if (myid == 0) then
*       print *, 'Node:',' intptr(6)',intptr(6), 'iblock',iblock
*       print '(4f20.16)', (coef(i),i=1,intptr(6))
*       print *, 'SUM', Sum(1:10)
*       print '(4F20.16)',SUm(1:20)
*         end if
      end if
      END

*-----------------------------------------------------------------------
*               U P D A T C
*-----------------------------------------------------------------------
*
*     Evaluate all coeffients of integrals.  This involves a
*     sequential pass through all the c.lst data.
*
*
      SUBROUTINE UPDATC_disk(iblock,wt,ncfg,maxev,jptr,ihh,ico,
     :         eigst_weight,tmp,max_col,nijcurr,last,isom_shift)

      IMPLICIT  double precision(A-H,O-Z)
      PARAMETER (NWD=70,NWC=20,NOD=220,NOFFD=800)
      parameter(MEIG=20,MTERM=20)
*
*     MPI stuff ***********************************************
*
        INCLUDE 'mpif.h'
        parameter (MAXPROC=100)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
****************************************************************

*
      POINTER (pkval,kval(1)),(pcoef,coef(1)),(pih,ih(1)),(pjh,jh(1)),
     :    (pcoeff,coeff(1)),(pnijptr,nijptr(1)),(pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PCOEF,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
*
      INTEGER       IN,OUT,PRI,OUC,OUF,OUH,OUD
      COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,iud,OUC,OUF,OUD,OUH,ISCW
      COMMON/LABEL/ EL(NWD),ATOM,TERM
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     : (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     : (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :           QMETH,QIEPTR

      COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
      dimension wt(ncfg*maxev), jptr(*), ico(*), ihh(*)
      character string*64

      double precision eigst_weight(meig,mterm)
      double precision isom_shift(meig,mterm)
      integer, allocatable, dimension(:) :: jptr_tmp,ico_tmp
      logical leigen(meig,mterm),lguess,last
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm)
      character term_bl(mterm)*3
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess
      double precision, allocatable, dimension(:) :: tmp_coef
      double precision, allocatable, dimension(:) :: tmp_tc
      dimension tmp(*)
*
   
*     .. put zeros in coef
      if ((iblock.eq.1)) then
        nint = intptr(6)
        coef(1:nint) = 0
      end if

      nijcurr = 1; 
      iih = 1 + myid;
      nij_tmp = 1; ncoef_tmp = 0;
      ncoef = 0
      jjh = 1 + myid 
      max_col = 1;
      coeff = 0;
      inptr = 0;

*  allocate temp_coef for computing isotope shift
      if (last) then
        itmpu = 0;
        do im1 = 1,maxev
          if (leigen(im1,iblock)) itmpu = itmpu + 1;
        end do
        allocate(tmp_tc(itmpu*intptr(6)));
        allocate(tmp_coef(itmpu*intptr(6)));
        tmp_coef(1:itmpu*intptr(6)) = 0.0D0;
        tmp_tc(1:itmpu*intptr(6)) = 0.0D0;
      end if

      do 
       read(39) num,(coeff(j),j=1,num),(inptr(j),j=1,num)
         do i=1,num
*           ..test for next non-zero element
            if (ncoef+i.gt.ico(nijcurr)) then
                nijcurr = nijcurr + 1
*               .. have we also changed column?
                if (nijcurr.gt.jptr(max_col)) then
                   jjh = jjh + nprocs;
                   max_col = max_col+1
                end if
             end if
             iih = ihh(nijcurr)
             im1 = 0;
             do j = 1,maxev
               ioffw = (j-1)*ncfg
               if (leigen(j,iblock)) then
                 wcoef = eigst_weight(j,iblock)
                 W = wcoef*wt(ioffw+iih)*wt(ioffw+jjh)
*             print *,coef(inptr(i)) 
                 T = W*coeff(i) 
                 IF (IIH .NE. JJH) T = T+T
                 coef(inptr(i)) = coef(inptr(i)) + T
                 if (last) then
                   W0 = wt(ioffw+iih)*wt(ioffw+jjh);
                   T0 = W0*coeff(i);
                   IF (IIH .NE. JJH) T0 = T0 + T0;
                   itmp_s = (inptr(i)) + im1*intptr(6);
                   tmp_coef(itmp_s) = tmp_coef(itmp_s) + T0;
                   im1 = im1 + 1;
                 end if
                end if
             end do
*	     if (inptr(i) .gt. intptr(5)) then
*         print '(2I4,3F15.10,I6,F15.10,i6)',
*     :  iih,jjh,wt(iih+ioffw),wt(jjh+ioffw),coeff(i),
*     :  inptr(i), coef(inptr(i)),nijcurr
*	end if
         end do
         ncoef = ncoef + num
         if (num .lt. ncodim) exit
      end do

      if (last) then
         call MPI_ALLREDUCE(tmp_coef,tmp_tc,itmpu*intptr(6),
     :          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
         if (myid == 0) then
           im2 = 0;
           do im1 = 1,maxev
             if (leigen(im1,iblock)) then
               istart = 1+intptr(6)*im2;
               call prprty(im1,iblock,tmp_tc(istart),intptr(6),
     :                         isom_shift(im1,iblock));
               im2 = im2 + 1;
             end if
           end do
         end if
         deallocate(tmp_coef);
         deallocate(tmp_tc);
      end if

!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if (iblock == nblock) then   
        !call gdsummpi(coef, intptr(6), tmp)
        call mpi_allr_dp(coef,intptr(6))
*
*  *****  DEFINE SUM(I)
*
         IBEGIN = INTPTR(5)+1
         IEND = INTPTR(6)
         DO I = IBEGIN,IEND
           call unpacki(6,i,kv,iel1,iel2,iel3,iel4)
           IF (IEL1 .EQ. IEL2) SUM(IEL1) = -2*COEF(I)
         end do 
         if (myid == 0) then
*       print *, 'Node:',' intptr(6)',intptr(6), 'iblock',iblock
*       print '(4f20.16)', (coef(i),i=1,intptr(6))
*       print *, 'SUM', Sum(1:10)
*       print '(4F20.16)',SUm(1:20)
         end if
      end if
      END

