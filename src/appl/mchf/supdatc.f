*-----------------------------------------------------------------------
*               U P D A T C 
*-----------------------------------------------------------------------
*
*     Evaluate all coeffients of integrals.  This involves a
*     sequential pass through all the c.lst data.
*
*     
      SUBROUTINE UPDATC_memory_all(iblock,wt,ncfg,maxev,jptr,ihh,ico,
     :           eigst_weight,max_col,nijcurr,cf_tot,last,isom_shift)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=94,NWC=20,NOD=220,NOFFD=800)
      parameter(MEIG=20,MTERM=20)
      POINTER (pkval,kval(1)),(pcoef,coef(1)),(pih,ih(1)),(pjh,jh(1)),
     :     (pcoeff,coeff(1)),(pnijptr,nijptr(1)),(pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PCOEF,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
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
      double precision eigst_weight(meig,mterm),isom_shift(MEIG,nblock)
      integer cf_tot(nblock)
      logical leigen(meig,mterm),lguess,last
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :   nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm),nze_max(mterm)
      character term_bl(mterm)*3
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess,nze_max
      double precision, allocatable, dimension(:) :: tmp_coef
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
        allocate(tmp_coef(itmpu*idim));      
        tmp_coef(1:itmpu*idim) = 0.0;
      end if
        
      nijcurr = 1; 
      iih = 1 ;
      nij_tmp = 1; 
      jjh = 1  ;
      max_col = 1;
      n_count_tmp = 0; 
      if (iblock == 1) ncoef = 0;
*        .. add contribution from this record to coef
*        .. make use of the fact that all coefficients are read in
*           increasing order of the nze elements of the matrix
        do i = 1, cf_tot(iblock);
*           ..test for next non-zero element
            n_count_tmp = ncoef+i
            if (i.gt.ico(nijcurr)) then
                nijcurr = nijcurr + 1
*               .. have we also changed column?
                if (nijcurr.gt.jptr(max_col)) then
                   jjh = jjh + 1;
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
                   itmp_s = (inptr(n_count_tmp))+idim*im1;
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
      
ctc
*      print *, 'ncoef=', ncoef
ctc

      if (last) then
         im2 = 0;
         do im1 = 1,maxev
           if (leigen(im1,iblock)) then
             istart = 1+idim*im2;
             call prprty(im1,iblock,tmp_coef(istart),
     :                         isom_shift(im1,iblock));
             im2 = im2 + 1;
           end if
         end do
         deallocate(tmp_coef); 
      end if
*
*  *****  DEFINE SUM(I)
*
         IBEGIN = INTPTR(5)+1
         IEND = INTPTR(6)

         DO I = IBEGIN,IEND
           call unpacki(6,i,kv,iel1,iel2,iel3,iel4)
           IF (IEL1 .EQ. IEL2) SUM(IEL1) = -2*COEF(I)
         end do 

*       print *, 'Node:',' intptr(6)',intptr(6), 'iblock',iblock
*       print '(4f20.16)', (coef(i),i=1,intptr(6))
*       print *, 'SUM', Sum(1:10)
*       print '(4F20.16)',SUm(1:20)
      END

****************************************************************************
      SUBROUTINE UPDATC_disk_ico(iblock,wt,ncfg,maxev,jptr,ihh,
     :        eigst_weight,max_col,nijcurr,last,isom_shift)
cgd  :        eigst_weight,max_col,nijcurr,cf_tot,last,isom_shift)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=94,NWC=20,NOD=220,NOFFD=800)
      parameter(MEIG=20,MTERM=20)
      POINTER (pkval,kval(1)),(pcoef,coef(1)),(pih,ih(1)),(pjh,jh(1)),
     :  (pcoeff,coeff(1)),(pnijptr,nijptr(1)),(pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PCOEF,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
      INTEGER       IN,OUT,PRI,OUC,OUF,OUH,OUD
      COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,iud,OUC,OUF,OUD,OUH,ISCW
      COMMON/LABEL/ EL(NWD),ATOM,TERM
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :  (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :  (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :           QMETH,QIEPTR

      COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
      dimension wt(ncfg*maxev), jptr(*), ihh(*)
      character string*64
      double precision eigst_weight(meig,mterm),isom_shift(meig,nblock)
      integer cf_tot(nblock)
      logical leigen(meig,mterm),lguess,last
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm),
     :        nze_max(mterm)
      character term_bl(mterm)*3
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess,nze_max
cgd   integer ico(nze_max(iblock))
      integer ico(ncodim)
      double precision, allocatable, dimension(:) :: tmp_coef

	  if ((iblock.eq.1)) then
        nint = intptr(6)
        coef(1:nint) = 0
      end if
      ico = 0.0
      nijcurr = 1; 
      nij_count = 1;
      iih = 1 ;
      nij_tmp = 1; ncoef_tmp = 0;
      ncoef = 0
      jjh = 1  
      max_col = 1;
!      coeff = 0;
!      inptr = 0;
      read (12) nico, (ico(jj), jj = 1,nico);
     
*  allocate temp_coef for computing isotope shift
      if (last) then
        itmpu = 0;
        do im1 = 1,maxev
          if (leigen(im1,iblock)) itmpu = itmpu + 1;
        end do
        allocate(tmp_coef(itmpu*idim));
        tmp_coef(1:itmpu*idim) = 0.0;
      end if

      do 
       read(39) num,(coeff(j),j=1,num),(inptr(j),j=1,num)
         do i=1,num
            if (ncoef+i.gt.ico(nij_count)) then
                nijcurr = nijcurr + 1
                nij_count = nij_count + 1;
                if (nijcurr.gt.jptr(max_col)) then
                   if (jjh < ncfg ) then
                      read(12) nico ,(ico(jj),jj=1,nico)
                      nij_count = 1
                   end if
                   jjh = jjh + 1;
                   max_col = max_col+1
                end if
             end if
             iih = ihh(nijcurr)
             im1 = 0;        
             do j = 1,maxev
                ioffw = (j-1)*ncfg
                if (leigen(j,iblock)) then
                   wcoef = eigst_weight(j,iblock)
                   W =  wcoef*wt(ioffw+iih)*wt(ioffw+jjh)
                   T = W*coeff(i) 
                   IF (IIH .NE. JJH) T = T+T
                   coef(inptr(i)) = coef(inptr(i)) + T
                   if (last) then
                     W0 = wt(ioffw+iih)*wt(ioffw+jjh);
                     T0 = W0*coeff(i);
                     IF (IIH .NE. JJH) T0 = T0 + T0;
                     itmp_s = (inptr(i))+idim*im1;
                     tmp_coef(itmp_s) = tmp_coef(itmp_s) + T0;
                     im1 = im1 + 1;
                   end if
                end if
             end do
         end do
         ncoef = ncoef + num
         if (num .lt. ncodim) exit
      end do
     
      if (last) then
         im2 = 0;
         do im1 = 1,maxev
           if (leigen(im1,iblock)) then
             istart = 1+idim*im2;
             call prprty(im1,iblock,tmp_coef(istart),
     :                         isom_shift(im1,iblock));
             im2 = im2 + 1;
           end if
         end do
         deallocate(tmp_coef);
      end if

      IBEGIN = INTPTR(5)+1
      IEND = INTPTR(6)
      DO I = IBEGIN,IEND
        call unpacki(6,i,kv,iel1,iel2,iel3,iel4)
        IF (IEL1 .EQ. IEL2) SUM(IEL1) = -2*COEF(I)
      end do 

      END

**********************************************************************
      SUBROUTINE UPDATC_disk_clst(iblock,wt,ncfg,maxev,jptr,ihh,ico,
     :          eigst_weight,max_col,nijcurr,cf_tot,last,isom_shift)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=94,NWC=20,NOD=220,NOFFD=800)
      parameter(MEIG=20,MTERM=20)
      POINTER (pkval,kval(1)),(pcoef,coef(1)),(pih,ih(1)),(pjh,jh(1)),
     :   (pcoeff,coeff(1)),(pnijptr,nijptr(1)),(pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PCOEF,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
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
      double precision isom_shift(MEIG,nblock);
      character string*64

      double precision eigst_weight(meig,mterm)
      integer cf_tot(nblock)
      logical leigen(meig,mterm),lguess,last
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm),
     :        nze_max(mterm)
      character term_bl(mterm)*3
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess,nze_max
      double precision, allocatable, dimension(:) :: tmp_coef

      if ((iblock.eq.1)) then
        nint = intptr(6)
        coef(1:nint) = 0
      end if
      nijcurr = 1; 
      iih = 1 ;
      nij_tmp = 1; ncoef_tmp = 0;
      ncoef = 0
      jjh = 1  
      max_col = 1;
      coeff = 0;
      inptr = 0;
     
*  allocate temp_coef for computing isotope shift
      if (last) then
        itmpu = 0;
        do im1 = 1,maxev
          if (leigen(im1,iblock)) itmpu = itmpu + 1;
        end do
        allocate(tmp_coef(itmpu*idim));
        tmp_coef(1:itmpu*idim) = 0.0; !print *, itmpu*idim;
      end if

      do 
       read(39) num,(coeff(j),j=1,num),(inptr(j),j=1,num)
         do i=1,num
            if (ncoef+i.gt.ico(nijcurr)) then
                nijcurr = nijcurr + 1
                if (Nijcurr.gt.jptr(max_col)) then
                   jjh = jjh + 1;
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
                  T = W*coeff(i) 
                  IF (IIH .NE. JJH) T = T+T
                  coef(inptr(i)) = coef(inptr(i)) + T
                   if (last) then
                     W0 = wt(ioffw+iih)*wt(ioffw+jjh);
                     T0 = W0*coeff(i);
                     IF (IIH .NE. JJH) T0 = T0 + T0;
                     itmp_s = (inptr(i))+idim*im1;
                     tmp_coef(itmp_s) = tmp_coef(itmp_s) + T0;
                     im1 = im1 + 1;
                   end if
                end if
             end do
         end do
         ncoef = ncoef + num
         if (num .lt. ncodim) exit
      end do

      if (last) then
         im2 = 0;
         do im1 = 1,maxev
           if (leigen(im1,iblock)) then
             istart = 1+idim*im2;
             call prprty(im1,iblock,tmp_coef(istart),
     :                         isom_shift(im1,iblock));
             im2 = im2 + 1;
           end if
         end do
         deallocate(tmp_coef);
      end if

      IBEGIN = INTPTR(5)+1
      IEND = INTPTR(6)
      DO I = IBEGIN,IEND
        call unpacki(6,i,kv,iel1,iel2,iel3,iel4)
        IF (IEL1 .EQ. IEL2) SUM(IEL1) = -2*COEF(I)
      end do 
      END

**********************************************************************
      SUBROUTINE UPDATC_disk_ih(iblock,wt,ncfg,maxev,jptr,ihh,
     :          eigst_weight,max_col,nijcurr,last,isom_shift)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=94,NWC=20,NOD=220,NOFFD=800)
      parameter(MEIG=20,MTERM=20)
      POINTER (pkval,kval(1)),(pcoef,coef(1)),(pih,ih(1)),(pjh,jh(1)),
     :  (pcoeff,coeff(1)),(pnijptr,nijptr(1)),(pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PCOEF,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
      INTEGER       IN,OUT,PRI,OUC,OUF,OUH,OUD
      COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,iud,OUC,OUF,OUD,OUH,ISCW
      COMMON/LABEL/ EL(NWD),ATOM,TERM
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :  (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :  (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :           QMETH,QIEPTR

      COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
      dimension wt(ncfg*maxev),jptr(*),ihh(*)
      double precision isom_shift(MEIG,nblock);
      character string*64
      double precision eigst_weight(meig,mterm)
      logical leigen(meig,mterm),lguess,last
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm),
     :        nze_max(mterm)
      character term_bl(mterm)*3
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess,nze_max
      integer ico(nze_max(iblock))
      character*2 ih_file
      double precision, allocatable, dimension(:) :: tmp_coef
   
      if ((iblock.eq.1)) then
        nint = intptr(6)
        coef(1:nint) = 0
      end if
      ico = 0.0
      nijcurr = 1; 
      nij_count = 1;
      iih = 1 ;
      nij_tmp = 1; ncoef_tmp = 0;
      ncoef = 0
      jjh = 1  
      max_col = 1;
      coeff = 0;
      inptr = 0;
*  allocate temp_coef for computing isotope shift
      if (last) then
        itmpu = 0;
        do im1 = 1,maxev
          if (leigen(im1,iblock)) itmpu = itmpu + 1;
        end do
        allocate(tmp_coef(itmpu*idim));
        tmp_coef(1:itmpu*idim) = 0.0;
      end if

      write(ih_file,'(I2.2)') iblock
      read (12) nico, (ico(jj), jj = 1,nico);
      open(unit=11,file='ih.'//ih_file//'.lst',status='old',
     :     form='unformatted');
      read(11) nze_ih, (ihh(jj),jj=1,nze_ih)
     
      do 
       read(39) num,(coeff(j),j=1,num),(inptr(j),j=1,num)
         do i=1,num
            if (ncoef+i.gt.ico(nij_count)) then
                nijcurr = nijcurr + 1
                nij_count = nij_count + 1;
                if (nijcurr.gt.jptr(max_col)) then
                   if (jjh < ncfg ) then
                      read(12) nico ,(ico(jj),jj=1,nico)
                      read(11) nze_ih,(ihh(jj),jj=1,nze_ih)
                      nij_count = 1
                   end if
                   jjh = jjh + 1;
                   max_col = max_col+1
                end if
             end if
             iih = ihh(nij_count)
             im1 = 0
             do j = 1,maxev
               ioffw = (j-1)*ncfg
               if (leigen(j,iblock)) then
                 wcoef = eigst_weight(j,iblock)
                 W = wcoef*wt(ioffw+iih)*wt(ioffw+jjh)
                 T = W*coeff(i) 
                 IF (IIH .NE. JJH) T = T+T
                 coef(inptr(i)) = coef(inptr(i)) + T
                 if (last) then
                   W0 = wt(ioffw+iih)*wt(ioffw+jjh);
                   T0 = W0*coeff(i);
                   IF (IIH .NE. JJH) T0 = T0 + T0;
                   itmp_s = (inptr(i))+idim*im1;
                   tmp_coef(itmp_s) = tmp_coef(itmp_s) + T0;
                   im1 = im1 + 1;
                 end if
               end if
             end do
         end do
         ncoef = ncoef + num
         if (num .lt. ncodim) exit
      end do

      if (last) then
         im2 = 0;
         do im1 = 1,maxev
           if (leigen(im1,iblock)) then
             istart = 1+idim*im2;
             call prprty(im1,iblock,tmp_coef(istart),
     :                         isom_shift(im1,iblock));
             im2 = im2 + 1;
           end if
         end do
         deallocate(tmp_coef);
      end if

      IBEGIN = INTPTR(5)+1
      IEND = INTPTR(6)
      DO I = IBEGIN,IEND
        call unpacki(6,i,kv,iel1,iel2,iel3,iel4)
        IF (IEL1 .EQ. IEL2) SUM(IEL1) = -2*COEF(I)
      end do 
      close(11)
      END

