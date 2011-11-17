*****:******************************************************************
*
* --- This SUBROUTINE allocates memory for the Interaction 
*     matrix
*****:******************************************************************

      Subroutine alcmat(skip,ncfg,nume,n_maxcol,n_hmx,lhmx_memory,
     :                  lhmx_disk,jj,njv)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NTERMD=31,NWD=60)
*        
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,iout,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk
      POINTER (iqhmx,hmx(1)),(iqjptr,jptr(1)),(iqimx,imx(1)),
     :        (iqeigvec,eigvec(1)),
     :        (iqkval,dummy1(1)),(iqcptr,dummy2(1)),
     :        (qtm,tm(1)),(qtp,tp(1)),(qdiag,hii(1)),(qiwork,iwork(1))
      COMMON /MX/iqhmx, iqjptr, iqimx, iqeigvec, iqkval, iqcptr,
     :           qtm,qtp,qdiag,qiwork
      logical skip, lhmx_disk, lhmx_memory, ltmp
      character*2 idstring
      integer n_maxcol, n_hmx
      common/hmx_nze/nze_t,nze_c
      POINTER (IQLSP,LSP(ncfg))
      COMMON /NCFGS/ IQLSP,index(ntermd),IJK,
     :               termsh(ntermd),nterm
      POINTER (qh,h(1)),(qjan,jan(1))

*

      myid = 0; nprocs = 1;
      write(idstring,'(I2.2)') jj
      lhmx_disk = .true.
      lhmx_memory = .true.
      idisk = 0;
      write (iscw,*)
      if (myid == 0) write (iscw,*) ':: Allocating memory for ',
     :                              'Block 2J = ',jj 

      call alloc(iqjptr,ncfg,4); 
      jptr(1:ncfg) = 0
      nh = ncfg*3;
      call alloc(qh,nh,8); 
      h(1:nh) = 0.0
      call alloc(qjan,nh,4); 
      jan(1:nh) = 0


      if (iLS .eq. 0) call lsp_comp(ncfg,lsp,jj)
      call jmem(ncfg,j,nze,ind_jj,nij,istart,njv,nze_c,nze_t)

      print*, 'in jmem: nze_c=', nze_c, 'nze_t = ', nze_t
!      lim = min(20+nume,ncfg)
      lim = min(2*nume+40,ncfg)
      iworksz = (2*ncfg+lim+nume+10)*lim + nume
      call alloc(iqeigvec,iworksz,8); eigvec(1:iworksz) = 0.0
      iiwsz =  6*lim + nume
      call alloc(qtm,ncfg,8); tm(1:ncfg) = 0.0
      call alloc(qtp,ncfg,8); tp(1:ncfg) = 0.0
      call alloc(qdiag,ncfg,8); hii(1:ncfg) = 0.0
      call alloc(qiwork,iiwsz,4); iwork(1:iiwsz) = 0.0

      if (qtm.eq.0.or.qtp.eq.0.or.qdiag.eq.0.or.qiwork.eq.0) then
	   write(ISCW,*) 'Insufficient memory for diagonalization ',
     :                   ' failed at ##1(qtm,qtp,qdiag,qiwork)'
      endif


      call p_alloc(iqhmx,n_hmx,8); hmx(1:n_hmx) = 0.0
      call p_alloc(iqimx,n_hmx,4); imx(1:n_hmx) = 0

!       iqhmx = 0; ! DEBUG code
      if (iqhmx.eq.0.or.iqimx.eq.0) then
        lhmx_memory = .false.;
          write(iscw,*) 'BLOCK 2J = ' , jj, 'Insufficient memory hmx= ',
     :    n_hmx, '...dvdson ON DISK'
      endif

      if (.not.lhmx_memory) then
        if (iqhmx .ne. 0) call p_dalloc(iqhmx,n_hmx)
        if (iqimx .ne. 0) call p_dalloc(iqimx,n_hmx)
      end if

      if (.not.lhmx_memory) then
         idisk = 1;
         call p_alloc(iqhmx,nze_c,8); !hmx(1:nze_c) = 0.0
         call p_alloc(iqimx,nze_c,4); !imx(1:nze_c) = 0
         if (iqhmx.eq.0.or.iqimx.eq.0) then
           lhmx_disk = .false.;
           if (iqhmx.ne.0) call p_dalloc(iqhmx,n_hmx)
           if (iqimx.ne.0) call p_dalloc(iqimx,n_hmx)
           if (myid == 0) then
             write(iscw,*) 'Insufficient memory for on-disk',
     :                     ' diagonalization: Exitting..'
           end if
         endif
      end if

         if(lhmx_memory) then
             write(iscw,*) ':: IN MEMORY: Block 2J = ',jj, 'with ',
     :                      n_hmx, ' matrix', ' elements'
         end if

      if (skip) then
          if (idisk .eq. 1) iouhm = iouhn
      else
          if (idisk .eq. 1)
     :         open(iouhm,file='mat.jj='//idstring,
     :              status='unknown',form='unformatted')
      end if

      write(0,*)
      write(0,'(A50,i15)') 'MEMORY AND DISK USE FOR J = ', JJ
      write(0,'(A50,i15)') '->MATRIX NONZERO ELEMENTS = ', nze_t 
      write(0,'(A50,i15)') '->MAX NONZERO ELEMENTS FOR A COLUMN = ', 
     :                       nze_c
      if (lhmx_memory) then
         write(0,'(A50)') ' DVDSON IN MEMORY'
      else 
         write(0,'(A50)') ' INSUFFICIENT MEMORY, DVDSON ON DISK'
         write(0,'(A50,A15)') ' MATRIX SAVED IN FILE ',
     :                           'mat.jj='//idstring 
      end if
      write(0,*)

      end subroutine alcmat



