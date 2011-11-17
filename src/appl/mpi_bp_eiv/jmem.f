      subroutine jmem(ncfg,j,nze,ind_jj,nij,istart,njv,nze_c,nze_t)


      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NTERMD=31,NWD=60)
*     MPI stuff ***********************************************
*
        INCLUDE 'mpif.h'
        parameter (MAXPROC=9999)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
****************************************************************
      LOGICAL  SKIP
      POINTER (qh,h(ncfg,3)),(qjan,jan(ncfg,3))
      COMMON /buffer/qh, qjan, nrow(3), iflag(3)
      POINTER (IQLSP,LSP(ncfg))
      COMMON /NCFGS/ IQLSP,index(ntermd),IJK,nterm,termsh(ntermd)
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,iout,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk
*
      INTEGER cptr
      POINTER (iqhmx,hmx(1)),(iqjptr,jptr(1)),(iqimx,imx(1)),
     :        (iqeigvec,eigvec(1)),(iqkval,kval(1)),(iqcptr,cptr(1)),
     :        (qtm,tm(1)),(qtp,tp(1)),(qdiag,hii(1)),(qiwork,iwork(1))
      COMMON /MX/iqhmx, iqjptr, iqimx, iqeigvec, iqkval, iqcptr,
     :           qtm,qtp,qdiag,qiwork
      INTEGER   q(8), ipos(3),iselec(NTERMD)
      LOGICAL iupper,hiend
      DATA iupper/.false./
      DATA crite,critc,critr,trhold,ortho/
     :       0.d-17,1.d-8,1.d-8,1.d-10,1.d-15/
!     :       0.d-17,1.d-12,1.d-20,1.d0,1.d0/
      pointer(pflsj,flsj(nterm,nterm,2,njv))
      double precision flsj;
      logical lhmx_memory,lhmx_disk
      integer ind_jj
*
      rewind(iouhn)
      rewind(iouhz)
      rewind(iouhs)
      mycol = 0;
      nze_t = 0;
      nze_c = 0;
      nij = 0
      do j = myid+1,ncfg,nprocs
         mycol = mycol + 1
         if (idisk .eq. 1) then
	    nij = 0
	 end if
	 if (mod(lsp(j),2) .eq. 0) then
*            .. csf does not contribute
	    nij = nij+1
	    jptr(mycol) = nij
	    if (j .eq. 1) then
	       read(iouhn) jb,m1,(h(i,1),i=1,m1),
     :                           (jan(i,1),i=1,m1)
	    else
	       read(iouhn) jb
	    end if
	    read(iouhz) jb
	    read(iouhs) jb
 	 else
            ipos(1:3) = 0
	    ntermj = mod(lsp(j),64)/2
	    read(iouhn) jb,m1,(h(i,1),i=1,m1),
     :                        (jan(i,1),i=1,m1)
	    read(iouhz) jb,m2,(h(i,2),i=1,m2),
     :                        (jan(i,2),i=1,m2)
            read(iouhs) jb,m3,(h(i,3),i=1,m3),
     :                        (jan(i,3),i=1,m3)
            call advance(1,m1,ipos,ncfg)
	    call advance(2,m2,ipos,ncfg)
	    call advance(3,m3,ipos,ncfg)
*           .. while there is more data
*             find which file has lowest row
            iipos = min(ipos(1),ipos(2),ipos(3))
            do while (iipos <= ncfg)
*             .. there is more data
	       nij = nij + 1
	       irow = min(nrow(1),nrow(2),nrow(3))
	       if (nrow(1) .eq. irow) then
	          call advance(1,m1,ipos,ncfg)
	       end if
	       if (nrow(2) .eq. irow) then
	          ntermi = mod(lsp(j),64)/2
	          call advance(2,m2,ipos,ncfg)
	       end if
	       if (nrow(3) .eq. irow) then
	          ntermi = mod(lsp(j),64)/2
	          call advance(3,m3,ipos,ncfg)
	       end if
               iipos = min(ipos(1),ipos(2),ipos(3));
            end do
	 end if
         nze_ctem = nij - nze_t
         if (nze_ctem > nze_c) nze_c = nze_ctem
         nze_t = nij;
      end do
      return
      end subroutine
