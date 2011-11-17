*****:******************************************************************
*
* --- This SUBROUTINE allocates memory for the Interaction 
*     matrix
*
*****:******************************************************************
*
      Subroutine alcmat(ncfg,nze,nume)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk
      POINTER (iqhmx,hmx(1)),(iqjptr,jptr(1)),(iqimx,imx(1)),
     :        (iqeigvec,eigvec(1)),
     :        (iqkval,dummy1(1)),(iqcptr,dummy2(1)),
     :        (qtm,tm(1)),(qtp,tp(1)),(qdiag,hii(1)),(qiwork,iwork(1))
      COMMON /MX/iqhmx, iqjptr, iqimx, iqeigvec, iqkval, iqcptr,
     :           qtm,qtp,qdiag,qiwork
*
!      lim = min(20,ncfg)
      lim = min(2*nume+40,ncfg)
      iworksz = (2*ncfg+lim+nume+10)*lim + nume
      call alloc(iqeigvec,iworksz,8)
      iiwsz =  6*lim + nume
      call alloc(qtm,ncfg,8)
      call alloc(qtp,ncfg,8)
      call alloc(qdiag,ncfg,8)
      call alloc(qiwork,iiwsz,4)
*
*     .. check if all pointers are OK
      if (qtm.eq.0 .or. qtp.eq.0 .or. qdiag.eq.0 .or. qiwork.eq.0)
     :      STOP 'Insufficient memory for diagonalization'
      call alloc(iqhmx,nze,8)
      call alloc(iqjptr,ncfg,4)
      call alloc(iqimx,nze,4)
*     idisk = 0
*     .. testing
*     call dalloc(iqhmx,nze)
*     iqhmx=0
      if (iqhmx.eq.0 .or. iqjptr.eq.0 .or. iqimx .eq.0) then
	write(iscw,*) 'Disk version used for diagonalization'
        if (iqhmx .ne. 0) call dalloc(iqhmx,nze)
        if (iqjptr .ne. 0) call dalloc(iqjptr,ncfg)
        if (iqimx .ne. 0) call dalloc(iqimx,nze)
	nze = ncfg
        call alloc(iqhmx,nze,8)
        call alloc(iqimx,nze,4)
	if (iqhmx.eq.0 .or. iqimx.eq.0) 
     :      STOP 'Insufficient memory for diagonalization'
	idisk = 1
      end if
      Print *, 'Alcmat allocations for idisk=',idisk, 'nze =',nze
      end
