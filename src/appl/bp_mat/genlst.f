*------------------------------------------------------------------------
*        G E N L S T
*------------------------------------------------------------------------
*
*     Three matrices are created -- LS (non-relativistic or with
*     relativistic shift), LSJ1 (Breit-Pauli operators other
*     than spin-spin), LSJ2 (spin-spin)
*
*     An LS calculation only uses the first whereas LSJ creates
*     J-dependent matrices from all three.
*
      SUBROUTINE GENLST(ncfg,jptr)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER (LSDIM=30000)
      POINTER (qcn,cn(lsdim)),(qinptr,inptr(lsdim)),
     :        (qnijptr,nijptr(lsdim)),(qjanc,janc(lsdim)),
     :        (qjbn,jbn(lsdim)),(qintptr,intptr(0:2*lmax+1,7)),
     :        (qpackn,ipackn(1)),(qlused,lused(1)),(qico,ico(ncfg)),
     :        (qvalue,value(1))
*      COMMON /buffer/qcn,qinptr,qpackn,qlused,qintptr,lmax,qnijptr,
*     :               qjan,qjbn,qico,qvalue
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,iout,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk

      COMMON /TABLE/qintptr,qpackn,qvalue,lmax
      POINTER  (qh,h(ncfg,3)), (qjan,jan(ncfg,3))
      DIMENSION jptr(*)
      DIMENSION idum(6),nrow(3),iflag(3)
 
        INTEGER maxorb,lmax,n

      call alloc(qh, 3*ncfg, 8)
!      h(1:ncfg,1:3) = 0.0
      call alloc(qjan, 3*ncfg,4); !jan(1:ncfg,1:3) = 0
      call alloc(qico,ncfg,4)   ; ico(1:ncfg) = 0
      call alloc(qjanc,ncfg, 4); janc(1:ncfg) = 0
      call alloc(qcn, lsdim, 8) ; 
      call alloc(qinptr,lsdim,4); 

*    
*     n  - counts coefficients processed from current record
*     nij- counts coefficients from all records so far
*     lij- counts matrix elements processed  for current column
*
*     TODO: IORBORB=1 is not correct,  because icase=9 (in 
*           mpi_breit gives a pointer inptr > intptr(0,4) and
*           then lcase = 2, while it has to stay 1 
*     May 11, 2000, git
* 

      myid = 0; nprocs = 1;
      n = 0
      nij = 1
      rewind(50)
      read (50) new,(cn(j),j=1,new),(inptr(j),j=1,new)
      jb_count = 1
      Do JB = 1+myid,ncfg,nprocs
         if(mod(jb,100).eq.0) write(0,'(A,I5)') '   jb = ',jb
         if(jb == ncfg) write(0,'(A,I5)') '   jb = ',jb
*        initialize arrays for a column
         nrow(1:3) = 0
         jan(1:ncfg,1:3) = 0
         h(1:ncfg,1:3) = 0.d0
         ico(1:ncfg) = 0
	 iflag(1:3)=0
         lij =  1
*       read info on matrix elements for the column
	 read (51) nih, (janc(i),i=1,nih)
	 read (52) nih, (ico(i),i=1,nih)
	 Do
           n = n+1;
	   if (n .gt. new) then
*           .. more data is needed
             read (50) new,(cn(j),j=1,new),(inptr(j),j=1,new)
	     n = 1
	   end if
*         .. determine case
           if (inptr(n) < 0) then
             lcase = 3
           else if (inptr(n) == 0) then
             print*,'inptr(',n,') can''t be 0'; lcase = 2
           else if (inptr(n).le.intptr(1,4))  then
	     lcase = 1
           else 
	     lcase = 2
           end if
*         .. determine matrix element
	   if ((nij > ico(lij))) then
*           .. we have new matrix element; set iflags to zero
	     lij = lij + 1
	     iflag(1:3) = 0
	   end if 
	   if (iflag(lcase) .eq. 0) then
	     nrow(lcase) = nrow(lcase) + 1
	     jan(nrow(lcase),lcase) = janc(lij)
	     iflag(lcase) = 1
           end if
*         .. add contribution to right matrix element
           v = value(abs(inptr(n)))
           h(nrow(lcase),lcase) = h(nrow(lcase),lcase) + cn(n)*v
           hh = h(nrow(lcase),lcase)
!        write(*,'(4i4,3F20.12)') n,jb,lcase,nrow(lcase),cn(n),v,hh
           nij = nij + 1
	   if (nij .eq. jptr(jb_count) + 1) exit
*          if (lij .eq. nih) exit
	 END DO
	 write(iouhn) jb,nrow(1),(h(i,1),i=1,nrow(1)),
     :                (jan(i,1),i=1,nrow(1))
	 write(iouhz) jb,nrow(2),(h(i,2),i=1,nrow(2)),
     :                (jan(i,2),i=1,nrow(2))
	 write(iouhs) jb,nrow(3),(h(i,3),i=1,nrow(3)),
     :                (jan(i,3),i=1,nrow(3))
!         write(*,'(A10)') 'iouhn'
!         write(*,'(2I5)') jb,nrow(1)
!         write(*,'(6F14.8)') (h(i,1),i=1,nrow(1))
!         write(*,'(6i14)')   (jan(i,1),i=1,nrow(1))
!         write(*,'(A10)') 'iouhz'
!         write(*,'(2I5)') jb,nrow(2)
!         write(*,'(6F14.8)') (h(i,2),i=1,nrow(2))
!         write(*,'(6i14)')   (jan(i,2),i=1,nrow(2))
!         write(*,'(A10)') 'iouhs'
!         write(*,'(2I5)') jb,nrow(3)
!         write(*,'(6F14.8)') (h(i,3),i=1,nrow(3))
!         write(*,'(6i14)')   (jan(i,3),i=1,nrow(3))
          jb_count = jb_count + 1
       END DO
       close(iouhn)
       close(iouhz)
       close(iouhs)
      END
