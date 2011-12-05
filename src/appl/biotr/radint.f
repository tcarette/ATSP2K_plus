*
*     ------------------------------------------------------------
*     R A D I N T
*     ------------------------------------------------------------
*
* --- Calculates the radial overlap and transition integrals
*
      SUBROUTINE radint(im)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=128,NOD=220)
*
      LOGICAL rel,vok
      CHARACTER el*3,im*1
      character*3 elras,elrasi,elrasf
      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON /DBG  /IBUGM
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON /NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7),el(nwd)
      POINTER(IQRL,RLINT(MAXORB,1)),(IQRV,RVINT(MAXORB,1)),
     :       (IQOV,OVLP(MAXORB,1))
      COMMON /RDINT/IQRL,IQRV,IQOV
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      COMMON /RAS/ ninac(11,2),nras1(11,10,2),nras2(11,2),
     :             nras3(11,10,2),NGAS1(2),NGAS3(2),
     :             nl(11,2),elras(2*nwd),elrasi(nwd),elrasf(nwd),
     :             itab(nwd,2),lmax(2)
*
* --- allocate the /RDINT/
*
      print*,' maxorb = ',maxorb
      if(ibugm.ne.0)print*,' iqrl    allocation: maxorb**2 = ',maxorb**2
      call alloc(iqrl,maxorb*maxorb,8)
      if(ibugm.ne.0)print*,' iqov    allocation: maxorb**2 = ',maxorb**2
      call alloc(iqov,maxorb*maxorb,8)
      if(ibugm.ne.0)print*,' iqrv    allocation: maxorb**2 = ',maxorb**2
      call alloc(iqrv,maxorb*maxorb,8)
      if(ibugm.ne.0) print*,' transition ',im,lam
!      do 1 i = 1,maxorb
!        do 2 j = 1,maxorb
!           if(itab(i,1).eq.0.or.itab(j,2).eq.0) go to 2
!           il = iabs(l(itab(i,1))-l(itab(j,2)))
!           ovlp(i,j) = quadr(itab(i,1),itab(j,2),0)
!           if (im.eq.'E') then
!             rlint(i,j) = quadr(itab(i,1),itab(j,2),lam)
!*            if (.not.rel) then
!               if (lam.eq.1.and.il.eq.1)
!     :              rvint(i,j) = grad(itab(i,1),itab(j,2))
!               if (lam.eq.2.and.(il.eq.0.or.il.eq.2))
!     :              rvint(i,j) = grad2(itab(i,1),itab(j,2))
!*            end if
!           else if (im.eq.'M') then
!             rlint(i,j) = quadr(itab(i,1),itab(j,2),lam-1)
!           end if
!    2   continue
!    1 continue
       return
       end
