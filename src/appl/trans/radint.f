*     ------------------------------------------------------------
*     R A D I N T
*     ------------------------------------------------------------
*
* --- Calculates the radial overlap and transition integrals
*
      SUBROUTINE radint(im)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=80,NOD=220)
*
      LOGICAL rel,vok
      CHARACTER el*3,im*1
      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON /NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7),el(nwd)
      POINTER(IQRL,RLINT(MAXORB,1)),(IQRV,RVINT(MAXORB,1)),
     :       (IQOV,OVLP(MAXORB,1))
      COMMON /RDINT/IQRL,IQRV,IQOV
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1))
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,MAXORB
      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
*
* --- allocate the /RDINT/
*
      print*,' maxorb = ',maxorb
      if(ibug1.ne.0)print*,' iqrl    allocation: maxorb**2 = ',maxorb**2
      call alloc(iqrl,maxorb*maxorb,8)
      if(ibug1.ne.0)print*,' iqov    allocation: maxorb**2 = ',maxorb**2
      call alloc(iqov,maxorb*maxorb,8)
      if(ibug1.ne.0)print*,' iqrv    allocation: maxorb**2 = ',maxorb**2
      call alloc(iqrv,maxorb*maxorb,8)
      if(ibug1.ne.0) print*,' transition ',im,lam
      do 1 i = 1,maxorb
        do 2 j = 1,maxorb
           il = iabs(l(i)-l(j))
Cdbg       print*,' in radint, i = ',i,' j = ',j
           ovlp(i,j) = quadr(i,j,0)
Cdbg  print*,' ovlp = ',ovlp(i,j)
           if (im.eq.'E') then
             rlint(i,j) = quadr(i,j,lam)
Cdbg  print*,' rlint = ',rlint(i,j)
             if (.not.rel) then
               if (lam.eq.1.and.il.eq.1) rvint(i,j) = grad(i,j)
               if (lam.eq.2.and.(il.eq.0.or.il.eq.2))
     :              rvint(i,j) = grad2(i,j)
             end if
           else if (im.eq.'M') then
             rlint(i,j) = quadr(i,j,lam-1)
Cdbg  print*,' rlint = ',rlint(i,j)
           end if
    2   continue
    1 continue
       return
       end
