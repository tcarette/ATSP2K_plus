*
*     ------------------------------------------------------------------
*       R E A D W
*     ------------------------------------------------------------------
*
      SUBROUTINE READW(el,rel,ec,nnel,skip,onlydvd,iounew)
*	WITH RESTARTING
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NWD=60,NOD=220)
	logical onlydvd
*     .. breit commons
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
*     .. radial commons
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,iout,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk
C GG
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MSOO,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
C GG
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON /NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      CHARACTER*3 el(NWD),el1
      CHARACTER*6 atom, term
      LOGICAL found(NWD),check,rel,skip
      DOUBLE PRECISION PT(NOD)
       
      myid = 0; nprocs = 1;

      call initr
*     copy from breit common to radial common
!      print *, 'l(:)', l(1:nwf)
!      print *, 'el(:)', el(1:nwf)
*
      found(1:nwf) = .false.
*
*  *****  READ THE RADIAL FUNCTIONS
*
      if (myid == 0) then
         IW = 1
  14     READ(iuw,END=13) ATOM,Term,EL1,M,ZZ,ETI,EKI,AZD,(PT(J),J=1,M)
         call eptr(el,el1,i,*999)
         if (.not. found(i)) then
            p(1:m,i) = pt(1:j)
	    az(i) = azd
	    max(i) = m
            p(m+1:no,i) = D0
	    found(i) = .true.
            IW = IW+1
         end if
         Z = ZZ
         IF (IW.LE.NWF) GO TO 14
   13    close(iuw)
*     check that all functions were found
         check = .true.
         do i = 1,nwf
	   check = check .and. found(i)
         end do
         if (.not. check) then
	    write(iscw,*) ' Not all radial functions were found'
	    write(iscw,'(1X,A3,2X,L1)') (el(i),found(i),i=1,nwf)
	    stop
         end if
      end if 

*     .. initialize R and Rydberg constant
      call initm(1)
*     .. compute energy of core
      call ecore(el,ec,rel)

      write(iouj,'(2X,A6,A,F5.1,A,I3,A,I6)' ) ATOM,'  Z = ',Z,
     :        '  NEL = ', NNEL, '   NCFG = ',NCFG
      return
  999 write(iscw,*) ' Electron',el1,' not found in electron list'
      write(iscw,*) ' List is',  el
      go to 14
      end
