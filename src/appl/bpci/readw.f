*
*     ------------------------------------------------------------------
*       R E A D W
*     ------------------------------------------------------------------
*
      SUBROUTINE READW(rel,ec,nnel,skip,onlydvd,iounew)
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
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,
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
      nwf = nclosd + maxorb
      call initr
*     copy from breit common to radial common
      do 1 i = 1,nclosd
	 l(i) = ljclsd(i)
	 write(el(i),'(A3)') iajcld(i)
	 found(i) = .false.
    1 continue
      do 2 i = nclosd+1,nwf
	 l(i) = ljcomp(i-nclosd)
	 n(i) = njcomp(i-nclosd)
	 write(el(i),'(A3)') iajcmp(i-nclosd)
	 found(i) = .false.
    2 continue
*
*
*  *****  READ THE RADIAL FUNCTIONS
*
      IW = 1
  14  READ(iuw,END=13) ATOM,Term,EL1,M,ZZ,ETI,EKI,AZD,(PT(J),J=1,M)
      call eptr(el,el1,i,*999)
      if (.not. found(i)) then
        Do 15 j= 1,m
 	  p(j,i) = pt(J)
  15    continue
	az(i) = azd
	max(i) = m
        DO 24 J = M+1,NO
   24     P(J,I) = D0
	found(i) = .true.
        IW = IW+1
      end if
      Z = ZZ
      IF (IW.LE.NWF) GO TO 14
   13 close(iuw)
*     check that all functions were found
      check = .true.
      do 16 i = 1,nwf
	 check = check .and. found(i)
  16  continue
      if (.not. check) then
	 write(iscw,*) ' Not all radial functions were found'
	 write(iscw,'(1X,A3,2X,L1)') (el(i),found(i),i=1,nwf)
	 stop
      end if
*     .. initialize R and Rydberg constant
      call initm(1)
*     .. compute energy of core
      call ecore(el,ec,rel)
	if (.not.onlydvd) close(unit=iounew, status='delete')
      if (skip) then 
	if (.not.onlydvd) iounew = ioul
	close(unit=iouj, status='delete')
      else
	if (.not.onlydvd) iounew = iouj
	close(unit=ioul, status='delete')
      end if
      write(iounew,'(2X,A6,A,F5.1,A,I3,A,I6)' ) ATOM,'  Z = ',Z,
     :        '  NEL = ', NNEL, '   NCFG = ',NCFG
      return
  999 write(iscw,*) ' Electron',el1,' not found in electron list'
      write(iscw,*) ' List is',  el
      go to 14
      end
