*     ------------------------------------------------------------------
*       R E A D W 2
*     ------------------------------------------------------------------
*
      SUBROUTINE readw2
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NWD=80,NOD=220)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCF
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,MAXORB
      COMMON /NOR/NCOM,NCLOSI,NCLOSF,NORBI,NORBF,IWAR
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      CHARACTER*3 el,el1
      CHARACTER*3 eli(nwd),elf(nwd)
      COMMON /NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7),el(nwd)
      CHARACTER*6 atom, term
      LOGICAL found(NWD,2),check
      DOUBLE PRECISION PT(NOD)
*
* --- allocate the /NEL/
*
      if(ibug1.ne.0) print*,' iqp     allocation: nod*maxorb = ',
     :               nod*maxorb
      call alloc(iqp,nod*maxorb,8)
      if(ibug1.ne.0) print*,' iqn     allocation: maxorb = ',maxorb
      call alloc(iqn,maxorb,4)
      if(ibug1.ne.0) print*,' iql     allocation: maxorb = ',maxorb
      call alloc(iql,maxorb,4)
      if(ibug1.ne.0) print*,' iqaz    allocation: maxorb = ',maxorb
      call alloc(iqaz,maxorb,8)
      if(ibug1.ne.0) print*,' iqmax   allocation: maxorb = ',maxorb
      call alloc(iqmax,maxorb,4)
*
* --- read the radial wavefunctions on units iuw(1) and iuw(2)
*     and sort them according the IAJCMP order
*
*
* --- transfer NDIMS/NON30 -> PARAM
*
      nwf = maxorb
      ncfg = ncf
      print*,' in readw2, nwf = ',nwf,' ncfg = ',ncfg
*
      do 1 i = 1,nwf
        l(i) = ljcomp(i)
        n(i) = njcomp(i)
        write(el(i),'(A3)') iajcmp(i)
        if (i.le.(ncom+norbi)) then
          write(eli(i),'(A3)') iajcmp(i)
          found(i,1) = .false.
        else
          k = i - ncom - norbi
          write(elf(k),'(A3)') iajcmp(i)
          found(k,2) = .false.
        end if
    1 continue
      print*,' el     = ',(el(i),i=1,nwf)
      print*,' eli    = ',(eli(i),i=1,(ncom+norbi))
      print*,' elf    = ',(elf(i),i=1,norbf)
*
*  *****  READ THE RADIAL FUNCTIONS
*  *****  Orthogonal case
*
      IFL = IUW(1)
      IW = 1
  314 READ(IFL,END=400)
     :   ATOM,Term,EL1,M,ZZ,ETI,EKI,AZD,(PT(J),J=1,M)
         call eptr(el,el1,i,*399)
      if (found(i,1)) go to 399
        Do 315 j= 1,m
            p(j,i) = pt(j)
  315   continue
	az(i) = azd
	max(i) = m
        DO 324 J = M+1,NO
  324       P(J,i) = D0
        found(i,1) = .true.
      Z = ZZ
      print*,' el1 = ',el1,' ifl = ',ifl,' iw = ', iw,
     :' component of P vector = ',i
      IW = IW+1
  399 IF(IW .LE. NCOM) GO TO 314
      go to 401
  400 IF (IFL .EQ. IUW(2)) GO TO 401
      IFL = IUW(2)
      GO TO 314
*
* ---  check that all functions were found
*
  401 close(iuw(1))
      close(iuw(2))
      check = .true.
      do 316 i = 1,ncom
	 check = check .and. found(i,1)
 316  continue
      if (.not. check) then
	 write(iscw,*) ' Not all radial functions were found'
      write(iscw,'(1X,A3,2X,L1)')(el(i),found(i,1),i=1,ncom)
	 stop
      end if
*
*     .. initialize R and Rydberg constant
      call initm2(1)
      return
      end
