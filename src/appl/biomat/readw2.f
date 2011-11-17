*
*     ------------------------------------------------------------------
*       R E A D W 2
*     ------------------------------------------------------------------
*
      SUBROUTINE readw2
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NWD=128,NOD=220)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCF
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      COMMON /NOR/NCOM,NCLOSI,NCLOSF,NORBI,NORBF,IWAR
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2),iut(2)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON /RAS/ ninac(11,2),nras1(11,10,2),nras2(11,2),
     :             nras3(11,10,2),NGAS1(2),NGAS3(2),
     :             nl(11,2),elras(2*nwd),elrasi(nwd),elrasf(nwd),
     :             itab(nwd,2),lmax(2)
      COMMON /NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7),el(nwd)
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      COMMON /FO/IWF(2),NCLOS(2),NOPEN(2),MAXNFO
      COMMON/DFG/ IBUGM
      CHARACTER*6 atom, term
      LOGICAL found(NWD,2),check
      CHARACTER*3 el,el1,elras,elrasi,elrasf,eltens
      COMMON /ELT/ELTENS(NWD)
      character*7 label(2)
      DOUBLE PRECISION PT(NOD)
      DATA LABEL/'Initial','Final  '/
*
* --- allocate the /NEL/
*
      if(ibugm.ne.0) print*,' iqp     allocation: nod*maxnfo = ',
     :               nod*maxnfo
      call alloc(iqp,nod*maxnfo,8)
      if(ibugm.ne.0) print*,' iqn     allocation: maxnfo = ',maxnfo
      call alloc(iqn,maxnfo,4)
      if(ibugm.ne.0) print*,' iql     allocation: maxnfo = ',maxnfo
      call alloc(iql,maxnfo,4)
      if(ibugm.ne.0) print*,' iqaz    allocation: maxnfo = ',maxnfo
      call alloc(iqaz,maxnfo,8)
      if(ibugm.ne.0) print*,' iqmax   allocation: maxnfo = ',maxnfo
      call alloc(iqmax,maxnfo,4)
*
* --- read the radial wavefunctions on units iuw(1) and iuw(2)
*     and sort them according the RAS order
*
*
* --- transfer NDIMS/NON30 -> PARAM
*
      nwf = maxnfo
      ncfg = ncf
      print*,' in readw2, nwf = ',nwf,' ncfg = ',ncfg
*
      do 1 i = 1,nwf
        if(i.le.iwf(1)) then
           elrasi(i) = elras(i)
           found(i,1) = .false.
        else
           k = i-iwf(1)
           elrasf(k) = elras(i)
           found(k,2) = .false.
        end if
    1 continue
      print*,' elrasi = ',(elrasi(i),i=1,iwf(1))
      print*,' elrasf = ',(elrasf(i),i=1,iwf(2))
*
*  *****  READ THE RADIAL FUNCTIONS
*
      IW = 1
  14  IF (IW .LE. (iwf(1))) THEN
          IK = 1
          READ(iuw(1)) ATOM,Term,EL1,M,ZZ,ETI,EKI,AZD,(PT(J),J=1,M)
          call eptr(elrasi,el1,i,*999)
      ELSE
          IK = 2
          READ(iuw(2)) ATOM,Term,EL1,M,ZZ,ETI,EKI,AZD,(PT(J),J=1,M)
          call eptr(elrasf,el1,i,*999)
      END IF
      if (.not. found(i,ik)) then
        if (ik.eq.1) then
          ij = i
        elseif (ik.eq.2) then
          ij = i+iwf(1)
        end if
        print*,label(ik),' electron ',el1,' P index = ',ij
        Do 15 j= 1,m
            p(j,ij) = pt(j)
   15   continue
        n(ij) = ichar(el1(2:2))-ichar('1')+1
        l(ij) = lval(el1(3:3))
	az(ij) = azd
	max(ij) = m
        DO 24 J = M+1,NO
   24       P(J,ij) = D0
        found(i,ik) = .true.
      end if
      Z = ZZ
      IW = IW+1
  999 IF (IW.LE.NWF) GO TO 14
   13 close(iuw(1))
      close(iuw(2))
*
* ---  check that all functions were found
*
      check = .true.
      do 16 i = 1,iwf(1)
	 check = check .and. found(i,1)
  16  continue
      do 17 i = 1,iwf(2)
	 check = check .and. found(i,2)
  17  continue
      if (.not. check) then
	 write(iscw,*) ' Not all radial functions were found'
         write(iscw,'(1X,A3,2X,L1)')(elrasi(i),found(i,1),i=1,iwf(1))
         write(iscw,'(1X,A3,2X,L1)')(elrasf(i),found(i,2),i=1,iwf(2))
	 stop
      end if
      do 20 i = 1,maxorb
        do 21 j=1,iwf(1)
          if(eltens(i).eq.elrasi(j)) itab(i,1)=j
  21    continue
        do 22 j=1,iwf(2)
          if(eltens(i).eq.elrasf(j)) itab(i,2)=j+iwf(1)
  22    continue
  20  continue
       print*
      print*
      print*,' Position of biorthonormal shells in P vector : '
      print*,' ---------------------------------------------- '
      print*,label(1),' :  ',(itab(i,1),i=1,maxorb)
      print*,label(2),' :  ',(itab(i,2),i=1,maxorb)
      print*
*     .. initialize R and Rydberg constant
      call initm2(1)
      return
      end
