* Last update: September, 1997 by G. Gaigalas.
*
*  Things to do: - clean the old nonorthogonalities.
*                - check a Breit-Pauli case (I never did!).
*                - test inactive shells.
*                - add a threshold for contributions to the line strength,
*                  after counter-transformation.
*  BE AWARE: you should use fractional parentage modules to get consistency
*            with the Gedeminas-adapted angular codes !!!
*
**************************************************************************
*
*  This program evaluates the gf values and transition probabilities 
*  of electric and magnetic transitions in the LS coupling scheme
*  (length and velocity forms for E1/E2) or LSJ Breit-Pauli intermediate
*  coupling scheme (length form only).  
*
*  M. R. Godefroid    Laboratoire de Chimie Physique Moleculaire
*                     Universite Libre de Bruxelles, Belgium
*  A. Hibbert         Department of Applied Mathematics
*                     Queen's University, Belfast, Northern Ireland
*  C. Froese Fischer  Department of Computer Science
*                     Vanderbilt University, Nashville, U.S.A.
*  P. Jonsson         Department of Physics,
*                     Lund Institute of Technology, Lund, Sweden
*  J. Olsen           Department of Theoretical Chemistry,
*                     Chemical Center, Lund, Sweden
*
*
*  References: 
*  ----------
*  1. A. Hibbert et al, Comput. Phys. Commun. 51(1988)285
*  2. C. Froese Fischer et al, Comput. Phys. Commun. 64(1991)486-500
*  3. C. Froese Fischer and M. Godefroid, Comput. Phys. Commun. 64(1991)
*     501-519
*  4. P.A. Malmqvist, Int.J. of Quantum Chemistry, XXX, 479-94 (1986)
*  5. J. Olsen, M.R. Godefroid, P. Jonsson, P.A. Malmqvist and 
*     C. Froese Fischer, Phys. Rev. A, submitted. (#AY5043)
*
**************************************************************************
*
      PROGRAM TRANS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      EXTERNAL INITT
      PARAMETER (NWD=128,NOD=220)
      LOGICAL rel,vok,Iprint
      CHARACTER*1 IM,YES,PP, m_or_c
      CHARACTER*24 CFILE(2),WFILE(2),JFILE(2),NAME(2),
     : LFILE(2),LINTF(2)
      character*24 tfile(2)*24,ttfile*64,trfile*64
      CHARACTER el*3,ligne*80
      character line*70,configi*64,configf*64
      character elc(8)*3,couple(15)*3
      character*3 elras,elrasi,elrasf
      integer q(8),idx(2),lsj(2)
      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON /DBG  /IBUGM
      COMMON /INFORM/ IREAD2,IWRITE2,IOUT,IREADF,ISC(7),ISCW2
      COMMON/INOUT/IREAD,IWRITE,ISCW,
     :         iuc(2),iuw(2),iul(2),iuj(2),iut(2)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON /NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7),el(nwd)
      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      POINTER (QIORTH,IORTH(1))
      COMMON /OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     1 ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     2 QIORTH
      COMMON /OVRINT/IOVEL(2),IOVER(2),IOVEP(2)
      COMMON /NTRM/NTERMS
      COMMON /NOR/NCOM,NCLOSI,NCLOSF,NORBI,NORBF,IWAR
      COMMON /RED/INDL(20),INDR(20),IORST(20),NCI1,KPL,KPR,NC,LIS(16),
     : JIST,JFST,NCONTR
      COMMON /STATE/QET1,QLBL1,QWT1,QJV1,QET2,QLBL2,QWT2,QJV2,NJV(2),
     : NVC(2),LGTH(2),NCF(2),QCFG1,QCFG2
      POINTER(QET1,ET1(1)),(QLBL1,LBL1(1)),(QWT1,WT1(1)),
     :       (QJV1,JV1(1)),
     :       (QET2,ET2(1)),(QLBL2,LBL2(1)),(QWT2,WT2(1)),
     :       (QJV2,JV2(1)),
     :       (QCFG1,CFG1(8,1)),(QCFG2,CFG2(8,1))
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,LCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /FO/IWF(2),NCLOS(2),NOPEN(2),MAXNFO
      COMMON /RAS/ ninac(11,2),nras1(11,10,2),nras2(11,2),
     :             nras3(11,10,2),NGAS1(2),NGAS3(2),
     :             nl(11,2),elras(2*nwd),elrasi(nwd),elrasf(nwd),
     :             itab(nwd,2),lmax(2)
    4 FORMAT(1H+,41X,2H_ /41X,3HI,J///5X,'COEFF      I    J     < NL||'
     :,A2,1H(,I1,8H)||N'L'>//)
    5 FORMAT(///20X,'========================'/
     :          20X,'  T R A N S B I O  99 '/
     :          20X,'========================'//)
    6 FORMAT(A1,I1)
    7 FORMAT(1H1,/,' ELECTRIC TRANSITION OF ORDER ',I3,/)
    8 FORMAT(1H1,/,' MAGNETIC TRANSITION OF ORDER ',I3,/)
    9 FORMAT('  Transition between files: '/2X,A24/2X,A24)
   10 FORMAT(1H ,41X,2H_ /19H < INITIAL STATE ||,A2,1H(,I1,75H)|| FINAL
     :STATE > = >  COEFF * WEIGHT(INITIAL,I) * WEIGHT(FINAL,J) * < NL||,
     :A2,1H(,I1,8H)||N'L'>)
CDBG88 FORMAT(////
CDBG :         ' IBUG1  =',I3,' (DEBUG FOR COEFF. OF FP IN CFP PROG.)'/
CDBG :         ' IBUG3  =',I3,' (DEBUG IN RECOUPLING PACKAGE)'/
CDBG :         ' NBUG6  =',I3,' (DEBUG IN TENSOR PACKAGE)'//)
*
*     Set machine dependent unit numbers:
*  IREAD  - The unit number of standard input
*  IWRITE - The unit number of standard output
*  ISCW   - The unit number of standard error (screen)
*  iuc(1) - The unit number for the <name>.c file for initial state
*  iuc(2) - The unit number for the <name>.c file for final   state
*  iuw(1) - The unit number for the <name>.w file for initial state
*  iuw(2) - The unit number for the <name>.w file for final   state
*  iul(1) - The unit number for the <name>.l file for initial state
*  iul(2) - The unit number for the <name>.l file for final   state
*  iuj(1) - The unit number for the <name>.j file for initial state
*  iuj(2) - The unit number for the <name>.j file for final   state
*  iuls   - The unit number for   'iul(1).iul(2).ls'    file (results)
*  iulsj  - The unit number for   'iuj(1).iuj(2).lsj'   file (results)
*
      iscw   = 0
      iscw2 = iscw
      iread  = 5
      iwrite = 6
      iwrite2 = iwrite
      iuc(1) = 1
      iuc(2) = 2
      iuw(1) = 3
      iuw(2) = 4
      iul(1) = 8
      iul(2) = 9
      iuj(1) = 10
      iuj(2) = 11
      iuls   = 13
      iulsj  = 14
C     iut(1) = 15
C     iut(2) = 16
*
      call initr
      call factrl(32)
      WRITE(IWRITE,5)
      rel = .false.
      IBUG1 = 1
*     IBUG1 = 1  ! Allows for output
      IBUG2 = 1
      IBUG3 = 1
      IBUGM = 1
      NBUG6 = 1
*
* --- determine debug information
CDBG  WRITE(ISCW,*)  ' Input IBUG1, IBUG3, NBUG6 (0/1) '
CDBG  READ (IREAD,*)  IBUG1,IBUG3,NBUG6
CDBG  WRITE(IWRITE,8) IBUG1,IBUG3,NBUG6
*
* --- get the configuration data for the states
    1 CONTINUE
CSUN  i = iargc()
CSUN  if (i .lt. 2) then
         WRITE(ISCW,*)  ' Name of Initial State'
         READ(IREAD,'(A)') name(1)
         WRITE(ISCW,*)  ' Name of Final State'
         READ(IREAD,'(A)') name(2)
CSUN  else
CSUN     call getarg(1,name(1))
CSUN     call getarg(2,name(2))
CSUN  end if
*
      do 2 i = 1,2
        j = index(name(i),' ')
        idx(i) = j
        if (j .eq. 1) then
          WRITE(ISCW,*) ' Names may not start with blanks'
          go to 1
        end if
        cfile(i) = name(i)(1:j-1)//'.c'
        wfile(i) = name(i)(1:j-1)//'.w'
        lfile(i) = name(i)(1:j-1)//'.l'
        jfile(i) = name(i)(1:j-1)//'.j'
        lintf(i) = name(i)(1:j-1)//'.t'
        tfile(i)(1:j-1) = name(i)(1:j-1)
    2 continue
      ttfile = tfile(1)(1:idx(1)-1)//'.'//tfile(2)(1:idx(2)-1)
      OPEN(UNIT=IUC(1),FILE=CFILE(1),STATUS='OLD')
      OPEN(UNIT=IUC(2),FILE=CFILE(2),STATUS='OLD')
Cdbg  OPEN(UNIT=iut(1),STATUS='scratch',FORM='FORMATTED')
Cdbg  OPEN(UNIT=iut(2),STATUS='scratch',FORM='FORMATTED')
*
* --- calculate the one-electron part for initial and final superpositions
      write(iscw,*) ' intermediate printing (y or n) ?  '
      read(iread,'(a1)') pp
      if ( pp .eq. 'y' .or. pp .eq. 'Y' ) then
         Iprint = .true.
         write(iscw,*) '  tolerance for printing ?  '
         read(iread,*) tol
      else
         Iprint = .false.
      endif

      write(iscw,*) ' Relativistic calculation ? (y/n) '
      READ(iread,'(A1)') YES
      if (YES.EQ.'Y'.OR.YES.EQ.'y') then
         rel = .true.
      else
         rel = .false.
      write(iscw,*) ' Input from an MCHF (M) or CI (C) calculation ? '
         read(iread,'(A1)') m_or_c
      end if

      WRITE(ISCW,*)
     :      ' Type of transition ? (E1, E2, M1, M2, .. or *) '
      READ (IREAD,6) IM, LAM

      write(iscw,*)
     :' One-electron coupling coefficients for initial state...'
      call nonh1(1)
      rewind (iuc(1))
      write(iscw,*)
     :' One-electron coupling coefficients for final state...'
      call nonh1(2)
      rewind (iuc(2))
*
      OPEN(UNIT=iuw(1),FILE=wfile(1),STATUS='OLD',FORM='UNFORMATTED')
      OPEN(UNIT=iuw(2),FILE=wfile(2),STATUS='OLD',FORM='UNFORMATTED')
*
      CALL CFGIN2(MCFG,KCFG)
*
      lcfg = ncfg
      ncf(1) = mcfg
      ncf(2) = kcfg
*
* --- specify the type of RAS calculation
      CALL RASIN
*
*     allocate memory for orthogonality
      if(norbi.ne.0.or.norbf.ne.0) then
         print*,' norbi or norbf is different from 0...'
         print*,' something wrong...'
      end if
      north = norbi*norbf
      print*,' north = ',north
      if (north.gt.0) then
      if(ibugm.ne.0) print*,' qiorth  allocation: north = ',north
         call alloc(qiorth,north,4)
         call orth
      end if
*
* --- allocate the /NEL/
* --- read the radial wavefunctions on units iuw(1) and iuw(2)
*     and sort them according the RAS order
      call readw2
      close(iuw(1))
      close(iuw(2))
*
* --- calculate the overlap matrix between initial and final
*     orbitals
      CALL BRKT
*
* --- open .l/.j and .ls/.lsj files according rel value
    3 j = index(ttfile,' ')
      if (j .eq. 1) then
        WRITE(ISCW,*) ' Names may not start with blanks'
        go to 3
      end if
*
      ici = 1
      if (rel) then
        OPEN(UNIT=iuj(1),FILE=jfile(1),STATUS='OLD')
        OPEN(UNIT=iuj(2),FILE=jfile(2),STATUS='OLD')
        trfile = ttfile(1:j-1)//'.lsj'
        OPEN(UNIT=iulsj,FILE=trfile,STATUS='UNKNOWN')
        write(iulsj,9) name(1),name(2)
      else
        if (m_or_c .eq. 'C' .or. m_or_c .eq. 'c') then
	   ici = 2
           OPEN(UNIT=iul(1),FILE=lfile(1),STATUS='OLD')
           OPEN(UNIT=iul(2),FILE=lfile(2),STATUS='OLD')
*          .. compute correct J values
	   nc = 1
	   do i =1,2
	   print *, 'j1qnrd', nc, j1qnrd(2*noccsh(nc)-1,nc)
	   jd = j1qnrd(2*noccsh(nc)-1,nc)
	   ja = mod(jd,64)
	   jd = jd/64
	   jb = mod(jd,64)
	   jc = jd/64
	   print *, 'seniority =', ja, '  2*L+1=',jb, ' 2S+1 =',jc
	   lsj(i) = jb+jc-2
	   nc = nc+ncf(1)
	   print *, 'nc,ncf', nc,ncf(1)
	   end do
        else
           ici = 0
        end if
        trfile = ttfile(1:j-1)//'.ls'
        OPEN(UNIT=iuls,FILE=trfile,STATUS='UNKNOWN')
        write(iuls,9) name(1),name(2)
      end if
      CALL EIGVEC(ICI,CONFIGI,CONFIGF,LSJ)
      print *, 'jv1,jv2', jv1(1),jv2(1)

      IF (IM .EQ. '*' .OR. IM .EQ. ' ' ) STOP ' END OF CASE'
      IF (IM .EQ. 'e') IM = 'E'
      IF (IM .EQ. 'm') IM = 'M'
      IF (IM .EQ. 'E') THEN
         IFL = 1
         WRITE(IWRITE,7) LAM
      ELSE
         IFL = 2
         WRITE(IWRITE,8) LAM
      ENDIF
      lam2 = lam + lam
      vok = .false.
*     .. let us compute velocity also for rel = .true.
*     if(.not.rel.and.ifl.eq.1.and.lam.le.2) vok = .true.
      if(ifl.eq.1.and.lam.le.2) vok = .true.
*
* --- allocate the /MULT/
      CALL ALMULT(NPAIR)
*
* --- Here,we go : call routine to perform tranformations 
*     of shells and CI coefficients so orbital bases become
*     biorthogonal
Cmrg  LBUF = 96*1027
      LBUF = 1
      MXL = MAX0(LMAX(1),LMAX(2))-1
*. Define TEST level for BIOTRN : 0 => silence, 1=> just a bit
Cmrg  NTESTB = 1
      NTESTB = 15
      IF(NTESTB.NE.0) THEN
        write(6,*) ' TRANS : MXL ', MXL
        write(6,*) ' I am going to call BIOTRA, wish me luck'
      END IF
      write(iscw,*) ' Call BIOTRN...'
      if(ibugm.ne.0) print*,' wt1(1) = ',wt1(1),' wt2(1) = ',wt2(1)
*
*. Here we are...
      CALL BIOTRN(P(1,1),NL(1,1),WT1(1),NCF(1),NVC(1),IDUMMY,
     &            P(1,IWF(1)+1),NL(1,2),WT2(1),NCF(2),NVC(2),IDUMMY,
     &            NOD,MXL,NINAC(1,1),
     &            LBUF,15,16,
     &            NTESTB)
*
*. Here we were...
      if(ibugm.ne.0) print*,' wt1(1) = ',wt1(1),' wt2(1) = ',wt2(1)
*
* --- calculate the rotated slopes at the origin
      call slope
      WRITE(6,'(A)') ' Test call of BRKT after BIOTRN'
      CALL BRKT
 1000 Continue
*
* --- calculate the overlap and radial one-electron transition integrals
*     velocity form only for E1/E2 using the non-relativistic option.
!      write(iscw,*) ' Calculation of radial integrals...'
      call RADINT(im)
*
* --- calculate the elements of vshell
      WRITE(IWRITE,10) IEM(IFL),LAM,IEM(IFL),LAM
      WRITE(IWRITE,4) IEM(IFL),LAM
      print*,' npair = ',npair
      NTERMS = 0
      IF(IFL .EQ. 2) IFL = 3
*     CALL CALCUL(NPAIR)
      CALL CALCUL(NPAIR,IPRINT,TOL)
      CALL PROBAB(ICI,IULS,IULSJ,NPAIR,CONFIGI,CONFIGF,IM,IPRINT)
      WRITE(ISCW,*)
     :      ' Type of transition ? (E1, E2, M1, M2, .. or *) '
      READ (IREAD,6) IM, LAM
      IF (IM .EQ. '*' .OR. IM .EQ. ' ' ) STOP ' END OF CASE'
      IF (IM .EQ. 'e') IM = 'E'
      IF (IM .EQ. 'm') IM = 'M'
      IF (IM .EQ. 'E') THEN
         IFL = 1
         WRITE(IWRITE,7) LAM
      ELSE
         IFL = 2
         WRITE(IWRITE,8) LAM
      ENDIF
      lam2 = lam + lam
      vok = .false.
*     .. let us compute velocity also for rel = .true.
*     if(.not.rel.and.ifl.eq.1.and.lam.le.2) vok = .true.
      if(ifl.eq.1.and.lam.le.2) vok = .true.
*
* --- allocate the /MULT/
      CALL ALMULT(NPAIR)
*
      Go to 1000
      END
