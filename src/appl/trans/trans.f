* Last modifications : September 14, 1997  
* tests to do : E2/M1
*               lhs = rhs 
**************************************************************************
*
*  This program evaluates the gf values and transition probabilities 
*  of electric and magnetic transitions in the LS coupling scheme
*  (length and velocity forms for E1/E2) or LSJ intermediate coupling 
*  scheme (length form only).  
*
*  M. Godefroid       Laboratoire de Chimie Physique Moleculaire CP160/09
*                     Universite Libre de Bruxelles, Belgium
*  A. Hibbert         Department of Applied Mathematics
*                     Queen's University, Belfast, Northern Ireland
*  C. Froese Fischer  Department of Computer Science
*                     Vanderbilt University, Nashville, U.S.A.
*  G. Gaigalas        Institute of Theoretical Physiscs
*                     and Astronomy, Vilnius, LITHUANIA
*
*  References: 
*  ----------
*  1. Hibbert et al,                Comput. Phys. Commun. 51(1988)285
*  2. Froese Fischer et al,         Comput. Phys. Commun. 64(1991)486-500
*  3. Froese Fischer and Godefroid, Comput. Phys. Commun. 64(1991)501-519
*
**************************************************************************
*
      PROGRAM TRANS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      EXTERNAL INITT
      PARAMETER (NWD=80,NOD=220)
      LOGICAL LORTH,LESS,rel,vok,IPRINT,ANGLES
      CHARACTER*1 IM,YES,PP,SYMBOL
      CHARACTER*24 CFILE(2),WFILE(2),JFILE(2),NAME(2),LFILE(2),tfile(2)
      character*64 configi,configf,ttfile,trfile,GFILE
      integer idx(2)
      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON /INFORM/ IREAD2,IWRITE2,IOUT,IREADF,ISC(7),ISCW2
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2)
      COMMON /INOUT2/IANG
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      POINTER (QIORTH,IORTH(1))
      COMMON /OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     1 ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     2 QIORTH
      COMMON /NTRM/NTERMS
      COMMON /NOR/NCOM,NCLOSI,NCLOSF,NORBI,NORBF,IWAR
      COMMON /RED/INDL(20),INDR(20),IORST(20),NCI1,KPL,KPR,NC,LIS(16),
     : JIST,JFST,NCONTR
      COMMON /STATE/QET1,QLBL1,QWT1,QJV1,QET2,QLBL2,QWT2,QJV2,NJV(2),
     : NVC(2),LGTH(2),NCF(2),QCFG1,QCFG2
      POINTER(QET1,ET1(1)),(QLBL1,LBL1(1)),(QWT1,WT1(1)),(QJV1,JV1(1)),
     :       (QET2,ET2(1)),(QLBL2,LBL2(1)),(QWT2,WT2(1)),(QJV2,JV2(1)),
     :       (QCFG1,CFG1(8,1)),(QCFG2,CFG2(8,1))
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,LCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
    9 FORMAT(///20X,'================='/
     :          20X,'  T R A N S   97  '/
     :          20X,'================='//)
    8 FORMAT(1H1,/,' ELECTRIC TRANSITION OF ORDER ',I3,/)
    7 FORMAT(1H1,/,' MAGNETIC TRANSITION OF ORDER ',I3,/)
    6 FORMAT(A1,I1)
CDBG88 FORMAT(////
CDBG :         ' IBUG1  =',I3,' (DEBUG FOR COEFF. OF FP IN CFP PROG.)'/
CDBG :         ' IBUG3  =',I3,' (DEBUG IN RECOUPLING PACKAGE)'/
CDBG :         ' NBUG6  =',I3,' (DEBUG IN TENSOR PACKAGE)'//)
  200 FORMAT(1H ,41X,2H_ /19H < INITIAL STATE ||,A2,1H(,I1,75H)|| FINAL
     1STATE > = >  COEFF * WEIGHT(INITIAL,I) * WEIGHT(FINAL,J) * < NL||,
     2A2,1H(,I1,8H)||N'L'>)
    4 FORMAT(1H+,41X,2H_ /41X,3HI,J///5X,'COEFF      I    J     < NL||'
     1,A2,1H(,I1,8H)||N'L'>//)
    5 FORMAT('  Transition between files: '/2x,a24/2x,a24)
*
*     Set machine dependent unit numbers:
*        IREAD  - The unit number of standard input
*        IWRITE - The unit number of standard output
*        ISCW   - The unit number of standard error (screen)
*        iuc(1) - The unit number for the <name>.c file for initial state
*        iuc(2) - The unit number for the <name>.c file for final   state
*        iuw(1) - The unit number for the <name>.w file for initial state
*        iuw(2) - The unit number for the <name>.w file for final   state
*        iul(1) - The unit number for the <name>.l file for initial state
*        iul(2) - The unit number for the <name>.l file for final   state
*        iuj(1) - The unit number for the <name>.j file for initial state
*        iuj(2) - The unit number for the <name>.j file for final   state
*        iuls   - The unit number for   'iul(1).iul(2).ls'    file (results)
*        iulsj  - The unit number for   'iuj(1).iuj(2).lsj'   file (results)
*        IANG
*
      iscw   =  0
      iread  =  5
      iwrite =  6
      iuc(1) =  1
      iuc(2) =  2
      iuw(1) =  3
      iuw(2) =  4
      iul(1) =  8
      iul(2) =  9
      iuj(1) = 10
      iuj(2) = 11
      iuls   = 13
      iulsj  = 14
      IANG   = 15
      call initr
      call factrl(32)
      WRITE(IWRITE,9)
      rel = .false.
      LORTH = .TRUE.
      IBUG1 = 1
      IBUG2 = 1
      IBUG3 = 1
      NBUG6 = 1
CDBG  WRITE(ISCW,*)  ' Input IBUG1, IBUG3, NBUG6 (0/1) '
CDBG  READ (IREAD,*)  IBUG1,IBUG3,NBUG6
CDBG  WRITE(IWRITE,88) IBUG1,IBUG3,NBUG6
  1   CONTINUE
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
        tfile(i)(1:j-1) = name(i)(1:j-1)
    2 continue
      ttfile = tfile(1)(1:idx(1)-1)//'.'//tfile(2)(1:idx(2)-1)
      write(iscw,*) ' intermediate printing (y or n) ?  '
      read(iread,'(a1)') pp
      if ( pp .eq. 'y' .or. pp .eq. 'Y' ) then
         IPRINT = .true.
         write(iscw,*) '  tolerance for printing ?  '
         read(iread,*) tol
      else
         IPRINT = .false.
      endif
      write(iscw,*) ' transitions only for E(initial) < ',
     :     'E(final) (y or n) ? '
      Read(iread,'(a1)') pp
      less = .false.
      if ( pp .eq. 'y' .or. pp .eq. 'Y' ) less = .true.
      OPEN(UNIT=iuw(1),FILE=wfile(1),STATUS='OLD',FORM='UNFORMATTED',
     :     action='READ')
      OPEN(UNIT=iuw(2),FILE=wfile(2),STATUS='OLD',FORM='UNFORMATTED',
     :     action='READ')
      CALL CFGIN2(MCFG,KCFG,LORTH,CFILE)
      lcfg = ncfg
      ncf(1) = mcfg
      ncf(2) = kcfg
*
*     allocate memory for orthogonality
      north = norbi*norbf
      if (north.gt.0) then
        if(ibug1.ne.0) print*,' qiorth  allocation: north = ',north
        call alloc(qiorth,north,4)
      CALL ORTH
      end if
*
* --- read the radial wavefunctions on units iuw(1) and iuw(2)
*     and sort them according the IAJCMP order
      CALL READW2
*
* --- specify the type of calculation
      write(iscw,*) ' Relativistic calculation ? (y/n) '
      READ(iread,'(A1)') YES
      if (YES.EQ.'Y'.OR.YES.EQ.'y') rel = .true.
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
        write(iulsj,5) name(1),name(2)
      else
        write(iscw,*) ' Input from an MCHF (M) or CI (C) calculation ? '
        read(iread,'(A1)') YES
        if (yes .eq. 'C' .or. yes .eq. 'c') then
           OPEN(UNIT=iul(1),FILE=lfile(1),STATUS='OLD')
           OPEN(UNIT=iul(2),FILE=lfile(2),STATUS='OLD')
        else
           ici = 0
        end if
        trfile = ttfile(1:j-1)//'.ls'
        OPEN(UNIT=iuls,FILE=trfile,STATUS='UNKNOWN')
        write(iuls,5) name(1),name(2)
      end if
*
* --- determine the number of different J-values (njv)
*               the total number of eigenvectors (nvc)
*               the total length of the wt vector to allocate (lgth)
* --- read in all the eigenvectors of the initial state
* --- read in all the eigenvectors of the final state
      CALL EIGVEC(ICI,CONFIGI,CONFIGF)
      WRITE(ISCW,*) ' Type of transition ? (E1, E2, M1, M2, .. or *) '
      READ (IREAD,6) IM, LAM
      IF (IM .EQ. '*' .OR. IM .EQ. ' ' ) STOP ' END OF CASE'
      CALL VALNUM(SYMBOL,LAM)
      IF (IM .EQ. 'e') IM = 'E'
      IF (IM .EQ. 'm') IM = 'M'
      gfile = ttfile(1:j-1)//'.'//IM(1:1)//SYMBOL(1:1)
      IF (IM .EQ. 'E') THEN
         IFL = 1
         WRITE(IWRITE,8) LAM
      ELSE
         IFL = 2
         WRITE(IWRITE,7) LAM
      ENDIF
      vok = .false.
      if(.not.rel.and.ifl.eq.1.and.lam.le.2) vok = .true.
*
      write(iscw,*) ' Use existing file for Angles ? (y/n) '
      READ(iread,'(A1)') YES
      if (YES.EQ.'Y'.OR.YES.EQ.'y') THEN
      ANGLES = .true.
        OPEN(UNIT=IANG,FILE=gfile,STATUS='OLD',form='unformatted' )
      ELSE
        ANGLES = .false.
        OPEN(UNIT=IANG,FILE=gfile,form='unformatted' )
      ENDIF
*
* --- allocate the /MULT/
* --- determine the number of (J,J') pairs satisfying selection rules
      CALL ALMULT(NPAIR)
*
* --- calculate the overlap and radial one-electron transition integrals
*     velocity form only for E1/E2 using the non-relativistic option.
      call radint(im)
*
      WRITE(IWRITE,200) IEM(IFL),LAM,IEM(IFL),LAM
      WRITE(IWRITE,4) IEM(IFL),LAM
      write(*,*) ' npair = ',npair
      NTERMS = 0
      IF(IFL .EQ. 2) IFL = 3
      IF(ANGLES) THEN
        CALL TAKE(NPAIR,GFILE)
      ELSE
        CALL CALCUL(NPAIR)
      ENDIF
      CLOSE(IANG)
        CALL PROBAB(ICI,IULS,IULSJ,NPAIR,CONFIGI,CONFIGF,LESS,IPRINT,IM)
      STOP
      END
