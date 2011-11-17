* ======================================================================
*
*     ISOTOPE SHIFT
*      
*     Written by -
*
*     G. GAIGALAS, INSTITUTE OF THEORETICAL PHYSICS
*        AND ASTRONOMY, VILNIUS, LITHUANIA
*
*     C. FROESE FISCHER, VANDERBILT UNIVERSITY
*        NASHVILLE, TN 37135 USA
*
*     FEBRUARY, 1994
*
* ======================================================================
*
      PROGRAM GISO
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWCD=20,LSDIM=30000, NWD=60, LCDIM=65536)
*
Cyy
      CHARACTER*1 ans
      CHARACTER*3 eliso1,eliso2,ellab
      COMMON/GRADINT/ggrad(nwd,nwd)
      POINTER(QWT,WT(1))
      COMMON/NAN/QWT
      COMMON/ISO/GIS,RIS,COR,CO,FIELCOR,FIEL
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/ELLABEL/ELLAB(NWD)
      COMMON/ANGCORE/CBANG(NWD,NWD)
      COMMON/INOUT/IN,IERR,IPRI,IUF
Cyy
      POINTER (qcn,cn(LSDIM,7)), (qpackn,ipackn(lsdim,7)),
     :        (qnijptr,nijptr(lsdim,7)),(qjan,jan(lsdim)),
     :        (qjbn,jbn(lsdim)) 
      COMMON/buffer/qcn, qpackn, qnijptr, qjan, qjbn
      COMMON/CLOSED/B1ELC(4),NCLOSD,IBK
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/DIMEN/KFL1,KFL2,KFL3,KFL4,KFL5,KFL6,KFL7,MXIHSH
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,
     : IALL,JSC(3),ISCW
      COMMON /DIAGNL/IDIAG,JA,JB
      POINTER  (qjptr, jptr(1))
      COMMON/fout/lcount(7),nrec(7),iflag,lij,nij,qjptr
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON/NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON/NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
CG
      COMMON /OPERAT/ICOLOM,ISOTOP
CG
      REAL TIME(2)
Cyy      Character*72 string
      Character*72 string,name(5)
*
      CHARACTER*16 INPUT
      EXTERNAL INITT
*
    1 FORMAT(//' IOUT =  FGR.LST (OUTPUT FILE)'/
     :         ' IBUG1  =',I3,' (DEBUG IN WEIGHTS OF 1-EL PART)'/
     :         ' IBUG2  =',I3,' (DEBUG IN WEIGHTS OF 2-EL PART)'/
     :         ' IBUG3  =',I3,' (DEBUG IN RECOUPLING PACKAGE)'//)
*
    2 FORMAT(/////10X,'   ==============================='/
     :            10X,'            G I S O    94 ',/
     :            10X,'   ==============================='////)
*
          i = iargc()
          if (i .eq. 0) then
             INPUT = 'clist'
          else
             call getarg(1,INPUT)
          end if
*
*
*     ...  THE FOLLOWING SECTION CONCERNS INPUT/OUTPUT AND MAY BE
*          SYSTEM DEPENDENT.  CHECK ALLOWED UNIT NUMBERS AND
*          FILE NAME CONVENTIONS - MODIFY, IF NECESSARY.
*
      IREAD=4
      IUF=14
      IN=5
      IERR=0
      IPRI=65
CTSS      ISCW = 5
          ISCW = 0
CHP    ISCW = 7
      IWRITE = 6
      IOUT=8
      WRITE(IWRITE,2)
      RIS=ZERO
      GIS=ZERO
*
Cyy
7     WRITE(IERR,*) 'Name of state'
      READ(IN,'(A)') NAME(1)
      K=INDEX(NAME(1),' ')
      IF (K.EQ.1) THEN
         WRITE(IERR,*) 'Names may not start with a blank'
         GO TO 7
      ELSE
         NAME(1)(K:K+1)='.c'
         NAME(2)=NAME(1)(1:K-1)//'.l'
         NAME(3)=NAME(1)(1:K-1)//'.g'
         NAME(4)=NAME(1)(1:K-1)//'.i'
         NAME(5)=NAME(1)(1:K-1)//'.w'
      ENDIF
30    WRITE(IERR,*)'Input from an MCHF (M) or CI (C) calculation ?'
      READ(IN,'(A)') ANS
      IF (ANS.EQ.'M'.OR.ANS.EQ.'m') THEN
         IMCHF=1
         NR=1
      ELSEIF (ANS.EQ.'C'.OR.ANS.EQ.'c') THEN
         IMCHF=2  
         OPEN(UNIT=63, FILE=NAME(2),STATUS='OLD')
         WRITE(0,*) 'Give the number of the cfg. in the CI expansion,'
         WRITE(0,*) 'for which the mass shift is to be calculated '
         READ(5,*) NR
      ELSE
         GO TO 30
      ENDIF
      OPEN(UNIT=IREAD,FILE=NAME(1),STATUS='OLD')
      OPEN(UNIT=64,FILE=NAME(3),STATUS='OLD')
      OPEN(UNIT=65,FILE=NAME(4),STATUS='UNKNOWN')
      OPEN(UNIT=14, FILE=NAME(5),STATUS='OLD',
     :     FORM='UNFORMATTED')
Cyy
Cyy      OPEN(UNIT=IREAD,FILE=INPUT,STATUS='UNKNOWN')
      OPEN(UNIT=IOUT, FILE='yint.lst',STATUS='UNKNOWN',
     :     FORM='unformatted')
CG
      ICOLOM=0
      ISOTOP=1
CG
*
*     ... END OF MACHINE DEPENDENT SECTION
*
      IBUG1 = 0
      IBUG2 = 0
      IBUG3 = 0
*     WRITE(ISCW, '(A)') ' IBUG1,IBUG2,IBUG3 ?  FORMAT(3(I1,1X)) '
*     READ(5,'(I1,1X,I1,1X,I1)') IBUG1,IBUG2,IBUG3
*     WRITE(IWRITE,1) IBUG1,IBUG2,IBUG3
*     WRITE(ISCW, '(A)') ' FULL PRINT-OUT ? (Y/N) '
      IFULL = 0
*     READ(5,'(A2)') ANS
*     IF ( ANS .EQ. 'Y' .OR. ANS .EQ. 'y') IFULL = 1
*
      NEW = 0
      NZER0 = 0
*
*  ---  Determine input data; non-orthogonal case
*
      CALL CFGN1
      call alloc(qjptr, ncfg,4)
*
*     WRITE(ISCW, '(/A)') ' ALL INTERACTIONS ? (Y/N) '
*     READ(5,'(A1)') ANS
*     IF (ANS .EQ. 'N' .OR. ANS .EQ. 'n' ) THEN
Cyy          WRITE(ISCW, *) ' NEW(ALL=0), N-ZERO(ALL=0)'
Cyy          READ(5,*) NEW, NZERO
*     END IF
      IF ( NEW .EQ. 0 ) NEW = NCFG
      IF (NZERO .EQ. 0) THEN
         NZERO = NCFG
         IFIRST = 0
      ELSE
         IFIRST = 1
      END IF
*
      IALL = 1
*     WRITE(ISCW, '(A)') ' MCHF or BREIT format ? (0 or 1) '
*     READ (5,'(I1)') IALL
*
*
Cyy
      write(*,*) 'Fore isowf'
      CALL ISOWF
      write(*,*) 'efter isowf'
* --- Read configuration weights

      REWIND(UNIT=IREAD)
      CALL READWT(MAXORB,IMCHF,NCFG,NR)
Cyy      do 933 k = 1,ncfg
Cyy         write(*,*) 'j,wt',k,wt(k)
Cyy933   continue

      MAX=MAXORB+NCLOSD
      read(64,193) etotal,number
      if (max.ne.number) then
         write(iscw,*) ' Number of radial orbitals differ'
         stop
      endif
      do 5898 i = 1,max
         write(IERR,'(a3)') ellab(i)
5898  continue
      do 5899 i = 1,max
         do 5900 j = 1,max
            read(64,194) eliso1,eliso2,grad
            CBANG(I,J)=ZERO
            n1 = 0
            n2 = 0
            do 5901 k = 1,max 
               if (eliso1.eq.ellab(k)) n1 = k
               if (eliso2.eq.ellab(k)) n2 = k
5901        continue
            if (n1.eq.0.or.n2.eq.0) then
               write(iscw,*) ' Radial orbitals does not match'
               stop
            endif
            ggrad(n1,n2) = grad
            write(IERR,194) ellab(n1),ellab(n2),ggrad(n1,n2)
5900     continue
5899  continue
193   format(f16.9,3x,i3)
194   format(a3,a3,1x,d25.18)
Cyy
*     SET FACTORIALS AND LOG OF FACTORIALS
*
      CALL INITA
      CALL INITISO
*
*     Initialize parameters for output
*
      do 5 i = 1,7
	lcount(i) = 0
	nrec(i) = 0
  5   continue
      lij = 0
      nij = 0
      ncol = 0
*     .. write out intial data about the problem
      write(iout) nclosd, maxorb, ncfg, lsdim, lcdim
      string =' '
      if (nclosd .gt. 0) write(string,'(24A3)') (iajcld(i),i=1,nclosd)
      write(iout) string
      do 10 i = 1, maxorb, 24
	m = min(maxorb,i+23)
	write(string,'(24A3)') (iajcmp(j),j=i,m)
	write(iout) string
  10  continue
*     .. allocate memory for buffered i/o
      call alloc(qcn,7*lsdim,8)
      call alloc(qpackn,lsdim*7,4)
      call alloc(qnijptr,lsdim*7,4)
      call alloc(qjan,lsdim,4)
      call alloc(qjbn,lsdim,4)
      do 200 jb = 1,ncfg
	CALL SHELLSJB(JB)
        CALL ANGMOMG(NEW,NZERO,IFIRST)
	ncol = ncol + 1
	jptr(ncol) = nij
  200 Continue
      CALL COREISO
      CALL COREOUT
      CALL PRISO
      CALL SPEC
      CALL PRFIELD
*     .. finish writing the jan, jbn information
      if (lij .eq. lsdim ) then
*       .. we just happen to have a full array
	write(0,*) ' Happened to have a full array'
        write(iout) lsdim, (jan(i),i=1,lsdim), (jbn(I),i=1,lsdim)
	lij = 1
	jan(lij) = 0
	jbn(lij) = 0
      end if
      write(iout) lij,(jan(i),i=1,lij), (jbn(i),i=1,lij)
      write(iout) ncol, (jptr(i),i=1,ncol)
*     .. clear out other arrays so memory can be used for sorting
      do 205 i = 1,7
	iou = 9+i
	n = lcount(i)
	if (n .ne. 0) then
	   if (nrec(i) .eq. 0) 
     :        OPEN(unit=iou,form='unformatted',status='scratch')
           write(iou) (cn(j,i),j=1,n),
     :                (ipackn(j,i),j=1,n),
     :                (nijptr(j,i),j=1,n)
	end if
	rewind(iou)
  205 continue
*     .. deallocate memory for buffered i/o
      nwf = maxorb
      if (nwf .gt. 1) call dalloc(qiorth,(nwf*(nwf-1))/2)
      call dalloc(qnoc,ncfg)
      call dalloc(qnelcsh,8*ncfg)
      call dalloc(qnocorb,8*ncfg)
      call dalloc(qj1,15*ncfg)
      call dalloc(qnjcomp,nwf)
      call dalloc(qljcomp,nwf)
      call dalloc(qljclsd,nwcd)
      call dalloc(qcn,7*lsdim)
      call dalloc(qpackn,lsdim*7)
      call dalloc(qnijptr,lsdim*7)
      call dalloc(qjan,lsdim)
      call dalloc(qjbn,lsdim)
      CLOSE(UNIT=IOUT)
      Nf = nrec(1)*lsdim + lcount(1)
      ng = nrec(2)*lsdim + lcount(2)
      nl = nrec(7)*lsdim + lcount(7)
      nr = (nrec(3)+nrec(4)+nrec(5)+nrec(6))*lsdim + 
     :      lcount(3)+lcount(4)+lcount(5)+lcount(6)
      ITOTAL = NF+NG+NR+NL
Cmrg  call wcfg(itotal,ntot,nij)
      call etime(time)
6     write(iscw,'(//A/A//A,F8.3,A//)') ' END OF CASE',' ===========',
     :      ' Total CPU time was ', TIME(1)/60, ' minutes'
Cyy
Cyy
      END
*
*     -----------------------------------------------------------------
*      I S O W F
*     -----------------------------------------------------------------
*
*     READ THE STARTING PARAMETERS AZ FROM THE   <NAME>.W   FILE.
*     AND SORT THEM ACCORDING TO THE ORDER IN ELLABEL
*
      SUBROUTINE ISOWF
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=60,NOD=220)
      CHARACTER EL1*3,ELLAB*3,ATOM*6,TERM*6
      DIMENSION P(NOD)
      COMMON/ELLABEL/ELLAB(NWD)
      COMMON/PAPIL/NCI,NWI,MAXIM
      COMMON/INOUT/IN,IERR,IPRI,IUF
      COMMON/PRA/Z
      POINTER(IQAZ,AZ(1))
      COMMON/NEL1/IQAZ
         write(*,*) 'In isowf,nwi,maxim',nwi,maxim
      CALL ALLOC(IQAZ,MAXIM,8)
      I=1
    1 READ(IUF,END=2) ATOM,TERM,EL1,M,ZZ,ETI,EKI,AZ(I),(P(J),J=1,M)
      DO 93 J=1,MAXIM
         K=0
         IF(EL1.EQ.ELLAB(J)) THEN
           AZ(I)=AZ(J)
           K=K+1
           I=I+1
         ENDIF
         IF(K.EQ.2) THEN
           WRITE(*,*) 'THERE ARE TWO IDENTICAL ORBITALS IN THE WFN'
           STOP
         ENDIF
93    CONTINUE
      Z=ZZ
      IF(I.LE.MAXIM) GO TO 1
    2 CLOSE(IUF)
      RETURN
      END
*
*     ------------------------------------------------------------
*      I N I T I S O
*     ------------------------------------------------------------
*
      SUBROUTINE INITISO
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=60,NOD=220)
      PARAMETER (NELEMENT=50)
      DIMENSION FZ(NELEMENT-9),MN(NELEMENT)
*
      COMMON/GISOCON/AME,RY,ISP,ISF
      COMMON/CLOSED/B1ELC(4),NCLOSD,IBK
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/INOUT/IN,IERR,IPRI,IUF
      COMMON/PRA/Z
      COMMON/RYDBERG/RYDBRG,ZMU,ENERGY
      COMMON/COREC/ZM1,ZM2,FZE,FZA,DR2
      CHARACTER YES
      DATA MN/1,4,7,9,11,12,14,16,
     : 19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,58,
     : 63,64,69,74,75,80,79,84,85,88,89,90,93,98,99,102,103,108,
     : 107,114,115,120/
      DATA FZ/163.5,198.4,239.6,283.9,332.7,386.1,444.2,507.8,
     : 576.4,650.7,731.2,817.9,911.2,1011.4,1119.3,1234.9,1359.3,
     : 1491.9,1634.4,1786.6,1948.7,2122.7,2308.2,2507.4,2718.5,
     : 2944.3,3187.5,3442.2,3720.3,4015.4,4326.6,4664.8,5022.4,
     : 5405.9,5811.9,6245.0,6706.8,7201.3,7728.7,8296.5,8889.3/
*
* ****   DETERMINE THE RYDBERG CONSTANT AND F(Z) PARAMETER 
*
      AME = 548.579903D-6
      RY = 109737.31534
      WRITE(IERR,'(/1X,A)') 'Default Rydberg constant (y/n)'
      READ(IN,'(A1)') YES
      IF(YES.EQ.'Y'.OR.YES.EQ.'y'.AND.(INT(Z).LE.NELEMENT)) THEN
        ZMU=DBLE(MN(INT(Z)))
      ELSE
        WRITE(IERR,'(1X,A)') 'Enter the mass of the atom'
        READ(IN,*) ZMU
      ENDIF
      RYDBRG=RY/(ONE+AME/ZMU)
        IF(INT(Z).LT.10) THEN
          FZE=1.566432
        ELSEIF(INT(Z).LE.NELEMENT) THEN
          FZE=FZ(INT(Z-9))
        ELSE
          WRITE(IERR,'(1X,A)') 'Enter the f(Z) value in MHz'
          READ(IN,*) FZE
        ENDIF
      ISP=0
      ISF=0
      IF(NCLOSD.NE.0) THEN
        ISP=1
        ISF=1
        WRITE(IERR,'(/1X,A)')
     :  'Contribution to the isotope mass shift (y/n)'
        READ(IN,'(A)') YES
        IF(YES.EQ.'Y'.OR.YES.EQ.'y') THEN
          ISP=2
          WRITE(IERR,'(A)')' Enter two masses'
          READ(IN,*) ZM1,ZM2
        ENDIF
        IF(INT(Z).LT.10) WRITE(IERR,'(A)') ' Warning Z below 10 !'
        WRITE(IERR,'(A)')
     : ' Contribution to the isotope field shift (y/n)'
        READ(IN,'(A)') YES
        IF(YES.EQ.'Y'.OR.YES.EQ.'y') THEN
          ISF=2
          WRITE(IERR,'(A)')
     :           ' Enter f(Z) (MHz/fm^2) and delta <r^2> (fm^2)'
          READ(IN,*) FZA,DR2
        ENDIF
      ENDIF
      RETURN
      END
*
*     -------------------------------------------------------------
*      C O R E I S O
*     -------------------------------------------------------------
*
*     CLOSED SHELL CONTRIBUTION
*
      SUBROUTINE COREISO
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=60,NOD=220)
      COMMON/ISO/GIS,RIS,COR,CO,FIELCOR,FIEL
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/GRADINT/GGRAD(NWD,NWD)
      COMMON/CLOSED/B1ELC(4),NCLOSD,IBK
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON/NON30/QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      POINTER(IQAZ,AZ(1))
      COMMON/NEL1/IQAZ
      GRADM2=ZERO
      CORER=ZERO
      FIELCOR=ZERO
      IF (NCLOSD.NE.0) THEN
        DO 1 I=1,NCLOSD
Cwwa      GRADM2=0.D0
          GRADM2=0.D0
   	  LA=LJCLSD(I)
          SUMI=DBLE(4*LA+2)
          IF(LJCLSD(I).EQ.0) THEN
            CFC=TWO*AZ(I)**2
            FIELCOR=FIELCOR+CFC
          ENDIF
          DO 2 J=1,I-1
  	    LB=LJCLSD(J)
	    IF(ITTK(2*LA,2*LB,2).NE.0) THEN
              GRADM1=ZERO
              SUMJ=DBLE(4*LB+2)
              C=CB(LA,LB,1)
              GRADM1=-C*GGRAD(I,J)**2
              GRADM2=SUMI*SUMJ*GRADM1
              COREG=COREG+GRADM2
	    ENDIF
   2      CONTINUE
   1    CONTINUE
      ENDIF
      COR=COREG
      RETURN
      END
*
*     -------------------------------------------------------------
*      C O R E O U T
*     -------------------------------------------------------------
*
*     CALCULATE THE MATRIX ELEMENTS BETEEN THE
*     OUTER ELECTRONS AND THE CORE ELECTRONS
*
      SUBROUTINE COREOUT
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=60,NOD=220)
      COMMON/ISO/GIS,RIS,COR,CO,FIELCOR,FIEL
      COMMON/GRADINT/GGRAD(NWD,NWD)
      COMMON/ANGCORE/CBANG(NWD,NWD)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/CLOSED/B1ELC(4),NCLOSD,IBK
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      POINTER(IQAZ,AZ(1))
      COMMON/NEL1/IQAZ
      GRADM2=ZERO
      CORE=ZERO
      CO=ZERO
      CFO=ZERO
      FIEL=ZERO
      DO 1 IAA=1,MAXORB
        LAA=LJCOMP(IAA)
	IAI=IAA+NCLOSD
        DO 2 IBB=IAA,MAXORB
          SUMI=-TWO*CBANG(IAA,IBB)
	  IF(DABS(SUMI).GT.EPS) THEN
	    IBI=IBB+NCLOSD
            LBB=LJCOMP(IBB)
	    IF(LBB.EQ.0) THEN
	      CFO=SUMI*AZ(IAI)*AZ(IBI)
	      FIEL=FIEL+CFO
	    ENDIF
            IF(NCLOSD.NE.0) THEN
              DO 3 J=1,NCLOSD
Cww             GRADM1=D0
                GRADM1=0.D0
	        LB=LJCLSD(J)
	        IF(ITTK(2*LAA,2*LB,2).NE.0) THEN
                    SUMJ=DBLE(4*LB+2)
                    C=CB(LAA,LB,1)
                    IF (IAA.EQ.IBB) THEN
                      GRADM1=-C*GGRAD(IAI,J)**2
                    ELSE
                      GRADM1=C*GGRAD(IAI,J)*GGRAD(J,IBI)
                    ENDIF
                    GRADM2=SUMI*SUMJ*GRADM1
                    CORE=CORE+GRADM2
                ENDIF
   3          CONTINUE
            ENDIF
          ENDIF
   2    CONTINUE
   1  CONTINUE
      CO=CORE
      RETURN
      END
*
*     -------------------------------------------------------------
*      S P E C
*     --------------------------------------------------------------
*
*     CALCULATE THE TOTAL CONTRIBUTION TO THE ISOTOPE SHIFT
* 
      SUBROUTINE SPEC
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=60,NOD=220)
*
      COMMON/GISOCON/AME,RY,ISP,ISF
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/ISO/GIS,RIS,COR,CO,FIELCOR,FIEL
      COMMON/COREC/ZM1,ZM2,FZE,FZA,DR2
      COMMON/INOUT/IN,IERR,IPRI,IUF
      COMMON/PRA/Z
      POINTER(IQAZ,AZ(1))
      COMMON/NEL1/IQAZ
      COMMON /RYDBERG/ RYDBRG,ZMU,ENERGY
      IF (ISP.EQ.0) RETURN
      RATIO=TWO*RY*AME/ZMU
      GRADM3= (GIS+RIS+COR+CO)
      TOTGRAD=RATIO*GRADM3
      WRITE(IPRI,108) ZMU
      WRITE(IPRI,109) 'Gradient Form       ',TOTGRAD
      DE1=-ENERGY*TWO*RYDBRG*AME/ZMU
      WRITE(IPRI,112)DE1
      IF(ISP.EQ.2) THEN
        WRITE(IPRI,114)' 3  ISOTOPE SHIFT FOR MASSES',ZM1,ZM2
        WRITE(IPRI,'(A)')'    a)  SPECIFIC MASS CONTRIBUTION'
        RATIO=(ZM2-ZM1)/((ZM1+AME)*(ZM2+AME))
        C=TWO*RY*AME
        DE1=ENERGY*TWO*RY*RATIO*AME
        RATIO=(ZM2-ZM1)/(ZM1*ZM2)
        TOTGRAD=C*GRADM3*RATIO
        WRITE(IPRI,109) 'Gradient Form       ',TOTGRAD
        WRITE(IPRI,116) DE1
      ENDIF
108   FORMAT(2(/)1X,'1  SPECIFIC MASS SHIFT CORRECTION FOR MASS'
     :,17X,F6.2)
109   FORMAT(6X,A20,26X,F14.8,'  cm-1')
110   FORMAT(6X,A20,26X,F14.4,'  cm-1')
112   FORMAT(/1X,'2  NORMAL MASS SHIFT CORRECTION  ',
     :18X,F14.8,'  cm-1')
114   FORMAT(2(/)A,23X,F6.2,' - ',F6.2,/)
116   FORMAT('    b)  NORMAL MASS CONTRIBUTION  ',
     :18X,F14.8,'  cm-1')
      RETURN
      END
*
*     -------------------------------------------------------------
*      P R I S O
*     -------------------------------------------------------------
*
      SUBROUTINE PRISO
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=60,NOD=220)
      CHARACTER*3 ELLAB,HEAD*15
      COMMON/ELLABEL/ELLAB(NWD)
      COMMON/INOUT/IN,IERR,IPRI,IUF
      COMMON/ISO/GIS,RIS,COR,CO,FIELCOR,FIEL
      COMMON/HE/HEAD
      COMMON/CLOSED/B1ELC(4),NCLOSD,IBK
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON/NON30/QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      COMMON/RYDBERG/RYDBRG,ZMU,ENERGY
      WRITE(IPRI,1)
      WRITE(IPRI,12) HEAD,ENERGY,RYDBRG
      WRITE(IPRI,10) NCLOSD,(ELLAB(I),I=1,NCLOSD)
      NCL=NCLOSD+1
      NAT=NCLOSD+MAXORB
      WRITE(IPRI,11) MAXORB,(ELLAB(I),I=NCL,NAT)
      WRITE(IPRI,13)
      WRITE(IPRI,3) COR
      WRITE(IPRI,2) GIS
      TINT=GIS+RIS
      WRITE(IPRI,5) TINT
      WRITE(IPRI,6) RIS
      WRITE(IPRI,4) CO
      TINT=GIS+RIS+COR+CO
C      TINTC=109737.31534*TINT
      TINTC=RYDBRG*TINT
      WRITE(IPRI,7) TINT,TINTC
    1 FORMAT(/8X,'............................................'//
     :        10X,'             ISOTOPE SHIFT     '/
     :        8X,'............................................'//)
   10 FORMAT(/10H There are ,I3,31H closed subshells common to all ,
     :  27H configurations as follows: //
     :  5X, 21(1X,A3))
   11 FORMAT(/10H There are,I3,21H orbitals as follows://
     : 5X,21(1X,A3):/5X,21(1X,A3)//)
   12 FORMAT(/4X,'   Atom  State   Energy (in a.u.)'4X,
     :'Rydberg constant was'/5X,A,F14.7,7X,F14.4,'cm-1'/)
   13 FORMAT(//4X,77('.')/
     ://20X,21('.'),//20X,'  MASS CONTRIBUTION  '/20X,21('.')/)
    2 FORMAT(/T50,F17.11,T70,'G integralas')
    3 FORMAT(/4X,'Total for Core',T30,F17.11)
    4 FORMAT(/4X,'Total for Core-Outer',T30,F17.11)
    5 FORMAT(4X,'Total for Integralas',T30,F17.11)	
    6 FORMAT(T50,F17.11,T70,'R integralas')
    7 FORMAT(/T30,52('.')/
     :	     /4X,'TOTAL CONTRIBUTION',T30,F17.11,T49,'(in a.u.)',
     :T62,F17.7,T81,'(in cm-1)')
      RETURN
      END
*
*     --------------------------------------------------------------
*      P R F I E L D 
*     --------------------------------------------------------------
*
*     CALCULATE THE FINITE VOLUME CORRECTION
*
*
      SUBROUTINE PRFIELD
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=60,NOD=220)
      PARAMETER (NELEMENT=50)
      COMMON/GISOCON/AME,RY,ISP,ISF
      COMMON/ISO/GIS,RIS,COR,CO,FIELCOR,FIEL
      COMMON/CLOSED/B1ELC(4),NCLOSD,IBK
      COMMON/INOUT/IN,IERR,IPRI,IUF
      COMMON/PRA/Z
*
      POINTER(IQAZ,AZ(1))
      COMMON/NEL1/IQAZ
      COMMON/COREC/ZM1,ZM2,FZE,FZA,DR2
      COMMON/RYDBERG/RYDBRG,ZMU,ENERGY
      DATA FVP/5.003461D-6/,FSP/8.339102D-6/
      WRITE(IPRI,101)
      WRITE(IPRI,103)'Total for core',FIELCOR
      WRITE(IPRI,103)'Total for outer electrons ',FIEL
      WRITE(IPRI,*)
      WRITE(IPRI,7)FIELCOR+FIEL
      A1=1.115*(ZMU**(1./3.))
      A2=2.151*(ZMU**(-1./3.))
      A3=1.742/ZMU
      REQ=A1+A2-A3
      REQ2=REQ*REQ
      DENS=FIELCOR+FIEL
      IF(INT(Z).LT.10) THEN
        FVC=FVP*DENS*FZE*REQ2*Z
      ELSE
        FVC=FVP*DENS*FZE*REQ2/Z
      ENDIF
      WRITE(IPRI,108) ZMU
      WRITE(IPRI,110) FVC
      WRITE(IPRI,111) REQ
Cww      IF(YF) THEN
      IF(ISF.NE.2) RETURN
        FSC=FSP*DENS*FZA*DR2/Z
        WRITE(IPRI,112) FSC
Cww      ENDIF
 101  FORMAT(///4X,77('.')///10X,46('.'),//10X,
     :'  ELECTRONIC CONTRIBUTION TO THE FIELD SHIFT  ',/10X,46('.')/)
 103  FORMAT(/4X,A,T47,F20.8)
    7 FORMAT(/T30,52('.')/
     :	     /4X,'TOTAL CONTRIBUTION',T30,F20.8,T50,'(in a.u.)')
 108  FORMAT(2(/)1X,'1  FINITE VOLUME CORRECTION FOR MASS'
     :,23X,F6.2)
 110  FORMAT(/6X,'Volume correction',27X,F16.8,'  cm-1')
 111  FORMAT(6X,'Req was',37X,F16.8,'  fm'/ )
 112  FORMAT(2(/)1X,'2  FIELD SHIFT CONTRIBUTION',22X,F16.8,'  cm-1')
      RETURN
      END
