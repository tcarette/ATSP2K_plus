*
*     ------------------------------------------------------------------
*	C F G I N 2
*     ------------------------------------------------------------------
*
      SUBROUTINE cfgin2(MCFG,KCFG)
*
* --- Read two sets of configurations and determine the orthogonality
* --- conditions between them
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=128)
      CHARACTER*3 EL(NWD), ELC(NWD), JAJCMP(NWD,3)*1
      CHARACTER INPUT(2)*24,HEADI*72,HEADF*72,HEADER*72
      CHARACTER*7 LABEL(2)
      character*3 elras,elrasi,elrasf
      DATA LABEL/'Initial','Final  '/
      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON /DBG  /IBUGM
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2),iut(2)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      COMMON /NOR/NCOM,NCLOSI,NCLOSF,NORBI,NORBF,IWAR
      COMMON /FO/IWF(2),NCLOS(2),NOPEN(2),MAXNFO
      COMMON /RAS/ ninac(11,2),nras1(11,10,2),nras2(11,2),
     :             nras3(11,10,2),NGAS1(2),NGAS3(2),
     :             nl(11,2),elras(2*nwd),elrasi(nwd),elrasf(nwd),
     :             itab(nwd,2),lmax(2)
*
    3 FORMAT(20(1X,A3))
    7 FORMAT(A30,I3,I4)
   22 FORMAT(// 7H STATE ,' (WITH',I5,' CONFIGURATIONS):'/1H ,31(1H-)/)
   23 FORMAT(/10H THERE ARE,I3,21H ORBITALS AS FOLLOWS://
     1 5X,21(1X,A3):/5X,21(1X,A3))
*
*
* --- ANALYZE INITIAL AND FINAL STATE DATA
*
      CALL ANALY2(MCFG,KCFG,EL)
*
      REWIND(UNIT=iuc(1))
      REWIND(UNIT=iuc(2))
*
      MAXORB = NCOM + NORBI + NORBF
      NCFG = MCFG + KCFG
*
* --- allocate the memory
*
      if(ibugm.ne.0) print*,' qiajcmp allocation: maxorb   = ',maxorb
      call alloc(qiajcmp,maxorb,4)
      if(ibugm.ne.0) print*,' qljcomp allocation: maxorb   = ',maxorb
      call alloc(qljcomp,maxorb,4)
      if(ibugm.ne.0) print*,' qnjcomp allocation: maxorb   = ',maxorb
      call alloc(qnjcomp,maxorb,4)
      if(ibugm.ne.0) print*,' qnoc    allocation: ncfg     = ',ncfg
      call alloc(qnoc,ncfg,4)
      if(ibugm.ne.0) print*,' qnelcsh allocation: 8*ncfg   = ',8*ncfg
      call alloc(qnelcsh,8*ncfg,4)
      if(ibugm.ne.0) print*,' qnocorb allocation: 8*ncfg   = ',8*ncfg
      call alloc(qnocorb,8*ncfg,4)
      if(ibugm.ne.0) print*,' qj1     allocation: 15*ncfg   = ',15*ncfg
      call alloc(qj1,15*ncfg,4)
*
* --- set up the electrons
*
      READ(EL,'(A3)') (IAJCMP(I),I=1,MAXORB)
      READ(EL,'(3A1)')((JAJCMP(I,J),J=1,3),I=1,MAXORB)
*
* --- set up of ljcomp
*
      DO 60 I = 1,MAXORB
      IF (JAJCMP(I,1) .EQ. ' ') THEN
         JAJCMP(I,1) = JAJCMP(I,2)
         JAJCMP(I,2) = JAJCMP(I,3)
         JAJCMP(I,3) = ' '
      ENDIF
      LJCOMP(I) = LVAL(JAJCMP(I,2))
      NJCOMP(I) = ICHAR(JAJCMP(I,1)) - ICHAR('1') + 1
   60 CONTINUE
*
* ---- check common closed shells
*
      IF (NCLOSI .NE. NCLOSF)
     :   STOP ' Common closed shells not the same in the two states'
*
      READ(iuc(1),7) HEADI,iclosdi,iwfi
      READ(iuc(2),7) HEADF,iclosdf,iwff
      HEADER = HEADI(1:30)//'=>'//HEADF(1:30)
*
* --- check closed shells further
*
      READ(iuc(1),3) (ELC(I),I=1,NCLOSI)
      READ(iuc(2),3) (EL(I),I=1,NCLOSF)
      DO 1 I = 1,NCLOSF
         J = 1
    2    IF (EL(I) .NE. ELC(J) ) THEN
            J = J+1
            IF (J .LE. NCLOSI) THEN
               GO TO 2
              ELSE
               STOP ' Common closed sub-shells not the same'
            END IF
         END IF
    1 CONTINUE
*
* --- Check if electron lists are present
*
      if (iwfi .gt. iclosdi) read(iuc(1),3) (elc(i),i=iclosdi+1,iwfi)
      if (iwff .gt. iclosdf) read(iuc(2),3) (el(i),i=iclosdf+1,iwff)
*
* --- get initial state configurations
*
      CALL GSTATE(1,MCFG)
      CALL GSTATE(MCFG+1,NCFG)
*
*  --- check the data
*
      CALL CFGTST(NCFG,QLJCOMP,QNOC,QNELCSH,QNOCORB,QJ1)
      REWIND(UNIT=iuc(1))
      REWIND(UNIT=iuc(2))
      RETURN
      END
