*
*     -------------------------------------------------------------
*       C A L C U L
*     -------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville           September 1997   *
*
      SUBROUTINE CALCUL(NPAIR,IPRINT,TOL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL rel,vok,iprint
      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2),iut(2)
      COMMON /MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :J1QN2(31,3),IJFUL(16)
      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      POINTER (QIORTH,IORTH(1))
c      COMMON /OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
c     1 ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
c     2 QIORTH
      COMMON /NOR/NCOM,NCLOSI,NCLOSF,NORBI,NORBF,IWAR
      COMMON /STATE/QET1,QLBL1,QWT1,QJV1,QET2,QLBL2,QWT2,QJV2,NJV(2),
     : NVC(2),LGTH(2),NCF(2),QCFG1,QCFG2
      POINTER(QET1,ET1(1)),(QLBL1,LBL1(1)),(QWT1,WT1(1)),(QJV1,JV1(1)),
     :       (QET2,ET2(1)),(QLBL2,LBL2(1)),(QWT2,WT2(1)),(QJV2,JV2(1)),
     :       (QCFG1,CFG1(8,1)),(QCFG2,CFG2(8,1))
      COMMON /MULT/QSL,QSV,QIL,QIR,QJVL,QJVR
      POINTER(QSL,SL(64)),(QSV,SV(64)),(QIL,IL(64)),(QIR,IR(64)),
     :       (QJVL,JVL(64)),(QJVR,JVR(64))
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      CHARACTER*64 configi,configf
      
    5 FORMAT(//,' JI =',I4,'  JF =',I4)
    6 FORMAT(/,1H ,'NORTH = ',I3,/
     1           ' JMU = ',I3,2X,'JNU = ',I3,/
     2           ' JMUP= ',I3,2X,'JNUP= ',I3,/
     3           ' NOVLPS = ',I3,/)
    
      print '(//A)', 'Entering calcul'
      IF (Iprint) then
        print *, 'Tolerance for print =', TOL
        print *,
     :     '  k  ', '  I1  ','  I2  ', 
     :      '  sl_cont  ', '  sv_cont  '
      END IF

!      Nbug6 = 1;
!      IBUG1=1
      KA=LAM
      KB=0
      IREZ=1
      IF(IFL.EQ.3)IREZ=2
      WRITE(ISCW,*) ' Double sum over CSF ...'
      DO 1 IIFLL=1,IREZ
        IF (IIFLL.EQ.2) THEN
          IFL=4
          KA=LAM-1
          KB=1
        ENDIF
        DO 2 JII=1,NCF(1)
          JI=JII
          IF(MOD(JI,100).EQ.0) WRITE(ISCW,*) '     JA = ',JI
          DO 3 JFI=1,NCF(2)
            JF=JFI
            if(npair.eq.0) then
              print*,' ji = ',ji,' jf = ',jf
              stop
            end if
            NCONTR = 0
            IF(NBUG6.NE.0) WRITE(IWRITE,5) JI,JF
            write(*,*) 
            write(*,*) ji,jf
            JFF = JF+NCF(1)
c            NOVLPS=0
c            JMUP=0
c            JNUP=0
c            JMU=0
c            JNU=0
*
            NC = 0
*
c            IF (NORTH .NE. 0) CALL NORTBP(JI,JFF)
c            IF (IWAR .EQ. 1) GO TO 3
c            IF(NBUG6.NE.0) 
c     :            WRITE(IWRITE,6) NORTH,JMU,JNU,JMUP,JNUP,NOVLPS
*
* --- set up the occupation and coupling arrays
*
            CALL SETUP(JI,JFF,LET)
            IF(LET.EQ.0) GO TO 3
*
* --- test selection rules delta(S)=0 for Ek and M1
*                          delta(L)=0 for M1
            N1=2*IHSH-1
            ITIK=1
            IF(ITTK(J1QN1(N1,2)-1,J1QN2(N1,2)-1,2*KA).EQ.0)ITIK=0
            IF(ITTK(J1QN1(N1,3)-1,J1QN2(N1,3)-1,2*KB).EQ.0)ITIK=0
            IF(ITIK.NE.0) THEN
              CALL NONTRANS(KA,KB,CL2,CV2)
!              print*,KA,KB,CL2,CV2, ' nontrans'
*
* --- calculate the contribution of <JI/ O /JF> to the line
*     strengths for the npair (J,J') found
              IF(DABS(CL2).GT.ZERO.OR.DABS(CV2).GT.ZERO) THEN
                do 4 k = 1,npair
                  kl=(il(k)-1)*ncf(1)+ji
                  kr=(ir(k)-1)*ncf(2)+jf
                  fww=fline(k) ! *wt1(kl)*wt2(kr)
                  sl(k) = sl(k) + cl2*fww
                  if(vok) sv(k) = sv(k) + cv2*fww
!                  if(iprint .and. (abs(cl2*fww) .gt. tol)) then
c                   Write(configi,'(8A8)') (cfg1(ii,il(k)),ii=1,8)
c                   WRITE(configf,'(8A8)') (cfg2(ii,ir(k)),ii=1,8)
                   print *,
     :                    k,kl,
     :                    kr,cl2*fww, cv2*fww
                  write(*,*) "   cl2=",Cl2,"cv2=",CV2, "fline=",fww
                  write(*,*) sl(k)
!                  end if
    4           continue
              ENDIF
            ENDIF
    3     CONTINUE
    2   CONTINUE
    1 CONTINUE
    
*              write(6,'(5F12.8)') SL
*               write(6,'(5F12.8)') sv

      print '(//A//)', 'Exiting calcul'

      RETURN
      END
