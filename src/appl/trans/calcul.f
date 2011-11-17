*
*     -------------------------------------------------------------
*       C A L C U L
*     -------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville           September 1997   *
*
      SUBROUTINE CALCUL(NPAIR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL rel,vok
      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2)
      COMMON /MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :J1QN2(31,3),IJFUL(16)
      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      COMMON /STATE/QET1,QLBL1,QWT1,QJV1,QET2,QLBL2,QWT2,QJV2,NJV(2),
     : NVC(2),LGTH(2),NCF(2),QCFG1,QCFG2
      POINTER(QET1,ET1(1)),(QLBL1,LBL1(1)),(QWT1,WT1(1)),(QJV1,JV1(1)),
     :       (QET2,ET2(1)),(QLBL2,LBL2(1)),(QWT2,WT2(1)),(QJV2,JV2(1)),
     :       (QCFG1,CFG1(8,1)),(QCFG2,CFG2(8,1))
      COMMON /MULT/QSL,QSV,QIL,QIR,QJVL,QJVR
      POINTER(QSL,SL(1)),(QSV,SV(1)),(QIL,IL(1)),(QIR,IR(1)),
     :       (QJVL,JVL(1)),(QJVR,JVR(1))
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
    5 FORMAT(//,' JI =',I4,'  JF =',I4)
      KA=LAM
      KB=0
      IREZ=1
      IF(IFL.EQ.3)IREZ=2
      DO 1 IIFLL=1,IREZ
        IF (IIFLL.EQ.2) THEN
          IFL=4
          KA=LAM-1
          KB=1
        ENDIF
        DO 2 JII=1,NCF(1)
	  JI=JII
          IF(MOD(JI,10).EQ.0) print*,' JA = ',JI
          DO 3 JFI=1,NCF(2)
	    JF=JFI
            if(npair.eq.0) then
              print*,' ji = ',ji,' jf = ',jf
              stop
            end if
            NCONTR = 0
            IF(NBUG6.NE.0) WRITE(IWRITE,5) JI,JF
            JFF = JF+NCF(1)
            CALL SETUP(JI,JFF)
*
* --- test selection rules delta(S)=0 for Ek and M1
*                          delta(L)=0 for M1
            N1=2*IHSH-1
            ITIK=1
            IF(ITTK(J1QN1(N1,2)-1,J1QN2(N1,2)-1,2*KA).EQ.0)ITIK=0
            IF(ITTK(J1QN1(N1,3)-1,J1QN2(N1,3)-1,2*KB).EQ.0)ITIK=0
	    IF(ITIK.NE.0) THEN
              CALL NONTRANS(KA,KB,CL2,CV2)
*
* --- calculate the contribution of <JI/ O /JF> to the line
*     strengths for the npair (J,J') found
              IF(DABS(CL2).GT.ZERO.OR.DABS(CV2).GT.ZERO) THEN
                do 4 k = 1,npair
                  kl=(il(k)-1)*ncf(1)+ji
                  kr=(ir(k)-1)*ncf(2)+jf
                  fww=fline(k)*wt1(kl)*wt2(kr)
                  sl(k) = sl(k) + cl2*fww
                  if(vok) sv(k) = sv(k) + cv2*fww
                  if(ibug1.ne.0) then
                    print*,' pair = ',k,
     :                   ' wt1 = ',wt1(kl),' wt2 = ',wt2(kr)
                    print*,'              sl(pair) = ',sl(k)
                  end if
    4           continue
              ENDIF
            ENDIF
    3     CONTINUE
    2   CONTINUE
    1 CONTINUE
      RETURN
      END
