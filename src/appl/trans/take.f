*
*     -------------------------------------------------------------
*      T A K E 
*     -------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville           September 1997   *
*
      SUBROUTINE TAKE(NPAIR,GFILE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL REL,VOK
      character GFILE*64
      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON /INOUT2/IANG
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1))
      POINTER(IQRL,RLINT(MAXORB,1)),(IQRV,RVINT(MAXORB,1)),
     :       (IQOV,OVLP(MAXORB,1))
      COMMON/RDINT/IQRL,IQRV,IQOV
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,MAXORB
      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      COMMON /NTRM/NTERMS
      COMMON /STATE/QET1,QLBL1,QWT1,QJV1,QET2,QLBL2,QWT2,QJV2,NJV(2),
     : NVC(2),LGTH(2),NCF(2),QCFG1,QCFG2
      POINTER(QET1,ET1(1)),(QLBL1,LBL1(1)),(QWT1,WT1(1)),(QJV1,JV1(1)),
     :       (QET2,ET2(1)),(QLBL2,LBL2(1)),(QWT2,WT2(1)),(QJV2,JV2(1)),
     :       (QCFG1,CFG1(8,1)),(QCFG2,CFG2(8,1))
      COMMON /MULT/QSL,QSV,QIL,QIR,QJVL,QJVR
      POINTER(QSL,SL(1)),(QSV,SV(1)),(QIL,IL(1)),(QIR,IR(1)),
     :       (QJVL,JVL(1)),(QJVR,JVR(1))
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      IFL1=IFL
      NTERMS=0
      CL=ZERO 
      CV=ZERO
      JI=0
      JF=0
    1 READ(IANG,END=3) IFL,JJJI,JJJF,LLB,LLD,CL1
      IF(JI.NE.JJJI) THEN
	IF(MOD(JJJI,10).EQ.0) PRINT*,' JA = ',JJJI
      ENDIF
      IF(JI.EQ.0.AND.JF.EQ.0) THEN
	IF(IFL.NE.IFL1) THEN
	  PRINT*,'  WRONG ANGLES FILE  ',GFILE
	  STOP
        ENDIF
	JI=JJJI
	JF=JJJF
	LB=LLB
	LD=LLD
        JFF = JF+NCF(1)
        CALL SETUP(JI,JFF)
        CL=CL1*RLINT(LB,LD)
        IF(VOK)CV=CL1*RVINT(LB,LD)
        CL2=CL
        CV2=CV
      ELSEIF(JJJI.EQ.JI.AND.JJJF.EQ.JF) THEN
	LB=LLB
	LD=LLD
        CL=CL1*RLINT(LB,LD)
        IF(VOK)CV=CL1*RVINT(LB,LD)
        CL2=CL2+CL
        CV2=CV2+CV
      ELSE
        do 2 k = 1,npair
          kl=(il(k)-1)*ncf(1)+ji
          kr=(ir(k)-1)*ncf(2)+jf
          fww=fline(k)*wt1(kl)*wt2(kr)
          sl(k) = sl(k) + cl2*fww
          if(vok) sv(k) = sv(k) + cv2*fww
          if(ibug1.ne.0) then
            print*,' pair = ',k,' wt1 = ',wt1(kl),' wt2 = ',wt2(kr)
            print*,'              sl(pair) = ',sl(k)
          end if
    2   continue
	JI=JJJI
	JF=JJJF
        JFF = JF+NCF(1)
	LB=LLB
	LD=LLD
        CALL SETUP(JI,JFF)
        CL=CL1*RLINT(LB,LD)
        IF(VOK)CV=CL1*RVINT(LB,LD)
        CL2=CL
        CV2=CV
      ENDIF
      NTERMS =NTERMS+1
      GO TO 1
    3 CONTINUE
      do 4 k = 1,npair
        kl=(il(k)-1)*ncf(1)+ji
        kr=(ir(k)-1)*ncf(2)+jf
        fww=fline(k)*wt1(kl)*wt2(kr)
        sl(k) = sl(k) + cl2*fww
        if(vok) sv(k) = sv(k) + cv2*fww
      if(ibug1.ne.0) then
        print*,' pair = ',k,' wt1 = ',wt1(kl),' wt2 = ',wt2(kr)
        print*,'              sl(pair) = ',sl(k)
      end if
    4 continue
      RETURN
      END
