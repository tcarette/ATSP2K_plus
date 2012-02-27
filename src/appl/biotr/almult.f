*
*     -------------------------------------------------------------
*      A L M U L T 
*     -------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville           September 1997   *
*
      SUBROUTINE ALMULT(NPAIR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAl REL,VOK
      COMMON /DBG/IBUGM
      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      COMMON /STATE/QET1,QLBL1,QWT1,QJV1,QET2,QLBL2,QWT2,QJV2,NJV(2),
     : NVC(2),LGTH(2),NCF(2),QCFG1,QCFG2
      POINTER(QET1,ET1(1)),(QLBL1,LBL1(1)),(QWT1,WT1(1)),(QJV1,JV1(1)),
     :       (QET2,ET2(1)),(QLBL2,LBL2(1)),(QWT2,WT2(1)),(QJV2,JV2(1)),
     :       (QCFG1,CFG1(8,1)),(QCFG2,CFG2(8,1))
      COMMON /MULT/QSL,QSV,QIL,QIR,QJVL,QJVR
      POINTER(QSL,SL(1)),(QSV,SV(1)),(QIL,IL(1)),(QIR,IR(1)),
     :       (QJVL,JVL(1)),(QJVR,JVR(1))
*
* --- allocate the /MULT/
*
      nmult = nvc(1)*nvc(2)
      if(ibugm.ne.0) print*,' nvc(1) = ',nvc(1),' nvc(2) = ',nvc(2),
     :' nmult = ',nmult
      if(ibugm.ne.0) print*,' qsl     allocation: nmult    = ',nmult
      call alloc(qsl,nmult,8)
      if(ibugm.ne.0) print*,' qsv     allocation: nmult    = ',nmult
      call alloc(qsv,nmult,8)
      if(ibugm.ne.0) print*,' qil     allocation: nmult    = ',nmult
      call alloc(qil,nmult,4)
      if(ibugm.ne.0) print*,' qir     allocation: nmult    = ',nmult
      call alloc(qir,nmult,4)
      if(ibugm.ne.0) print*,' qjvl    allocation: nmult    = ',nmult
      call alloc(qjvl,nmult,4)
      if(ibugm.ne.0) print*,' qjvr    allocation: nmult    = ',nmult
      call alloc(qjvr,nmult,4)
*
* --- determine the number of (J,J') pairs satisfying selection rules
*
      write(*,*) "in almult"
      npair = 0
      do 1 i = 1,nvc(1)
        do 2 j = 1,nvc(2)
          write(*,*) jv1(i),jv2(j),LAM+LAM
          atst = tritst(jv1(i),jv2(j),LAM+LAM)
          if (atst .eq. 0.0) then
            npair = npair + 1
            il(npair) = i
            ir(npair) = j
            jvl(npair) = jv1(i)
            jvr(npair) = jv2(j)
            sl(npair) = 0.0
            sv(npair) = 0.0
            print*,' npair = ',npair,' jvl = ', jvl(npair),
     :             ' jvr = ',jvr(npair)
          end if
    2   continue
    1 continue
      RETURN
      END
