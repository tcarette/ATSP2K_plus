*
*     -------------------------------------------------------------
*      S A V E N O N
*     -------------------------------------------------------------
*                                                                  *
*     THIS SUBROUTINE FOR            B I O R T O G                 *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville           September 1997   *
*
      SUBROUTINE SAVENON(I,A,KL,LA,LB,LC,LD,JA,JB,IPTR)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL REL,VOK
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      POINTER(IQRL,RLINT(MAXORB,1)),(IQRV,RVINT(MAXORB,1)),
     :       (IQOV,OVLP(MAXORB,1))
      COMMON/RDINT/IQRL,IQRV,IQOV
      COMMON/EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      COMMON/NTRM/NTERMS
      COMMON/SAVECOM/LRHO,LSIG,CL2,CV2
      COMMON/SIGNF/SIGNFA
      IF(I.EQ.200) THEN
        CL=ZERO 
        CV=ZERO
        NTERMS=NTERMS+1
        RMET=RMETR(LRHO,LSIG)
        CL=A*RMET
        IF(VOK)CV=CL ! *RVINT(LB,LD)
!        CL=CL*RLINT(LB,LD)
        CL2=CL2+CL
        CV2=CV2+CV
      ELSE
        A=A*SIGNFA
        CALL SAVELS(I,A,KL,LA,LB,LC,LD,JA,JB,IPTR)
      ENDIF
      RETURN
      END
