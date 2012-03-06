*
*     -------------------------------------------------------------
*      S A V E N O N 
*     -------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville           September 1997   *
*
      SUBROUTINE SAVENON(I,A,KL,LA,LB,LC,LD,JJI,JJF,IPTR)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL REL,VOK
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1))
!      POINTER(IQRL,RLINT(MAXORB,1)),(IQRV,RVINT(MAXORB,1)),
!     :       (IQOV,OVLP(MAXORB,1))
!      COMMON/RDINT/IQRL,IQRV,IQOV
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,MAXORB
      COMMON/EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      COMMON/NTRM/NTERMS
      COMMON/SAVECOM/LRHO,LSIG,CL2,CV2
      COMMON /INOUT2/IANG
      CL=ZERO 
      CV=ZERO
      NTERMS=NTERMS+1
      RMET=RMETR(LRHO,LSIG)
      CL=A*RMET
*
      WRITE(IANG) IFL,JI,JF,LB,LD,CL
*
      IF(VOK)CV=CL ! *RVINT(LB,LD)
      CL=CL ! *RLINT(LB,LD)
      CL2=CL2+CL
      CV2=CV2+CV
      RETURN
      END
