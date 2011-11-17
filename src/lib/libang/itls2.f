*
*     -------------------------------------------------------------
*      I T L S 2
*     -------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      FUNCTION ITLS2(IK,ID,KG1,KG2,BD,IBT,BT,ITP,ITG,IQ)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/RIBOLS/IMPTLS(40),IMGTLS(40),IMPNLS(40),IMGNLS(40)
      COMMON/RIBOLSF/IMPTLSF(8),IMGTLSF(8),IMPNLSF(8),IMGNLSF(8)
      COMMON/RIBOF/IMPTF(238),IMGTF(238),IMPNF(238),IMGNF(238)
      COMMON/RIBOLS3/IMPTLS3(90),IMGTLS3(90),IMPNLS3(90),IMGNLS3(90)
      DIMENSION ID(7),IK(7),IBT(7),BT(3),BD(3)
      ITLS2=0
      IF(ID(3).GT.9) RETURN
      IF(ITTK(ID(5),IK(5),KG1).EQ.0)RETURN
      IF(ITTK(ID(6),IK(6),KG2).EQ.0)RETURN
      ITK=IK(1)
      ITD=ID(1)
      IF(ID(3).LT.3) THEN
        ITP1=IMPNLS(ITK)
        ITP=IMPTLS(ITD)
        IF(ITP1.NE.ITP)RETURN
        ITG1=IMGNLS(ITK)
        ITG=IMGTLS(ITD)
      ELSEIF(ID(3).EQ.3) THEN
	IF(ITK.GT.300) THEN
	  IF(ITD.LT.300) CALL MES(52)
          IF(ID(4).GT.2) CALL MES(12)
          IF(IK(4).GT.2) CALL MES(12)
	  ITK=ITK-300
	  ITD=ITD-300
          ITP1=IMPNLSF(ITK)
          ITP=IMPTLSF(ITD)
          IF(ITP1.NE.ITP)RETURN
          ITG1=IMGNLSF(ITK)
          ITG=IMGTLSF(ITD)
	ELSE
	  IF(ITD.GT.300) CALL MES(52)
          ITP1=IMPNF(ITK)
          ITP=IMPTF(ITD)
          IF(ITP1.NE.ITP)RETURN
          ITG1=IMGNF(ITK)
          ITG=IMGTF(ITD)
	ENDIF
      ELSE
        IF(ID(4).GT.2) CALL MES(12)
        IF(IK(4).GT.2) CALL MES(12)
        ITP1=IMPNLS3(ITK)
        ITP=IMPTLS3(ITD)
        IF(ITP1.NE.ITP)RETURN
        ITG1=IMGNLS3(ITK)
        ITG=IMGTLS3(ITD)
      ENDIF
      IF(ITG1.NE.ITG)RETURN
      ITLS2=1
      IBT(2)=ID(2)
      IBT(3)=ID(3)
      IBT(4)=ID(4)+IQ
      BT(3)=BD(3)+HALF*DBLE(IQ)
      RETURN
      END
