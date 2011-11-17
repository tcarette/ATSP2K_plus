*
*     -------------------------------------------------------------
*      I T L S
*     -------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      FUNCTION ITLS(IK,ID,KG,KGG,BD,IBT,BT,KG1,KGG1,ITP,ITG,IQ)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/RIBOLS/IMPTLS(40),IMGTLS(40),IMPNLS(40),IMGNLS(40)
      COMMON/RIBOLSF/IMPTLSF(8),IMGTLSF(8),IMPNLSF(8),IMGNLSF(8)
      COMMON/RIBOF/IMPTF(238),IMGTF(238),IMPNF(238),IMGNF(238)
      COMMON/RIBOLS3/IMPTLS3(90),IMGTLS3(90),IMPNLS3(90),IMGNLS3(90)
      DIMENSION ID(7),IK(7),IBT(7),BT(3),BD(3)
      ITLS=0
      IF(ID(3).GT.9) RETURN
      KG1=2*KG
      KGG1=2*KGG
      IF(ITTK(ID(5),IK(5),KG1).EQ.0)RETURN
      IF(ITTK(ID(6),IK(6),KGG1).EQ.0)RETURN
      ITK=IK(1)
      ITD=ID(1)
      IF(ID(3).LT.3) THEN
        ITP1=IMPTLS(ITK)
        ITP=IMPTLS(ITD)
        IF(ITP1.NE.ITP)RETURN
        ITG1=IMGTLS(ITK)
        ITG=IMGTLS(ITD)
      ELSEIF(ID(3).EQ.3) THEN
	IF(ITK.GT.300) THEN
	  IF(ITD.LT.300) CALL MES(51)
          IF(ID(4).GT.2) CALL MES(11)
          IF(IK(4).GT.2) CALL MES(11)
	  ITK=ITK-300
	  ITD=ITD-300
          ITP1=IMPTLSF(ITK)
          ITP=IMPTLSF(ITD)
          IF(ITP1.NE.ITP)RETURN
          ITG1=IMGTLSF(ITK)
          ITG=IMGTLSF(ITD)
	ELSE
	  IF(ITD.GT.300) CALL MES(51)
          ITP1=IMPTF(ITK)
          ITP=IMPTF(ITD)
          IF(ITP1.NE.ITP)RETURN
          ITG1=IMGTF(ITK)
          ITG=IMGTF(ITD)
	ENDIF
      ELSE
        IF(ID(4).GT.2) CALL MES(11)
        IF(IK(4).GT.2) CALL MES(11)
        ITP1=IMPTLS3(ITK)
        ITP=IMPTLS3(ITD)
        IF(ITP1.NE.ITP)RETURN
        ITG1=IMGTLS3(ITK)
        ITG=IMGTLS3(ITD)
      ENDIF
      IF(ITG1.NE.ITG)RETURN
      ITLS=1
      IBT(2)=ID(2)
      IBT(3)=ID(3)
      IBT(4)=ID(4)+IQ
      BT(3)=BD(3)+HALF*DBLE(IQ)
      RETURN
      END
