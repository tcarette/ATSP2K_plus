*
*     -------------------------------------------------------------
*      R U M T
*     -------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE RUMT(KNT,LL,LQ,LS,L)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/MT/MT(40)
      COMMON/SKMT2/MTF(8),MT3(90)
      IF(LL.LT.3) THEN
        KT=MT(KNT)
      ELSEIF(LL.EQ.3) THEN
	IF(KNT.GT.300) THEN
	  KNTMIN=KNT-300
          KT=MTF(KNTMIN)
        ELSE
          CALL RUMT67(KNT,NR,LQ,LS,L)
          RETURN
        ENDIF
      ELSE
        KT=MT3(KNT)
      ENDIF
      LQ=JTHN(KT,3,100)
      LS=JTHN(KT,2,100)
      L=JTHN(KT,1,100)
      RETURN
      END
