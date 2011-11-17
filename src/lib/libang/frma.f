*     ..........................................................   :
*                                                                  :
*          Block                                                   :
*                  Standard Quantities  -  S Q L S F               :
*                         Part Two                                 :
*                                                                  :
*     Written by G. Gaigalas,                                      :
*                Institute of Theoretical Physics and Astronomy    :
*                Vilnius,  Lithuania                               :
*                                                       May 1995   :
*                                                                  :
*     ..........................................................   :
*
*
*     -----------------------------------------------------------------
*      F R M A
*     -----------------------------------------------------------------
*
      SUBROUTINE FRMA(J1,LQ,LL,LS,J2,LQS,LLS,LSS,COEF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/RIBOF/IMPTF(238),IMGTF(238),IMPNF(238),IMGNF(238)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IV(78)
      DATA IV/1,49,147,1,539,49,33,1,11,1176,12936,
     :2646,4116,181104,296352,543312,7392,
     :2,162,22176,20790,66528,3234,
     :0,0,0,0,0,
     :0,0,0,0,0,
     :0,0,0,0,0,0,
     :0,0,0,0,0,0,0,
     :980628,1,2,196,1176,51744,7392,25872,37044,5488,
     :4116,1992144,121968,2646,1008,1,
     :0,0,0,0,0,0,0,
     :12603360,491040,11319,142688,7546,
     :0,0,0,0/
      ISKA=0
      JI1=IMPTF(J1)
      JI2=IMPNF(J2)
      IF(JI1.NE.JI2) RETURN
      IF(J1.LT.J2) THEN
        JJ1=J1
        JJ2=J2
        IQ=LQ
        IL=LL
        IS=LS
        IQS=LQS
        ILS=LLS
        ISS=LSS
      ELSE
        JJ1=J2
        JJ2=J1
        IQ=LQS
        IL=LLS
        IS=LSS
        IQS=LQ
        ILS=LL
        ISS=LS
      ENDIF
      IF(JJ1.LT.12) THEN
        CALL FRMA01(JJ1,JJ2,ISKA)
        IVAR=IV(JJ1)
      ELSEIF(JJ1.LT.18) THEN
        CALL FRMA02(JJ1,JJ2,ISKA)
        IVAR=IV(JJ1)
      ELSEIF(JJ1.LT.24) THEN
        CALL FRMA03(JJ1,JJ2,ISKA)
        IVAR=IV(JJ1)
      ELSEIF(JJ1.LT.29) THEN
        CALL FRMA04(JJ1,JJ2,ISKA,IV(JJ1))
        IVAR=IV(JJ1)
      ELSEIF(JJ1.LT.34) THEN
        CALL FRMA05(JJ1,JJ2,ISKA,IV(JJ1))
        IVAR=IV(JJ1)
      ELSEIF(JJ1.LT.40) THEN
        CALL FRMA06(JJ1,JJ2,ISKA,IV(JJ1))
        IVAR=IV(JJ1)
      ELSEIF(JJ1.LT.47) THEN
        CALL FRMA07(JJ1,JJ2,ISKA,IV(JJ1))
        IVAR=IV(JJ1)
      ELSEIF(JJ1.LT.57) THEN
        CALL FRMA08(JJ1,JJ2,ISKA)
        IVAR=IV(JJ1)
      ELSEIF(JJ1.LT.63) THEN
        CALL FRMA09(JJ1,JJ2,ISKA)
        IVAR=IV(JJ1)
      ELSEIF(JJ1.LT.70) THEN
        CALL FRMA10(JJ1,JJ2,ISKA,IV(JJ1))
        IVAR=IV(JJ1)
      ELSEIF(JJ1.LT.75) THEN
        CALL FRMA11(JJ1,JJ2,ISKA)
        IVAR=IV(JJ1)
      ELSEIF(JJ1.LT.79) THEN
        CALL FRMA12(JJ1,JJ2,ISKA,IV(JJ1))
        IVAR=IV(JJ1)
      ELSEIF(JJ1.LT.86) THEN
        CALL FRMA13(JJ1,JJ2,ISKA,IVAR)
      ELSEIF(JJ1.LT.91) THEN
        CALL FRMA14(JJ1,JJ2,ISKA,IVAR)
      ELSEIF(JJ1.LT.96) THEN
        CALL FRMA15(JJ1,JJ2,ISKA,IVAR)
      ELSEIF(JJ1.LT.101) THEN
        CALL FRMA16(JJ1,JJ2,ISKA,IVAR)
      ELSEIF(JJ1.LT.107) THEN
        CALL FRMA17(JJ1,JJ2,ISKA,IVAR)
      ELSEIF(JJ1.LT.113) THEN
        CALL FRMA18(JJ1,JJ2,ISKA,IVAR)
      ELSEIF(JJ1.LT.120) THEN
        CALL FRMA19(JJ1,JJ2,ISKA,IVAR)
      ENDIF
      IF(ISKA.GT.0) THEN
        COEF=DSQRT(DBLE(ISKA)/DBLE(IVAR))
      ELSEIF(ISKA.EQ.0) THEN
        COEF=ZERO
        RETURN
      ELSE
        COEF=-DSQRT(-DBLE(ISKA)/DBLE(IVAR))
      ENDIF
      IF(J1.GT.J2) THEN
        IE1=LQ+LL+LS-LQS-LLS-LSS+2*3
        IF((IE1/4)*4.NE.IE1)COEF=-COEF
      ENDIF
      RETURN
      END
