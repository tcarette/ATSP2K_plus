*
*     -----------------------------------------------------------------
*     N U M T E R
*     -----------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      FUNCTION NUMTER(I2Q,I2S,I2L,L,NK,ND)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      DIMENSION LP(3),LG(3),LP3(6),LG3(6)
      COMMON/MT/MT(40)
      COMMON/SKMT2/MTF(8),MT3(90)
      COMMON /MT67/ M76(238)
      EXTERNAL TRMF,TERMLS,TRMF15
      DATA LP/1,3,9/,LG/2,8,40/,
     *LP3/1,11,23,37,53,71/,LG3/10,22,36,52,70,90/
      NUMTER=0
      IF(L.GT.9)RETURN
      IN=(I2Q*100+I2S)*100+I2L
      IF(L.LT.3) THEN
        IL=L+1
        JP=LP(IL)
        JG=LG(IL)
        J=JP
    2   IF(IN-MT(J))3,1,3
    3   J=J+1
        IF(J.GT.JG)RETURN
        GO TO 2
    1   NUMTER=J
      ELSEIF(L.EQ.3) THEN
	IF(MAX0(NK,ND).LT.3) THEN
          JP=1
          JG=8
          J=JP
    6     IF(IN-MTF(J))4,5,4
    4     J=J+1
          IF(J.GT.JG)RETURN
          GO TO 6
    5     NUMTER=J+300
	ELSE
          IF((I2Q/2)*2.EQ.I2Q) THEN
            JP=1
            JG=119
          ELSE
            JP=120
            JG=238
          ENDIF
          J=JP
   12     MF=JTHN(M76(J),1,1000000)
          IF(IN-MF)13,11,13
   13     J=J+1
          IF(J.GT.JG)RETURN
          GO TO 12
   11     NUMTER=J
        ENDIF
      ELSE
        IL=L-3
        JP=LP3(IL)
        JG=LG3(IL)
        J=JP
   22   IF(IN-MT3(J))23,21,23
   23   J=J+1
        IF(J.GT.JG)RETURN
        GO TO 22
   21   NUMTER=J
      ENDIF
      RETURN
      END
