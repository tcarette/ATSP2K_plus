*
*     -------------------------------------------------------------
*       N U M V A L
*     -------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      INTEGER FUNCTION NUMVAL(SYMBOL)
      CHARACTER*1 SYMBOL
      CHARACTER*13 SET
      DATA SET/'0123456789 aA'/
*
      LOCATE = INDEX(SET,SYMBOL)
      IF ( LOCATE .LE. 10) THEN
        NUMVAL = LOCATE - 1
      ELSEIF( LOCATE .EQ. 11) THEN
        NUMVAL = 0
      ELSE
        NUMVAL = 10
      ENDIF
      RETURN
      END
