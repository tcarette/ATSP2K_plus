*     -------------------------------------------------------------
*       L V A L
*     -------------------------------------------------------------
*
*     Modified by Gediminas Gaigalas,                September 1997
*
*
      INTEGER FUNCTION LVAL(SYMBOL)
      CHARACTER*1 SYMBOL
      CHARACTER*26 SET
      DATA SET/'spdfghiklmnoqSPDFGHIKLMNOQ'/
*
      LOCATE = INDEX(SET,SYMBOL)
      IF ( LOCATE .LE. 13) THEN
            LVAL = LOCATE - 1
         ELSE
            LVAL = LOCATE - 14
      ENDIF
      RETURN
      END
