*     -----------------------------------------------------------------
*       L V A L
*     -----------------------------------------------------------------
*
*
      INTEGER FUNCTION LVAL(SYMBOL) 
      CHARACTER*1 SYMBOL 
      CHARACTER*22 SET
      DATA SET/'spdfghiklmnSPDFGHIKLMN'/ 
*
      LOCATE = INDEX(SET,SYMBOL) 
      IF ( LOCATE .LE. 11) THEN 
            LVAL = LOCATE - 1 
         ELSE 
            LVAL = LOCATE - 12 
      ENDIF 
      RETURN 
      END 
