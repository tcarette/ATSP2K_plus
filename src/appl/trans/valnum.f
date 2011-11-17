*     -------------------------------------------------------------
*       V A L N U M
*     -------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville           September 1997   *
*
      SUBROUTINE VALNUM(SYMBOL,I)
      CHARACTER*1 SET(9),SYMBOL
      DATA SET/'1','2','3','4','5','6','7','8','9'/
      SYMBOL = SET(I)
      RETURN
      END
