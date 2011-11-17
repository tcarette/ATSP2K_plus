*
*     -------------------------------------------------------------
*      C A L C U L A T I O N
*     -------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville            February 1994   * 
*
      LOGICAL FUNCTION CALCULATION(KL)
      COMMON /OPERAT/ ICOLOM,ISOTOP,IORBORB
      CALCULATION=.FALSE.
      IF(ICOLOM.EQ.1) THEN
        CALCULATION=.TRUE.
      ELSEIF(IORBORB.EQ.1) THEN
        CALCULATION=.TRUE.
      ELSEIF(ISOTOP.EQ.1) THEN
        IF(KL.EQ.1) CALCULATION=.TRUE.
      ENDIF
      RETURN
      END
