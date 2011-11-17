*
*     ------------------------------------------------------------------
*	F I N D
*     ------------------------------------------------------------------
*
      CHARACTER*3 FUNCTION FIND (I,OF,EL)
*
*  ---  THIS ROUTINE FINDS ELECTRONS IN ONE OF THREE LISTS
*
      PARAMETER (NWD=30)
      CHARACTER*3 OF(NWD,2),EL(NWD)
*
      IF ( I .LE. (NWD)) THEN
         FIND = EL(I)
        ELSE IF ( I .LE. (NWD)*2 ) THEN
         FIND = OF(I-(NWD),1)
        ELSE
         FIND = OF(I-2*(NWD),2)
      END IF
      RETURN
      END
