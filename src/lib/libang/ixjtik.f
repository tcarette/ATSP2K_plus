*
*     -------------------------------------------------------------
*      I X J T I K
*     -------------------------------------------------------------
*
*                                                                  *
*     CHESKED TRIANGULAR CONDITIONS FOR 6j COEFFICIENT             *
*                                                                  *
*     | I/2  J/2  K/2 |            IXJTIK=1 - IF NOT SATISFY       *
*     | L/2  M/2  N/2 |            IXJTIK=0 - IN OVER CASES        *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                          December 1993   *
*
      FUNCTION IXJTIK(I,J,K,L,M,N)
      IXJTIK=0
      IF(ITTK(I,J,K).EQ.0)RETURN
      IF(ITTK(I,M,N).EQ.0)RETURN
      IF(ITTK(L,J,N).EQ.0)RETURN
      IF(ITTK(L,M,K).EQ.0)RETURN
      IXJTIK=1
      RETURN
      END
