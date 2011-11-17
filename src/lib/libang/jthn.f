*
*     -------------------------------------------------------------
*      J T H N
*     -------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      FUNCTION JTHN(K,N,I)
      IF(N.EQ.1) THEN
        JTHN=MOD(K,I)
      ELSEIF(N.EQ.2) THEN
        JTHN=MOD(K/I,I)
      ELSEIF(N.EQ.3) THEN
        JTHN=MOD(K/(I*I),I)
      ELSEIF(N.EQ.4) THEN
        JTHN=MOD(K/(I*I*I),I)
      ELSEIF(N.EQ.5) THEN
        JTHN=MOD(K/(I*I*I*I),I)
      ELSE
	WRITE(6,'(A)') ' ERROR IN JTHN '
	STOP
      ENDIF
      RETURN
      END
