*
*     -------------------------------------------------------------
*      S L S
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
*                                                    (qls)         *
*     REDUCED MATRIX ELEMENT             (nl QLS::: a  :::nl QLS)  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE SLS(L,IT,LQ,LL,LS,ITS,LQS,LLS,LSS,S)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IC(5,2)
      DATA IC/3,5,6,8,2,9,24,25,40,5/
      IF(L.GT.3) THEN
        IF(L.GT.9) RETURN
        IF(LSS.EQ.1) THEN
          CALL SUBLS(IT,ITS,L,S)
        ELSE
          CALL SUBLS(ITS,IT,L,S)
          IE1=LQ+LL+LS-LQS-LLS-LSS+2*L
          IF((IE1/4)*4.NE.IE1)S=-S
        ENDIF
      ELSEIF(L.EQ.3) THEN
	IF(IT.GT.300) THEN
	  IF(ITS.LT.300) THEN
	    WRITE(6,'(A)') ' ERROR IN    S L S '
	    STOP
          ENDIF
          IF(LSS.EQ.1) THEN
            CALL SUBLS(IT,ITS,L,S)
          ELSE
            CALL SUBLS(ITS,IT,L,S)
            IE1=LQ+LL+LS-LQS-LLS-LSS+2*L
            IF((IE1/4)*4.NE.IE1)S=-S
          ENDIF
	ELSE
          CALL FRMA(IT,LQ,LL,LS,ITS,LQS,LLS,LSS,S)
	ENDIF 
      ELSE
        CALL RMEALS(L,IT,LQ,LL,LS,ITS,LQS,LLS,LSS,S)
      ENDIF
      RETURN
      END
