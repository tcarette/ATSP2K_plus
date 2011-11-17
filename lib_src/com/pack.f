* 
* 
*     ------------------------------------------------------------------ 
*	P A C K 
*     ------------------------------------------------------------------ 
*  Subroutine written by Bin LIU 
* 
*  Rules for encoding 
*  1. All blanks deleted 
*  2. If Qi=1, omit Qi 
*  3. If Qi=1 or Qi>=4l+1, omit ALFAi 
*  4. If i=1 or (Qi=4l+2 and i<>m), insert '.'; else _BETAi. 
* 
      SUBROUTINE PACK (M, EL, Q, COUPLE, STR) 
 
	  INTEGER FULL,Q(5),CONST 
	  CHARACTER*3 EL(5),COUPLE(9),CH3
	  CHARACTER   CH1
	  CHARACTER*66 STR
* 
*  FULL   :  4l+2 
*  CONST  :   constant for converting lowercase to uppercase 
*  CH*    :  temporary variables 
* 
	  CONST = ICHAR('a') - ICHAR('A') 
	  STR=' ' 
	  J = 0
* 
*   -----  begin to encode  ----- 
* 
	  DO 100 I=1,M 
	    N = Q(I)
	    IF (N .EQ. 0) GO TO 100
	    K = 3
	    IF (EL(I)(3:3) .EQ. ' ') K = 2
	    IF (EL(I)(1:1) .EQ. ' ') THEN
		EL(I)=EL(I)(2:3)//' ' 
		K = 2
	    END IF
	    CH1=EL(I)(2:2) 
	    IF ((CH1.GE.'A') .AND. (CH1.LE.'Z')) 
     :            EL(I)(2:2)=CHAR(ICHAR(CH1)+CONST) 
	    FULL=4*LVAL(CH1)+2 
* 
*  -----  convert Qi into character  ----- 
* 
	    WRITE(CH3,'(I2)') Q(I) 
	    STR=STR(1:J)//EL(I)(1:K)
	    J= J + K
* 
*  -----  If Qi<>1, add Qi 
*           If Qi<4l+1, add TERMi for the shell ----- 
* 
	    IF (N .NE. 1) THEN 
 		IF (N .GT. 9) THEN 
	          STR=STR(1:J)//'('//CH3(1:2)//')' 
		  J = J + 4
	        ELSE 
	          STR=STR(1:J)//'('//CH3(2:2)//')' 
		  J = J + 3
	        ENDIF 
	        IF (N .LT. FULL-1 .AND. M .NE. 1) THEN 
		    CH3=COUPLE(I) 
		    CH1= CH3(2:2) 
		    IF (CH1.GE.'a' .AND. CH1.LE.'z') 
     :			CH3(2:2) = CHAR(ICHAR(CH1)-CONST) 
	            STR=STR(1:J)//CH3 
	            J= J + 3
	        ENDIF 
	    ENDIF 
* 
*  -----  If i=1 or Qi=4l+2 and i<>m, 
*           insert '.'; else _RESULTANTi.  ----- 
* 
50          IF ((I.NE.1 .AND. N.NE.FULL )
     :           .OR. I.EQ.M  ) THEN 
		CH3=COUPLE(M+I-1) 
		CH1 = CH3(2:2) 
		IF (CH1.GE.'a' .AND. CH1.LE.'z') 
     :			CH3(2:2) = CHAR(ICHAR(CH1)-CONST) 
		K = 2 
		IF (M .EQ. 1) K = 3 
	        STR=STR(1:J)//'_'//CH3(1:K) 
		J = J + K + 1
	    ENDIF 
	    IF (I .NE. M .AND. N.NE.0 ) THEN 
		J = J + 1
	        STR(J:J)='.' 
	    ENDIF 
100	  CONTINUE 
*
*>>>>>    Because of a compiler error on the SUN, the following is
*         needed to have the string printed correctly.
	  STR = STR(1:J)
	  RETURN 
	END 
