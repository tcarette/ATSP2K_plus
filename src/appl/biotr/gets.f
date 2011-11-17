*
*     ------------------------------------------------------------------
*	G E T S    
*     ------------------------------------------------------------------
*
      SUBROUTINE GETS(S,NSHLI,NSHLF)
*
* Obtain overlap matrix between all I shells and
* all F shells
*
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      double precision SCR(NSHLI*NSHLF);
      double precision S(NSHLI,NSHLF)
       
      
      !S = reshape(SCR,(/NSHLI,NSHLF/))
*
      DO 100 II = 1, NSHLI
       DO 50 IF = 1, NSHLF
        S(II,IF) = QUADR(II,IF+NSHLI,0)
   50  CONTINUE
  100 CONTINUE
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
       WRITE(6,*) ' S matrix from GETS' 
       CALL WRTMAT(S,NSHLI,NSHLF,NSHLI,NSHLF)
      END IF
*
      RETURN
      END
