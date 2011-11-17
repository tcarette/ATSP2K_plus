*
*     ------------------------------------------------------------------
*       R M E T R
*     ------------------------------------------------------------------
*            
* --- evaluates the angular part of the one-electron transition reduced 
*     matrix element. See equations (4) and (5) of paper II.
*
      DOUBLE PRECISION FUNCTION RMETR(L1,L2)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL REL,VOK
      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      IF(IFL.EQ.1) THEN
*       electric multipole transitions
        RMETR=RME(L1,L2,LAM)
	IF(MOD(L2-L1+LAM,4).NE.0) RMETR=-RMETR
      ELSE
*       magnetic multipole transitions
        RMETR=RME(L1,L2,LAM-1)
        IF(DABS(RMETR).GT.EPS) THEN
	  IF(MOD(L2-L1+LAM-1,4).NE.0) RMETR=-RMETR
          IF(IFL.EQ.3) THEN
            LL2=L2+L2
            L3=LAM+LAM
	    CALL SIXJ(2,LL2,LL2,L1+L1,L3,L3-2,1,W)
            RMETR=-RMETR*W*SQRT(DBLE(L2*(L2+1)*(LL2+1)*(L3+1)))
          ELSE
            RMETR=RMETR*SQRT(HALF*THREE)
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
