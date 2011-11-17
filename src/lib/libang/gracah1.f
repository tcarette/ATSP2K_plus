*
*     -----------------------------------------------------------------
*      G R A C A H 1
*     -----------------------------------------------------------------
*
      SUBROUTINE GRACAH1(I,J,K,L,M,N,RAC) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
* 
*      SUBROUTINE TO CALCULATE RACAH COEFFICIENTS. 
*      THE ARGUMENTS I,J,K,L,M,N SHOULD BE TWICE THEIR ACTUAL VALUE. 
*      WRITTEN BY N. S. SCOTT 
*      Modified by C. Froese Fischer, March 11, 1988 to use
*         table look-up
*      and
*                  G. Gaigalas, September 14, 1997
*
      COMMON/FACT/GAM(100) 
      DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/ 
* 
      J1 = I+J+M 
      J2 = K+L+M 
      J3 = I+K+N 
      J4 = J+L+N 
      IF (MOD(J1,2) .EQ. 0  .AND.  MOD(J2,2) .EQ. 0   .AND. 
     :    MOD(J3,2) .EQ. 0  .AND.  MOD(J4,2) .EQ. 0 )  THEN 
          J1 = J1/2 
          J2 = J2/2 
          J3 = J3/2 
          J4 = J4/2 
          IF (MAX(I,J,M) .LE. J1 .AND.  MAX(K,L,M) .LE. J2  .AND. 
     :        MAX(I,K,N) .LE. J3 .AND.  MAX(J,L,N) .LE. J4  )  THEN 
              J5 = (I+J+K+L)/2 
              J6 = (I+L+M+N)/2 
              J7 = (J+K+M+N)/2 
              NUMIN = MAX(J1, J2, J3, J4) + 1 
              NUMAX = MIN(J5, J6, J7)     + 1 
              RAC = ONE 
              ICOUNT = 0 
              DO 10 KK = NUMIN+1,NUMAX 
                KI = NUMAX - ICOUNT 
                RAC = ONE - (RAC*(KI*(J5-KI+2)*(J6-KI+2)*(J7-KI+2)))/ 
     :                   ((KI-1-J1)*(KI-1-J2)*(KI-1-J3)*(KI-1-J4)) 
                ICOUNT = ICOUNT+1 
  10          CONTINUE 
              RAC = RAC*EXP(
     :              (GAM(NUMIN+1) - GAM(NUMIN-J1) - GAM(NUMIN-J2) - 
     :               GAM(NUMIN-J3) - GAM(NUMIN-J4) - GAM(J5+2-NUMIN)- 
     :               GAM(J6+2-NUMIN)-GAM(J7+2-NUMIN)) + 
     :              (GAM(J1+1-I)+GAM(J1+1-J)+GAM(J1+1-M)-GAM(J1+2) + 
     :               GAM(J2+1-K)+GAM(J2+1-L)+GAM(J2+1-M)-GAM(J2+2) + 
     :               GAM(J3+1-I)+GAM(J3+1-K)+GAM(J3+1-N)-GAM(J3+2) + 
     :               GAM(J4+1-J)+GAM(J4+1-L)+GAM(J4+1-N)-GAM(J4+2))/TWO) 
	      IF (MOD(J5+NUMIN,2) .EQ. 0) RAC = -RAC
	  ELSE
	      RAC = ZERO
	  END IF
      ELSE 
         RAC = ZERO 
      END IF 
      RETURN 
      END
