* 
*     ------------------------------------------------------------------ 
*              H W F 
*     ------------------------------------------------------------------ 
* 
*       Returns the value of an unnormalized (nl) hydrogenic function 
*   with nuclear charge ZZ and radius r. 
* 
* 
      DOUBLE PRECISION FUNCTION HWF(N,L,ZZ,R) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID 
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS 
      K = N-L-1 
      P = D1 
      A = D1 
      B = K 
      C = N+ L 
      X = -D2*ZZ*R/N 
* 
*  *****  TEST IF UNDERFLOW MAY OCCUR, IF SO SET HWF = 0 
* 
      IF ( X .LT. -150.D0 ) GO TO 5 
      IF (K) 1,2,3 
3     DO 4 I = 1,K 
      P = D1 + A/B*P/C*X 
      A = A + D1 
      B = B - D1 
4     C = C - D1 
2     HWF = P*DEXP(X/D2)*(-X)**(L+1) 
      RETURN 
1     WRITE(6,7) N,L,ZZ,R 
7     FORMAT(51H FORBIDDEN COMBINATION OF N AND L IN HWF SUBPROGRAM/ 
     :    4H N = ,I4,6H   L = ,I4,6H   Z = ,F6.1,6H   R = ,F8.4) 
      STOP 
5     HWF = D0 
      RETURN 
      END 
