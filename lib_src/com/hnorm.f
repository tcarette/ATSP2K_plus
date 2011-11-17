*     ------------------------------------------------------------------ 
*              H N O R M 
*     ------------------------------------------------------------------ 
* 
*       Returns the value of the normalization constant for an (nl) 
*   hydrogenic function with nuclear charge ZZ. 
* 
* 
      DOUBLE PRECISION FUNCTION HNORM(N,L,ZZ) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID 
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS 
      M = L + L + 1 
      A = N + L 
      B = M 
      T = A 
      D = B 
      M = M - 1 
      IF (M .EQ. 0) GO TO 2 
      DO 1 I = 1,M 
      A = A - D1 
      B = B - D1 
      T = T*A 
1     D = D*B 
2     HNORM = DSQRT(ZZ*T)/( N*D) 
      RETURN 
      END 
