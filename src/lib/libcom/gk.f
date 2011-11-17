*
*     ------------------------------------------------------------------
*              G K
*     ------------------------------------------------------------------
*                             k
*       Returns the value of G (i,j).
*
*
      DOUBLE PRECISION FUNCTION GK(I,J,K,REL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL REL
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      CALL YKF(I,J,K,REL)
      GK = QUADS(I,J,1)
      IF (MASS .GT. 0) THEN
         IF (MASS .EQ. 1) THEN
            IF (K .EQ. 1) GK = GK + RMASS*GRAD(I,J)**2
         ELSE
            GK = GK*(D1 + RMASS/D2)
            IF (K .EQ. 1) GK = GK + Z*RMASS*QUADR(I,J,1)*QUADR(J,I,-2)
         END IF
      END IF
      RETURN
      END
