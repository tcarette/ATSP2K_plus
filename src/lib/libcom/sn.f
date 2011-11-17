*
*     ------------------------------------------------------------------
*              S N
*     ------------------------------------------------------------------
*
*                                      3              k
*       Evaluates the integral of (1/r)  P (r) P (r) Z (i, j; r)  with
*                                         i     j
*   respect to r.
*
      DOUBLE PRECISION FUNCTION SN(I,J,II,JJ,K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      CALL ZK(J,JJ,K)
      SN = QUADS(I,II,3)*FINE
      RETURN
      END
