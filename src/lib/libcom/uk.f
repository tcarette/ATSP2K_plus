C----------------------------------------------------------------------
C        U K
C---------------------------------------------------------------------


      DOUBLE PRECISION FUNCTION UK(I,J,II,JJ,K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
c
      CALL DZK(J,JJ,K)
      UK1 = QUADS(I,II,2)
c
      CALL YKK(J,JJ,K,1)
      UK2 = QUADS(I,II,2)
c
      UK = (UK1 - (k+2.)/(k+k+1.)*UK2) * fine  
      RETURN
      END
