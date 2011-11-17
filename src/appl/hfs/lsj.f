*-----------------------------------------------------------------------
*     L S J
*-----------------------------------------------------------------------
*
*     This subroutine determines 2L, 2S, 2Jmax, and 2Jmin for the
*     configuration with the order number (index) NR.

      SUBROUTINE LSJ(QNOC,QJ1,NR,LL,SS,JJMAX,JJMIN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER SS
      POINTER(QNOC,NOCCSH(1)),(QJ1,J1QNRD(15,1))
      MI=2*NOCCSH(NR)-1
      MQ=J1QNRD(MI,NR)/64
      LL=MOD(MQ,64)-1
      SS=MQ/64-1
      JJMIN=IABS(LL-SS)
      JJMAX=LL+SS
      RETURN
      END
