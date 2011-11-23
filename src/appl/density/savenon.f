*
*     -------------------------------------------------------------
*      S A V E N O N 
*     -------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville           September 1997   *
*
      SUBROUTINE SAVENON(I,A,KL,LA,LB,LC,LD,JJI,JJF,IPTR)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/HYPER/VHY(20),NHY,IRHY(20),ISHY(20),JRHY(20),JSHY(20),
     1 IMUHY(20),IMUPHY(20),JMUHY(20),JMUPHY(20),
     2 INUHY(20),INUPHY(20),JNUHY(20),JNUPHY(20)
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      COMMON/PAPIL/ IIRHO,IISIG
      NHY=NHY+1
      IRHY(NHY)=IJFUL(IIRHO)
      ISHY(NHY)=IJFUL(IISIG)
      VHY(NHY)=A
      RETURN
      END
