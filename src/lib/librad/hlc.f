*
*    --------------------------------------------------------------
*         H L C
*    --------------------------------------------------------------
*
      DOUBLE PRECISION FUNCTION HLC(EL,I,J,REL)
*
* *** COMPUTES HL(I,J) MODIFIED BY THE INTERACTIONS WITH THE CLOSED
*     SHELL
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MSOO,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
*
      CHARACTER EL(*)*3
      LOGICAL REL
      HLC = HL(EL,I,J,REL)
*
      LI=L(i)
      DO 10 IP = 1,NCLOSD
       LP=L(IP)
       SUMIP = 4*LP+2
       T = RK(I,IP,J,IP,0,REL)            ! direct interaction
       K1 = IABS(LI-LP)
       K2 = LI+LP
       DO 20 K = K1,K2,2
        CB=ZCB(LI,K,LP)/2
*                                            ! exchange interaction
           T = T - CB*RK(I,IP,IP,J,K,REL)
*                                            ! o-o - interaction
       if(MSOO.GT.2.and.LI.gt.0.and.LP.gt.0) then
         CN = CB * (LI+LP+k+2)*(LI+LP-k)*(LI-LP+k+1)*(LP-LI+k+1)
         CN = CN / ((k+1)*(k+2))
         if(abs(CN).gt.0.0)
     :     T = T + CN*(SN(I,IP,IP,J,K) + SN(IP,I,J,IP,K))
       end if
*
   20 Continue    ! over k
*
       HLC = HLC - 2*SUMIP*T
*
   10 Continue     ! over NCLOSD
      Return
      End
