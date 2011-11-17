*    ------------------------------------------------------------------
*     E C O R E
*    ------------------------------------------------------------------
*
*
      SUBROUTINE ECORE(EL,EC,REL)
*
* *** COMPUTES THE ENERGY OF THE COMMON CLOSED SHELLS
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MSOO,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
      CHARACTER EL(*)*3
      LOGICAL REL
*
      EC = 0.D0
      DO 10 I = 1,NCLOSD
       SUMI = 4*L(I)+2
       TI = RK(I,I,I,I,0,REL)
       DO 20 K = 2,2*L(I),2
        CA = ZCB(L(i),K,L(i)) * (2*L(I)+1)/(4*L(i)+1)
        TI = TI - CA*RK(I,I,I,I,K,REL)
   20  Continue
*                                                 o-o contribution
       IF(MSOO.GT.2) THEN
        if(L(I).eq.1) TI = TI +  8.D0/5.D0 *SN(I,I,I,I,0)
        if(L(I).eq.2) TI = TI + 56.D0/21.D0*SN(I,I,I,I,0)
     :                        + 16.D0/21.D0*SN(I,I,I,I,2)
        if(L(I).eq.3) TI = TI + 48.D0/13.D0*SN(I,I,I,I,0)
     :                        + 16.D0/13.D0*SN(I,I,I,I,2)
     :                        + 80.D0/143.D0*SN(I,I,I,I,4)
       end if
*
       EC = EC + SUMI*((SUMI-1)*TI - HL(EL,I,I,REL))/D2
*
       DO 30 J = 1,I-1
         SUMJ = 4*L(J)+2
         TIJ = RK(I,J,I,J,0,REL)
         DO 40 K=IABS(L(I)-L(J)),L(I)+L(J),2
            CB = ZCB(L(i),K,L(j))/2
            TIJ = TIJ -CB*RK(I,J,J,I,K,REL)
*
         if(MSOO.GT.2.and.L(i).gt.0.and.L(j).gt.0) then
          CN = CB * (L(i)+L(j)+k+2)*(L(i)+L(j)-k)
          CN = CN * (L(i)-L(j)+k+1)*(L(j)-L(i)+k+1)
          CN = CN / ((k+1)*(k+2))
          if(abs(CN).gt.0.0) TIJ = TIJ + CN*2*SN(I,J,J,I,K)
         end if
*
  40     Continue
*
       EC = EC + SUMI*SUMJ*TIJ
*
  30   Continue
  10  Continue
*
      END

