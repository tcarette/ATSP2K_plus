*     ------------------------------------------------------------------
*              T K
*     ------------------------------------------------------------------
*

      DOUBLE PRECISION FUNCTION  TK(I,II,J,JJ,K)
      PARAMETER (NOD=220)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
c
      Call DZK(II,JJ,k)
      CALL YKK(II,JJ,K,1)
c
      DEN = L(I)+ L(J)+ 2
      FACT = L(J)
      D =  FACT*P(3,I)*P(3,J)*YK(3)/DEN
c
      MX = MIN0(MAX(I),MAX(J)) - 1
      Do M =3,MX
       S = (-P(M+2,J) + D8*(P(M+1,J)-P(M-1,J)) + P(M-2,J))/(D6*H)
       YK(M) = D5*P(M,I)*(S - P(M,J))*YK(M)
      End do
c
      S1=D0
      S2=D0
      Do M = 4,MX,2
       S1 = S1 + YK(m)
       S2 = S2 + YK(m+1)
      End do
c
      TK = D + H1*(S2 + D2*S1 + D5*YK(3))
      TK = TK / (k+k+1) * FINE
      RETURN
      END
