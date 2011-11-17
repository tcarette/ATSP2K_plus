*     ------------------------------------------------------------------
*            Y K K
*     ------------------------------------------------------------------
*
*       Stores in YK-array the values of the YK integral
*
      SUBROUTINE YKK(I,J,k,kk)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
c
      B = EH**k
      A = B*EH**kk
      AA = A*A
      A = D4*A
      C = 2*k+kk
      HH = C*H3
      F2 = YK(ND)*B
      F1 = F2*B
      DO 9 MM = 3,NO
      M = NO -MM+1
      F3 =YK(M)
      YK(M) = YK(M+2)*AA + HH*(F3 +A*F2 + AA*F1)
      F1 = F2
9     F2 = F3
      RETURN
      END
