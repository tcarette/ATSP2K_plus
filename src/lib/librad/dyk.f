*
*     ------------------------------------------------------------------
*              D Y K
*     ------------------------------------------------------------------
*
*       Stores in YK the values of the integral of
*              k
*       P (s/r) (dP /ds - P /s) integrated over the interval (0,r)
*        i         j       j
*
*   which enter into the spin-orbit calculation.
*
*
      SUBROUTINE DYK(I,J,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
      DIMENSION F(NOD),G(NOD)
*
      DEN = L(I)+ L(J)+ 2 + K
      FACT = L(J)
      C = FACT/DEN
      DO 1 JJ =1,2
1     YK(JJ) = C*P(JJ,I)*P(JJ,J)*R(JJ)
      A = EH**K
      AA = A*A
      A = D4*A
      C = D1/(D6*H)
      HH = D5*H3
      F(1) = D2*DEN*YK(1)
      F(2) = D2*DEN*YK(2)
      DO 8 M =3,ND
        G(M) = C*(-P(M+2,J) + D8*(P(M+1,J)-P(M-1,J)) + P(M-2,J))
8     CONTINUE
      DO 9 M = 3,ND
        F(M) = P(M,I)*(G(M) - P(M,J))*R(M)
9     CONTINUE
      DO 10 M = 3,ND
        G(M) = HH*(AA*F(M-2) + A*F(M-1) + F(M))
10    CONTINUE
      DO 11 M = 3,ND
        YK(M) = YK(M-2)*AA + G(M)
11    CONTINUE
      A = A*(EH)**3
      AA = A*A/D16
      C = 2*K+3
      HH = C*H3
      YK(NO) = YK(ND)
      YK(NO-1) = YK(ND)
      DO 12 M = NO-2,1,-1
        G(M) = HH*(YK(M) + A*YK(M+1) + AA*YK(M+2))
12    CONTINUE
      DO 13 M = NO-2,1,-1
        YK(M) = YK(M+2)*AA + G(M)
13    CONTINUE
      RETURN
      END
