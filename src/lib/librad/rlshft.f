*
*     ------------------------------------------------------------------
*              S H I F T
*     ------------------------------------------------------------------
*
*
*       Computes the mass velocity  and one-body   Darwin term
*   corrections for the relativistic shift in the energy of the electron
*   including non-diagonal corrections
*
*
      DOUBLE PRECISION  FUNCTION RLSHFT(I1,I2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      PARAMETER(NOD=220)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
*
*  *****  FORM DD  -L(L+1)/RR|P(I)>
*
      FL = L(I1)
      C = (FL+D5)**2
      LL = L(I1) + 1
      L2 = 2*L(I1) + 1
      L3 = 2*L(I1) + 3
      ZZ = Z*Z
      HH = 180.D0*H*H
      MX = MAX0(MAX(I1), MAX(I2))
      DO 1 J = 2,MX
1     YK(J) = -D1/RR(J)
*
*  *****  FORM THE INTEGRAND
*
      I = I1
      A1 = D0
      B1 = D0
      DO 2 KK = 1,2
      B2 = B1
      B1 = YK(4)
      YY = (P(3,I)+P(1,I) - D2*P(2,I))/(H*H) - C*P(2,I)
      YK(2) = YY*YK(2)
      YK(3) = YK(3)*((-(P(5,I)+P(1,I)) + D16*(P(4,I)+P(2,I))
     :        -D30*P(3,I))/(D12*H*H) - C*P(3,I))
      MM = MAX(I) - 3
      DO 11 K =  4,MM
      YY = D2*(P(K+3,I)+P(K-3,I))
      YY = YY - 27.D0*(P(K+2,I)+P(K-2,I))
      YY = YY +  270.D0*(P(K+1,I)+P(K-1,I)) - 490.D0*P(K,I)
      YY = YY/HH - C*P(K,I)
      YK(K) = YY*YK(K)
11    CONTINUE
      B1 = (YK(4)/(B1*(D2*Z*P(4,I)*R(4))) + D1)/R(4)
      MM = MM + 1
      YK(MM) = YK(MM)*((-(P(MM+2,I)+P(MM-2,I))+D16*(P(MM+1,I)+P(MM-1,I))
     :      -D30*P(MM,I))/(D12*H*H) - C*P(MM,I))
      MM = MM + 1
      YK(MM) = YK(MM)*((P(MM+1,I) + P(MM-1,I) - D2*P(MM,I))/(H*H)
     :         - C*P(MM,I))
      A2 = A1
      A1 = (P(1,I)/(AZ(I)*R(1)**L(I)*R2(1)) - D1 +Z*R(1)/LL)/RR(1)
      I = I2
2     CONTINUE
*
*  ***** DETERMINE CONTRIBUTION FROM NEAR THE NUCLEUS
*
      A = (Z/LL - L2*(B1 + B2)/D2)/LL
      B = (L2*B1*B2 - D2*(A1 + A2) + (Z/LL**2)*(D2*Z*(D1 + D1/LL)
     :     - L2*(B1 + B2)))/L3
      RELSH = -P(4,I1)*P(4,I2)*(D1 + A*R(4) + B*RR(4))*D4*ZZ/L2
      RELSH = RELSH/H1 - D5*YK(4)
      RELSH2 = D0
*
*  *****  INTEGRATE
*
      DO 50 J = 5,MX,2
      RELSH2 = RELSH2 + YK(J)
50    RELSH = RELSH  + YK(J-1)
      RELSH  = (RELSH + D2*RELSH2)*H1
      IF ( L(I1) .EQ. 0 ) RELSH = RELSH + Z*AZ(I1)*AZ(I2)
      RLSHFT = RELSH*D5*FINE
      RETURN
      END
