!     ------------------------------------------------------------------
!              D Z K
!     ------------------------------------------------------------------
!
!       Stores in YK the values of the integral of
!              k
!       P (s/r) (dP /ds - P /s) integrated over the interval (0,r)
!        i         j       j
!
!   which enter into the spin-orbit calculation.
!
!
      SUBROUTINE DZK(I, J, K) 
      use param_C
      use radial_C
      use nel_C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!...Translated by Pacific-Sierra Research 77to90  4.3E  23:30:36  11/15/01  
!...Switches:                     
!...Translated by Pacific-Sierra Research 77to90  4.3E  23:28:47  11/15/01
!...Switches:
!      PARAMETER(NOD=220)
!      COMMON/PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
!     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
!      COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
!      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
!     :       (IQMAX,MAX_(1))
!      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
!
      LI = L(I) + 1 
      LJ = L(J) 
      ZR = Z*R(4) 
      BJ = (P(4,J)/(AZ(J)*R2(4)*R(4)**LJ)-D1+ZR/(LJ+1))/ZR**2 
!
      DEN = LI + LJ + K + 1 
      DI = D1/LI 
      IF (LJ .GT. 0) DJ = D1/LJ 
      IF (LJ .EQ. 0) DJ = D2*BJ 
      FACT = (DI + DJ)/(DEN + D1) 
      ZR = Z*R(1) 
      F1 = AZ(J)*R(1)**LJ*(LJ - ZR + BJ*(LJ + 2)*ZR**2) 
      F1 = R(1)*R2(1)*P(1,I)*F1 
      YK(1) = F1*(D1 + ZR*FACT)/DEN 
      ZR = Z*R(2) 
      F2 = AZ(J)*R(2)**LJ*(LJ - ZR + BJ*(LJ + 2)*ZR**2) 
      F2 = R(2)*R2(2)*P(2,I)*F2 
      YK(2) = F2*(D1 + ZR*FACT)/DEN 
!
      A = EH**K 
      AA = A*A 
      A = D4*A 
      DO M = 3, ND 
         F3 = (((-P(M+2,J)))+D8*(P(M+1,J)-P(M-1,J))+P(M-2,J))/(D6*H) 
         F3 = D5*P(M,I)*(F3 - P(M,J))*R(M) 
         YK(M) = YK(M-2)*AA + H3*(F3 + A*F2 + AA*F1) 
         F1 = F2 
         F2 = F3 
      END DO 
      RETURN  
      END SUBROUTINE DZK 
 
