!     ------------------------------------------------------------------
!              T K
!     ------------------------------------------------------------------
!
 
      DOUBLE PRECISION FUNCTION TK (I, II, J, JJ, K) 
      use param_C
      use radial_C
      use nel_C
!
!      PARAMETER (NOD=220)
!
!      PARAMETER (NOD=220)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!...Translated by Pacific-Sierra Research 77to90  4.3E  00:45:17  11/16/01  
!...Switches:                     
!      COMMON/PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
!     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
!      COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
!      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
!     :       (IQMAX,MAX_(1))
!      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
!
      CALL DZK (II, JJ, K) 
      CALL YKK (II, JJ, K, 1) 
!
      DEN = L(I) + L(J) + 2 
      FACT = L(J) 
      D = FACT*P(3,I)*P(3,J)*YK(3)/DEN 
!
      MX = MIN0(MAX_(I),MAX_(J)) - 1 
      DO M = 3, MX 
         S = ((-P(M+2,J))+D8*(P(M+1,J)-P(M-1,J))+P(M-2,J))/(D6*H) 
         YK(M) = D5*P(M,I)*(S - P(M,J))*YK(M) 
      END DO 
!
      S1 = D0 
      S2 = D0 
      DO M = 4, MX, 2 
         S1 = S1 + YK(M) 
         S2 = S2 + YK(M+1) 
      END DO 
!
      TK = D + H1*(S2 + D2*S1 + D5*YK(3)) 
      TK = TK/(K + K + 1)*FINE 
      RETURN  
      END FUNCTION TK 
