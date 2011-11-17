!     ------------------------------------------------------------------
!            Y K K
!     ------------------------------------------------------------------
!
!       Stores in YK-array the values of the YK integral
!
      SUBROUTINE YKK(I, J, K, KK) 
      use param_C
      use radial_C
      use nel_C
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!...Translated by Pacific-Sierra Research 77to90  4.3E  00:48:14  11/16/01  
!...Switches:                     
!      PARAMETER(NOD=220)
!      COMMON/PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
!     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
!      COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
!      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
!     :       (IQMAX,MAX_(1))
!      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
!
      B = EH**K 
      A = B*EH**KK 
      AA = A*A 
      A = D4*A 
      C = 2*K + KK 
      HH = C*H3 
      F2 = YK(ND)*B 
      F1 = F2*B 
      DO MM = 3, NO 
         M = NO - MM + 1 
         F3 = YK(M) 
         YK(M) = YK(M+2)*AA + HH*(F3 + A*F2 + AA*F1) 
         F1 = F2 
         F2 = F3 
      END DO 
      RETURN  
      END SUBROUTINE YKK 
