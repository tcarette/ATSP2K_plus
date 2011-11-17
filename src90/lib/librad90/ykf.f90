!
!     ------------------------------------------------------------------
!              Y K F
!     ------------------------------------------------------------------
!
!               k
!       Stores Y (i, j; r) in the array YK
!
!
      SUBROUTINE YKF(I, J, K, REL) 
      use param_C
      use radial_C
      use nel_C
!
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!...Translated by Pacific-Sierra Research 77to90  4.3E  00:46:39  11/16/01  
!...Switches:                     
!      PARAMETER(NOD=220)
!      COMMON/PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
!     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
!      COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
!      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
!     :       (IQMAX,MAX_(1))
!      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
      DIMENSION :: F(NOD) 
      LOGICAL :: REL 
!
      CALL ZK (I, J, K) 
      A = EH**(K + 1) 
      C = 2*K + 1 
      A2 = A*A 
      H90 = C*H3/D30 
      A3 = A2*A*H90 
      AI = H90/A 
      AN = 114.D0*A*H90 
      A34 = 34.D0*H90 
      DO M = ND - 1, 2, -1 
         F(M) = AN*YK(M+1) + A34*(YK(M)+A2*YK(M+2)) - YK(M-1)*AI - YK(M+3)*A3 
      END DO 
      F(1) = C*H3*(YK(1)+D4*A*YK(2)+A2*YK(3)) 
      DO M = ND - 1, 1, -1 
         YK(M) = YK(M+2)*A2 + F(M) 
      END DO 
      IF (.NOT.REL) RETURN  
      MM = MAX0(MAX_(I),MAX_(J)) 
      C = C*FINE 
      DO M = 1, MM 
         YK(M) = YK(M) + C*P(M,I)*P(M,J) 
      END DO 
      RETURN  
      END SUBROUTINE YKF 
