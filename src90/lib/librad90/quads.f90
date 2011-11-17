!
!     ------------------------------------------------------------------
!              Q U A D S
!     ------------------------------------------------------------------
!
!
!                                       kk
!       Evaluates the integral of  (1/r)   YK(r) P (r) P (r)  with
!                                                 i     j
!   respect to r.
!
!
      DOUBLE PRECISION FUNCTION QUADS (I, J, KK) 
      use param_C
      use radial_C
      use nel_C
!
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!...Translated by Pacific-Sierra Research 77to90  4.3E  00:41:52  11/16/01  
!...Switches:                     
!      PARAMETER(NOD=220)
!      COMMON/PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
!     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
!      COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
!      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
!     :       (IQMAX,MAX_(1))
!      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
      DIMENSION :: G(NOD) 
!
      DEN = L(I) + L(J) + 3 
      K = 2 - KK 
      CD = D1 + Z*R(1)*(DEN - D1)/((DEN + D1)*((L(I)+1)*(L(J)+1))) 
      D = YK(1)*P(1,I)*P(1,J)*R(1)**K*(CD/(DEN*H1) + D5) 
      DD = D0 
      MX = MIN0(MAX_(I),MAX_(J)) 
      DO M = 2, MX 
         G(M) = YK(M)*P(M,I)*P(M,J)*R(M)**K 
      END DO 
      DO M = 2, MX - 1, 2 
         DD = DD + G(M) 
         D = D + G(M+1) 
      END DO 
      QUADS = H1*(D + D2*DD) 
      RETURN  
      END FUNCTION QUADS 
