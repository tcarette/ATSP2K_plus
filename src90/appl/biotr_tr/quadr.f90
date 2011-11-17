!
!     ------------------------------------------------------------------
!              Q U A D R
!     ------------------------------------------------------------------
!
!                                   kk
!       Evaluates the integral of  r   P (r) P (r) with respect to r
!                                       i     j
!
      DOUBLE PRECISION FUNCTION QUADR (I, J, KK) 
      use param_C
      use radial_C
      use nel_C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!...Translated by Pacific-Sierra Research 77to90  4.3E  00:41:09  11/16/01  
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
      K = KK + 2 
      LI = L(I) 
      LJ = L(J) 
      DEN = LI + LJ + 1 + K 
      ZR = Z*R(4) 
      BI = (P(4,I)/(AZ(I)*R2(4)*R(4)**LI)-D1+ZR/(LI+1))/ZR**2 
      BJ = (P(4,J)/(AZ(J)*R2(4)*R(4)**LJ)-D1+ZR/(LJ+1))/ZR**2 
      ALPHA = (D1/(LI + 1) + D1/(LJ + 1))/(DEN + D1) 
      ZR = Z*R(1) 
      BETA = (DEN + D1)*ALPHA**2 - D2*(BI + BJ + D1/((LI + 1)*(LJ + 1)))/(DEN&
          + D2) 
      D = P(1,I)*P(1,J)*R(1)**K*(((BETA*ZR + ALPHA)*ZR + D1)/(DEN*H1) + D5) 
      DD = D0 
      M = MIN0(MAX_(I),MAX_(J)) 
      DO JJ = 1, M 
         G(JJ) = P(JJ,I)*P(JJ,J)*R(JJ)**K 
      END DO 
      DO JJ = 3, M, 2 
         D = D + G(JJ) 
      END DO 
      DO JJ = 2, M, 2 
         DD = DD + G(JJ) 
      END DO 
      QUADR = H1*(D + D2*DD) 
      RETURN  
      END FUNCTION QUADR 
