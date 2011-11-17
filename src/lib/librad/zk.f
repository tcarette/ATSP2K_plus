*
*     ------------------------------------------------------------------
*              Z K
*     ------------------------------------------------------------------
*
*               k
*       Stores Z (i, j; r) in the array YK.
*
*
      SUBROUTINE ZK(I,J,K)
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
      DEN = L(I) + L(J) + 3+ K
      FACT = (D1/(L(I)+1) + D1/(L(J)+1))/(DEN + D1)
      A = EH**K
      A2 = A*A
      H90 = H/90.D0
      A3 = A2*A*H90
      AI = H90/A
      AN = 114.D0*A*H90
      A34 = 34.D0*H90
      DO 10 M = 1,NO
      F(M) = RR(M)*P(M,I)*P(M,J)
   10 CONTINUE
      YK(1) = F(1)*(D1 + Z*R(1)*FACT)/DEN
      YK(2) = F(2)*(D1 + Z*R(2)*FACT)/DEN
      YK(3) = F(3)*(D1 + Z*R(3)*FACT)/DEN
*      YK(3) = YK(1)*A2 + H3*(F3 + D4*A*F2 + A2*F1)
      DO 9 M = 5,NO
       G(M) = AN*F(M-2) + A34*(F(M-1)+A2*F(M-3))-F(M)*AI-F(M-4)*A3
    9 CONTINUE
      DO 8 M = 5,NO
        YK(M-1) = YK(M-3)*A2 + G(M)
    8 CONTINUE
      YK(NO) = A*YK(NO-1)
      IF (IABS(I-J)  +  IABS(K) .NE. 0) GO TO 2
*
*  *****  FOR Y0(I,I) SET THE LIMIT TO 1 AND REMOVE OSCILLATIONS
*  *****  INTRODUCED BY THE USE OF SIMPSON'S RULE
*
      M1 = (NO/2)*2 - 1
      M2 = M1 - 1
      C1 = D1 - YK(M1)
      C2 = D1 - YK(M2)
      DO 3 M = 1,NO,2
        YK(M) = YK(M) + C1
3     CONTINUE
      DO 4 M = 2,NO,2
        YK(M) = YK(M) + C2
4     CONTINUE
2     RETURN
      END
