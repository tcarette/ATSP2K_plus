*
*     ------------------------------------------------------------
*        Z Z
*     ------------------------------------------------------------
*
      Double precision FUNCTION ZZ(i1,i2,i3,i4,K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Parameter (NOD=220)
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
      C= 2*k*(k+1)
      ZZ = - C * (TK(I1,I2,I3,I4,K+1) - TK(I1,I2,I3,I4,K-1))
*
      C1= L(i1)*(L(i1)+1) - L(i3)*(L(i3)+1) - k*(k+1)
      ZZ = ZZ - C1 * (UK(I1,I2,I3,I4,K+1) - UK(I1,I2,I3,I4,K-1))
*
      C2= L(i2)*(L(i2)+1) - L(i4)*(L(i4)+1) - k*(k+1)
      ZZ = ZZ - C2 * (UK(I2,I1,I4,I3,K+1) - UK(I2,I1,I4,I3,K-1))
*
      C= C1*C2/2
      C1= C*(k-2)/k/(k+k-1)
      C2= C*(k+3)/(k+1)/(k+k+3)
      ZZ = ZZ - C1 * (SN(I1,I2,I3,I4,K-2) + SN(I2,I1,I4,I3,K-2) )
     :        + C2 * (SN(I1,I2,I3,I4,K)   + SN(I2,I1,I4,I3,K) )
      RETURN
      END
