!
!     ------------------------------------------------------------
!        Z Z
!     ------------------------------------------------------------
!
      REAL(KIND(0.0D0)) FUNCTION ZZ (I1, I2, I3, I4, K) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE1=>DOUBLE 
      use nel_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  00:56:55  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE tk_I 
      USE uk_I 
      USE sn_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I1 
      INTEGER  :: I2 
      INTEGER  :: I3 
      INTEGER  :: I4 
      INTEGER  :: K 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NWD = 128 
      INTEGER, PARAMETER :: NOD = 220 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE1) :: C, C1, C2 
!-----------------------------------------------
 
!      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
!     :       (IQMAX,MAX_(1))
!      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
      C = 2*K*(K + 1) 
      ZZ = -C*(TK(I1,I2,I3,I4,K + 1) - TK(I1,I2,I3,I4,K - 1)) 
!
      C1 = L(I1)*(L(I1)+1) - L(I3)*(L(I3)+1) - K*(K + 1) 
      ZZ = ZZ - C1*(UK(I1,I2,I3,I4,K + 1) - UK(I1,I2,I3,I4,K - 1)) 
!
      C2 = L(I2)*(L(I2)+1) - L(I4)*(L(I4)+1) - K*(K + 1) 
      ZZ = ZZ - C2*(UK(I2,I1,I4,I3,K + 1) - UK(I2,I1,I4,I3,K - 1)) 
!
      C = C1*C2/2 
      C1 = C*(K - 2)/K/(K + K - 1) 
      C2 = C*(K + 3)/(K + 1)/(K + K + 3) 
      ZZ = ZZ - C1*(SN(I1,I2,I3,I4,K - 2) + SN(I2,I1,I4,I3,K - 2)) + C2*(SN(I1,&
         I2,I3,I4,K) + SN(I2,I1,I4,I3,K)) 
      RETURN  
      END FUNCTION ZZ 
