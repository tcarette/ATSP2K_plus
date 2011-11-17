!
!     ---------------------------------------------------------------
!        R K
!     ---------------------------------------------------------------
!
!      MSOO     Mass     OO
!
!       0       No       No
!       1       Yes      No
!       2       Yes      No
!       3       No       Yes
!       4       Yes      Yes
!       5       Yes      Yes
!
      DOUBLE PRECISION FUNCTION RK (I1, I2, I3, I4, K, REL) 
      use param_C
      use radial_C
      use nel_C
!
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!...Translated by Pacific-Sierra Research 77to90  4.3E  00:43:00  11/16/01  
!...Switches:                     
!      Parameter (NOD=220)
      LOGICAL :: REL 
!      COMMON/PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MSOO,NCFG,IB,IC,ID
!     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
!      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
!     :       (IQMAX,MAX_(1))
!      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
 
      CALL YKF (I1, I3, K, REL) 
      RK = QUADS(I2,I4,1) 
      IF (MSOO .EQ. 0) RETURN  
      IF (MSOO.EQ.1 .OR. MSOO.EQ.4) THEN 
         IF (K .EQ. 1) RK = RK - RMASS*GRAD(I1,I3)*GRAD(I2,I4) 
      ELSE IF (MSOO.EQ.2 .OR. MSOO.EQ.5) THEN 
         RK = RK*(D1 + RMASS/D2) 
         IF (K .EQ. 1) RK = RK + Z*RMASS/D2*(QUADR(I1,I3,1)*QUADR(I2,I4,-2) + &
            QUADR(I1,I3,-2)*QUADR(I2,I4,1)) 
      ENDIF 
! ! o-o  interaction
      IF (K.EQ.0 .OR. MSOO.LT.3) RETURN  
      IF (I1.EQ.I2 .AND. I1.EQ.I3 .AND. I1.EQ.I4) RETURN  
      RK = RK + (ZZ(I1,I2,I3,I4,K) + ZZ(I2,I1,I4,I3,K) + ZZ(I3,I4,I1,I2,K) + ZZ&
         (I4,I3,I2,I1,K))/4.D0 
      RETURN  
      END FUNCTION RK 
