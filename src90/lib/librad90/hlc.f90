!
!    --------------------------------------------------------------
!         H L C
!    --------------------------------------------------------------
!
      DOUBLE PRECISION FUNCTION HLC (EL, I, J, REL) 
      use param_C
      use radial_C
      use nel_C, EL_=>EL
!
! *** COMPUTES HL(I,J) MODIFIED BY THE INTERACTIONS WITH THE CLOSED
!     SHELL
!
!
! *** COMPUTES HL(I,J) MODIFIED BY THE INTERACTIONS WITH THE CLOSED
!     SHELL
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!...Translated by Pacific-Sierra Research 77to90  4.3E  00:39:53  11/16/01  
!...Switches:                     
!      PARAMETER(NOD=220)
!      COMMON/PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MSOO,NCFG,IB,IC,ID
!     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
!      COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
!      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
!     :       (IQMAX,MAX_(1))
!      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
!
      CHARACTER*1 :: EL(*)*3 
      LOGICAL :: REL 
      HLC = HL(EL,I,J,REL) 
!
      LI = L(I) 
      DO IP = 1, NCLOSD 
         LP = L(IP) 
         SUMIP = 4*LP + 2 
         T = RK(I,IP,J,IP,0,REL)                 ! direct interaction 
         K1 = IABS(LI - LP) 
         K2 = LI + LP 
         DO K = K1, K2, 2 
            CB = ZCB(LI,K,LP)/2 
!                                            ! exchange interaction
            T = T - CB*RK(I,IP,IP,J,K,REL) 
!                                            ! o-o - interaction
            IF (MSOO.GT.2 .AND. LI.GT.0 .AND. LP.GT.0) THEN 
               CN = CB*(LI + LP + K + 2)*(LI + LP - K)*(LI - LP + K + 1)*(LP - &
                  LI + K + 1) 
               CN = CN/((K + 1)*(K + 2)) 
               IF (ABS(CN) .GT. 0.0) T = T + CN*(SN(I,IP,IP,J,K) + SN(IP,I,J,IP&
                  ,K)) 
            ENDIF 
!
         END DO 
!
         HLC = HLC - 2*SUMIP*T 
!
      END DO 
      RETURN  
      END FUNCTION HLC 
