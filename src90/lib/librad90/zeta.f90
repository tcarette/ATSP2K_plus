!
!     -----------------------------------------------------------------
!       Z E T A
!     -----------------------------------------------------------------
!
!
      DOUBLE PRECISION FUNCTION ZETA (I1, I2) 
      use blume_C
      use param_C
      use radial_C
      use nel_C
!
!
!  ***** COMPUTES THE NUCLEAR SPIN-ORBIT PARAMETER AND THE
!        CORRECTIONS FOR THE COMMON CLOSED SHELLS
!
!
!
!  ***** COMPUTES THE NUCLEAR SPIN-ORBIT PARAMETER AND THE
!        CORRECTIONS FOR THE COMMON CLOSED SHELLS
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!...Translated by Pacific-Sierra Research 77to90  4.3E  00:50:57  11/16/01  
!...Switches:                     
!      PARAMETER(NOD=220)
!      COMMON/BLUME/COEFN2(4),COEFNK(4),COEFVK(4)
!      COMMON/PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
!     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
!      COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
!      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
!     :       (IQMAX,MAX_(1))
!      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
!
      ZETA = FINE*Z*QUADR(I1,I2,-3) 
      LB = L(I1) 
      DO I = 1, NCLOSD 
         LA = L(I) 
         ZETA = ZETA - (4*LA + 2)*SN(I1,I,I2,I,0) 
         CALL BWINT (LA, LB) 
         KE1 = 2 
         IF (LA .NE. LB) KE1 = IABS(LA - LB) 
         IP = 0 
         DO K = KE1, LA + LB, 2 
            IP = IP + 1 
            ZETA = ZETA + COEFN2(IP)*SN(I1,I,I,I2,K - 2) + COEFNK(IP)*SN(I,I1,&
               I2,I,K) + COEFVK(IP)*(VK(I1,I,I,I2,K - 1) - VK(I,I1,I2,I,K - 1)) 
         END DO 
      END DO 
      RETURN  
      END FUNCTION ZETA 
