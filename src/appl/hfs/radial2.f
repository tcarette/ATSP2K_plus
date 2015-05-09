
*--------------------------------------------------------------------
*     R A D I A L 2
*--------------------------------------------------------------------
*
      SUBROUTINE RADIAL2(AZSQR1,maxorb)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,NWD=94)
      DIMENSION AZSQR1(NWD,NWD)
      COMMON/ADATA/P(NOD,NWD),R(NOD),RR(NOD),R2(NOD),
     :       ZED,AZ(NWD),L(NWD),MAX(NWD)
      COMMON/HYPER/VHY(20),NHY,IRHY(20),ISHY(20),JRHY(20),JSHY(20),
     1 IMUHY(20),IMUPHY(20),JMUHY(20),JMUPHY(20),
     2 INUHY(20),INUPHY(20),JNUHY(20),JNUPHY(20)
      COMMON/OVRINT/IOVEL(2),IOVER(2),IOVEP(2)

      DO 10 J = 1,maxorb
         DO 20 K = 1,maxorb
            AZSQR1(J,K) = AZ(J)*AZ(K)
20       CONTINUE
10    CONTINUE
      RETURN
      END
