* 
*     ------------------------------------------------------------------ 
*              H L 
*     ------------------------------------------------------------------ 
* 
*       Returns the value of <i|L|j>, using a special formula to 
*  preserve symmetry. 
* 
      DOUBLE PRECISION FUNCTION HL(EL,I,J,REL) 
      PARAMETER (NOD=220,NWD=30,NWD2=2*NWD) 
      PARAMETER(IWRITE=6)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      CHARACTER EL(*)*3 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD2),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD2),L(NWD2),MAX(NWD2),N(NWD2)
      LOGICAL REL
      IF (IABS(L(I)-L(J)) .EQ. 0) GO TO 3 
      WRITE(IWRITE,4) EL(I),L(I),EL(J),L(J) 
4     FORMAT(10X,'UNALLOWED L VALUES OCCURRED IN HL SUBROUTINE'/ 
     :   2(10X,A3,' HAS L = ',I3)) 
3     LI = L(I) 
      C = 2*LI + 1 
      A1 = -D2/(C*(LI+1)) 
      A2 = A1/((C+D2)*(LI+1)) 
      A3 = A2/((LI+2)*(LI+1)) 
      ZR = Z*R(1) 
      HL = H*C*P(1,I)*P(1,J)*(D1+ZR*(A1+ZR*(A2+ZR*A3))) 
      MM = MIN0(MAX(I)+3,MAX(J)+3,ND-1) 
      K = 2 
      C = D4/D3 
      DI1 = P(K+1,I) - P(K-1,I) 
      DI2 = P(K+1,I) - D2*P(K,I) + P(K-1,I) 
      DJ1 = P(K+1,J) - P(K-1,J) 
      DJ2 = P(K+1,J) - D2*P(K,J) + P(K-1,J) 
      HL = HL + DI1*DJ1 + C*DI2*DJ2 
      DO 1 K = 4,MM,2 
      DI1 = P(K+1,I) - P(K-1,I) 
      DI2 = P(K+1,I) - D2*P(K,I) + P(K-1,I) 
      DI4 = P(K+2,I) - D4*(P(K+1,I)+P(K-1,I)) + D6*P(K,I) +P(K-2,I) 
      DI3 = P(K+2,I) - P(K-2,I) - D2*DI1 
      DI5 = P(K+3,I)-P(K-3,I) - D4*(P(K+2,I)-P(K-2,I)) 
     :   + 5.D0*(P(K+1,I)-P(K-1,I)) 
      DI6 = P(K+3,I)+P(K-3,I) - D6*(P(K+2,I)+P(K-2,I)) 
     :   + 15.D0*(P(K+1,I)+P(K-1,I)) - 20.D0*P(K,I) 
      DJ1 = P(K+1,J) - P(K-1,J) 
      DJ2 = P(K+1,J) - D2*P(K,J) + P(K-1,J) 
      DJ4 = P(K+2,J) - D4*(P(K+1,J)+P(K-1,J)) + D6*P(K,J) +P(K-2,J) 
      DJ3 = P(K+2,J) - P(K-2,J) - D2*DJ1 
      DJ5 = P(K+3,J)-P(K-3,J) - D4*(P(K+2,J)-P(K-2,J)) 
     :   + 5.D0*(P(K+1,J)-P(K-1,J)) 
      DJ6 = P(K+3,J)+P(K-3,J) - D6*(P(K+2,J)+P(K-2,J)) 
     :   + 15.D0*(P(K+1,J)+P(K-1,J)) - 20.D0*P(K,J) 
1     HL = HL + DI1*DJ1 + C*DI2*DJ2 + (DI3*DJ3 + DI2*DJ4+DI4*DJ2)/45.D0 
     :  -(DI3*DJ5+DI5*DJ3)/252.D0 - (DI2*DJ6+DI6*DJ2-1.1*DI4*DJ4)/378.D0 
      TZ = Z + Z 
      C = (LI + D5)**2 
      HL2 = D5*(TZ*R(1) - C)*P(1,I)*P(1,J) 
      DO 2 K = 2,MM,2 
2     HL2 = HL2 + D2*(TZ*R(K) - C)*P(K,I)*P(K,J) 
     :   + (TZ*R(K+1) - C)*P(K+1,I)*P(K+1,J) 
      HL = -HL/(D2*H) + HL2*H1 
      IF (REL) HL=HL-D2*RLSHFT(I,J)
      RETURN 
      END 
