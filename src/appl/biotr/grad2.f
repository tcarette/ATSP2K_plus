*
*     ------------------------------------------------------------------- 
* *** GRAD2 
*     ------------------------------------------------------------------
* 
* *** THE GRAD2 FUNCTION SUBPROGRAM COMPUTES THE FOLLOWING DIRECTLY 
* ***        <P(I)[R.D + F(DELTA)[ P(J)>  
* 
      DOUBLE PRECISION FUNCTION grad2(I,J)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*3 el
      PARAMETER (NWD=128,NOD=220)
      PARAMETER (IWRITE=6)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON /NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7),el(nwd)
      DIMENSION Q(NOD)
      JJ=I
      II=J
      DO1K=1,NO
1     Q(K)=P(K,JJ)*R(K)
      LI=L(I)
      LJ=L(J)
      IL = IABS(LI - LJ)
      IF (IL .NE. 0 .AND. IL .NE. 2) GO TO 100
      A1=(LI+LJ+D2)/((LI+D1)*(LJ+D1))
      A2=((LJ+D5)*(LJ+D1)+(LJ+D1+D5)*(LI+D1))/((LI+D1)*(LJ+D1))
      A=A1-A2*(LI+LJ+D3)/((LI+LJ+D4)*(LJ+D5))
      FACT=(LJ+D5)/(LI+LJ+D3)
      G=R(1)**2*P(1,I)*P(1,J)*FACT*(1.+A*Z*R(1))
      MM=MIN0(MAX(I)+1,MAX(J)+1,ND)
      K=2
      F1=D5*(P(K+1,II)-P(K-1,II))
      F2=P(K+1,II)-D2*P(K,II)+P(K-1,II)
      G0=Q(K)*R(K)
      G1=D5*(Q(K+1)*R(K+1)-Q(K-1)*R(K-1))
      G2=Q(K+1)*R(K+1)-D2*Q(K)*R(K)+Q(K-1)*R(K-1)
      G=G+D2*F1*G0+(D2*F2*G1+F1*G2)/D3
      DO2K=4,MM,2
      F1=D5*(P(K+1,II)-P(K-1,II))
      F2=P(K+1,II)-D2*P(K,II)+P(K-1,II)
      F3=D5*(P(K+1,II)-P(K-2,II))-D2*F1
      F4=P(K+2,II)+P(K-2,II)-D4*(P(K+1,II)+P(K-1,II))
     1   +D6*P(K,II)
      G0=Q(K)*R(K)
      G1=D5*(Q(K+1)*R(K+1)-Q(K-1)*R(K-1))
      G2=Q(K+1)*R(K+1)-D2*Q(K)*R(K)+Q(K-1)*R(K-1)
      G3=D5*(Q(K+2)*R(K+2)-Q(K-2)*R(K-2))-D2*G1
      G4=Q(K+2)*R(K+2)+Q(K-2)*R(K-2)-D4*(Q(K+1)*R(K+1)
     1   +Q(K-1)*R(K-1))+D6*Q(K)*R(K)
      G=G+D2*F1*G0+(D2*F2*G1+F1*G2)/D3
     1  -(F1*G4-F4*G1+D4*(F2*G3-F3*G2))/90.E0
2     CONTINUE
      U=QUADR(JJ,II,0)
      G=G-D5*U
      DELTA=LJ-LI
      IF(DELTA)10,11,12
10    GRAD2=G-(LI-D2)*U
      RETURN
11    GRAD2=G+(D1+D5)*U
      RETURN
12    GRAD2=G+(LI+D3)*U
      RETURN
100   WRITE(6,101) I,J
101   FORMAT(5X,'L(I)-L(J) NOT=0,2 FOR I = ',I2,' AND J = ',I2)
      STOP
      END
