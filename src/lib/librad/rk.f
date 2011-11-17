*
*     ---------------------------------------------------------------
*        R K
*     ---------------------------------------------------------------
*
*      MSOO     Mass     OO
*
*       0       No       No
*       1       Yes      No
*       2       Yes      No
*       3       No       Yes
*       4       Yes      Yes
*       5       Yes      Yes
*
      Double precision FUNCTION RK(i1,i2,i3,i4,K,REL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      Parameter (NOD=220)
      LOGICAL REL
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MSOO,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
      CALL YKF(i1,i3,K,REL)
      RK = QUADS(i2,i4,1)
      IF(MSOO .EQ. 0) RETURN
      IF(MSOO .EQ. 1 .OR. MSOO.EQ. 4) THEN
        IF (K .EQ. 1) RK = RK - RMASS*GRAD(i1,i3)*GRAD(i2,i4)
      ELSEIF(MSOO .EQ. 2 .OR. MSOO.EQ.5) THEN
        RK = RK*(D1 + RMASS/D2)
            IF (K .EQ. 1) RK = RK + Z*RMASS/D2*(
     :    QUADR(i1,i3,1)*QUADR(i2,i4,-2)+QUADR(i1,i3,-2)*QUADR(i2,i4,1))
         END IF
c ! o-o  interaction
      if(k.eq.0.or.MSOO.LT.3) Return             
      if(i1.eq.i2.and.i1.eq.i3.and.i1.eq.i4) Return
      RK=RK+(ZZ(i1,i2,i3,i4,K)+ZZ(i2,i1,i4,i3,K)+
     :       ZZ(i3,i4,i1,i2,K)+ZZ(i4,i3,i2,i1,K))/4.D0
      RETURN
      END
