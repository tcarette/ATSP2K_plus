*            
*     --------------------------------------------------------------
*       F L I N E
*     --------------------------------------------------------------
*            
*     evaluates the line factor for the (J,J') pair.
*                                                                  *
*     Modified by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville           September 1997   *
*
      DOUBLE PRECISION FUNCTION FLINE(K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL REL,VOK
      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      COMMON /MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :J1QN2(31,3),IJFUL(16)
      COMMON /MULT/QSL,QSV,QIL,QIR,QJVL,QJVR
      POINTER(QSL,SL(1)),(QSV,SV(1)),(QIL,IL(1)),(QIR,IR(1)),
     :       (QJVL,JVL(1)),(QJVR,JVR(1))
      COMMON /LSJ/LL1,LL2,IS1,IS2
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      IVL=JVL(K)
      IVR=JVR(K)
      I2HSH=2*IHSH-1
      LL1=J1QN1(I2HSH,2)-1
      LL2=J1QN2(I2HSH,2)-1
      IS1=J1QN1(I2HSH,3)-1
      IS2=J1QN2(I2HSH,3)-1
      LAM2=LAM+LAM
      IF(IFL.NE.4) THEN
        CALL SIXJ(LL1,IS1,IVL,IVR,LAM2,LL2,1,F)
      IF(MOD(LL1+IS1+IVR+LAM2,4).NE.0) F=-F
        IF(IFL.EQ.3) F=F/DBLE(LAM+1)
      FLINE=F
      ELSE
        L4=LAM2-2
        CALL NINELS(IVL,LL1,IS1,IVR,LL2,IS2,LAM2,L4,2,1,IN,F)
      IF(IN.EQ.1) THEN
          CALL NINELS(IVL,LL1,IS1,IVR,LL2,IS2,LAM2,L4,2,0,IN,F)
          FLINE=F*SQRT(DBLE(LAM2+1))
      ELSE
        FLINE=ZERO
      ENDIF
      ENDIF
      IF(MOD(LL1+IS1+LL2+IS2-IVL-IVR,4).NE.0) FLINE=-FLINE
      RETURN
      END
