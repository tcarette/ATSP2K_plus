*
*     -------------------------------------------------------------
*      R W L S
*     -------------------------------------------------------------
*
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
*                                              (k1 k2 k3)          *
*     REDUCED MATRIX ELEMENT       (nl QLS::: W        :::nl QLS)  *
*                                                                  *
*     Written and extended by G. Gaigalas,                         *
*     Vilnius,  Lithuania                          December 1993   *
*     Bruxelles,  Belgium                          December 1995   *
*
      SUBROUTINE RWLS(K1,K2,K3,L,J1,J2,W)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
CGGf                                with    f-shell   ***  beginning
      COMMON/RIBOF/IMPTF(238),IMGTF(238),IMPNF(238),IMGNF(238)
CGGf                                with    f-shell   ***        end
      CALL RUMT(J1,L,L2QJ1,L2SJ1,L2LJ1)
      CALL RUMT(J2,L,L2QJ2,L2SJ2,L2LJ2)
      W=ZERO
      KK1=K1*2
      IF(ITTK(L2QJ1,L2QJ2,KK1).EQ.0)RETURN
      KK2=K2*2
      IF(ITTK(L2LJ1,L2LJ2,KK2).EQ.0)RETURN
      KK3=K3*2
      IF(ITTK(L2SJ1,L2SJ2,KK3).EQ.0)RETURN
      QJ1=DBLE(L2QJ1)/TWO
      LJ1=L2LJ1/2
      SJ1=DBLE(L2SJ1)/TWO
      IF(K2-1)11,14,2
   11 IF(K1.EQ.1)GO TO 13
      IF(K3.EQ.1)GO TO 12
      D=-(TWO*DBLE(L)+ONE)
      GO TO 15
   12 D=-TWO*SQRT(SJ1*(SJ1+ONE))
      GO TO 15
   13 D=-TWO*SQRT(QJ1*(QJ1+ONE))
      GO TO 15
   14 IF(K1.NE.0.OR.K3.NE.0) GO TO 2
      SAKNIS=DBLE(3*LJ1*(LJ1+1))
      SAKNI2=DBLE(L*(L+1))
      D=-SQRT(SAKNIS/SAKNI2)
   15 IF(J1.EQ.J2) GO TO 16
      RETURN
   16 SAKNIS=DBLE((L2QJ1+1)*(L2LJ1+1)*(L2SJ1+1))
      SAKNI2=DBLE(2*L+1)
      W=D*SQRT(SAKNIS/SAKNI2)
      RETURN
    2 CONTINUE
      IF(L.EQ.1) THEN
        CALL RMEWPLS(J1,J2,K1,K2,K3,W)
      ELSEIF(L.EQ.2) THEN
        CALL RMEWDLS(J1,J2,K1,K2,K3,W)
CGGf                                with    f-shell   ***  beginning
      ELSEIF(L.EQ.3) THEN
        IP=IMPNF(J2)
        IG=IMGNF(J2)
        S=ZERO
        DO 7 I=IP,IG
          CALL RUMT(I,L,L2QI,L2SI,L2LI)
          CALL SLS(L,J1,L2QJ1,L2LJ1,L2SJ1,I,L2QI,L2LI,L2SI,S1)
          IF(ABS(S1).LE.EPS)GO TO 7
          CALL SLS(L,I,L2QI,L2LI,L2SI,J2,L2QJ2,L2LJ2,L2SJ2,S2)
          S1=S1*S2
          IF(ABS(S1).LE.EPS)GO TO 7
          CALL SIXJ(2*L,2*L,KK2,L2LJ2,L2LJ1,L2LI,1,S2)
          S1=S1*S2
          CALL SIXJ(1,1,KK1,L2QJ2,L2QJ1,L2QI,1,S2)
          S1=S1*S2
          CALL SIXJ(1,1,KK3,L2SJ2,L2SJ1,L2SI,1,S2)
          S1=S1*S2
    7   S=S+S1
        SS=DBLE((KK1+1)*(KK2+1)*(KK3+1))
        W=S*SQRT(SS)
        LAS=L2QJ1+L2LJ1+L2SJ1+L2QJ2+L2SJ2+L2LJ2+KK1+KK2+KK3
        IF(MOD(LAS,4).NE.0)W=-W
CGGf                                with    f-shell   ***        end
      ENDIF
      RETURN
      END
