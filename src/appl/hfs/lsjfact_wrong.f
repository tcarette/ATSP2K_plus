*     ----------------------------------------------------------------
*     L S J F A C T
*     ----------------------------------------------------------------
*
*     This subroutine calculates the LSJ dependent angular factors
*     for the orbital,spindipolar,contact and quadrupol contribution and
*     saves them in vectors ORBF(K),DIPF(K),CONTF(K),QUADF(K).
*     K is a variable to order the different combinations of J and J'.
*     Note on the phase; MCHF_ASP assumes S + L coupling while the
*     decoupling formulas uses the L + S coupling of Cowan.
*     The needed (-1)**(L+S-J) faktors for l.h.s and r.h.s have to be
*     included for consistency.

      SUBROUTINE LSJFACT(JI,JF,NJQ,IDIAG,ndens)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER SS1,SS2,ALPHJ(24)
      COMMON/COEF/ORBF(20),DIPF(20),CONTF(20),QUADF(20),J(20),JP(20),
     :JJMAX,JJMIN,JJ1MAX,JJ1MIN,JJ2MAX,JJ2MIN,LL1,SS1,LL2,SS2,VOLF(20)
*                       DATA STATEMENTS
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DATA ALPHJ/4H  0 ,4H 1/2,4H  1 ,4H 3/2,4H  2 ,4H 5/2,4H  3 ,4H 7/2
     :          ,4H  4 ,4H 9/2,4H  5 ,4H11/2,4H  6 ,4H13/2,4H  7 ,4H15/2
     :          ,4H  8 ,4H17/2,4H  9 ,4H19/2,4H  10,4H21/2,4H  11,4H23/2
     :/
      DATA GS/2.002319D0/
      NJQ=0
Cwww
      NPH2=(LL2+SS2-LL1-SS1)/2
      PHASE=DFLOAT((-1)**NPH2)
Cwww
      DO 10 JJ=JJMIN,JJMAX,2
         IF (JJ.EQ.0) GO TO 10
         JJP=JJ
         NJQ=NJQ+1

*        J'=J

         IF (((JJ.GE.JJ1MIN).AND.(JJ.LE.JJ1MAX)).AND.
     :   ((JJ.GE.JJ2MIN).AND.(JJ.LE.JJ2MAX))) THEN
	    CALL SIXJ1(LL1,JJ,SS1,JJ,LL2,1,W61)
	    IF(MOD(LL1+SS1+2+JJ,4).NE.0) W61=-W61
	    CALL NINELS(LL1,SS1,JJ,LL2,SS2,JJ,4,2,2,1,INGG,S9J)
	    IF(INGG.EQ.0) THEN
	      S9J=ZERO
	    ELSE
	      CALL NINE12(LL1,SS1,JJ,LL2,SS2,4,S9J)
	    ENDIF
	    CALL SIXJ1(SS1,JJ,LL1,JJ,SS2,1,W62)
	    IF(MOD(LL1+SS1+2+JJ,4).NE.0)W62=-W62
            CALL SIXJ2(LL1,JJ,SS1,JJ,LL2,1,W63)
	    IF(MOD(LL1+SS1+4+JJ,4).NE.0)W63=-W63
            IF(NDENS.EQ.1) THEN
	      CALL SIXJ(LL1,SS1,JJ,JJ,0,LL2,1,W64)
	      IF(MOD(LL1+SS1+JJ,4).NE.0)W64=-W64
	    ENDIF
            SQRT1=DSQRT(DFLOAT((JJ+1)*4)/DFLOAT(JJ*(JJ+2)))
            SQRT2=DSQRT(DFLOAT(JJ*(JJ+1)*(JJ-1))/DFLOAT((JJ+2)*(JJ+3)))

            if (ndens.eq.1) sqrt3 = dsqrt(jj + 1.d0)

            NPH1=(SS2-SS1)/2
            ORBF(NJQ)=SQRT1*W61*2.D0*PHASE
            DIPF(NJQ)=(-1.D0)*DSQRT(30.D0)*SQRT1*S9J*GS*PHASE
         CONTF(NJQ)=DFLOAT((-1)**(NPH1))*SQRT1*W62*2.D0*GS*PHASE/3.D0

            if (ndens.eq.1) volf(njq) = sqrt3*w64*phase

            QUADF(NJQ)=(-1.D0)*SQRT2*W63*2.D0*PHASE
         ELSE
            ORBF(NJQ)=0.D0
            DIPF(NJQ)=0.D0
            CONTF(NJQ)=0.D0

            if (ndens.eq.1) volf(njq) = 0.d0

            QUADF(NJQ)=0.D0
         ENDIF

         J(NJQ)=ALPHJ(JJ+1)
         JP(NJQ)=ALPHJ(JJP+1)

*        J'=J-1

Cwww
Cwww     Observe that A(LS;J,J-1)=(-1)A(SL;J,J-1) therefore 
Cwww     we multiply with PHASE and not PHASE*(-1)
Cwww
         IF (JJP.EQ.JJMIN.OR.IDIAG.EQ.1) GO TO 10
         JJP=JJ-2
         NJQ=NJQ+1

         if (ndens.eq.1) volf(njq) = 0.d0

         IF (((JJ.GE.JJ1MIN).AND.(JJ.LE.JJ1MAX)).AND.
     :   ((JJP.GE.JJ2MIN).AND.(JJP.LE.JJ2MAX))) THEN
	    CALL SIXJ1(LL1,JJ,SS1,JJP,LL2,1,W61)
	    IF(MOD(LL1+SS1+2+JJP,4).NE.0)W61=-W61
	    CALL NINELS(LL1,SS1,JJ,LL2,SS2,JJP,4,2,2,1,INGG,S9J1)
	    IF(INGG.EQ.0) THEN
	      S9J1=ZERO
	    ELSE
	      CALL NINE13(LL1,SS1,JJ,LL2,SS2,4,S9J1)
	    ENDIF
	    CALL SIXJ1(SS1,JJ,LL1,JJP,SS2,1,W62)
	    IF(MOD(LL1+SS1+2+JJP,4).NE.0)W62=-W62
	    CALL SIXJ2(LL1,JJ,SS1,JJP,LL2,1,W63)
	    IF(MOD(LL1+SS1+4+JJP,4).NE.0)W63=-W63
            SQRT1=DSQRT(DFLOAT(2)/DFLOAT(JJ))
            SQRT2=DSQRT(DFLOAT(JJ*(JJ-2))/DFLOAT((JJ+2)))
            NPH1=(SS2-SS1+2)/2
            ORBF(NJQ)=SQRT1*W61*2.D0*PHASE
            DIPF(NJQ)=(-1.D0)*DSQRT(30.D0)*SQRT1*S9J1*GS*PHASE
            CONTF(NJQ)=DFLOAT((-1)**(NPH1))*SQRT1*W62*2.D0*GS*PHASE
     :      /3.D0
            QUADF(NJQ)=(-1.D0)*SQRT2*W63*PHASE/(DSQRT(2.D0)*2.D0)
         ELSE
             ORBF(NJQ)=0.D0
             DIPF(NJQ)=0.D0
             CONTF(NJQ)=0.D0
             QUADF(NJQ)=0.D0
         ENDIF

         J(NJQ)=ALPHJ(JJ+1)
         JP(NJQ)=ALPHJ(JJP+1)

*        J'=J-2

         IF (JJP.EQ.JJMIN) GO TO 10
         JJP=JJ-4
         NJQ=NJQ+1

         if (ndens.eq.1) volf(njq) = 0.d0

         IF (((JJ.GE.JJ1MIN).AND.(JJ.LE.JJ1MAX)).AND.
     :   ((JJP.GE.JJ2MIN).AND.(JJP.LE.JJ2MAX))) THEN
	    CALL SIXJ2(LL1,JJ,SS1,JJP,LL2,1,W61)
	    IF(MOD(LL1+SS1+4+JJP,4).NE.0)W61=-W61
            SQRT2=DSQRT(DFLOAT((JJ-2)*JJ*(JJ-1)))
            QUADF(NJQ)=(-1.D0)*SQRT2*W61*PHASE/8.D0
         ELSE
            ORBF(NJQ)=0.D0
            DIPF(NJQ)=0.D0
            CONTF(NJQ)=0.D0
         ENDIF

40       J(NJQ)=ALPHJ(JJ+1)
         JP(NJQ)=ALPHJ(JJP+1)
10    CONTINUE
      RETURN
      END
