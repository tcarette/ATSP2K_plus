*-----------------------------------------------------------------------
*     M U L T W T
*-----------------------------------------------------------------------
*
*     The J dependent weights for cfg JI and JF are  multiplied and 
*     are saved in WTJIJF(NJQ). NJQ is a parameter that orders the 
*     different J and J' combinations. This order agrees with the order 
*     in LSJFACT.

      SUBROUTINE MULTWT(JI,JF,IMCHF,IDIAG,WTJIJF)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION WTJIJF(20)
      INTEGER SS1,SS2
      POINTER(QWT,WT(20,1))
      COMMON /NAN/ QWT
      COMMON/COEF/ORBF(20),DIPF(20),CONTF(20),QUADF(20),J(20),JP(20),
     :JJMAX,JJMIN,JJ1MAX,JJ1MIN,JJ2MAX,JJ2MIN,LL1,SS1,LL2,SS2,volf(20)
      NJQ=0
      DO 10 JJ=JJMIN,JJMAX,2
Cper     IF (JJ.EQ.0) GO TO 10

*        J'=J

         JJP=JJ
         NJQ=NJQ+1
         IF (((JJ.GE.JJ1MIN).AND.(JJ.LE.JJ1MAX)).AND.
     :   ((JJ.GE.JJ2MIN).AND.(JJ.LE.JJ2MAX))) THEN
            IF (IMCHF.LT.3) THEN
               K1=1
            ELSE
               K1=(JJMAX-JJ)/2+1
            ENDIF
            WTJIJF(NJQ)=WT(K1,JI)*WT(K1,JF)
ctc i d
*            write(*,*) K1, JJ, JJMAX, JI, JF
*            write(*,*) JJ1MAX, JJ1MIN, JJ2MAX, JJ2MIN
*            write(*,*) WT(K1,JI), WT(K1,JF), WTJIJF(NJQ)
ctc f d
         ELSE
            WTJIJF(NJQ)=0.D0
ctc i d
*            write(*,*) K1, JJ, JI, JF
*            write(*,*) JJ1MAX, JJ1MIN, JJ2MAX, JJ2MIN
*            write(*,*) 'hellow'
ctc f d
         ENDIF

*        J'=J-1

         IF (JJP.EQ.JJMIN.OR.IDIAG.EQ.1) GO TO 10
         JJP=JJ-2
         NJQ=NJQ+1
         IF (((JJ.GE.JJ1MIN).AND.(JJ.LE.JJ1MAX)).AND.
     :   ((JJP.GE.JJ2MIN).AND.(JJP.LE.JJ2MAX))) THEN
            IF (IMCHF.LT.3) THEN
               K1=1
               K2=1
            ELSE
               K1=(JJMAX-JJ)/2+1
               K2=(JJMAX-JJP)/2+1
            ENDIF
            WTJIJF(NJQ)=WT(K1,JI)*WT(K2,JF)
         ELSE
            WTJIJF(NJQ)=0.D0
         ENDIF

*        J'=J-2

         IF (JJP.EQ.JJMIN) GO TO 10
         JJP=JJ-4
         NJQ=NJQ+1
         IF (((JJ.GE.JJ1MIN).AND.(JJ.LE.JJ1MAX)).AND.
     :   ((JJP.GE.JJ2MIN).AND.(JJP.LE.JJ2MAX))) THEN
            IF (IMCHF.LT.3) THEN
               K1=1
               K2=1
            ELSE
               K1=(JJMAX-JJ)/2+1
               K2=(JJMAX-JJP)/2+1
            ENDIF
            WTJIJF(NJQ)=WT(K1,JI)*WT(K2,JF)
         ELSE
            WTJIJF(NJQ)=0.D0
         ENDIF
10    CONTINUE
      RETURN
      END
