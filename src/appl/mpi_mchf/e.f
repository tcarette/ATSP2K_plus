        DOUBLE PRECISION FUNCTION E(I,J)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
        PARAMETER (NWD=94,NOD=220,NOFFD=800)
        POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :     (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :     (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
        COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :      QMETH,QIEPTR
*
        COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
*
        IBEGIN = 1
        IF (I .GT. 1) IBEGIN = IEPTR(I-1) + 1
        IEND = IEPTR(I)
        E = 0.D0
        DO 10 II = IBEGIN,IEND
           IF (IJE(II) .EQ. J) THEN
              E = EIJ(II)
              RETURN
           END IF
 10     CONTINUE
        END
