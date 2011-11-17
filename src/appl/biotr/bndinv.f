*
*     ------------------------------------------------------------------
*	B N D I N V
*     ------------------------------------------------------------------
*
        SUBROUTINE BNDINV(A,EL,N,DETERM,EPSIL,ITEST,NSIZE)
C
C       DOUBLE PRECISION MATRIX INVERSION SUBROUTINE
C       FROM "DLYTAP".
C
C*      DOUBLE PRECISION E,F
C*      DOUBLE PRECISION A,EL,D,DSQRT,C,S,DETERP
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(NSIZE,1),EL(NSIZE,1)
        IF(N.LT.2)GO TO 140
        ISL2=0
        K000FX=2
        IF(ISL2.EQ.0)INDSNL=2
        IF(ISL2.EQ.1)INDSNL=1
C       CALL SLITET(2,INDSNL)
C       CALL OVERFL(K000FX)
C       CALL DVCHK(K000FX)
C
C       SET EL = IDENTITY MATRIX
        DO 30 I=1,N
        DO 10 J=1,N
 10     EL(I,J)=0.0D0
 30     EL(I,I)=1.0D0
C
C       TRIANGULARIZE A, FORM EL
C
        N1=N-1
        M=2
        DO 50 J=1,N1
        DO 45 I=M,N
        IF(A(I,J).EQ.0.0D0)GO TO 45
        D=DSQRT(A(J,J)*A(J,J)+A(I,J)*A(I,J))
        C=A(J,J)/D
        S=A(I,J)/D
 38     DO 39 K=J,N
        D=C*A(J,K)+S*A(I,K)
        A(I,K)=C*A(I,K)-S*A(J,K)
        A(J,K)=D
 39     CONTINUE
        DO 40 K=1,N
        D=C*EL(J,K)+S*EL(I,K)
        EL(I,K)=C*EL(I,K)-S*EL(J,K)
        EL(J,K)=D
 40     CONTINUE
 45     CONTINUE
 50     M=M+1
C       CALL OVERFL(K000FX)
C       GO TO (140,51),K000FX
C
C       CALCULATE THE DETERMINANT
 51     DETERP=A(1,1)
        DO 52 I=2,N
 52     DETERP=DETERP*A(I,I)
        DETERM=DETERP
C       CALL OVERFL(K000FX)
C       GO TO (140,520,520),K000FX
C
C       IS MATRIX SINGULAR
 520    F=A(1,1)
        E=A(1,1)
        DO 58 I=2,N
        IF(DABS(F).LT.DABS(A(I,I)))F=A(I,I)
        IF(DABS(E).GT.DABS(A(I,I)))E=A(I,I)
 58     CONTINUE
        EPSILP=EPSIL
        IF(EPSILP.LE.0)EPSILP=1.0E-8
        RAT=E/F
        IF(ABS(RAT).LT.EPSILP)GO TO 130
C
C       INVERT TRIANGULAR MATRIX
        J=N
        DO 100 J1=1,N
C       CALL SLITE(2)
        I=J
        ISL2=1
        DO 90 I1=1,J
C       CALL SLITET(2,K000FX)
        IF(ISL2.EQ.0)K000FX=2
        IF(ISL2.EQ.1)K000FX=1
        IF(ISL2.EQ.1)ISL2=0
        GO TO (70,75),K000FX
 70     A(I,J)=1.0D0/A(I,I)
        GO TO 90
 75     KS=I+1
        D=0.0D0
        DO 80 K=KS,J
 80     D=D+A(I,K)*A(K,J)
        A(I,J)=-D/A(I,I)
 90     I=I-1
 100    J=J-1
C       CALL OVERFL(K000FX)
C       GO TO (140,103,103),K000FX

C103    CALL DVCHK(K000FX)
C       GO TO (140,105),K000FX
C
C       PREMULTIPLY EL BY INVERTED TRIANGULAR MATRIX
 105    M=1
        DO 120 I=1,N
        DO 118 J=1,N
        D=0.0D0
        DO 107 K=M,N
 107    D=D+A(I,K)*EL(K,J)
        EL(I,J)=D
 118    CONTINUE
 120    M=M+1
C       CALL OVERFL(K000FX)
C       GO TO (140,123,123),K000FX
C
C       RECOPY EL TO A
 123    DO 124 I=1,N
        DO 124 J=1,N
 124    A(I,J)=EL(I,J)
        ITEST=0
C126    IF(INDSNL.EQ.1)CALL SLITE(2)
 126    IF(INDSNL.EQ.1)ISL2=1
        RETURN
C
 130    ITEST=1
        GO TO 126
 140    ITEST=-1
        GO TO 126
        END
