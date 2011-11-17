*
*     -------------------------------------------------------------
*        B I S E C T
*     -------------------------------------------------------------
*
      SUBROUTINE BISECT(N,EPS1,D,E,E2,LB,UB,MM,M,W,IND,IERR,RV4,RV5) 
* 
      INTEGER I,J,K,L,M,N,P,Q,R,S,II,MM,M1,M2,TAG,IERR,ISTURM 
      DOUBLE PRECISION D(N),E(N),E2(N),W(MM),RV4(N),RV5(N) 
      DOUBLE PRECISION U,V,LB,T1,T2,UB,XU,X0,X1,EPS1,MACHEP 
      DOUBLE PRECISION DABS,DMAX1,DMIN1,DFLOAT 
      INTEGER IND(MM) 
* 
*     THIS SUBROUTINE IS A TRANSLATION OF THE BISECTION TECHNIQUE 
*     IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON. 
*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971). 
* 
*     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL 
*     SYMMETRIC MATRIX WHICH LIE IN A SPECIFIED INTERVAL, 
*     USING BISECTION. 
* 
*     ON INPUT: 
* 
*        N IS THE ORDER OF THE MATRIX; 
* 
*        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED 
*          EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE, 
*          IT IS RESET FOR EACH SUBMATRIX TO A DEFAULT VALUE, 
*          NAMELY, MINUS THE PRODUCT OF THE RELATIVE MACHINE 
*          PRECISION AND THE 1-NORM OF THE SUBMATRIX; 
* 
*        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX; 
* 
*        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX 
*          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY; 
* 
*        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E. 
*          E2(1) IS ARBITRARY; 
* 
*        LB AND UB DEFINE THE INTERVAL TO BE SEARCHED FOR EIGENVALUES. 
*          IF LB IS NOT LESS THAN UB, NO EIGENVALUES WILL BE FOUND; 
* 
*        MM SHOULD BE SET TO AN UPPER BOUND FOR THE NUMBER OF 
*          EIGENVALUES IN THE INTERVAL.  WARNING: IF MORE THAN 
*          MM EIGENVALUES ARE DETERMINED TO LIE IN THE INTERVAL, 
*          AN ERROR RETURN IS MADE WITH NO EIGENVALUES FOUND. 
* 
*     ON OUTPUT: 
* 
*        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS 
*          (LAST) DEFAULT VALUE; 
* 
*        D AND E ARE UNALTERED; 
* 
*        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED 
*          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE 
*          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES. 
*          E2(1) IS ALSO SET TO ZERO; 
* 
*        M IS THE NUMBER OF EIGENVALUES DETERMINED TO LIE IN (LB,UB); 
* 
*        W CONTAINS THE M EIGENVALUES IN ASCENDING ORDER; 
* 
*        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES 
*          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W -- 
*          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM 
*          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.; 
* 
*        IERR IS SET TO 
*          ZERO       FOR NORMAL RETURN, 
*          3*N+1      IF M EXCEEDS MM; 
* 
*        RV4 AND RV5 ARE TEMPORARY STORAGE ARRAYS. 
* 
*     THE ALGOL PROCEDURE STURMCNT CONTAINED IN TRISTURM 
*     APPEARS IN BISECT IN-LINE. 
* 
*     NOTE THAT SUBROUTINE TQL1 OR IMTQL1 IS GENERALLY FASTER THAN 
*     BISECT, IF MORE THAN N/4 EIGENVALUES ARE TO BE FOUND. 
* 
*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW, 
*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY 
* 
*     ------------------------------------------------------------------ 
* 
*     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING 
*                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC. 
*                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC 
*                ON S360 :::::::::: 
      DATA MACHEP/1.D-12/ 
* 
      IERR = 0 
      TAG = 0 
      T1 = LB 
      T2 = UB 
*     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ENTRIES :::::::::: 
      DO 40 I = 1, N 
         IF (I .EQ. 1) GO TO 20 
         IF (DABS(E(I)) .GT. MACHEP * (DABS(D(I)) + DABS(D(I-1)))) 
     :      GO TO 40 
   20    E2(I) = 0.0D0 
   40 CONTINUE 
*     :::::::::: DETERMINE THE NUMBER OF EIGENVALUES 
*                IN THE INTERVAL :::::::::: 
      P = 1 
      Q = N 
      X1 = UB 
      ISTURM = 1 
      GO TO 320 
   60 M = S 
      X1 = LB 
      ISTURM = 2 
      GO TO 320 
   80 M = M - S 
      IF (M .GT. MM) GO TO 980 
      Q = 0 
      R = 0 
*     :::::::::: ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING 
*                INTERVAL BY THE GERSCHGORIN BOUNDS :::::::::: 
  100 IF (R .EQ. M) GO TO 1001 
      TAG = TAG + 1 
      P = Q + 1 
      XU = D(P) 
      X0 = D(P) 
      U = 0.0D0 
* 
      DO 120 Q = P, N 
         X1 = U 
         U = 0.0D0 
         V = 0.0D0 
         IF (Q .EQ. N) GO TO 110 
         U = DABS(E(Q+1)) 
         V = E2(Q+1) 
  110    XU = DMIN1(D(Q)-(X1+U),XU) 
         X0 = DMAX1(D(Q)+(X1+U),X0) 
         IF (V .EQ. 0.0D0) GO TO 140 
  120 CONTINUE 
* 
  140 X1 = DMAX1(DABS(XU),DABS(X0)) * MACHEP 
      IF (EPS1 .LE. 0.0D0) EPS1 = -X1 
      IF (P .NE. Q) GO TO 180 
*     :::::::::: CHECK FOR ISOLATED ROOT WITHIN INTERVAL :::::::::: 
      IF (T1 .GT. D(P) .OR. D(P) .GE. T2) GO TO 940 
      M1 = P 
      M2 = P 
      RV5(P) = D(P) 
      GO TO 900 
  180 X1 = X1 * DFLOAT(Q-P+1) 
      LB = DMAX1(T1,XU-X1) 
      UB = DMIN1(T2,X0+X1) 
      X1 = LB 
      ISTURM = 3 
      GO TO 320 
  200 M1 = S + 1 
      X1 = UB 
      ISTURM = 4 
      GO TO 320 
  220 M2 = S 
      IF (M1 .GT. M2) GO TO 940 
*     :::::::::: FIND ROOTS BY BISECTION :::::::::: 
      X0 = UB 
      ISTURM = 5 
* 
      DO 240 I = M1, M2 
         RV5(I) = UB 
         RV4(I) = LB 
  240 CONTINUE 
*     :::::::::: LOOP FOR K-TH EIGENVALUE 
*                FOR K=M2 STEP -1 UNTIL M1 DO -- 
*                (-DO- NOT USED TO LEGALIZE COMPUTED-GO-TO) :::::::::: 
      K = M2 
  250    XU = LB 
*     :::::::::: FOR I=K STEP -1 UNTIL M1 DO -- :::::::::: 
         DO 260 II = M1, K 
            I = M1 + K - II 
            IF (XU .GE. RV4(I)) GO TO 260 
            XU = RV4(I) 
            GO TO 280 
  260    CONTINUE 
* 
  280    IF (X0 .GT. RV5(K)) X0 = RV5(K) 
*     :::::::::: NEXT BISECTION STEP :::::::::: 
  300    X1 = (XU + X0) * 0.5D0 
         IF ((X0 - XU) .LE. (2.0D0 * MACHEP * 
     :      (DABS(XU) + DABS(X0)) + DABS(EPS1))) GO TO 420 
*     :::::::::: IN-LINE PROCEDURE FOR STURM SEQUENCE :::::::::: 
  320    S = P - 1 
         U = 1.0D0 
* 
         DO 340 I = P, Q 
            IF (U .NE. 0.0D0) GO TO 325 
            V = DABS(E(I)) / MACHEP 
            GO TO 330 
  325       V = E2(I) / U 
  330       U = D(I) - X1 - V 
            IF (U .LT. 0.0D0) S = S + 1 
  340    CONTINUE 
* 
         GO TO (60,80,200,220,360), ISTURM 
*     :::::::::: REFINE INTERVALS :::::::::: 
  360    IF (S .GE. K) GO TO 400 
         XU = X1 
         IF (S .GE. M1) GO TO 380 
         RV4(M1) = X1 
         GO TO 300 
  380    RV4(S+1) = X1 
         IF (RV5(S) .GT. X1) RV5(S) = X1 
         GO TO 300 
  400    X0 = X1 
         GO TO 300 
*     :::::::::: K-TH EIGENVALUE FOUND :::::::::: 
  420    RV5(K) = X1 
      K = K - 1 
      IF (K .GE. M1) GO TO 250 
*     :::::::::: ORDER EIGENVALUES TAGGED WITH THEIR 
*                SUBMATRIX ASSOCIATIONS :::::::::: 
  900 S = R 
      R = R + M2 - M1 + 1 
      J = 1 
      K = M1 
* 
      DO 920 L = 1, R 
         IF (J .GT. S) GO TO 910 
         IF (K .GT. M2) GO TO 940 
         IF (RV5(K) .GE. W(L)) GO TO 915 
* 
         DO 905 II = J, S 
            I = L + S - II 
            W(I+1) = W(I) 
            IND(I+1) = IND(I) 
  905    CONTINUE 
* 
  910    W(L) = RV5(K) 
         IND(L) = TAG 
         K = K + 1 
         GO TO 920 
  915    J = J + 1 
  920 CONTINUE 
* 
  940 IF (Q .LT. N) GO TO 100 
      GO TO 1001 
*     :::::::::: SET ERROR -- UNDERESTIMATE OF NUMBER OF 
*                EIGENVALUES IN INTERVAL :::::::::: 
  980 IERR = 3 * N + 1 
 1001 LB = T1 
      UB = T2 
      RETURN 
*     :::::::::: LAST CARD OF BISECT :::::::::: 
      END 
