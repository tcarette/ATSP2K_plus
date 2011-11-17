*
*     --------------------------------------------------------------
*       T R E D 1
*     --------------------------------------------------------------
*
*
      SUBROUTINE TRED1(NM,N,A,D,E,E2) 
* 
      INTEGER I,J,K,L,N,II,NM,JP1 
      DOUBLE PRECISION A(NM,N),D(N),E(N),E2(N) 
      DOUBLE PRECISION F,G,H,SCALE 
      DOUBLE PRECISION DSQRT,DABS,DSIGN 
* 
*     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1, 
*     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON. 
*     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971). 
* 
*     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX 
*     TO A SYMMETRIC TRIDIAGONAL MATRIX USING 
*     ORTHOGONAL SIMILARITY TRANSFORMATIONS. 
* 
*     ON INPUT: 
* 
*        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL 
*          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM 
*          DIMENSION STATEMENT; 
* 
*        N IS THE ORDER OF THE MATRIX; 
* 
*        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE 
*          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED. 
* 
*     ON OUTPUT: 
* 
*        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS- 
*          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER 
*          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED; 
* 
*        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX; 
* 
*        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL 
*          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO; 
* 
*        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E. 
*          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED. 
* 
*     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW, 
*     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY 
* 
*     ------------------------------------------------------------------ 
* 
      DO 100 I = 1, N 
  100 D(I) = A(I,I) 
*     :::::::::: FOR I=N STEP -1 UNTIL 1 DO -- :::::::::: 
      DO 300 II = 1, N 
         I = N + 1 - II 
         L = I - 1 
         H = 0.0D0 
         SCALE = 0.0D0 
         IF (L .LT. 1) GO TO 130 
*     :::::::::: SCALE ROW (ALGOL TOL THEN NOT NEEDED) :::::::::: 
         DO 120 K = 1, L 
  120    SCALE = SCALE + DABS(A(I,K)) 
* 
         IF (SCALE .NE. 0.0D0) GO TO 140 
  130    E(I) = 0.0D0 
         E2(I) = 0.0D0 
         GO TO 290 
* 
  140    DO 150 K = 1, L 
            A(I,K) = A(I,K) / SCALE 
            H = H + A(I,K) * A(I,K) 
  150    CONTINUE 
* 
         E2(I) = SCALE * SCALE * H 
         F = A(I,L) 
         G = -DSIGN(DSQRT(H),F) 
         E(I) = SCALE * G 
         H = H - F * G 
         A(I,L) = F - G 
         IF (L .EQ. 1) GO TO 270 
         F = 0.0D0 
* 
         DO 240 J = 1, L 
            G = 0.0D0 
*     :::::::::: FORM ELEMENT OF A*U :::::::::: 
            DO 180 K = 1, J 
  180       G = G + A(J,K) * A(I,K) 
* 
            JP1 = J + 1 
            IF (L .LT. JP1) GO TO 220 
* 
            DO 200 K = JP1, L 
  200       G = G + A(K,J) * A(I,K) 
*     :::::::::: FORM ELEMENT OF P :::::::::: 
  220       E(J) = G / H 
            F = F + E(J) * A(I,J) 
  240    CONTINUE 
* 
         H = F / (H + H) 
*     :::::::::: FORM REDUCED A :::::::::: 
         DO 260 J = 1, L 
            F = A(I,J) 
            G = E(J) - H * F 
            E(J) = G 
* 
            DO 260 K = 1, J 
               A(J,K) = A(J,K) - F * E(K) - G * A(I,K) 
  260    CONTINUE 
* 
  270    DO 280 K = 1, L 
  280    A(I,K) = SCALE * A(I,K) 
* 
  290    H = D(I) 
         D(I) = A(I,I) 
         A(I,I) = H 
  300 CONTINUE 
* 
      RETURN 
*     :::::::::: LAST CARD OF TRED1 :::::::::: 
      END 
