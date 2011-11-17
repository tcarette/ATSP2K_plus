!
!     -------------------------------------------------------------
!        B I S E C T
!     -------------------------------------------------------------
!
      SUBROUTINE BISECT(N, EPS1, D, E, E2, LB, UB, MM, M, W, IND, IERR, RV4, &
         RV5) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: MM 
      INTEGER , INTENT(INOUT) :: M 
      INTEGER , INTENT(OUT) :: IERR 
      REAL(DOUBLE) , INTENT(INOUT) :: EPS1 
      REAL(DOUBLE) , INTENT(INOUT) :: LB 
      REAL(DOUBLE) , INTENT(INOUT) :: UB 
      INTEGER , INTENT(INOUT) :: IND(MM) 
      REAL(DOUBLE) , INTENT(IN) :: D(N) 
      REAL(DOUBLE) , INTENT(IN) :: E(N) 
      REAL(DOUBLE) , INTENT(INOUT) :: E2(N) 
      REAL(DOUBLE) , INTENT(INOUT) :: W(MM) 
      REAL(DOUBLE) , INTENT(INOUT) :: RV4(N) 
      REAL(DOUBLE) , INTENT(INOUT) :: RV5(N) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, K, L, P, Q, R, S, II, M1, M2, TAG, ISTURM 
      REAL(DOUBLE) :: U, V, T1, T2, XU, X0, X1, MACHEP 
!-----------------------------------------------
!
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE BISECTION TECHNIQUE
!     IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
!
!     THIS SUBROUTINE FINDS THOSE EIGENVALUES OF A TRIDIAGONAL
!     SYMMETRIC MATRIX WHICH LIE IN A SPECIFIED INTERVAL,
!     USING BISECTION.
!
!     ON INPUT:
!
!        N IS THE ORDER OF THE MATRIX;
!
!        EPS1 IS AN ABSOLUTE ERROR TOLERANCE FOR THE COMPUTED
!          EIGENVALUES.  IF THE INPUT EPS1 IS NON-POSITIVE,
!          IT IS RESET FOR EACH SUBMATRIX TO A DEFAULT VALUE,
!          NAMELY, MINUS THE PRODUCT OF THE RELATIVE MACHINE
!          PRECISION AND THE 1-NORM OF THE SUBMATRIX;
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
!          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;
!
!        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
!          E2(1) IS ARBITRARY;
!
!        LB AND UB DEFINE THE INTERVAL TO BE SEARCHED FOR EIGENVALUES.
!          IF LB IS NOT LESS THAN UB, NO EIGENVALUES WILL BE FOUND;
!
!        MM SHOULD BE SET TO AN UPPER BOUND FOR THE NUMBER OF
!          EIGENVALUES IN THE INTERVAL.  WARNING: IF MORE THAN
!          MM EIGENVALUES ARE DETERMINED TO LIE IN THE INTERVAL,
!          AN ERROR RETURN IS MADE WITH NO EIGENVALUES FOUND.
!
!     ON OUTPUT:
!
!        EPS1 IS UNALTERED UNLESS IT HAS BEEN RESET TO ITS
!          (LAST) DEFAULT VALUE;
!
!        D AND E ARE UNALTERED;
!
!        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
!          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
!          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
!          E2(1) IS ALSO SET TO ZERO;
!
!        M IS THE NUMBER OF EIGENVALUES DETERMINED TO LIE IN (LB,UB);
!
!        W CONTAINS THE M EIGENVALUES IN ASCENDING ORDER;
!
!        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
!          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
!          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
!          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.;
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          3*N+1      IF M EXCEEDS MM;
!
!        RV4 AND RV5 ARE TEMPORARY STORAGE ARRAYS.
!
!     THE ALGOL PROCEDURE STURMCNT CONTAINED IN TRISTURM
!     APPEARS IN BISECT IN-LINE.
!
!     NOTE THAT SUBROUTINE TQL1 OR IMTQL1 IS GENERALLY FASTER THAN
!     BISECT, IF MORE THAN N/4 EIGENVALUES ARE TO BE FOUND.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
!     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
!                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
!                ON S360 ::::::::::
      DATA MACHEP/ 1.D-12/  
!
      IERR = 0 
      TAG = 0 
      T1 = LB 
      T2 = UB 
!     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ENTRIES ::::::::::
      DO I = 1, N 
         IF (I == 1) GO TO 20 
         IF (DABS(E(I)) > MACHEP*(DABS(D(I))+DABS(D(I-1)))) CYCLE  
   20    CONTINUE 
         E2(I) = 0.0D0 
      END DO 
!     :::::::::: DETERMINE THE NUMBER OF EIGENVALUES
!                IN THE INTERVAL ::::::::::
      P = 1 
      Q = N 
      X1 = UB 
      ISTURM = 1 
      GO TO 320 
   60 CONTINUE 
      M = S 
      X1 = LB 
      ISTURM = 2 
      GO TO 320 
   80 CONTINUE 
      M = M - S 
      IF (M > MM) GO TO 980 
      Q = 0 
      R = 0 
!     :::::::::: ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
!                INTERVAL BY THE GERSCHGORIN BOUNDS ::::::::::
  100 CONTINUE 
      IF (R == M) GO TO 1001 
      TAG = TAG + 1 
      P = Q + 1 
      XU = D(P) 
      X0 = D(P) 
      U = 0.0D0 
!
      DO Q = P, N 
         X1 = U 
         U = 0.0D0 
         V = 0.0D0 
         IF (Q /= N) THEN 
            U = DABS(E(Q+1)) 
            V = E2(Q+1) 
         ENDIF 
         XU = DMIN1(D(Q)-(X1+U),XU) 
         X0 = DMAX1(D(Q)+(X1+U),X0) 
         IF (V /= 0.0D0) CYCLE  
         EXIT  
      END DO 
!
      X1 = DMAX1(DABS(XU),DABS(X0))*MACHEP 
      IF (EPS1 <= 0.0D0) EPS1 = -X1 
      IF (P == Q) THEN 
!     :::::::::: CHECK FOR ISOLATED ROOT WITHIN INTERVAL ::::::::::
         IF (T1>D(P) .OR. D(P)>=T2) GO TO 940 
         M1 = P 
         M2 = P 
         RV5(P) = D(P) 
         GO TO 900 
      ENDIF 
      X1 = X1*DFLOAT(Q - P + 1) 
      LB = DMAX1(T1,XU - X1) 
      UB = DMIN1(T2,X0 + X1) 
      X1 = LB 
      ISTURM = 3 
      GO TO 320 
  200 CONTINUE 
      M1 = S + 1 
      X1 = UB 
      ISTURM = 4 
      GO TO 320 
  220 CONTINUE 
      M2 = S 
      IF (M1 > M2) GO TO 940 
!     :::::::::: FIND ROOTS BY BISECTION ::::::::::
      X0 = UB 
      ISTURM = 5 
!
      RV5(M1:M2) = UB 
      RV4(M1:M2) = LB 
!     :::::::::: LOOP FOR K-TH EIGENVALUE
!                FOR K=M2 STEP -1 UNTIL M1 DO --
!                (-DO- NOT USED TO LEGALIZE COMPUTED-GO-TO) ::::::::::
      K = M2 
  250 CONTINUE 
      XU = LB 
!     :::::::::: FOR I=K STEP -1 UNTIL M1 DO -- ::::::::::
      DO II = M1, K 
         I = M1 + K - II 
         IF (XU >= RV4(I)) CYCLE  
         XU = RV4(I) 
         EXIT  
      END DO 
!
      X0 = MIN(RV5(K),X0) 
!     :::::::::: NEXT BISECTION STEP ::::::::::
  300 CONTINUE 
      X1 = (XU + X0)*0.5D0 
      IF (X0 - XU <= 2.0D0*MACHEP*(DABS(XU) + DABS(X0)) + DABS(EPS1)) GO TO 420 
!     :::::::::: IN-LINE PROCEDURE FOR STURM SEQUENCE ::::::::::
  320 CONTINUE 
      S = P - 1 
      U = 1.0D0 
!
      DO I = P, Q 
         IF (U == 0.0D0) THEN 
            V = DABS(E(I))/MACHEP 
         ELSE 
            V = E2(I)/U 
         ENDIF 
         U = D(I) - X1 - V 
         IF (U >= 0.0D0) CYCLE  
         S = S + 1 
      END DO 
!
      GO TO (60,80,200,220,360) ISTURM 
!     :::::::::: REFINE INTERVALS ::::::::::
  360 CONTINUE 
      IF (S >= K) GO TO 400 
      XU = X1 
      IF (S >= M1) GO TO 380 
      RV4(M1) = X1 
      GO TO 300 
  380 CONTINUE 
      RV4(S+1) = X1 
      RV5(S) = MIN(X1,RV5(S)) 
      GO TO 300 
  400 CONTINUE 
      X0 = X1 
      GO TO 300 
!     :::::::::: K-TH EIGENVALUE FOUND ::::::::::
  420 CONTINUE 
      RV5(K) = X1 
      K = K - 1 
      IF (K >= M1) GO TO 250 
!     :::::::::: ORDER EIGENVALUES TAGGED WITH THEIR
!                SUBMATRIX ASSOCIATIONS ::::::::::
  900 CONTINUE 
      S = R 
      R = R + M2 - M1 + 1 
      J = 1 
      K = M1 
!
      DO L = 1, R 
         IF (J <= S) THEN 
            IF (K > M2) EXIT  
            IF (RV5(K) >= W(L)) GO TO 915 
!
            W(L+S-J+1:L+1:(-1)) = W(L+S-J:L:(-1)) 
            IND(L+S-J+1:L+1:(-1)) = IND(L+S-J:L:(-1)) 
         ENDIF 
!
         W(L) = RV5(K) 
         IND(L) = TAG 
         K = K + 1 
         CYCLE  
  915    CONTINUE 
         J = J + 1 
      END DO 
!
  940 CONTINUE 
      IF (Q < N) GO TO 100 
      GO TO 1001 
!     :::::::::: SET ERROR -- UNDERESTIMATE OF NUMBER OF
!                EIGENVALUES IN INTERVAL ::::::::::
  980 CONTINUE 
      IERR = 3*N + 1 
 1001 CONTINUE 
      LB = T1 
      UB = T2 
      RETURN  
!     :::::::::: LAST CARD OF BISECT ::::::::::
      END SUBROUTINE BISECT 
