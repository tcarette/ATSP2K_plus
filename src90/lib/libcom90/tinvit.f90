!
!     ---------------------------------------------------------------
!        T I N V I T
!     ---------------------------------------------------------------
!
!
      SUBROUTINE TINVIT(NM, N, D, E, E2, M, W, IND, Z, IERR, RV1, RV2, RV3, RV4&
         , RV6) 
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
      INTEGER , INTENT(IN) :: NM 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: M 
      INTEGER , INTENT(OUT) :: IERR 
      INTEGER , INTENT(IN) :: IND(M) 
      REAL(DOUBLE) , INTENT(IN) :: D(N) 
      REAL(DOUBLE) , INTENT(IN) :: E(N) 
      REAL(DOUBLE) , INTENT(IN) :: E2(N) 
      REAL(DOUBLE) , INTENT(IN) :: W(M) 
      REAL(DOUBLE) , INTENT(INOUT) :: Z(NM,M) 
      REAL(DOUBLE) , INTENT(INOUT) :: RV1(N) 
      REAL(DOUBLE) , INTENT(INOUT) :: RV2(N) 
      REAL(DOUBLE) , INTENT(INOUT) :: RV3(N) 
      REAL(DOUBLE) , INTENT(INOUT) :: RV4(N) 
      REAL(DOUBLE) , INTENT(INOUT) :: RV6(N) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, P, Q, R, S, II, IP, JJ, ITS, TAG, GROUP 
      REAL(DOUBLE) :: U, V, UK, XU, X0, X1, EPS2, EPS3, EPS4, NORM, ORDER, &
         MACHEP 
!-----------------------------------------------
!
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE INVERSE ITERATION TECH-
!     NIQUE IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
!
!     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
!     SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
!     USING INVERSE ITERATION.
!
!     ON INPUT:
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT;
!
!        N IS THE ORDER OF THE MATRIX;
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
!          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;
!
!        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,
!          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
!          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
!          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM
!          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN
!          0.0D0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0D0
!          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,
!          TRIDIB, OR  IMTQLV  HAS BEEN USED TO FIND THE EIGENVALUES,
!          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE;
!
!        M IS THE NUMBER OF SPECIFIED EIGENVALUES;
!
!        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER;
!
!        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
!          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
!          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
!          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.
!
!     ON OUTPUT:
!
!        ALL INPUT ARRAYS ARE UNALTERED;
!
!        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.
!          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO;
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
!                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS;
!
!        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
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
      IF (M /= 0) THEN 
         TAG = 0 
         ORDER = 1.0D0 - E2(1) 
         Q = 0 
!     :::::::::: ESTABLISH AND PROCESS NEXT SUBMATRIX ::::::::::
  100    CONTINUE 
         P = Q + 1 
!
         DO Q = P, N 
            IF (Q == N) EXIT  
            IF (E2(Q+1) /= 0.0D0) CYCLE  
            EXIT  
         END DO 
!     :::::::::: FIND VECTORS BY INVERSE ITERATION ::::::::::
         TAG = TAG + 1 
         S = 0 
!
         DO R = 1, M 
            IF (IND(R) /= TAG) CYCLE  
            ITS = 1 
            X1 = W(R) 
            IF (S /= 0) GO TO 510 
!     :::::::::: CHECK FOR ISOLATED ROOT ::::::::::
            XU = 1.0D0 
            IF (P == Q) THEN 
               RV6(P) = 1.0D0 
               GO TO 870 
            ENDIF 
            NORM = DABS(D(P)) 
            IP = P + 1 
!
            DO I = IP, Q 
               NORM = NORM + DABS(D(I)) + DABS(E(I)) 
            END DO 
!     :::::::::: EPS2 IS THE CRITERION FOR GROUPING,
!                EPS3 REPLACES ZERO PIVOTS AND EQUAL
!                ROOTS ARE MODIFIED BY EPS3,
!                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ::::::::::
            EPS2 = 1.0D-3*NORM 
            EPS3 = MACHEP*NORM 
            UK = DFLOAT(Q - P + 1) 
            EPS4 = UK*EPS3 
            UK = EPS4/DSQRT(UK) 
            S = P 
  505       CONTINUE 
            GROUP = 0 
            GO TO 520 
!     :::::::::: LOOK FOR CLOSE OR COINCIDENT ROOTS ::::::::::
  510       CONTINUE 
            IF (DABS(X1 - X0) >= EPS2) GO TO 505 
            GROUP = GROUP + 1 
            IF (ORDER*(X1 - X0) <= 0.0D0) X1 = X0 + ORDER*EPS3 
!     :::::::::: ELIMINATION WITH INTERCHANGES AND
!                INITIALIZATION OF VECTOR ::::::::::
  520       CONTINUE 
            V = 0.0D0 
!
            DO I = P, Q 
               RV6(I) = UK 
               IF (I /= P) THEN 
                  IF (DABS(E(I)) >= DABS(U)) THEN 
!     :::::::::: WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
!                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY ::::::::::
                     XU = U/E(I) 
                     RV4(I) = XU 
                     RV1(I-1) = E(I) 
                     RV2(I-1) = D(I) - X1 
                     RV3(I-1) = 0.0D0 
                     IF (I /= Q) RV3(I-1) = E(I+1) 
                     U = V - XU*RV2(I-1) 
                     V = -XU*RV3(I-1) 
                     CYCLE  
                  ENDIF 
                  XU = E(I)/U 
                  RV4(I) = XU 
                  RV1(I-1) = U 
                  RV2(I-1) = V 
                  RV3(I-1) = 0.0D0 
               ENDIF 
               U = D(I) - X1 - XU*V 
               IF (I == Q) CYCLE  
               V = E(I+1) 
            END DO 
!
            IF (U == 0.0D0) U = EPS3 
            RV1(Q) = U 
            RV2(Q) = 0.0D0 
            RV3(Q) = 0.0D0 
!     :::::::::: BACK SUBSTITUTION
!                FOR I=Q STEP -1 UNTIL P DO -- ::::::::::
  600       CONTINUE 
            DO II = 1, Q - P + 1 
               RV6(Q+1-II) = (RV6(Q+1-II)-U*RV2(Q+1-II)-V*RV3(Q+1-II))/RV1(Q+1-&
                  II) 
               V = U 
               U = RV6(Q+1-II) 
            END DO 
!     :::::::::: ORTHOGONALIZE WITH RESPECT TO PREVIOUS
!                MEMBERS OF GROUP ::::::::::
            IF (GROUP /= 0) THEN 
               J = R 
!
               DO JJ = 1, GROUP 
  630             CONTINUE 
                  J = J - 1 
                  DO WHILE(IND(J) /= TAG) 
                     J = J - 1 
                  END DO 
                  XU = 0.0D0 
!
                  XU = SUM(RV6(P:Q)*Z(P:Q,J)) 
!
                  RV6(P:Q) = RV6(P:Q) - XU*Z(P:Q,J) 
!
               END DO 
            ENDIF 
!
            NORM = 0.0D0 
!
            DO I = P, Q 
               NORM = NORM + DABS(RV6(I)) 
            END DO 
!
            IF (NORM < 1.0D0) THEN 
!     :::::::::: FORWARD SUBSTITUTION ::::::::::
               IF (ITS == 5) GO TO 830 
               IF (NORM == 0.0D0) THEN 
                  RV6(S) = EPS4 
                  S = S + 1 
                  IF (S > Q) S = P 
               ELSE 
                  XU = EPS4/NORM 
!
                  RV6(P:Q) = RV6(P:Q)*XU 
               ENDIF 
!     :::::::::: ELIMINATION OPERATIONS ON NEXT VECTOR
!                ITERATE ::::::::::
               DO I = IP, Q 
                  U = RV6(I) 
!     :::::::::: IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
!                WAS PERFORMED EARLIER IN THE
!                TRIANGULARIZATION PROCESS ::::::::::
                  IF (RV1(I-1) == E(I)) THEN 
                     U = RV6(I-1) 
                     RV6(I-1) = RV6(I) 
                  ENDIF 
                  RV6(I) = U - RV4(I)*RV6(I-1) 
               END DO 
!
               ITS = ITS + 1 
               GO TO 600 
!     :::::::::: SET ERROR -- NON-CONVERGED EIGENVECTOR ::::::::::
  830          CONTINUE 
               IERR = -R 
               XU = 0.0D0 
            ELSE 
!     :::::::::: NORMALIZE SO THAT SUM OF SQUARES IS
!                1 AND EXPAND TO FULL ORDER ::::::::::
               U = 0.0D0 
!
               U = SUM(RV6(P:Q)**2) 
!
               XU = 1.0D0/DSQRT(U) 
!
            ENDIF 
  870       CONTINUE 
            Z(:N,R) = 0.0D0 
!
            Z(P:Q,R) = RV6(P:Q)*XU 
!
            X0 = X1 
         END DO 
!
         IF (Q < N) GO TO 100 
      ENDIF 
      RETURN  
!     :::::::::: LAST CARD OF TINVIT ::::::::::
      END SUBROUTINE TINVIT 
