!
!     --------------------------------------------------------------
!       T R E D 1
!     --------------------------------------------------------------
!
!
      SUBROUTINE TRED1(NM, N, A, D, E, E2) 
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
      REAL(DOUBLE) , INTENT(INOUT) :: A(NM,N) 
      REAL(DOUBLE) , INTENT(INOUT) :: D(N) 
      REAL(DOUBLE) , INTENT(INOUT) :: E(N) 
      REAL(DOUBLE) , INTENT(OUT) :: E2(N) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, K, L, II, JP1 
      REAL(DOUBLE) :: F, G, H, SCALE 
!-----------------------------------------------
!
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED1,
!     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX
!     TO A SYMMETRIC TRIDIAGONAL MATRIX USING
!     ORTHOGONAL SIMILARITY TRANSFORMATIONS.
!
!     ON INPUT:
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT;
!
!        N IS THE ORDER OF THE MATRIX;
!
!        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE
!          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
!
!     ON OUTPUT:
!
!        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
!          FORMATIONS USED IN THE REDUCTION IN ITS STRICT LOWER
!          TRIANGLE.  THE FULL UPPER TRIANGLE OF A IS UNALTERED;
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX;
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
!          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO;
!
!        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
!          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
      DO I = 1, N 
         D(I) = A(I,I) 
      END DO 
!     :::::::::: FOR I=N STEP -1 UNTIL 1 DO -- ::::::::::
      DO II = 1, N 
         I = N + 1 - II 
         L = I - 1 
         H = 0.0D0 
         SCALE = 0.0D0 
         IF (L < 1) GO TO 130 
!     :::::::::: SCALE ROW (ALGOL TOL THEN NOT NEEDED) ::::::::::
         DO K = 1, L 
            SCALE = SCALE + DABS(A(I,K)) 
         END DO 
!
         IF (SCALE /= 0.0D0) GO TO 140 
  130    CONTINUE 
         E(I) = 0.0D0 
         E2(I) = 0.0D0 
         GO TO 290 
!
  140    CONTINUE 
         A(I,:L) = A(I,:L)/SCALE 
         H = SUM(A(I,:L)*A(I,:L)) 
!
         E2(I) = SCALE*SCALE*H 
         F = A(I,L) 
         G = -DSIGN(DSQRT(H),F) 
         E(I) = SCALE*G 
         H = H - F*G 
         A(I,L) = F - G 
         IF (L /= 1) THEN 
            F = 0.0D0 
!
            DO J = 1, L 
               G = 0.0D0 
!     :::::::::: FORM ELEMENT OF A*U ::::::::::
               G = SUM(A(J,:J)*A(I,:J)) 
!
               JP1 = J + 1 
               IF (L >= JP1) THEN 
!
                  G = G + SUM(A(JP1:L,J)*A(I,JP1:L)) 
               ENDIF 
!     :::::::::: FORM ELEMENT OF P ::::::::::
               E(J) = G/H 
               F = F + E(J)*A(I,J) 
            END DO 
!
            H = F/(H + H) 
!     :::::::::: FORM REDUCED A ::::::::::
            DO J = 1, L 
               F = A(I,J) 
               G = E(J) - H*F 
               E(J) = G 
!
               A(J,:J) = A(J,:J) - F*E(:J) - G*A(I,:J) 
            END DO 
         ENDIF 
!
         A(I,:L) = SCALE*A(I,:L) 
!
  290    CONTINUE 
         H = D(I) 
         D(I) = A(I,I) 
         A(I,I) = H 
      END DO 
!
      RETURN  
!     :::::::::: LAST CARD OF TRED1 ::::::::::
      END SUBROUTINE TRED1 
