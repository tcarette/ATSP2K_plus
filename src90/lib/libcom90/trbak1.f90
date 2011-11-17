!
!     --------------------------------------------------------------
!        T R B A K 1
!     --------------------------------------------------------------
!
!
      SUBROUTINE TRBAK1(NM, N, A, E, M, Z) 
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
      REAL(DOUBLE) , INTENT(IN) :: A(NM,N) 
      REAL(DOUBLE) , INTENT(IN) :: E(N) 
      REAL(DOUBLE) , INTENT(INOUT) :: Z(NM,M) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, K, L 
      REAL(DOUBLE) :: S 
!-----------------------------------------------
!
!
!     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK1,
!     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
!     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
!     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED1.
!
!     ON INPUT:
!
!        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
!          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
!          DIMENSION STATEMENT;
!
!        N IS THE ORDER OF THE MATRIX;
!
!        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANS-
!          FORMATIONS USED IN THE REDUCTION BY  TRED1
!          IN ITS STRICT LOWER TRIANGLE;
!
!        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
!          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;
!
!        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED;
!
!        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
!          IN ITS FIRST M COLUMNS.
!
!     ON OUTPUT:
!
!        Z CONTAINS THE TRANSFORMED EIGENVECTORS
!          IN ITS FIRST M COLUMNS.
!
!     NOTE THAT TRBAK1 PRESERVES VECTOR EUCLIDEAN NORMS.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ------------------------------------------------------------------
!
      IF (M /= 0) THEN 
         IF (N /= 1) THEN 
!
            DO I = 2, N 
               L = I - 1 
               IF (E(I) == 0.0D0) CYCLE  
!
               DO J = 1, M 
!
                  S = SUM(A(I,:L)*Z(:L,J)) 
!     :::::::::: DIVISOR BELOW IS NEGATIVE OF H FORMED IN TRED1.
!                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ::::::::::
                  S = (S/A(I,L))/E(I) 
!
                  Z(:L,J) = Z(:L,J) + S*A(I,:L) 
!
               END DO 
!
            END DO 
!
         ENDIF 
      ENDIF 
      RETURN  
!     :::::::::: LAST CARD OF TRBAK1 ::::::::::
      END SUBROUTINE TRBAK1 
