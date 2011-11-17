      SUBROUTINE DPBDI(ABD, LDA, N, M, DET) 
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
      INTEGER , INTENT(IN) :: LDA 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: M 
      REAL(DOUBLE) , INTENT(IN) :: ABD(LDA,1) 
      REAL(DOUBLE) , INTENT(INOUT) :: DET(2) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
      REAL(DOUBLE) :: S 
!-----------------------------------------------
!
!     DPBDI COMPUTES THE DETERMINANT
!     OF A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE BAND MATRIX
!     USING THE FACTORS COMPUTED BY DPBCO OR DPBFA.
!     IF THE INVERSE IS NEEDED, USE DPBSL  N  TIMES.
!
!     ON ENTRY
!
!        ABD     DOUBLE PRECISION(LDA, N)
!                THE OUTPUT FROM DPBCO OR DPBFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        M       INTEGER
!                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!
!     ON RETURN
!
!        DET     DOUBLE PRECISION(2)
!                DETERMINANT OF ORIGINAL MATRIX IN THE FORM
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                WITH  1.0 .LE. DET(1) .LT. 10.0
!                OR  DET(1) .EQ. 0.0 .
!
!     LINPACK.  THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!
!     INTERNAL VARIABLES
!
!
!     COMPUTE DETERMINANT
!
      DET(1) = 1.0D0 
      DET(2) = 0.0D0 
      S = 10.0D0 
      DO I = 1, N 
         DET(1) = ABD(M+1,I)**2*DET(1) 
!     ...EXIT
         IF (DET(1) == 0.0D0) EXIT  
   10    CONTINUE 
         IF (DET(1) >= 1.0D0) GO TO 20 
         DET(1) = S*DET(1) 
         DET(2) = DET(2) - 1.0D0 
         GO TO 10 
   20    CONTINUE 
   30    CONTINUE 
         IF (DET(1) < S) CYCLE  
         DET(1) = DET(1)/S 
         DET(2) = DET(2) + 1.0D0 
         GO TO 30 
      END DO 
      RETURN  
      END SUBROUTINE DPBDI 
