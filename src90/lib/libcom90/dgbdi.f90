      SUBROUTINE DGBDI(ABD, LDA, N, ML, MU, IPVT, DET) 
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
      INTEGER , INTENT(IN) :: ML 
      INTEGER , INTENT(IN) :: MU 
      INTEGER , INTENT(IN) :: IPVT(1) 
      REAL(DOUBLE) , INTENT(IN) :: ABD(LDA,1) 
      REAL(DOUBLE) , INTENT(INOUT) :: DET(2) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, M 
      REAL(DOUBLE) :: TEN 
!-----------------------------------------------
!
!     DGBDI COMPUTES THE DETERMINANT OF A BAND MATRIX
!     USING THE FACTORS COMPUTED BY DGBCO OR DGBFA.
!     IF THE INVERSE IS NEEDED, USE DGBSL  N  TIMES.
!
!     ON ENTRY
!
!        ABD     DOUBLE PRECISION(LDA, N)
!                THE OUTPUT FROM DGBCO OR DGBFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  ABD .
!
!        N       INTEGER
!                THE ORDER OF THE ORIGINAL MATRIX.
!
!        ML      INTEGER
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
!
!        MU      INTEGER
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
!
!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM DGBCO OR DGBFA.
!
!     ON RETURN
!
!        DET     DOUBLE PRECISION(2)
!                DETERMINANT OF ORIGINAL MATRIX.
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0
!                OR  DET(1) = 0.0 .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     FORTRAN DABS
!
!     INTERNAL VARIABLES
!
!
!
      M = ML + MU + 1 
      DET(1) = 1.0D0 
      DET(2) = 0.0D0 
      TEN = 10.0D0 
      DO I = 1, N 
         IF (IPVT(I) /= I) DET(1) = -DET(1) 
         DET(1) = ABD(M,I)*DET(1) 
!     ...EXIT
         IF (DET(1) == 0.0D0) EXIT  
   10    CONTINUE 
         IF (DABS(DET(1)) >= 1.0D0) GO TO 20 
         DET(1) = TEN*DET(1) 
         DET(2) = DET(2) - 1.0D0 
         GO TO 10 
   20    CONTINUE 
   30    CONTINUE 
         IF (DABS(DET(1)) < TEN) CYCLE  
         DET(1) = DET(1)/TEN 
         DET(2) = DET(2) + 1.0D0 
         GO TO 30 
      END DO 
      RETURN  
      END SUBROUTINE DGBDI 
