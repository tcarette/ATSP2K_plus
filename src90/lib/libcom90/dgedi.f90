      SUBROUTINE DGEDI(A, LDA, N, IPVT, DET, WORK, JOB) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dscal_I 
      USE daxpy_I 
      USE dswap_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: LDA 
      INTEGER  :: N 
      INTEGER , INTENT(IN) :: JOB 
      INTEGER , INTENT(IN) :: IPVT(1) 
      REAL(DOUBLE)  :: A(LDA,1) 
      REAL(DOUBLE) , INTENT(INOUT) :: DET(2) 
      REAL(DOUBLE) , INTENT(INOUT) :: WORK(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, K, KB, KP1, L, NM1 
      REAL(DOUBLE) :: T, TEN 
!-----------------------------------------------
!
!     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
!     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
!
!     ON ENTRY
!
!        A       DOUBLE PRECISION(LDA, N)
!                THE OUTPUT FROM DGECO OR DGEFA.
!
!        LDA     INTEGER
!                THE LEADING DIMENSION OF THE ARRAY  A .
!
!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .
!
!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM DGECO OR DGEFA.
!
!        WORK    DOUBLE PRECISION(N)
!                WORK VECTOR.  CONTENTS DESTROYED.
!
!        JOB     INTEGER
!                = 11   BOTH DETERMINANT AND INVERSE.
!                = 01   INVERSE ONLY.
!                = 10   DETERMINANT ONLY.
!
!     ON RETURN
!
!        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE UNCHANGED.
!
!        DET     DOUBLE PRECISION(2)
!                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
!                OTHERWISE NOT REFERENCED.
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0
!                OR  DET(1) .EQ. 0.0 .
!
!     ERROR CONDITION
!
!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
!        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
!        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
!        AND IF DGECO HAS SET RCOND .GT. 0.0 OR DGEFA HAS SET
!        INFO .EQ. 0 .
!
!     LINPACK. THIS VERSION DATED 08/14/78 .
!     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
!
!     SUBROUTINES AND FUNCTIONS
!
!     BLAS DAXPY,DSCAL,DSWAP
!     FORTRAN DABS,MOD
!
!     INTERNAL VARIABLES
!
!
!
!     COMPUTE DETERMINANT
!
      IF (JOB/10 /= 0) THEN 
         DET(1) = 1.0D0 
         DET(2) = 0.0D0 
         TEN = 10.0D0 
         DO I = 1, N 
            IF (IPVT(I) /= I) DET(1) = -DET(1) 
            DET(1) = A(I,I)*DET(1) 
!        ...EXIT
            IF (DET(1) == 0.0D0) EXIT  
   10       CONTINUE 
            IF (DABS(DET(1)) >= 1.0D0) GO TO 20 
            DET(1) = TEN*DET(1) 
            DET(2) = DET(2) - 1.0D0 
            GO TO 10 
   20       CONTINUE 
   30       CONTINUE 
            IF (DABS(DET(1)) < TEN) CYCLE  
            DET(1) = DET(1)/TEN 
            DET(2) = DET(2) + 1.0D0 
            GO TO 30 
         END DO 
      ENDIF 
      IF (MOD(JOB,10) /= 0) THEN 
         DO K = 1, N 
            A(K,K) = 1.0D0/A(K,K) 
            T = -A(K,K) 
            CALL DSCAL (K - 1, T, A(1,K), 1) 
            KP1 = K + 1 
            IF (N < KP1) CYCLE  
            DO J = KP1, N 
               T = A(K,J) 
               A(K,J) = 0.0D0 
               CALL DAXPY (K, T, A(1,K), 1, A(1,J), 1) 
            END DO 
         END DO 
!
!        FORM INVERSE(U)*INVERSE(L)
!
         NM1 = N - 1 
         IF (NM1 >= 1) THEN 
            DO KB = 1, NM1 
               K = N - KB 
               KP1 = K + 1 
               WORK(KP1:N) = A(KP1:N,K) 
               A(KP1:N,K) = 0.0D0 
               DO J = KP1, N 
                  T = WORK(J) 
                  CALL DAXPY (N, T, A(1,J), 1, A(1,K), 1) 
               END DO 
               L = IPVT(K) 
               IF (L == K) CYCLE  
               CALL DSWAP (N, A(1,K), 1, A(1,L), 1) 
            END DO 
         ENDIF 
      ENDIF 
      RETURN  
      END SUBROUTINE DGEDI 
