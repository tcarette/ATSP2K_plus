!
!
!     ------------------------------------------------------------------
!     I N V M A T
!     ------------------------------------------------------------------
!
      SUBROUTINE INVMAT(A, B, MATDIM, NDIM) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
! FIND INVERSE OF MATRIX A
! INPUT :
!        A : MATRIX TO BE INVERTED
!        B : SCRATCH ARRAY
!        MATDIM : PHYSICAL DIMENSION OF MATRICES
!        NDIM :   DIMENSION OF SUBMATRIX TO BE INVERTED
!
! OUTPUT : A : INVERSE MATRIX ( ORIGINAL MATRIX THUS DESTROYED )
! WARNINGS ARE ISSUED IN CASE OF CONVERGENCE PROBLEMS )
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:43:44  11/20/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE bndinv_I 
      USE wrtmat_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: MATDIM 
      INTEGER  :: NDIM 
      REAL(DOUBLE)  :: A(MATDIM,MATDIM) 
      REAL(DOUBLE)  :: B(MATDIM,MATDIM) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ITEST, NTEST 
      REAL(DOUBLE) :: DETERM, EPSIL 
!-----------------------------------------------
!
      ITEST = 0 
      IF (NDIM == 1) THEN 
         IF (A(1,1) /= 0.0D0) THEN 
            A(1,1) = 1.0D0/A(1,1) 
         ELSE 
            ITEST = 1 
         ENDIF 
      ELSE 
         DETERM = 0.0D0 
         EPSIL = 0.0D0 
         CALL BNDINV (A, B, NDIM, DETERM, EPSIL, ITEST, MATDIM) 
      ENDIF 
!
      IF (ITEST /= 0) WRITE (6, '(A,I3)') ' INVERSION PROBLEM NUMBER..', ITEST 
      NTEST = 1 
      IF (NTEST /= 0) THEN 
         WRITE (6, *) ' INVERTED MATRIX ' 
         CALL WRTMAT (A, NDIM, NDIM, MATDIM, MATDIM) 
      ENDIF 
!
      RETURN  
      END SUBROUTINE INVMAT 
