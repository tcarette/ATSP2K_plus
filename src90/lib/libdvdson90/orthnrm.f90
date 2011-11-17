!=======================================================================
      SUBROUTINE ORTHNRM(N, LIM, ORTHO, KPASS, NNCV, SCRA1, BASIS, RESTART) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!=======================================================================
!
!       It orthogonalizes the new NNCV basis vectors starting from the
!       kpass+1, to the previous vectors of the basis and to themselves.
!       A Gram-Schmidt method is followed after which the residuals
!       should be orthogonal to the BASIS. Because of machine arithmetic
!       errors this orthogonality may be lost, and a reorthogonalization
!       procedure is adopted whenever orthogonality loss is above a
!       ORTHO. If after some reorthogonalizations the procedure does not
!       converge to orthogonality, the basis is collapsed to the
!       current eigenvector approximations.
!
!       Subroutines called:
!       DAXPY, DDOT, DSCAL
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:22  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ddot_I 
      USE daxpy_I 
      USE dcopy_I 
      USE dscal_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N 
      INTEGER , INTENT(IN) :: LIM 
      INTEGER , INTENT(IN) :: KPASS 
      INTEGER , INTENT(INOUT) :: NNCV 
      REAL(DOUBLE) , INTENT(IN) :: ORTHO 
      LOGICAL , INTENT(OUT) :: RESTART 
      REAL(DOUBLE) , INTENT(INOUT) :: SCRA1(N) 
      REAL(DOUBLE)  :: BASIS(N*LIM) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ICUR, IV, IBSTART, I 
      REAL(DOUBLE) :: DPREV, DCUR 
!-----------------------------------------------
!-----------------------------------------------------------------------
!   on entry
!   --------
!   N           The order of the matrix A
!   LIM         The limit on the size of the expanding Basis
!   ORTHO       The orthogonality threshold
!   KPASS       The number of basis vectors already in Basis
!   NNCV        The number of new vectors in the basis
!   SCRA1       Scratch vector of size N
!   BASIS       the expanding basis having kpass vectors
!
!   on exit
!   -------
!   BASIS       the new basis orthonormalized
!   RESTART     Logical, if true the algoritm will collapse BASIS.
!-----------------------------------------------------------------------
!
! ORTHOGONALIZATION
!
      RESTART = .FALSE. 
      ICUR = KPASS*N + 1 
!
!       .. do iv=1,nncv
      IV = 1 
   30 CONTINUE 
      DPREV = 1.D+7 
    5 CONTINUE 
      DCUR = 0.D0 
      IBSTART = 1 
      DO I = 1, KPASS + IV - 1 
         SCRA1(I) = DDOT(N,BASIS(IBSTART),1,BASIS(ICUR),1) 
         DCUR = MAX(DCUR,ABS(SCRA1(I))) 
         IBSTART = IBSTART + N 
      END DO 
      IBSTART = 1 
      DO I = 1, KPASS + IV - 1 
         CALL DAXPY (N, (-SCRA1(I)),BASIS(IBSTART), 1, BASIS(ICUR), 1) 
         IBSTART = IBSTART + N 
      END DO 
 
      IF (DCUR >= ORTHO) THEN 
         IF (DCUR > DPREV) THEN 
            RESTART = .TRUE. 
!                ..Adjust the number of added vectors.
            NNCV = IV - 1 
            RETURN  
         ELSE 
            DPREV = DCUR 
            GO TO 5 
         ENDIF 
      ENDIF 
!
! NORMALIZATION
!
      SCRA1(1) = DDOT(N,BASIS(ICUR),1,BASIS(ICUR),1) 
      SCRA1(1) = SQRT(SCRA1(1)) 
      IF (SCRA1(1) < 1D-14) THEN 
         CALL DCOPY (N, BASIS(N*(NNCV-1)+1), 1, BASIS(ICUR), 1) 
         NNCV = NNCV - 1 
      ELSE 
         CALL DSCAL (N, 1/SCRA1(1), BASIS(ICUR), 1) 
         ICUR = ICUR + N 
         IV = IV + 1 
      ENDIF 
      IF (IV <= NNCV) GO TO 30 
 
      RETURN  
      END SUBROUTINE ORTHNRM 
