!=======================================================================
      SUBROUTINE NEWVEC(N, NUME, LIM, MBLOCK, KPASS, CRITR, ORTHO, NNCV, INCV, &
         DIAG, SVEC, EIGVAL, AB, BASIS, ICV, RESTART, DONE) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!=======================================================================
!
!       Called by: DVDRVR
!
!       It calculates the new expansion vectors of the basis.
!       For each one of the vectors in INCV starting with the largest
!       megnitude one, calculate its residual Ri= DCi-liBCi and check
!       the ||Ri|| for convergence. If it is converged do not add it
!       but look for the immediate larger coefficient and its vector.
!       The above procedure continues until MBLOCK vectors have been
!       added to the basis, or the upper limit has been encountered.
!       Thus only  the required MBLOCK residuals are computed. Then,
!       calculate the first order correction on the added residuals
!       Ri(j) = Ri(j)/(li-Ajj) and orthonormalizes the new vectors
!       to the basis and to themselves.
!
!       Subroutines called:
!       ORTHNRM, DDOT, DGEMV
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:22  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dgemv_I 
      USE dnrm2_I 
      USE orthnrm_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N 
      INTEGER , INTENT(IN) :: NUME 
      INTEGER  :: LIM 
      INTEGER , INTENT(IN) :: MBLOCK 
      INTEGER  :: KPASS 
      INTEGER  :: NNCV 
      REAL(DOUBLE) , INTENT(IN) :: CRITR 
      REAL(DOUBLE)  :: ORTHO 
      LOGICAL  :: RESTART 
      LOGICAL , INTENT(OUT) :: DONE 
      INTEGER , INTENT(INOUT) :: INCV(NUME) 
      INTEGER , INTENT(OUT) :: ICV(NUME) 
      REAL(DOUBLE) , INTENT(IN) :: DIAG(N) 
      REAL(DOUBLE)  :: SVEC(LIM*NUME) 
      REAL(DOUBLE) , INTENT(IN) :: EIGVAL(LIM) 
      REAL(DOUBLE)  :: AB(N*LIM) 
      REAL(DOUBLE)  :: BASIS(N*LIM) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NEWSTART, NADDED, ICVC, LIMADD, ICUR, I, INDX, IROW 
      REAL(DOUBLE) :: SS, DG 
!-----------------------------------------------
!-----------------------------------------------------------------------
!   on entry
!   --------
!   N           The order of the matrix A
!   NUME        The largest index of the eigenvalues wanted.
!   LIM         The limit on the size of the expanding Basis
!   MBLOCK      Maximum number of vectora to enter the basis
!   KPASS       the current dimension of the expanding basis
!   CRITR       Convergence threshold for residuals
!   ORTHO       Orthogonality threshold to be passed to ORTHNRM
!   NNCV        Number of Non ConVerged pairs (MBLOCK will be targeted)
!   INCV        Index to the corresponding SVEC columns of these pairs.
!   DIAG        Array of size N with the diagonal elements of A
!   SVEC,EIGVAL Arrays holding the eigenvectors and eigenvalues of S
!   AB          Array with the vectors D=AB
!   BASIS       the expanding basis having kpass vectors
!   ICV         Index of converged eigenpairs (ICV(i)=1 <=>i converged)
 
!   on exit
!   -------
!   NNCV        The number of vectors finally added to the basis.
!   BASIS       The new basis incorporating the new vectors at the end
!   ICV         Index of converged eigenpairs (updated)
!   DONE        logical, if covergance has been reached.
!   RESTART     logical, if because of extreme loss of orthogonality
!               the Basis should be collapsed to current approximations.
!-----------------------------------------------------------------------
      DONE = .FALSE. 
      NEWSTART = KPASS*N + 1 
      NADDED = 0 
      ICVC = 0 
      LIMADD = MIN(LIM,MBLOCK + KPASS) 
      ICUR = NEWSTART 
!
! Compute RESIDUALS for the MBLOCK of the NNCV not converged vectors.
!
      DO I = 1, NNCV 
         INDX = INCV(I) 
!          ..Compute  Newv=BASIS*Svec_indx , then
!          ..Compute  Newv=AB*Svec_indx - eigval*Newv and then
!          ..compute the norm of the residual of Newv
         CALL DGEMV ('N', N, KPASS, 1.D0, BASIS, N, SVEC((INDX-1)*KPASS+1), 1, &
            0.D0, BASIS(ICUR), 1) 
         CALL DGEMV ('N', N, KPASS, 1.D0, AB, N, SVEC((INDX-1)*KPASS+1), 1, (-&
            EIGVAL(INDX)),BASIS(ICUR), 1) 
         SS = DNRM2(N,BASIS(ICUR),1) 
!
!          ..Check for convergence of this residual
!
         IF (SS < CRITR) THEN 
!             ..Converged,do not add. Go for next non converged one
            ICVC = ICVC + 1 
            ICV(INDX) = 1 
            IF (ICVC < NNCV) CYCLE  
!             ..All have converged.
            DONE = .TRUE. 
            RETURN  
         ELSE 
!             ..Not converged. Add it in the basis
            NADDED = NADDED + 1 
            INCV(NADDED) = INDX 
            IF (NADDED + KPASS == LIMADD) EXIT  
!             ..More to be added in the block
            ICUR = ICUR + N 
         ENDIF 
      END DO 
 
      NNCV = NADDED 
!
! Diagonal preconditioning: newvect(i)=newvect(i)/(l-Aii)
! If (l-Aii) is very small (or zero) divide by 10.D-6
!
      ICUR = NEWSTART - 1 
      DO I = 1, NNCV 
         DO IROW = 1, N 
            DG = EIGVAL(INCV(I)) - DIAG(IROW) 
            IF (ABS(DG) > 1.D-13) THEN 
               BASIS(ICUR+IROW) = BASIS(ICUR+IROW)/DG 
            ELSE 
               BASIS(ICUR+IROW) = BASIS(ICUR+IROW)/1.D-13 
            ENDIF 
         END DO 
         ICUR = ICUR + N 
      END DO 
!
! ORTHONORMALIZATION
!
      CALL ORTHNRM (N, LIM, ORTHO, KPASS, NNCV, AB(NEWSTART), BASIS, RESTART) 
 
      RETURN  
      END SUBROUTINE NEWVEC 
