 
!=======================================================================
      SUBROUTINE OVFLOW(NUME, LIM, S, SVEC, EIGVAL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!=======================================================================
!       Called by: DVDRVR
!       Called when the upper limit (LIM) has been reached for the basis
!       expansion. The new S is computed as S'(i,j)=l(i)delta(i,j) where
!       l(i) eigenvalues, and delta of Kronecker, i,j=1,NUME. The new
!       eigenvectors of the small matrix are the unit vectors.
!
!       Subroutines called:
!       DCOPY, DINIT
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:22  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
       
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NUME 
      INTEGER , INTENT(IN) :: LIM 
      REAL(DOUBLE)  :: S((LIM*(LIM + 1))/2) 
      REAL(DOUBLE)  :: SVEC(LIM*NUME) 
      REAL(DOUBLE) , INTENT(IN) :: EIGVAL(LIM) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IND, ICUR, I 
!-----------------------------------------------
!-----------------------------------------------------------------------
!   on entry
!   -------
!   NUME        The largest index of eigenvalues wanted.
!   SVEC        the kpass eigenvectors of the smaller system solved
!   EIGVAL      the eigenvalues of this small system
!   on exit
!   -------
!   S           The new small matrix to be solved.
!-----------------------------------------------------------------------
!
! calculation of the new upper S=diag(l1,...,l_NUME) and
! its matrix Svec of eigenvectors (e1,...,e_NUME)
!
      CALL DINIT ((NUME*(NUME + 1))/2, 0.D0, S, 1) 
      CALL DINIT (NUME*NUME, 0.D0, SVEC, 1) 
      IND = 0 
      ICUR = 0 
      DO I = 1, NUME 
         S(IND+I) = EIGVAL(I) 
         SVEC(ICUR+I) = 1 
         ICUR = ICUR + NUME 
         IND = IND + I 
      END DO 
 
      RETURN  
      END SUBROUTINE OVFLOW 
