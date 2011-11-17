!=======================================================================
      LOGICAL FUNCTION TSTSEL (KPASS, NUME, NEIG, ISELEC, SVEC, EIGVAL, ICV, &
         CRITE, CRITC, ROWLAST, IND, OLDVAL, NNCV, INCV) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!=======================================================================
!
!       Called by: DVDRVR
!
!       It first checks if the wanted eigenvalues have reached
!       convergence and updates OLDVAL. Second, for each wanted and non
!       converged eigenvector, it finds the largest absolute coefficient
!       of the NNCV last added vectors (from SVEC) and if not coverged,
!       places it in ROWLAST. IND has the corresponding indices.
!       Third, it sorts ROWLAST in decreasing order and places the
!       corresponding indices in the array INCV. The index array INCV
!       and the number of unconverged pairs NNCV, are passed to DVDRVR.
!       Later in NEWVEC only the first MBLOCK of NNCV pairs will be
!       targeted, since if ROWLAST(i) > ROWLAST(j)
!       then approximately RESIDUAL(i) > RESIDUAL(j)
!
!       Subroutines called
!       IDAMAX
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:22  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE idamax_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: KPASS 
      INTEGER , INTENT(IN) :: NUME 
      INTEGER , INTENT(IN) :: NEIG 
      INTEGER , INTENT(INOUT) :: NNCV 
      REAL(DOUBLE) , INTENT(IN) :: CRITE 
      REAL(DOUBLE) , INTENT(IN) :: CRITC 
      INTEGER , INTENT(IN) :: ISELEC(NEIG) 
      INTEGER , INTENT(INOUT) :: ICV(NUME) 
      INTEGER , INTENT(INOUT) :: IND(NEIG) 
      INTEGER , INTENT(OUT) :: INCV(NEIG) 
      REAL(DOUBLE) , INTENT(IN) :: SVEC(KPASS*NUME) 
      REAL(DOUBLE) , INTENT(IN) :: EIGVAL(NUME) 
      REAL(DOUBLE)  :: ROWLAST(NEIG) 
      REAL(DOUBLE) , INTENT(IN) :: OLDVAL(NUME) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NNCE, I, IVAL, ICNT, ICUR, L, INDX, ITEMP 
      REAL(DOUBLE) :: TMAX, TEMP 
      LOGICAL :: DONE 
!-----------------------------------------------
!-----------------------------------------------------------------------
!
!   on entry
!   -------
!   KPASS       current dimension of the expanding Basis
!   NUME        Largest index of the wanted eigenvalues.
!   NEIG        number of wanted eigenvalues of original matrix
!   ISELEC      index array of the wanted eigenvalues.
!   SVEC        the eigenvectors of the small system
!   EIGVAL      The NUME lowest eigenvalues of the small problem
!   ICV         Index of converged eigenpairs.ICV(i)=1 iff eigenpair i
!               has converged, and ICV(i)=0 if eigenpair i has not.
!   CRITE,CRITC Convergence thresholds for eigenvalues and coefficients
!   ROWLAST     scratch array, keeping the largest absolute coefficient
!               of the NNCV last rows of Svec.
!   IND         scratch array, temporary keeping the indices of Rowlast
!   OLDVAL      The previous iteration's eigenvalues.
!
!   on exit
!   -------
!   NNCV         Number of non converged eigenvectors (to be targeted)
!   INCV         Index to these columns in decreasing order of magnitude
!   TSTSEL       true if convergence has been reached
!
!-----------------------------------------------------------------------
 
      DONE = .FALSE. 
!
! Test all wanted eigenvalues for convergence under CRITE
!
      NNCE = 0 
      DO I = 1, NEIG 
         IVAL = ISELEC(I) 
         IF (ABS(OLDVAL(IVAL)-EIGVAL(IVAL)) < CRITE) CYCLE  
         NNCE = NNCE + 1 
      END DO 
      IF (NNCE == 0) THEN 
         TSTSEL = .TRUE. 
         RETURN  
      ENDIF 
!
! Find the maximum element of the last NNCV coefficients of unconverged
! eigenvectors. For those unconverged coefficients, put their indices
! to IND and find their number NNCV
!
      ICNT = 0 
      DO I = 1, NEIG 
         IF (ICV(ISELEC(I)) /= 0) CYCLE  
!             ..Find coefficient and test for convergence
         ICUR = KPASS*ISELEC(I) 
         TMAX = ABS(SVEC(ICUR)) 
         DO L = 1, NNCV - 1 
            TMAX = MAX(TMAX,ABS(SVEC(ICUR-L))) 
         END DO 
         IF (TMAX < CRITC) THEN 
!                ..this  coefficient converged
            ICV(ISELEC(I)) = 1 
         ELSE 
!                ..Not converged. Add it to the list.
            ICNT = ICNT + 1 
            IND(ICNT) = ISELEC(I) 
            ROWLAST(ICNT) = TMAX 
         ENDIF 
      END DO 
 
      NNCV = ICNT 
      IF (NNCV == 0) DONE = .TRUE. 
!
! Sort the ROWLAST elements interchanging their indices as well
!
      DO I = 1, NNCV 
         INDX = IDAMAX(NNCV - I + 1,ROWLAST(I),1) 
         INCV(I) = IND(INDX+I-1) 
 
         TEMP = ROWLAST(INDX+I-1) 
         ROWLAST(INDX+I-1) = ROWLAST(I) 
         ROWLAST(I) = TEMP 
         ITEMP = IND(INDX+I-1) 
         IND(INDX+I-1) = IND(I) 
         IND(I) = ITEMP 
      END DO 
 
      TSTSEL = DONE 
      RETURN  
      END FUNCTION TSTSEL 
