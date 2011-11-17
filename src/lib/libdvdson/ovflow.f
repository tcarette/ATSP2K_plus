
*=======================================================================
        SUBROUTINE OVFLOW(NUME,LIM,S,SVEC,EIGVAL)
*=======================================================================
*       Called by: DVDRVR
*       Called when the upper limit (LIM) has been reached for the basis
*       expansion. The new S is computed as S'(i,j)=l(i)delta(i,j) where
*       l(i) eigenvalues, and delta of Kronecker, i,j=1,NUME. The new
*       eigenvectors of the small matrix are the unit vectors.
*
*       Subroutines called:
*       DCOPY, DINIT
*-----------------------------------------------------------------------
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION SVEC(LIM*NUME),S((LIM*(LIM+1))/2)
        DIMENSION EIGVAL(LIM)
*-----------------------------------------------------------------------
*   on entry
*   -------
*   NUME        The largest index of eigenvalues wanted.
*   SVEC        the kpass eigenvectors of the smaller system solved
*   EIGVAL      the eigenvalues of this small system
*   on exit
*   -------
*   S           The new small matrix to be solved.
*-----------------------------------------------------------------------
*
* calculation of the new upper S=diag(l1,...,l_NUME) and
* its matrix Svec of eigenvectors (e1,...,e_NUME)
*
        CALL DINIT((NUME*(NUME+1))/2,0.d0,S,1)
        CALL DINIT(NUME*NUME,0.d0,SVEC,1)
        IND=0
        ICUR=0
        DO 10 I=1,NUME
           S(IND+I)=EIGVAL(I)
           SVEC(ICUR+I)=1
           ICUR=ICUR+NUME
   10      IND=IND+I

        RETURN
        END
