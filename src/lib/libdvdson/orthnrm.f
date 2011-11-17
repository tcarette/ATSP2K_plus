*=======================================================================
        SUBROUTINE ORTHNRM(N,LIM,ORTHO,KPASS,NNCV,SCRA1,
     :                     BASIS,RESTART)
*=======================================================================
*
*       It orthogonalizes the new NNCV basis vectors starting from the
*       kpass+1, to the previous vectors of the basis and to themselves.
*       A Gram-Schmidt method is followed after which the residuals
*       should be orthogonal to the BASIS. Because of machine arithmetic
*       errors this orthogonality may be lost, and a reorthogonalization
*       procedure is adopted whenever orthogonality loss is above a 
*       ORTHO. If after some reorthogonalizations the procedure does not
*       converge to orthogonality, the basis is collapsed to the 
*       current eigenvector approximations.
*
*       Subroutines called:
*       DAXPY, DDOT, DSCAL
*-----------------------------------------------------------------------
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION BASIS(N*LIM)
        DIMENSION SCRA1(N)
        LOGICAL RESTART
*-----------------------------------------------------------------------
*   on entry
*   --------
*   N           The order of the matrix A
*   LIM         The limit on the size of the expanding Basis
*   ORTHO       The orthogonality threshold
*   KPASS       The number of basis vectors already in Basis
*   NNCV        The number of new vectors in the basis
*   SCRA1       Scratch vector of size N
*   BASIS       the expanding basis having kpass vectors
*
*   on exit
*   -------
*   BASIS       the new basis orthonormalized
*   RESTART     Logical, if true the algoritm will collapse BASIS.
*-----------------------------------------------------------------------
*
* ORTHOGONALIZATION
*
        RESTART=.false.
        ICUR=KPASS*N+1
*
*       .. do iv=1,nncv 
        IV = 1
 30     CONTINUE

           DPREV=1.D+7
 5         DCUR=0.D0
           IBSTART=1
           DO 10 I=1,KPASS+IV-1
              SCRA1(I)=DDOT(N,BASIS(IBSTART),1,BASIS(ICUR),1)
              DCUR=MAX(DCUR,ABS(SCRA1(I)))
              IBSTART=IBSTART+N
 10        CONTINUE
           IBSTART=1
           DO 20 I=1,KPASS+IV-1
              CALL DAXPY(N,-SCRA1(I),BASIS(IBSTART),1,BASIS(ICUR),1)
              IBSTART=IBSTART+N
 20        CONTINUE

           IF (DCUR.GE.ORTHO) THEN
              IF (DCUR.GT.DPREV) THEN
                 RESTART=.true.
*                ..Adjust the number of added vectors.
                 NNCV=IV-1
                 RETURN
              ELSE
                 DPREV=DCUR
                 GOTO 5
              ENDIF
           ENDIF
*
* NORMALIZATION
*
           SCRA1(1)=DDOT(N,BASIS(ICUR),1,BASIS(ICUR),1)
           SCRA1(1)=SQRT(SCRA1(1))
           IF (SCRA1(1).LT.1D-14) THEN
              CALL DCOPY(N,BASIS( N*(NNCV-1)+1),1,BASIS(ICUR),1)
              NNCV=NNCV-1
           ELSE
              CALL DSCAL(N,1/SCRA1(1),BASIS(ICUR),1)
              ICUR=ICUR+N
              IV = IV +1
           ENDIF
        IF (IV.LE.NNCV) GOTO 30

        RETURN
        END
