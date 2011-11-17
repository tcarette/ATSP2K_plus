*=======================================================================
        SUBROUTINE SETUP(A,INDROW,INDCOL,IUPPER,NZER,TM,TP,
     :                   N,LIM,NUME,HIEND,DIAG,MINELEM,
     :                   BASIS,AB,S,NIV)
*=======================================================================
*       Subroutine for setting up (i) the initial BASIS if not provided,
*       (ii) the product of the matrix A with the Basis into matrix AB,
*       and (iii) the small matrix S=B^TAB. If no initial estimates are
*       available, the BASIS =(e_i1,e_i2,...,e_iNUME), where i1,i2,...,
*       iNUME are the indices of the NUME lowest diagonal elements, and
*       e_i the i-th unit vector. (ii) and (iii) are handled by ADDABS.
*-----------------------------------------------------------------------
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        LOGICAL IUPPER
        DIMENSION A(NZER),INDROW(NZER),INDCOL(N)
        DIMENSION TM(N),TP(N)
        DIMENSION DIAG(N),BASIS(N*LIM),AB(N*LIM)
        DIMENSION S(LIM*(LIM+1)/2),MINELEM(LIM)
*-----------------------------------------------------------------------
*   on entry
*   --------
*   OP          The block matrix-vector operation, passed to ADDABS
*   N           the order of the matrix A
*   LIM         The limit on the size of the expanding Basis
*   NUME        Largest index of the wanted eigenvalues.
*   HIEND       Logical. True only if the highest eigenpairs are needed.
*   DIAG        Array of size N with the diagonal elements of A
*   MINELEM     Array keeping the indices of the NUME lowest diagonals.
*
*   on exit
*   -------
*   BASIS       The starting basis.
*   AB, S       The starting D=AB, and small matrix S=B^TAB
*   NIV         The starting dimension of BASIS.
*-----------------------------------------------------------------------

        IF ((NIV.GT.LIM).OR.(NIV.LT.NUME)) THEN
*
*          ..Initial estimates are not available. Give as estimates unit
*          ..vectors corresponding to the NUME minimum diagonal elements
*          ..First find the indices of these NUME elements (in MINELEM).
*          ..Array AB is used temporarily as a scratch array.
*
           CALL DINIT(N,-1.D0,AB,1)
           DO 10 I=1,NUME
*             ..imin= the first not gotten elem( NUME<=N )
              DO 20 J=1,N
 20              IF (AB(J).LT.0) GOTO 30
 30           IMIN=J
              DO 40 J=IMIN+1,N
 40              IF ((AB(J).LT.0).AND.
     :              (DIAG(J).LT.DIAG(IMIN))) IMIN=J
              MINELEM(I)=IMIN
              AB(IMIN)=1.D0
 10        CONTINUE
*
*          ..Build the Basis. B_i=e_(MINELEM(i))
*
           CALL DINIT(N*LIM,0.D0,BASIS,1)
           DO 50 J=1,NUME
              I=(J-1)*N+MINELEM(J)
              BASIS(I)=1
 50        CONTINUE

           NIV=NUME
        ENDIF
*
* Find the matrix AB by matrix-vector multiplies, as well as the
* small matrix S = B^TAB.
*
        KPASS=0
        CALL ADDABS(A,INDROW,INDCOL,IUPPER,NZER,TM,TP,
     :              N,LIM,HIEND,KPASS,NIV,BASIS,AB,S)

        RETURN
        END
