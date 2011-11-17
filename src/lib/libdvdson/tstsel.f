*=======================================================================
        LOGICAL FUNCTION TSTSEL(KPASS,NUME,NEIG,ISELEC,SVEC,EIGVAL,ICV,
     :                       CRITE,CRITC,ROWLAST,IND,OLDVAL,NNCV,INCV)
*=======================================================================
*
*       Called by: DVDRVR

*       It first checks if the wanted eigenvalues have reached
*       convergence and updates OLDVAL. Second, for each wanted and non
*       converged eigenvector, it finds the largest absolute coefficient
*       of the NNCV last added vectors (from SVEC) and if not coverged,
*       places it in ROWLAST. IND has the corresponding indices.
*       Third, it sorts ROWLAST in decreasing order and places the
*       corresponding indices in the array INCV. The index array INCV
*       and the number of unconverged pairs NNCV, are passed to DVDRVR.
*       Later in NEWVEC only the first MBLOCK of NNCV pairs will be
*       targeted, since if ROWLAST(i) > ROWLAST(j)
*       then approximately RESIDUAL(i) > RESIDUAL(j)
*
*       Subroutines called
*       IDAMAX
*-----------------------------------------------------------------------
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        LOGICAL DONE
        DIMENSION SVEC(KPASS*NUME),EIGVAL(NUME)
        DIMENSION ICV(NUME)
        DIMENSION ROWLAST(NEIG),IND(NEIG),OLDVAL(NUME)
        DIMENSION INCV(NEIG),ISELEC(NEIG)
*-----------------------------------------------------------------------
*
*   on entry
*   -------
*   KPASS       current dimension of the expanding Basis
*   NUME        Largest index of the wanted eigenvalues.
*   NEIG        number of wanted eigenvalues of original matrix
*   ISELEC      index array of the wanted eigenvalues.
*   SVEC        the eigenvectors of the small system
*   EIGVAL      The NUME lowest eigenvalues of the small problem
*   ICV         Index of converged eigenpairs.ICV(i)=1 iff eigenpair i
*               has converged, and ICV(i)=0 if eigenpair i has not.
*   CRITE,CRITC Convergence thresholds for eigenvalues and coefficients
*   ROWLAST     scratch array, keeping the largest absolute coefficient
*               of the NNCV last rows of Svec.
*   IND         scratch array, temporary keeping the indices of Rowlast
*   OLDVAL      The previous iteration's eigenvalues.
*
*   on exit
*   -------
*   NNCV         Number of non converged eigenvectors (to be targeted)
*   INCV         Index to these columns in decreasing order of magnitude
*   TSTSEL       true if convergence has been reached
*
*-----------------------------------------------------------------------

        DONE=.False.
*
* Test all wanted eigenvalues for convergence under CRITE
*
        NNCE=0
        DO 10 I=1,NEIG
           IVAL=ISELEC(I)
 10        IF (ABS(OLDVAL(IVAL)-EIGVAL(IVAL)).GE.CRITE) NNCE=NNCE+1
        IF (NNCE.EQ.0) THEN
           TSTSEL=.TRUE.
           RETURN
        ENDIF
*
* Find the maximum element of the last NNCV coefficients of unconverged
* eigenvectors. For those unconverged coefficients, put their indices
* to IND and find their number NNCV
*
        ICNT=0
        DO 30 I=1,NEIG
           IF (ICV(ISELEC(I)).EQ.0) THEN
*             ..Find coefficient and test for convergence
              ICUR=KPASS*ISELEC(I)
              TMAX=ABS( SVEC(ICUR) )
              DO 20 L=1,NNCV-1
 20              TMAX=MAX( TMAX, ABS(SVEC(ICUR-L)) )
              IF (TMAX.LT.CRITC) THEN
*                ..this  coefficient converged
                 ICV(ISELEC(I))=1
              ELSE
*                ..Not converged. Add it to the list.
                 ICNT=ICNT+1
                 IND(ICNT)=ISELEC(I)
                 ROWLAST(ICNT)=TMAX
              ENDIF
           ENDIF
 30     CONTINUE

        NNCV=ICNT
        IF (NNCV.EQ.0) DONE=.TRUE.
*
* Sort the ROWLAST elements interchanging their indices as well
*
        DO 40 I=1,NNCV
           INDX=IDAMAX(NNCV-I+1,ROWLAST(I),1)
           INCV(I)=IND(INDX+I-1)

           TEMP=ROWLAST(INDX+I-1)
           ROWLAST(INDX+I-1)=ROWLAST(I)
           ROWLAST(I)=TEMP
           ITEMP=IND(INDX+I-1)
           IND(INDX+I-1)=IND(I)
           IND(I)=ITEMP
 40     CONTINUE

        TSTSEL=DONE
        RETURN
        END
