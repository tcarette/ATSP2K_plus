*
*     ------------------------------------------------------------------
*	P A M T M T
*     ------------------------------------------------------------------
*
      SUBROUTINE PAMTMT(X,T,WORK,NORB)
*
* GENERATE PER AKE'S T MATRIX FROM A
* ORBITAL ROTATION MATRIX X
*
* T IS OBTAINED AS A STRICTLY LOWER TRIANGULAR
* MATRIX TL AND AN UPPER TRIANGULAR MATRIX TU
*
*         TL = 1 - L
*         TU = U ** -1
*
* WHERE L AND U ARISES FROM A LU DECOMPOSITION OF
* X :
*         X = L * U
* WITH L BEING A LOWER TRIANGULAR MATRIX WITH UNIT ON THE
* DIAGONAL AND U IS AN UPPER TRIANGULAR MATRIX
*
* JEPPE OLSEN OCTOBER 1988
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(NORB,NORB),T(NORB,NORB)
      DIMENSION WORK(*)
* DIMENSION OF WORK : NORB ** 2 + NORB*(NORB+1) / 2
*
      NTEST = 0
*. Allocate local memory
      KLFREE = 1
C     KLL = KFLREE
      KLL = KLFREE
      KLFREE = KLL + NORB*(NORB+1)/2
      KLU = KLFREE
      KLFREE = KLU + NORB ** 2
*.LU factorize X
      CALL LULU(X,WORK(KLL),WORK(KLU),NORB)
*.Expand U to full matrix
      CALL SETVEC(T,0.0D0,NORB ** 2 )
      DO 10 I = 1,NORB
      DO 10 J = I,NORB
        T(I,J) = WORK(KLU-1+J*(J-1)/2+I)
   10 CONTINUE
      IF ( NTEST .GE. 10 ) THEN
        WRITE(6,*) ' MATRIX TO BE INVERTED '
        CALL WRTMAT(T,NORB,NORB,NORB,NORB)
      END IF
*.Invert U
      CALL INVMAT(T,WORK(KLU),NORB,NORB)
      IF ( NTEST .GE. 10 ) THEN
        WRITE(6,*) ' INVERTED MATRIX '
        CALL WRTMAT(T,NORB,NORB,NORB,NORB)
      END IF
*.Subtract L
      DO 20 I = 1, NORB
      DO 20 J = 1,I-1
       T(I,J)= - WORK(KLL-1+I*(I-1)/2+J)
   20 CONTINUE
*
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' INPUT X MATRIX '
        CALL WRTMAT(X,NORB,NORB,NORB,NORB)
        WRITE(6,*) ' T MATRIX '
        CALL WRTMAT(T,NORB,NORB,NORB,NORB)
      END IF
*
      RETURN
      END
