*
*     ------------------------------------------------------------------
*	T I I N I  
*     ------------------------------------------------------------------
*
      SUBROUTINE TIINI(CIIN,NCSF,NCIV,I,L,CONST,LBUF,LU,
     &                 CIOUT,BUF,IBUF,NTESTG)
*
* Calculate the action of the  operator
*     Const ** E(li,li) on a set of vectors
*
* Written in Bruxelles May 1993 for the BIOTRN program
*
* =====
* Input
* =====
* CIIN : Input CI vectors
* NCSCF : Length of CI expansion
* NCIV  : Number of CI vectors
* I     : Shell number
* L     : Lvalue of excitations
* NSHL  : Number of shells of this L
* CONST : The constant
* LBUF  : Number of one-electron excitations fetched
*       : per call
* LU    : Unit number for accessing Racah coefficients
*
* ======
* Output
* ======
*
* CIOUT : List of output CI vectors
*
* =======
* Scratch
* =======
*
* BUF  : LBUF
* IBUF : 4*LBUF
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*. Input
      DIMENSION CIIN(NCSF,NCIV)
*. Output
      DIMENSION CIOUT(NCSF,NCIV)
*. Scratch
      DIMENSION BUF(LBUF),IBUF(4,LBUF)
*
      NTESTL = 00
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.10) WRITE(6,*) ' Entering TIINI'
      CALL SETVEC(CIOUT,0.0D0,NCSF*NCIV)
*.
*.     Loop over batches of coupling coefficients
*.     involving shell I
       IBATCH=0
 900   CONTINUE
         IBATCH = IBATCH + 1
*
         IF(IBATCH.EQ.1) THEN
           IFIRST = 1
         ELSE
           IFIRST = IFIRST + NFOUND
         END IF
*. Get next lbuf coupling coefficients involving I,
*. Starting with element IFIRST
*
      JLEI = 2
* =1 => j smaller then i
* =2 => j Equal   to   i
* =3 => j smaller then i
      CALL GTRAC1(I,L,IFIRST,NFOUND,BUF,IBUF,LBUF,JLEI,LU,NTESTG)
*. Expected return from GTRC1
*     NFOUND : Number of coefficients obtained
*     NFOUND = 0 indicates all coefficients involving this i
*                has been found
*     BUF    : actual RACAH coefficient <CSF(L)!E(jl,il)!CSF(R)>
*     IBUF(1,*) j of E(jl,il)
*     IBUF(2,*) i of E(jl,il)
*     IBUF(3,*) R
*     IBUF(4,*) L
*
        DO 800 IELMNT = 1, NFOUND
          IVAL = BUF(IELMNT)
          CONSTN = CONST ** IVAL
          ILEFT = IBUF(3,IELMNT)
          DO 700 IVEC=1, NCIV
            CIOUT(ILEFT,IVEC) = CONSTN * CIIN(ILEFT,IVEC)
  700     CONTINUE
  800   CONTINUE
*
      IF(NFOUND.GT.0) GOTO 900
*
*. The previous provided us with all
*  terms with nonvanishing occupation.
*  For terms with vanishing occupation of il,
*  just copy coefficients, since (x) ** 0 = 1
*
      DO 300 IVEC = 1, NCIV
       DO 200 IELMNT = 1, NCSF
        IF(CIOUT(IELMNT,IVEC).EQ.0.0D0)
     &     CIOUT(IELMNT,IVEC) = CIIN(IELMNT,IVEC)
  200  CONTINUE
  300 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Input and output vectors from TIINI'
        CALL WRTMAT(CIIN,NCSF,NCIV,NCSF,NCIV)
        WRITE(6,*)
        CALL WRTMAT(CIOUT,NCSF,NCIV,NCSF,NCIV)
        WRITE(6,*)
      END IF
*
      IF(NTEST.GE.10) WRITE(6,*) ' Leaving  TIINI'
      RETURN
      END
