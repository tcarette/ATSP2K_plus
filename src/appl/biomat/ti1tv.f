*
*     ------------------------------------------------------------------
*	T I 1 T V  
*     ------------------------------------------------------------------
*
      SUBROUTINE TI1TV(CIIN,NCSF,NCIV,I,L,T,NSHL,JLEI,LBUF,LU,
     &                 CIOUT,BUF,IBUF,NTESTG)
* ciout=scr
* Calculate the action of the  operator
*     Sum( j) T(jl) E(jl,il) on a set of vectors
*
* The summation can be restrited by the use of JLEI
* JLEI = 0 => no restriction
* JLEI = 1 => only j<i terms included
* JLEI = 2 => only j=i terms included
* JLEI = 3 => only j<i and j=i terms included
*
* Written in Bruxelles May 1993 for the BIOTRN program
* Correction : December 1993
*
* =====
* Input
* =====
* CIIN : Input CI vectors
* NCSCF : Length of CI expansion
* NCIV  : Number of CI vectors
* I     : Shell number
* L     : L value of shells of this excitation
* NSHL  : Number of shells with this L
* T     : The integrals T(*,il)
* LBUF  : Number of one-electron excitations fetched
*       : per call
* LU    : Unit number for Racah coefficients
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
      DIMENSION CIIN(NCSF,NCIV),T(NSHL)
*. Output
      DIMENSION CIOUT(NCSF,NCIV)
*. Scratch
      DIMENSION BUF(LBUF),IBUF(4,LBUF)
*
      NTESTL = 00
      NTEST = MAX(NTESTL,NTESTG)
*
      IF(NTEST.GE.10) WRITE(6,*) ' Entering TI1TV'
      CALL SETVEC(CIOUT,0.0D0,NCIV*NCSF)
*.    Loop over batches of coupling coefficients
*.    involving shell I
      IBATCH=0
 900  CONTINUE
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
         CALL GTRAC1
     &   (I,L,IFIRST,NFOUND,BUF,IBUF,LBUF,JLEI,LU,NTESTG )
*. Expected return from GTRC1
*        NFOUND : Number of coefficients obtained
*        NFOUND = 0 indicates all coefficients involving this i
*                   has been found
*        BUF    : actual RACAH coefficient <CSF(L)!E(jl,il)!CSF(R)>
*        IBUF(1,*) j of E(jl,il)
*        IBUF(2,*) i of E(jl,il)
*        IBUF(3,*) R
*        IBUF(4,*) L
*
         DO 800 IELMNT = 1, NFOUND
          RACAH = BUF(IELMNT)
          J = IBUF(1,IELMNT)
          X = T(J)*RACAH
Cmrg      ILEFT = IBUF(4,IELMNT)
Cmrg      IRIGHT = IBUF(3,IELMNT)
          ILEFT = IBUF(3,IELMNT)
          IRIGHT = IBUF(4,IELMNT)
          DO 100 IVEC =1, NCIV
           CIOUT(ILEFT,IVEC) =  CIOUT(ILEFT,IVEC) + X*CIIN(IRIGHT,IVEC)
  100     CONTINUE
  800   CONTINUE
*
      IF(NFOUND.GT.0) GOTO 900
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Input and output vectors from TI1TV'
        CALL WRTMAT(CIIN,NCSF,NCIV,NCSF,NCIV)
        WRITE(6,*)
        CALL WRTMAT(CIOUT,NCSF,NCIV,NCSF,NCIV)
        WRITE(6,*)
      END IF
*
      IF(NTEST.GE.10) WRITE(6,*) ' LEAVING TI1TV'
      RETURN
      END
