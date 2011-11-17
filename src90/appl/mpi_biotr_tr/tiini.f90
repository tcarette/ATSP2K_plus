!
!     ------------------------------------------------------------------
!     T I I N I
!     ------------------------------------------------------------------
!
      SUBROUTINE TIINI(CIIN, NCSF, NCIV, I, L, CONST, LBUF, LU, CIOUT, BUF, &
         IBUF, NTESTG) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!
! Calculate the action of the  operator
!     Const ** E(li,li) on a set of vectors
!
! Written in Bruxelles May 1993 for the BIOTRN program
!
! =====
! Input
! =====
! CIIN : Input CI vectors
! NCSCF : Length of CI expansion
! NCIV  : Number of CI vectors
! I     : Shell number
! L     : Lvalue of excitations
! NSHL  : Number of shells of this L
! CONST : The constant
! LBUF  : Number of one-electron excitations fetched
!       : per call
! LU    : Unit number for accessing Racah coefficients
!
! ======
! Output
! ======
!
! CIOUT : List of output CI vectors
!
! =======
! Scratch
! =======
!
! BUF  : LBUF
! IBUF : 4*LBUF
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:14:32  11/20/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE setvec_I 
      USE wrtmat_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NCSF 
      INTEGER  :: NCIV 
      INTEGER  :: I 
      INTEGER  :: L 
      INTEGER  :: LBUF 
      INTEGER  :: LU 
      INTEGER  :: NTESTG 
      REAL(DOUBLE) , INTENT(IN) :: CONST 
      INTEGER  :: IBUF(4,LBUF) 
      REAL(DOUBLE)  :: CIIN(NCSF,NCIV) 
      REAL(DOUBLE)  :: CIOUT(NCSF,NCIV) 
      REAL(DOUBLE)  :: BUF(LBUF) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NTESTL, NTEST, IBATCH, IFIRST, NFOUND, JLEI, IELMNT, IVAL, &
         ILEFT, IVEC 
      REAL(DOUBLE) :: CONSTN 
!-----------------------------------------------
!. Input
!. Output
!. Scratch
!
      NTESTL = 0 
      NTEST = MAX(NTESTL,NTESTG) 
      IF (NTEST >= 10) WRITE (6, *) ' Entering TIINI' 
      CALL SETVEC (CIOUT, 0.0D0, NCSF*NCIV) 
!.
!.     Loop over batches of coupling coefficients
!.     involving shell I
      IBATCH = 0 
  900 CONTINUE 
      IBATCH = IBATCH + 1 
!
      IF (IBATCH == 1) THEN 
         IFIRST = 1 
      ELSE 
         IFIRST = IFIRST + NFOUND 
      ENDIF 
!. Get next lbuf coupling coefficients involving I,
!. Starting with element IFIRST
!
      JLEI = 2 
! =1 => j smaller then i
! =2 => j Equal   to   i
! =3 => j smaller then i
      CALL GTRAC1 (I, L, IFIRST, NFOUND, BUF, IBUF, LBUF, JLEI, LU, NTESTG) 
!. Expected return from GTRC1
!     NFOUND : Number of coefficients obtained
!     NFOUND = 0 indicates all coefficients involving this i
!                has been found
!     BUF    : actual RACAH coefficient <CSF(L)!E(jl,il)!CSF(R)>
!     IBUF(1,*) j of E(jl,il)
!     IBUF(2,*) i of E(jl,il)
!     IBUF(3,*) R
!     IBUF(4,*) L
!
      DO IELMNT = 1, NFOUND 
         IVAL = BUF(IELMNT) 
         CONSTN = CONST**IVAL 
         ILEFT = IBUF(3,IELMNT) 
         CIOUT(ILEFT,:NCIV) = CONSTN*CIIN(ILEFT,:NCIV) 
      END DO 
!
      IF (NFOUND > 0) GO TO 900 
!
!. The previous provided us with all
!  terms with nonvanishing occupation.
!  For terms with vanishing occupation of il,
!  just copy coefficients, since (x) ** 0 = 1
!
      WHERE (CIOUT(:NCSF,:NCIV) == 0.0D0)  
         CIOUT(:NCSF,:NCIV) = CIIN(:NCSF,:NCIV) 
      END WHERE 
!
      IF (NTEST >= 100) THEN 
         WRITE (6, *) 
         WRITE (6, *) ' Input and output vectors from TIINI' 
         CALL WRTMAT (CIIN, NCSF, NCIV, NCSF, NCIV) 
         WRITE (6, *) 
         CALL WRTMAT (CIOUT, NCSF, NCIV, NCSF, NCIV) 
         WRITE (6, *) 
      ENDIF 
!
      IF (NTEST >= 10) WRITE (6, *) ' Leaving  TIINI' 
      RETURN  
      END SUBROUTINE TIINI 
