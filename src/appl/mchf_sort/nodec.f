*
*     ------------------------------------------------------------------
*    3-21      N O D E C
*     ------------------------------------------------------------------
*
*      Counts the number of nodes of the function PDE(j) in the range
*   j = 40,...,M-10.   The node counting procedure counts the local max
*   and min values.   Only nodes between sufficiently large max and
*   min values are counted.
*
*
      FUNCTION NODEC(M)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
        PARAMETER (NWD=60,NOD=220,NOFFD=800)
        COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
*
*   ***** FIND MAX|PDE(J)|
*
       MM = M - 10
      DM = 0.D0
      DO 1 J = 40,MM
1     DM = DMAX1(DM, DABS(PDE(J)))
*
*   *****  COUNT THE NUMBER OF LOCAL MAX OR MIN'S
*
      NCC = 0
      SIGN = 0.D0
      DIFF1 = PDE(40) - PDE(39)
      DO 2 J = 40, MM
      DIFF2 = PDE(J+1) - PDE(J)
      IF (DIFF2*DIFF1 .GT. 0.D0 .OR. DIFF1 .EQ. 0.D0) GO TO 2
*
*   *****  A MAX OR MIN HAS BEEN FOUND.   TEST IF IT IS
*          SUFFICIENTLY LARGE
*
      IF ( DABS(PDE(J))/DM .LT. 0.05D0 ) GO TO 2
*
*   ***** CHECK IF THIS IS THE FIRST SIGNIFICANT MAXIMUM
*
      IF (SIGN .NE. 0.D0 ) GO TO 4
      M = J
      GO TO 3
*
*   ***** IF NOT THE FIRST, TEST WHETHER A SIGN CHANGE HAS
*         OCCURRED SINCE THE LAST SIGNIFICANT MAX OR MIN
*
4     IF (PDE(J)*SIGN .GT. 0.D0 ) GO TO 2
      NCC = NCC + 1
*
*   ***** RESET FOR THE NEXT NODE
*
3     SIGN = PDE(J)
2     DIFF1 = DIFF2
      NODEC = NCC
      RETURN
      END
