*
*     ------------------------------------------------------------------
*	B W I N T
*     ------------------------------------------------------------------
*
      SUBROUTINE BWINT(LC,LO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/BLUME/COEFN2(4),COEFNK(4),COEFVK(4)
*
* ... LC IS THE L-VALUE OF THE FILLED SUBSHELL, LO IS THE L-VALUE
*     OF THE PARTIALLY-FILLED SUBSHELL.
*
      IF(LC.LE.3.AND.LO.LE.4) GO TO 1
      PRINT 100, LC,LO
  100 FORMAT (37H INCORRECT CALLING OF BWINT WITH LC =,I2,6H, LO =,I2)
    1 LC1 = LC + 1
      GO TO (10,20,30,40), LC1
   10 GO TO (11,12,13,14), LO
*
* ... S-P
*
   11 COEFNK(1) = 1.D0
      COEFN2(1) = -2.D0
      COEFVK(1) = 1.D0
      RETURN
*
* ... S-D
*
   12 COEFNK(1) = 6.D0/5.D0
      COEFN2(1) = -9.D0/5.D0
      COEFVK(1) = 3.D0/5.D0
      RETURN
*
* ... S-F
*
   13 COEFNK(1) = 9.D0/7.D0
      COEFN2(1) = -12.D0/7.D0
      COEFVK(1) = 3.D0/7.D0
      RETURN
*
* ... S-G
*
   14 COEFNK(1) = 4.D0/3.D0
      COEFN2(1) = -5.D0/3.D0
      COEFVK(1) = 1.D0/3.D0
      RETURN
   20 GO TO (21,22,23,24), LO
*
* ... P-P
*
   21 COEFNK(1) = 0.D0
      COEFN2(1) = 3.D0
      COEFVK(1) = 9.D0/5.D0
      RETURN
*
* ... P-D
*
   22 COEFNK(1) = 3.D0/7.D0
      COEFNK(2) = 36.D0/35.D0
      COEFN2(1) = -12.D0/5.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = 3.D0/5.D0
      COEFVK(2) = 36.D0/35.D0
      RETURN
*
* ... P-F
*
   23 COEFNK(1) = 1.D0/7.D0
      COEFNK(2) = 10.D0/7.D0
      COEFN2(1) = -18.D0/7.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = 18.D0/35.D0
      COEFVK(2) = 5.D0/7.D0
      RETURN
*
* ... P-G
*
   24 COEFNK(1) = 5.D0/77.D0
      COEFNK(2) = 18.D0/11.D0
      COEFN2(1) = -18.D0/7.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = 3.D0/7.D0
      COEFVK(2) = 6.D0/11.D0
      RETURN
   30 GO TO (31,32,33,34), LO
*
* ... D-P
*
   31 COEFNK(1) = 59.D0/7.D0
      COEFNK(2) = -18.D0/7.D0
      COEFN2(1) = -4.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = -1.D0
      COEFVK(2) = 18.D0/7.D0
      RETURN
*
* ... D-D
*
   32 COEFNK(1) = 6.D0/7.D0
      COEFNK(2) = 0.D0
      COEFN2(1) = 3.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = 3.D0/7.D0
      COEFVK(2) = 10.D0/7.D0
      RETURN
*
* ... D-F
*
   33 COEFNK(1) = 9.D0/7.D0
      COEFNK(2) = -13.D0/77.D0
      COEFNK(3) = 75.D0/77.D0
      COEFN2(1) = -18.D0/7.D0
      COEFN2(2) = 0.D0
      COEFN2(3) = 0.D0
      COEFVK(1) = 3.D0/7.D0
      COEFVK(2) = 3.D0/7.D0
      COEFVK(3) = 75.D0/77.D0
      RETURN
*
* ... D-G
*
   34 COEFNK(1) = 741.D0/693.D0
      COEFNK(2) = -215.D0/429.D0
      COEFNK(3) = 210.D0/143.D0
      COEFN2(1) = -3.D0
      COEFN2(2) = 0.D0
      COEFN2(3) = 0.D0
      COEFVK(1) = 3.D0/7.D0
      COEFVK(2) = 255.D0/693.D0
      COEFVK(3) = 105.D0/143.D0
      RETURN
   40 GO TO (41,42,43,44), LO
*
* ... F-P
*
   41 COEFNK(1) = 52.D0/3.D0
      COEFNK(2) = -20.D0/3.D0
      COEFN2(1) = -9.D0
      COEFN2(2) = 0.D0
      COEFVK(1) = -9.D0/5.D0
      COEFVK(2) = 10.D0/3.D0
      RETURN
*
* ... F-D
*
   42 COEFNK(1) = 5.D0
      COEFNK(2) = 142.D0/55.D0
      COEFNK(3) = -20.D0/11.D0
      COEFN2(1) = -18.D0/5.D0
      COEFN2(2) = 0.D0
      COEFN2(3) = 0.D0
      COEFVK(1) = -3.D0/5.D0
      COEFVK(2) = 2.D0/5.D0
      COEFVK(3) = 20.D0/11.D0
      RETURN
*
* ... F-F
*
   43 COEFNK(1) = 1.D0
      COEFNK(2) = 5.D0/11.D0
      COEFNK(3) = 0.D0
      COEFN2(1) = 3.D0
      COEFN2(2) = 0.D0
      COEFN2(3) = 0.D0
      COEFVK(1) = 1.D0/5.D0
      COEFVK(2) = 5.D0/11.D0
      COEFVK(3) = 175.D0/143.D0
      RETURN
*
* ... F-G
*
   44 COEFNK(1) = 53.D0/33.D0
      COEFNK(2) = 57.D0/143.D0
      COEFNK(3) = -115.D0/429.D0
      COEFNK(4) = 392.D0/429.D0
      COEFN2(1) = -8.D0/3.D0
      COEFN2(2) = 0.D0
      COEFN2(3) = 0.D0
      COEFN2(4) = 0.D0
      COEFVK(1) = 1.D0/3.D0
      COEFVK(2) = 3.D0/11.D0
      COEFVK(3) = 57.D0/143.D0
      COEFVK(4) = 392.D0/429.D0
      RETURN
      END
