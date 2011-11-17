*
*     ------------------------------------------------------------------
*
*       T R I T S T
*     ------------------------------------------------------------------
*
      DOUBLE PRECISION FUNCTION TRITST(L,M,N)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
*
*      IF  TRITST=1.0   THE TRIANGLE RELATION IS NOT SATISFIED
*      IF  TRITST=0.0   THE TRIANGLE RELATION IS SATISFIED
*
      LMN=IABS(L-M)
      LM=L+M
      IF(N-LMN) 1,2,2
    2 IF(LM-N) 1,3,3
    3 TRITST=0.D0
      RETURN
    1 TRITST=1.D0
      RETURN
      END
