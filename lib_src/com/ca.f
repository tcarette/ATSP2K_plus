*
*     ------------------------------------------------------------------
*          C A
*     ------------------------------------------------------------------
*
*
      DOUBLE PRECISION FUNCTION CA(L,K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /EAV/CCA(10),CCB(35)
*
      IF (L .LE. 4) THEN
         CA = CCA((L*(L-1) + K)/2)
       ELSE
         CA = RME(L,L,K)**2
      END IF
      END
