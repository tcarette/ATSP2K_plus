*
*     -------------------------------------------------------------------
*          C A
*     ------------------------------------------------------------------
*
*
      DOUBLE PRECISION FUNCTION CA(L,K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /EAV/CCA(10),CCB(35)
*
      IF (k .eq. 0) then
	CA = 1.d0
      else if (L .LE. 4) THEN
        CA = CCA((L*(L-1) + K)/2)
      else
        CA = RME(L,L,K)**2/((2*L+1)*(4*L+1))
      end if
      END
