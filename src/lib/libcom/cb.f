*
*     -----------------------------------------------------------------
*          C B
*     -----------------------------------------------------------------
*
*
      DOUBLE PRECISION FUNCTION CB(L,LP,K)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /EAV/CCA(10),CCB(35)
      INTEGER ICBPTR(0:4)
      DATA ICBPTR/1,6,14,23,31/
*
             IF (L .LE. LP) THEN
                 L1 = L
                 L2 = LP
              ELSE
                 L1 = LP
                 L2 = L
             END IF
             IF ( L2 .LE. 4) THEN
                CB = CCB(ICBPTR(L1)+(K+L1-L2)/2+(L1+1)*(L2-L1))
               ELSE
                CB = RME(L,LP,K)**2/(2*(2*L+1)*(2*LP+1))
             END IF
      END
