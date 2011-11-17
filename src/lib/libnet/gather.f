      subroutine gather(n,a,b,index)
C
C     This subroutine collects array elements accessed via an
C         integer pointer to contiguous storage.
C
      INTEGER N,INDEX(1)
      REAL A(1),B(1)

      DO 10 I = 1,N
			A(I) = B(INDEX(I))
   10 CONTINUE

      RETURN
      END
