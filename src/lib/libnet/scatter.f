      subroutine scatter(n,a,index,b)
C
C     This subroutine disperses array elements accessed via an
C         integer pointer from contiguous storage to the appropriate
C         location.
C
      INTEGER N,INDEX(1)
      REAL A(1),B(1)

      DO 10 I = 1,N
			A(INDEX(I)) = B(I)
   10 CONTINUE

      RETURN
      END
