      subroutine  dfill(n,alpha,dx,incx)
c
c     fills a vector, x, with the scalar value of alpha
c
      double precision dx(1), alpha
      integer i,incx,n
c
      if(n.le.0)return
      if ( incx .eq. 1) then
        do 10 i = 1,n
 	  dx(i) = alpha
   10   continue
      else
	jx = 1
	do 20 i = 1,n
	  dx(jx) = alpha
	  jx = jx + incx
   20   continue
      end if
      return
      end
