      subroutine  ifill(n,ialpha,ix,incx)
c
c     fills a vector, x, with the scalar value of alpha
c
      integer i,ialpha,ix(1),incx,n
c
      if(n.le.0)return
      if ( incx .eq. 1) then
        do 10 i = 1,n
 	  ix(i) = ialpha
   10   continue
      else
	jx = 1
	do 20 i = 1,n
	  ix(jx) = ialpha
	  jx = jx + incx
   20   continue
      end if
      return
      end
