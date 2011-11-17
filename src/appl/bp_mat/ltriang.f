*
*------------------------------------------------------------------------
*        L T R I A N G
*------------------------------------------------------------------------
*
	LOGICAL FUNCTION ltriang(k,li,lj)
	INTEGER k,li,lj
*
*       Return a value of TRUE if k,li,lj form a triangle
*
	IF ( mod(k+li+lj,2) .ne. 0) then
	   ltriang = .false.
	else
	  IF (abs(li-lj) .gt. k ) then
	    ltriang = .false.
	  else if (li+lj .lt. k) then
	    ltriang = .false.
	  else
	    ltriang = .true.
	  end if
	end if
	end
