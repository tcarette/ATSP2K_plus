*
*     ------------------------------------------------------------------
*	R S C H E C K
*     ------------------------------------------------------------------
*
*. check if the order found by nonh is consistent with the built-in 
*  ras order: we want to use the (j.leq.i) restriction adopted in savels.
*
      subroutine rscheck(nq,lq,lmx,mx)
      implicit double precision(a-h,o-z)
      parameter (NWD=128)
      dimension nq(nwd),lq(nwd)
*
      do 1 l = 0,lmx
        nst = 0
	do 2 i = 1,mx
	  if(lq(i) .ne. l) go to 2
	  if(nq(i) .lt. nst) then
	     print*,' Problem: orbital order found by nonh inconsistent'
             print*,'          with the built-in ras order in savels!'
	     print*,' Sorry about that. '
	     stop
          end if
	  print*,' l = ',l,' i = ',i,' lq(i) = ',lq(i),' nq(i) = ',nq(i)
          nst = nq(i)
    2   continue
    1 continue
      print*
      print*,' rscheck report: Congratulations: RAS test satisfied.'
      print*
      return
      end
