*-----------------------------------------------------------------------
*        Q S O R T
*-----------------------------------------------------------------------

*     The method use to sort the data is quick sort with a pivot value
*     that is the larger value of the first 2 different value from the
*     the sublist to be sorted.
*     This sorting method used a stack to maintain the unsorted section,
*     and sorting will be finished when the stack is empty.
      subroutine qsort(n,key,pt,stack,nstack,ierr)
      integer top, from, to, pivot, left, right
      integer key(*),pt(*),stack(*)
*
*     Set the initial pointer values
*
      do 10 i=1,n
10    pt(i)=i
*
*     Initialize the stack and error indicator
*
      top=1
      stack(top)=1
      ierr = 0
*
*     Repeat Until the Stack is empty
*
100   continue
*
*     Determine the next section
*
      if (top .ne. 0) then
        from=stack(top)
        if (top.ne.1) then
          to=stack(top-1) - 1
        else
          to=n
        endif
        top=top-1
*
*        Find the position k of the pivot value that partitions
*        the current section. Return a value of k=0 when there is
*        no disinct value.
*
        if (from .eq. to) then
          k=0
        else
          k=0
          ismall=key( pt(from) )
          do 210 i=from+1,to
            if (key( pt(i) ) .ne. ismall) then
              k=i
              goto 200
            endif
210       continue
        endif
200     continue
        if (k.ne.0) then
          if( ismall .gt. key(pt(k)) ) then
            k=from
          endif
        endif
        if (k .ne. 0)then
*
*         Rearange the section of the keys such that all values
*         smaller than the pivot value are stored on the left and
*         larger values on the right.
*
          pivot=key(pt(k))
          left=from
          right=to
300       continue
310       if ( key(pt(left)) .lt. pivot ) then
            left = left+1
            goto 310
          endif
320       if ( key(pt(right)) .ge. pivot ) then
            right = right-1
            goto 320
          endif
          if (left .lt. right) then
            i=pt(left)
            pt(left)=pt(right)
            pt(right)=i
            goto 300
          endif
          kp=left
          if (top+2 .le. nstack) then
            stack(top+1)=kp
            stack(top+2)=from
            top=top+2
          else
            ierr =1
            return
          endif
        endif
        goto 100
      endif
*
*     Keys are sorted
*
      return
      end
