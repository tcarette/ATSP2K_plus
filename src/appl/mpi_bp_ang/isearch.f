*
*-----------------------------------------------------------------------
*               I S E A R CH
*-----------------------------------------------------------------------
*
*     Given the type of integral, the packed integral,
*     find its position in the list ipackn.  An open
*     addressing scheme will be used.
*
      INTEGER FUNCTION isearch(icase,int,qpackn,qintptr,lmax)
      POINTER (qpackn,ipackn(1)),(qintptr,intptr(0:2*lmax+1,7))
      
      LOGICAL found
*
*     Find k
*       Note that Sk (icase=8) is the same integral as Nk (icase=6)
*       The packed value for Nk is k = K+1 (K, the real k)
*       jcase introduced to avoid side-effects

      if (icase .eq. 8) then
         jcase = 6
      else if (icase .eq. 9) then
         jcase = 6
      else
         jcase = icase
      end if
      k = int
      if (jcase.le.2 .or. jcase.eq.4 .or. jcase.eq.5) then
        k = k/64/64
      else if (jcase .eq. 3 .or. jcase .gt. 5) then
        k = k/64/64/64/64
      end if

*     Find range of integral to search

      if (k .eq. 0) then
        if (jcase .eq. 1) then
          left = 0
        else 
          left = intptr(2*lmax+1,jcase-1)
        end if
      else
        left = intptr(k-1,jcase)
      end if
      iright = intptr(k,jcase)+1
*
*     .. begin searching the appropriate range

      found = .false.
10    mid = (left + iright)/2    
      if (left .lt. mid .and. .not. found) then
C        write(*,*) ipackn(mid)
        if (int .eq. ipackn(mid)) then
          isearch = mid
          found = .true.
        else if (int .gt. ipackn(mid)) then
          left = mid
        else
          iright = mid
        end if
        go to 10
      end if
      if (.not. found) then
        write(6,*) ' The integral', int,' for icase = ',icase,
     :             ' not found in the list'
	write(6,*) 'left,iright', left,iright
	write(6,*) 'k', k
        stop
      end if

      end

