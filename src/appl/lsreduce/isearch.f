*-----------------------------------------------------------------------
*               I S E A R CH
*-----------------------------------------------------------------------
*
*     Given the type of integral, the packed integral,
*     find its position in the list ipackn.  An open
*     addressing scheme will be used.
*
      INTEGER FUNCTION isearch(icase,int,qpackn,qlused,qintptr,lmax)
      POINTER (qpackn,ipackn(1)),(qintptr,intptr(0:2*lmax,4)),
     :        (qlused,lused(1))
      
      LOGICAL found
      LOGICAL lused
*     Find k
      k = int
      if (icase.le.2 .or. icase.eq.4 .or. icase.eq.5) then
ctc  NWD limited to 63
ctc        k = k/64/64
        k = k/71/71
ctc
      else if (icase .eq. 3 .or. icase .gt. 5) then
ctc  NWD limited to 63
ctc        k = k/64/64/64/64
        k = k/71/71/71/71
ctc
      end if
*     Find range of integral to search
      if (k .eq. 0) then
        if (icase .eq. 1) then
          left = 0
        else 
          left = intptr(2*lmax,icase-1)
        end if
      else
        left = intptr(k-1,icase)
      end if
      iright = intptr(k,icase)+1
*
*     .. begin searching the appropriate range
      found = .false.
10    mid = (left + iright)/2    
      if (left .lt. mid .and. .not. found) then
        if (int .eq. ipackn(mid)) then
          isearch = mid
          found = .true.
          lused(isearch) = .true.
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
        stop
      end if
      end






