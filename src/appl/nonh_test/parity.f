*-----------------------------------------------------------------------
*     P A R I T Y       written by Thomas Caretten, Stockholm, Nov 2011
*-----------------------------------------------------------------------
      character*1 function parity(line)

      character*64 line

      character*1 lstr(8)
      integer occ(8)

    5 FORMAT(8(3X,A1,1X,I2,1X))

      N=len_trim(line)/8
      READ(LINE,5)(lstr(J),occ(J),J=1,N)
      ip = 0
      Do iip = 1,n
        iq   = occ(iip)
        ip = ip + lval(lstr(iip))*iq
      end do
      if ((ip/2)*2 .eq. ip) then
        parity ='e'
      else
        parity ='o'
      end if

      end function
