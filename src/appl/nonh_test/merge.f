*-----------------------------------------------------------------------
*     M E R G E         written by Per Jonsson, Lund, April 1992
*
*             Modified by Thomas Carette, Stockholm, Nov 2011
*
*-----------------------------------------------------------------------
*

      subroutine merge(lmax,term)

      implicit double precision(a-h,o-z)



      character*64 line1,line2
      common/inform/iread,iwrite,iout,isc0,isc1,isc2,isc3,
     : iall,jsc(3),iscw
      character*3 term
      character*1 parity

9     format(a)

      ireadp=iread+1
      ireadr=iread+2

*
*  Determine where to start reading the CSFs in the multireference list
*
      nbegin1 = 0
      do
         read(ireadr,9) line1
         if (line1(5:5).eq.'(') exit; 
         nbegin1 = nbegin1 + 1
      end do
*
*  Determine where to start reading the CSFs in the clist
*
      nbegin2 = 0
      do 
         read(ireadp,9) line1
         if (line1(5:5).eq.'(') exit;
         nbegin2 = nbegin2 + 1
      end do
*
*  Read the header information in clist and write it to the new list (unit=16)
*
      rewind (ireadp)

      do i = 1,nbegin2
         read(ireadp,9) line1
         k = 64
 70      if (line1(k:k) .eq. ' ') then
            k = k-1
            if (k .gt. 1) go to 70
         end if
         write(iread,9) line1(1:k)
      end do
*
*  Goto right position and read the CSFs in the multiref
*  If the final LS coupling is the one considered,
*  write them as the first CSFs in the new list

      rewind(ireadr)
      do k = 1,nbegin1
         read(ireadr,9) line1
      end do

      do 
         read(ireadr,9,end=25) line1
         if (line1(1:1).ne.'*') then
           read(ireadr,9) line2
           if(parity(line1).eq.term(3:3))then
              k2 = LEN_TRIM(line2);
              if ( line2(k2-1:k2) .eq. term(1:2))then
                k1 = LEN_TRIM(line1);
  
                write(iread,9) line1(1:k1)
                write(iread,9) line2(1:k2)         
              endif
            endif
         end if 
      end do
25    continue

*
*  Read the CSFs in clist (parent list).
*  Couple the extra/excited electron.
*  Write the result to the new list.
*

      call coupex(lmax,term)

      write(iread,9) '*'
      return
      end
