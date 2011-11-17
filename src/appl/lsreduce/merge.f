*-----------------------------------------------------------------------
*     M E R G E         written by Per Jonsson, Lund, April 1992
*-----------------------------------------------------------------------
*

      subroutine merge(nmr)

      implicit double precision(a-h,o-z)
      character*64 line1,line2,cmr1(10000),cmr2(10000)
      common/inform/iread,iwrite,iout,isc0,isc1,isc2,isc3,
     : iall,jsc(3),iscw
      logical endloop

9     format(a)
*
*  Determine where to start reading the CSFs in the multireference list
*
      nbegin1 = 0
      do
         read(15,9) line1
         if (line1(5:5).eq.'(') exit; 
         nbegin1 = nbegin1 + 1
      end do
*
*  Determine where to start reading the CSFs in the clist
*
      nbegin2 = 0
      do 
         read(17,9) line1
         if (line1(5:5).eq.'(') exit;
         nbegin2 = nbegin2 + 1
      end do
*
*  Read the header information in clist and write it to the new list (unit=16)
*
      rewind (17)

      do i = 1,nbegin2
         read(17,9) line1
         k = 64
 70      if (line1(k:k) .eq. ' ') then
            k = k-1
            if (k .gt. 1) go to 70
         end if
         write(iread,9) line1(1:k)
      end do
*
*  Goto right position and read the CSFs in the multiref
*  and write them as the first CSFs in the new list

      rewind(15)
      do k = 1,nbegin1
         read(15,9) line1
      end do

      ncfg1 = 1
      do 
         read(15,9,end=25) cmr1(ncfg1)
         if (cmr1(ncfg1)(1:1).ne.'*') then
            read(15,9) cmr2(ncfg1)
            ncfg1 = ncfg1 + 1
         end if 
      end do
25    continue

      nmr = ncfg1 - 1
      do i = 1,ncfg1 - 1
         line1 = cmr1(i)
         line2 = cmr2(i)
        k = LEN_TRIM(line1);
        write(iread,9) line1(1:k)
        k = LEN_TRIM(line2);
        write(iread,9) line2(1:k)         
      end do
*
*  Read the CSFs in clist. Check if they are equal to any CSF in the
*  multireference. If not write them to the new list.
*
      nfound = 0
50    read(17,9,end=60) line1
      if (line1(1:1).ne.'*') then
         read(17,9) line2
         do i = 1,ncfg1
!            print*, 'cmr1(i):::',cmr1(i),':::';
!            print*, 'line1  :::',line1,':::';
!            print*, 'cmr2(i):::',cmr2(i),':::';
!            print*, 'line2  :::',line2,':::';
            if (cmr1(i).eq.line1.and.cmr2(i).eq.line2) nfound = 1
!            if (cmr1(i).eq.line1) nfound = 1 
         end do
         if (nfound.ne.1) then
            k = 64
 78         if (line1(k:k) .eq. ' ') then
               k = k-1
               if (k .gt. 1) go to 78
            end if
            write(iread,9) line1(1:k) 
            k = 64
 79         if (line1(k:k) .eq. ' ') then
               k = k-1
               if (k .gt. 1) go to 79
            end if
            write(iread,9) line2(1:k) 
         endif
         nfound = 0
      goto 50
      endif

60    continue

      write(iread,9) '*'
      return
      end
