*
*     ------------------------------------------------------------------
*       I N T C F G
*     ------------------------------------------------------------------
*
*     This routine writes the clist file containing only the CSFs
*     that interacts with the multireference set.
*
      subroutine intcfg
      implicit double precision(a-h,o-z)

      pointer(qmrsdint,mrsdint(1))
      common/mrsdcom/qmrsdint
      common/inform/iread,iwrite,iout,isc0,isc1,isc2,isc3,
     : iall,jsc(3),iscw
      
      character line1*74
      character line2*74
*
      rewind(iread)
      open(unit=25,file='clist.out',status='unknown')
*
      nbegin = 0
5     read(iread,'(a)') line1
      if (line1(5:5).ne.'(') then
         nbegin = nbegin + 1
      goto 5
      endif
      
      rewind(iread)

      do i = 1,nbegin
         read(iread,'(a)') line1
         k = 74
 76      if (line1(k:k) .eq. ' ') then
       	    k = k-1
            if (k .gt. 1) go to 76
         end if
         write(25,'(a)') line1(1:k)
      end do

      i = 0
      j = 0

*  Read a CSF with the coupling information and check if it should be saved

 77   read(iread,'(a)',end=99) line1
      if (line1(1:1).eq.'*') then
         write (25,'(A1)') '*'
         goto 99
      end if
      read(iread,'(a)') line2
      i = i + 1

*  Save CSF if mrsdint > 0

      if (mrsdint(i).gt.0) then
         j = j + 1
         k = 74
 78      if (line1(k:k) .eq. ' ') then
            k = k-1
            if (k .gt. 1) go to 78
         end if
         write(25,'(a)') line1(1:k)
         k = 74
 79      if (line2(k:k) .eq. ' ') then
            k = k-1
            if (k .gt. 1) go to 79
         end if
         write(25,'(a)') line2(1:k)
      end if
      go to 77
 
99    endfile 25
      
      close(unit=25)
      write(iscw,*)
      write(iscw,*) ' Number of CSFs in the reduced list  ',j
      return
      end
