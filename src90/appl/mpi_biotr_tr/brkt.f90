!
!     ------------------------------------------------------------------
!     B R K T    (for "bra-ket")
!     ------------------------------------------------------------------
!
      subroutine brkt 
      use fo_C
      use ras_C
      use nel_C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:50:09  11/18/01  
!...Switches:                     

      dimension :: braket(nwd,nwd) 
      do i = 1, iwf(1) 
         do j = 1, iwf(2) 
            braket(i,j) = quadr(i,j + iwf(1),0) 
            if (l(i) .eq. l(j+iwf(1))) write (6, *) '<', elrasi(i), '|', elrasf&
               (j), '> = ', braket(i,j) 
         end do 
      end do 
      return  
      end subroutine brkt 
