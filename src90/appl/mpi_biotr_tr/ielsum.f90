!
!     ------------------------------------------------------------------
!     I E L S U M
!     ------------------------------------------------------------------
!
      INTEGER FUNCTION IELSUM (IVEC, NELMNT) 
!
! Sum elements of integer array
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  17:00:23  11/18/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NELMNT 
      INTEGER , INTENT(IN) :: IVEC(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ISUM, IEL 
!-----------------------------------------------
!
      ISUM = 0 
!?    WRITE(6,*) ' IELSUM, NELMNT ', NELMNT
      ISUM = SUM(IVEC(:NELMNT)) 
!?     write(6,*) ' IELSUM IEL IVEC ISUM '
!?     write(6,*) IEL,IVEC(IEL),ISUM
!
      IELSUM = ISUM 
!
      RETURN  
      END FUNCTION IELSUM 
