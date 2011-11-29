!**********************************************************************
!
! Interface module for nonh (ATSP2K)
!
! Thomas Carette, Nov 2011
!
!**********************************************************************

module mdata

!**********************************************************************
!
!     orbital variables
!

      INTEGER, PARAMETER :: NWD=70
      INTEGER :: nclosd,nowf,nwf,lmax
      CHARACTER EL(nwd)*3
      INTEGER, ALLOCATABLE :: ossn(:),ossl(:)
 

!**********************************************************************
!
!    CSF list variables
!

!         Block and global variables

      INTEGER nterm, nblock
      integer, ALLOCATABLE :: nze_bl(:), ncfg_bl(:), niv_bl(:), &
                              nze_max(:),ncn_bl(:)
      integer nze_tot,ncfg_tot,ncn_tot,nze_col_max
      character, ALLOCATABLE :: term_bl(:)*3

!         CSFs variables

      integer, allocatable :: jptr(:),jan(:),ico(:)

!**********************************************************************
!
!     integral variables
!

      LOGICAL, ALLOCATABLE :: lused(:)
      INTEGER, ALLOCATABLE :: ipackn(:),inptr(:)
      integer nint,intptr(4)

!     Integral coefficients

      REAL(kind(1d0)), allocatable :: cn(:)

!**********************************************************************
!
!     interface variables
!

      integer lsdim

!**********************************************************************



contains

      INCLUDE 'data.f90'
      INCLUDE 'spintgrl.f90'

end module mdata
