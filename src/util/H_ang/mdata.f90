!**********************************************************************
!
! Interface module for nonh (ATSP2K)
!
! Thomas Carette, Nov 2011
!
!     Notes:
!       - the Hang type is to fit for an optimal code runing for
!         big matrices(?)
!
!**********************************************************************

module mdata

!**********************************************************************
!
!    types
!

  type iang 
    real(kind(1d0))c ! angular int numerical value
                     ! = coefficient of radial int
    integer inptr    ! address of the int in the int list
    integer icase    ! 1:F, 2:G, 3:R, 4:I
                     !
                     ! radial integrals characteristics
!    integer ipack    ! packed code for the int (ATSP2K:ipackn,kval)
    integer k
    integer l1
    integer r1
    integer l2
    integer r2
  end type iang

  type Hang
    integer ni
    type(iang), allocatable :: int(:)
  end type Hang

  type bha
    type(hang), allocatable :: h(:,:)
  end type bha

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
      type(bha), allocatable :: tot_hang(:)

!**********************************************************************
!
!     interface variables
!

      integer lsdim

!**********************************************************************



contains

      INCLUDE 'data.f90'
      INCLUDE 'spintgrl.f90'

      character*72 function intfmt(ic)
        integer ic
        character*9 s
        character*30 st1
        character*50 st2
        character*28 st3

        s='a3,",",a3'
        select case (ic)
          case (1)
            st1='"F",i2,"(",'//s//',")"'
            intfmt=st1
          case (2)
            st1='"G",i2,"(",'//s//',")"'
            intfmt=st1
          case (3)
            st2='"R",i2,"(",'//s//',";",'//s//',")"'
            intfmt=st2
          case (4)
            st3='"I  (",'//s//',")"'
            intfmt=st3
        end select
      end function intfmt

      INCLUDE 'unsparse.f90'

end module mdata
