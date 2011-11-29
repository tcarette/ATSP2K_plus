!
!     ------------------------------------------------------------------
!    3-3       D A T A
!
!                from ATSP2K [CPC 176 (2007) 559]
!                modfied November 2011
!     ------------------------------------------------------------------
!
!       Data concerning the number of configurations (NCFG), the number
!   and type of electrons in each  configuration,  as  well  as  data
!   associated with the energy expression are read and stored.
!
!
      SUBROUTINE DATA

      use safeio

      IMPLICIT NONE
!
      CHARACTER*3   STRING*64
      INTEGER fp,ierr
      INTEGER nb,J,i
! external
      INTEGER lval
!

      fp=safe_open('cfg.h',ST='O',FO='F',AC='R')

      rewind(fp) 
      read(fp,'(i4,A64)') nclosd, string
      read(fp,'(18(1x,A3))') (el(i), i=1,nclosd)
      read(fp,'(i4,A64)') nowf, string
      nwf = nowf + nclosd
      IF(NWF .GT. NWD) THEN
         write(0,'(A,I4)') 'NWF too large for this code: MAX = ',NWD
         stop
      END IF
      read(fp,'(18(1x,A3))') (el(i), i=nclosd+1,nwf)
      read(fp,'(I3,2I8,3X)') nblock,nint,lsdim

!     allocate CSF list variables

      allocate(term_bl(nblock))
      allocate(ncfg_bl(nblock))
      allocate(NZE_bl(nblock))
      allocate(ncn_bl(nblock))
      allocate(nze_max(nblock))

!    read cfg.h 

      ncfg_tot=0
      nze_tot=0
      ncn_tot=0
      do nb = 1,nblock
        READ(fp,'(2X,A3,I6,I10,I15,I10)')  &
          term_bl(nb),ncfg_bl(nb),NZE_bl(nb),ncn_bl(nb),nze_max(nb)

          term_bl(nb) = trim(adjustl(term_bl(nb)))

          ncfg_tot=ncfg_tot+ncfg_bl(nb)
          nze_tot=nze_tot+nze_bl(nb)
          ncn_tot=ncn_tot+ncn_bl(nb)
      end do

      if(safe_close(fp)/=0)call crash

!
!  *****  DETERMINE NCLOSD SHELLS
!

      allocate(ossN(nwf))
      allocate(ossL(nwf))

      do i=1,nclosD
        J = 3
        IF (EL(I)(1:1) .NE. ' ') J = 2
        ossL(I) = LVAL(EL(I)(J:J))
        ossN(I) = ICHAR(EL(I)(J-1:J-1)) - ICHAR('1') + 1
      end do

!
!  *****  DETERMINE THE OTHER ELECTRONS AND ORDER
!
19    FORMAT(/' There are ',I3,' orbitals as follows:'/(1X,18(1X,A3)))
      WRITE(0,19) NWF,(EL(J),J=1,NWF)
!

      do i = nclosd+1,nwf
        J = 2
        IF (EL(I)(1:1) .EQ. ' ') J = 3
        ossL(I) = LVAL(EL(I)(J:J))
        ossN(I) = ICHAR(EL(I)(J-1:J-1)) - ICHAR('1') + 1
      end do

!  .. find maximum l-value
      lmax = 0
      do i=1,nwf
         lmax = max (lmax,ossl(i))
      end do


!
!  .. INITIALIZE radial functions
!
!
! ..  initialize eigenvectors, if available
!
!
! ..  initialize integral data
!
      intptr(1:4) = 0
      call spintgrl
!
      RETURN
      END subroutine DATA
