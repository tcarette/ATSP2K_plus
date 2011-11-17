************************************************************************
*                                                                      *
      SUBROUTINE DALLOC (PTR, LENGTH)
*                                                                      *
*   This  subprogram deallocates the  memory locations that are asso-  *
*   ciated  with the pointer PTR.  LENGTH retained from earlier        *
*   requirements
*                                                                      *
*   Original code by Charlotte F. Fischer.                             *
*                                                                      *
*   This revision by Farid A. Parpia.     Last revision: 17 Sep 1992   *
*   Updated by Anders Ynnerman            Last revision: 26 Jan 1994   *
*   Updated for IBM and DEC by Bieron     Last revision: 05 Apr 1995   *
*                                                                      *
************************************************************************
*
c
cbieron DEC pointer
c
      POINTER (PTR,ptrdummy)
      INTEGER LENGTH

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      CALL FREE (PTR)                    ! SUN & DEC
!      CALL FREE (%val(PTR))              ! IBM
!      CALL DISCLAIM (%val(PTR))          ! IBM
!      IERR = MZDALLOC(PTR,NSIZE)         ! TSS, Needs nsize 
!      IABORT = 0
!       CALL HPDEALLC(PTR,IERR,IABORT)     ! CRAY 
!      ptr = 0;
!      IF (IERR.NE.0) THEN
!         PRINT *, 'DEALLOCATION ERROR ',IERR
!      END IF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Choose from above

        CALL FREE (PTR)                    ! SUN & DEC & SGI

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
