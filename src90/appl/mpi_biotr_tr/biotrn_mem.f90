!
!     ------------------------------------------------------------------
!     B I O T R N_M E M : compute sizes of arrays
!     ------------------------------------------------------------------
!
      SUBROUTINE BIOTRN_MEM(NLI, NCSFI, NCII, NLF, NCSFF, NCIF, NGRID, MXL, &
         LBUF, NTESTG, KFREE, KLPI, KLPF, KLPIF, KLSTOT, KLSIF, KLSIFI, KLCI, &
         KLCF, KLSCR, KBUF, KIBUF, KLCISC, KLCIF, NLIMX, NLFMX, NLIFMX, NLTI, &
         NLTF) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE INOUT_C 
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:11:50  11/20/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ifnmnx_I 
      USE ielsum_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NCSFI 
      INTEGER , INTENT(IN) :: NCII 
      INTEGER , INTENT(IN) :: NCSFF 
      INTEGER , INTENT(IN) :: NCIF 
      INTEGER , INTENT(IN) :: NGRID 
      INTEGER , INTENT(IN) :: MXL 
      INTEGER , INTENT(IN) :: LBUF 
      INTEGER , INTENT(IN) :: NTESTG 
      INTEGER , INTENT(OUT) :: KFREE 
      INTEGER , INTENT(OUT) :: KLPI 
      INTEGER , INTENT(OUT) :: KLPF 
      INTEGER , INTENT(OUT) :: KLPIF 
      INTEGER , INTENT(OUT) :: KLSTOT 
      INTEGER , INTENT(OUT) :: KLSIF 
      INTEGER , INTENT(OUT) :: KLSIFI 
      INTEGER , INTENT(OUT) :: KLCI 
      INTEGER , INTENT(OUT) :: KLCF 
      INTEGER , INTENT(OUT) :: KLSCR 
      INTEGER , INTENT(OUT) :: KBUF 
      INTEGER , INTENT(OUT) :: KIBUF 
      INTEGER , INTENT(OUT) :: KLCISC 
      INTEGER , INTENT(OUT) :: KLCIF 
      INTEGER , INTENT(OUT) :: NLIMX 
      INTEGER , INTENT(OUT) :: NLFMX 
      INTEGER , INTENT(OUT) :: NLIFMX 
      INTEGER , INTENT(OUT) :: NLTI 
      INTEGER , INTENT(OUT) :: NLTF 
      INTEGER  :: NLI(MXL) 
      INTEGER  :: NLF(MXL) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(MXL) :: NINSHL 
      INTEGER :: NTESTL, NTEST, ILI, ILF 
      REAL(DOUBLE), DIMENSION(NCSFI,NCII) :: CII 
      REAL(DOUBLE), DIMENSION(NCSFF,NCIF) :: CIF 
      REAL(DOUBLE), ALLOCATABLE, TARGET, DIMENSION(:) :: SCR 
!-----------------------------------------------
!
!. Scratch
!     !DIMENSION SCR(*)
!     !pointer(PSCR,SCR(1))
!
      NTESTL = 0 
      NTEST = MAX(NTESTL,NTESTG) 
      IF (NTEST >= 1) THEN 
         WRITE (6, *) 
         WRITE (6, *) '                   *************************' 
         WRITE (6, *) '                   *   Entering BIOTRN     *' 
         WRITE (6, *) '                   * May it bring you luck *' 
         WRITE (6, *) '                   *   and good numbers    *' 
         WRITE (6, *) '                   *************************' 
         WRITE (6, *) 
      ENDIF 
!
!. Scratch should at least be of length ??
      ILI = 1 
      ILF = 1 
!. Largest number of shells of a given symmetry
      NLIMX = IFNMNX(NLI,MXL + 1,1) 
      NLFMX = IFNMNX(NLF,MXL + 1,1) 
      NLIFMX = MAX(NLIMX,NLFMX) 
!mrg  IF(NTEST.GE.10) THEN
!        write(6,*) ' nlimx,nlfmx nlifmx ',
!     &               nlimx,nlfmx,nlifmx
!mrg  END IF
!. Total numner of shells
      NLTI = IELSUM(NLI,MXL + 1) 
      NLTF = IELSUM(NLF,MXL + 1) 
!mrg  IF(NTEST.GE.10)
!mrg &write(6,*) ' NLTI NLTF', NLTI,NLTF
!      write(6,*) ' NLTI NLTF', NLTI,NLTF
!
! Scratch space for orbital rotations
      KFREE = 1 
!
      KLPI = KFREE 
      KFREE = KFREE + NLIMX*NGRID 
      WRITE (6, *) ' In biotrn: KLPI    = ', KLPI 
      WRITE (6, *) '            KFREE   = ', KFREE 
!
      KLPF = KFREE 
      KFREE = KFREE + NLFMX*NGRID 
      WRITE (6, *) ' In biotrn: KLPF    = ', KLPF 
      WRITE (6, *) '            KFREE   = ', KFREE 
!
      KLPIF = KFREE 
      KFREE = KFREE + NLIFMX*NGRID 
      WRITE (6, *) ' In biotrn: KLPIF   = ', KLPIF 
      WRITE (6, *) '            KFREE   = ', KFREE 
!. Total overlap matrix
      KLSTOT = KFREE 
      KFREE = KFREE + NLTI*NLTF 
      WRITE (6, *) ' In biotrn: KLSTOT  = ', KLSTOT 
      WRITE (6, *) '            KFREE   = ', KFREE 
!
      KLSIF = KFREE 
      KFREE = KFREE + NLIFMX**2 
      WRITE (6, *) ' In biotrn: KLSIF   = ', KLSIF 
      WRITE (6, *) '            KFREE   = ', KFREE 
!
      KLSIFI = KFREE 
      KFREE = KFREE + NLIFMX**2 
      WRITE (6, *) ' In biotrn: KLSIFI  = ', KLSIFI 
      WRITE (6, *) '            KFREE   = ', KFREE 
!
      KLCI = KFREE 
      KFREE = KFREE + NLIFMX**2 
      WRITE (6, *) ' In biotrn: KLCI    = ', KLCI 
      WRITE (6, *) '            KFREE   = ', KFREE 
!
      KLCF = KFREE 
      KFREE = KFREE + NLIFMX**2 
      WRITE (6, *) ' In biotrn: KLCF    = ', KLCF 
      WRITE (6, *) '            KFREE   = ', KFREE 
!
      KLSCR = KFREE 
      KFREE = KFREE + NLIFMX**2 + NLIFMX*(NLIFMX + 1) 
      WRITE (6, *) ' In biotrn: KLSCR   = ', KLSCR 
      WRITE (6, *) '            KFREE   = ', KFREE 
!. For CI transformation
      KBUF = KFREE 
      KFREE = KFREE + LBUF 
      WRITE (6, *) ' In biotrn: KBUF    = ', KBUF 
      WRITE (6, *) '            KFREE   = ', KFREE 
!
      KIBUF = KFREE 
      KFREE = KFREE + 4*LBUF 
      WRITE (6, *) ' In biotrn: KIBUF   = ', KIBUF 
      WRITE (6, *) '            KFREE   = ', KFREE 
!
      KLCISC = KFREE 
      KFREE = KFREE + MAX(NCSFI*NCII,NCSFF*NCIF) 
      WRITE (6, *) ' In biotrn: KLCISC  = ', KLCISC 
      WRITE (6, *) '            KFREE   = ', KFREE 
!
      KLCIF = KFREE 
      KFREE = KFREE + MAX(NCSFI*NCII,NCSFF*NCIF) 
      WRITE (6, *) ' In biotrn: KLCIF   = ', KLCIF 
      WRITE (6, *) '            KFREE   = ', KFREE 
      WRITE (6, *) '             NCSFI  =   ', NCSFI 
      WRITE (6, *) '             NCII   =   ', NCII 
      WRITE (6, *) '             NCSFF  =   ', NCSFF 
      WRITE (6, *) '             NCIF   =   ', NCIF 
      WRITE (6, *) '         =>  FREE =     ', KFREE 
!. Here a scratch space of length KFREE-1 could be allocated
!. Check length of scratch
!
      RETURN  
      END SUBROUTINE BIOTRN_MEM 
