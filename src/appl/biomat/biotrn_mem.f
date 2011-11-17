*
*     ------------------------------------------------------------------
*	B I O T R N_M E M : compute sizes of arrays
*     ------------------------------------------------------------------
*
      SUBROUTINE BIOTRN_MEM(NLI,NCSFI,NCII,NLF,NCSFF,NCIF,
     : NGRID,MXL,LBUF,NTESTG,KFREE,KLPI,KLPF,KLPIF,KLSTOT,KLSIF,KLSIFI,
     : KLCI,KLCF,KLSCR,KBUF,KIBUF,KLCISC,KLCIF,NLIMX,NLFMX,NLIFMX,
     : NLTI,NLTF)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2),iut(2)
      DIMENSION NLI(MXL),CII(NCSFI,NCII)
      DIMENSION NLF(MXL),CIF(NCSFF,NCIF)
      DIMENSION NINSHL(MXL)
*. Scratch
      !DIMENSION SCR(*)
      !pointer(PSCR,SCR(1))
      double precision, dimension(:), allocatable,target :: scr
*
      NTESTL = 0    
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.1) THEN
        WRITE(6,*)
        WRITE(6,*) '                   *************************'
        WRITE(6,*) '                   *   Entering BIOTRN     *'
        WRITE(6,*) '                   * May it bring you luck *'
        WRITE(6,*) '                   *   and good numbers    *'
        WRITE(6,*) '                   *************************'
        WRITE(6,*)
      END IF
*
*. Scratch should at least be of length ??
      ILI = 1
      ILF = 1
*. Largest number of shells of a given symmetry
      NLIMX = IFNMNX(NLI,MXL+1,1)
      NLFMX = IFNMNX(NLF,MXL+1,1)
      NLIFMX = MAX(NLIMX,NLFMX)
Cmrg  IF(NTEST.GE.10) THEN
C        write(6,*) ' nlimx,nlfmx nlifmx ',
C     &               nlimx,nlfmx,nlifmx
Cmrg  END IF
*. Total numner of shells
      NLTI = IELSUM(NLI, MXL + 1)
      NLTF = IELSUM(NLF, MXL + 1)
Cmrg  IF(NTEST.GE.10)
Cmrg &write(6,*) ' NLTI NLTF', NLTI,NLTF
C      write(6,*) ' NLTI NLTF', NLTI,NLTF
*
* Scratch space for orbital rotations
      KFREE = 1
*
      KLPI = KFREE
      KFREE = KFREE + NLIMX * NGRID
      print*,' In biotrn: KLPI    = ',KLPI
      print*,'            KFREE   = ',KFREE
*
      KLPF = KFREE
      KFREE = KFREE + NLFMX * NGRID
      print*,' In biotrn: KLPF    = ',KLPF
      print*,'            KFREE   = ',KFREE
*
      KLPIF = KFREE
      KFREE = KFREE + NLIFMX * NGRID
      print*,' In biotrn: KLPIF   = ',KLPIF
      print*,'            KFREE   = ',KFREE
*. Total overlap matrix
      KLSTOT = KFREE
      KFREE = KFREE + NLTI*NLTF
      print*,' In biotrn: KLSTOT  = ',KLSTOT
      print*,'            KFREE   = ',KFREE
*
      KLSIF = KFREE
      KFREE = KFREE + NLIFMX ** 2
      print*,' In biotrn: KLSIF   = ',KLSIF
      print*,'            KFREE   = ',KFREE
*
      KLSIFI = KFREE
      KFREE = KFREE + NLIFMX ** 2
      print*,' In biotrn: KLSIFI  = ',KLSIFI
      print*,'            KFREE   = ',KFREE
*
      KLCI = KFREE
      KFREE = KFREE + NLIFMX ** 2
      print*,' In biotrn: KLCI    = ',KLCI
      print*,'            KFREE   = ',KFREE
*
      KLCF = KFREE
      KFREE = KFREE + NLIFMX ** 2
      print*,' In biotrn: KLCF    = ',KLCF
      print*,'            KFREE   = ',KFREE
*
      KLSCR = KFREE
      KFREE = KFREE + NLIFMX ** 2 + NLIFMX*(NLIFMX+1)
      print*,' In biotrn: KLSCR   = ',KLSCR
      print*,'            KFREE   = ',KFREE
*. For CI transformation
      KBUF = KFREE
      KFREE = KFREE + LBUF
      print*,' In biotrn: KBUF    = ',KBUF
      print*,'            KFREE   = ',KFREE
*
      KIBUF = KFREE
      KFREE = KFREE + 4 * LBUF
      print*,' In biotrn: KIBUF   = ',KIBUF
      print*,'            KFREE   = ',KFREE
*
      KLCISC = KFREE
      KFREE = KFREE + MAX(NCSFI*NCII,NCSFF*NCIF)
      print*,' In biotrn: KLCISC  = ',KLCISC
      print*,'            KFREE   = ',KFREE
*
      KLCIF  = KFREE
      KFREE = KFREE + MAX(NCSFI*NCII,NCSFF*NCIF)
      print*,' In biotrn: KLCIF   = ',KLCIF
      print*,'            KFREE   = ',KFREE
      print*,'             NCSFI  =   ',ncsfi
      print*,'             NCII   =   ',ncii
      print*,'             NCSFF  =   ',ncsff
      print*,'             NCIF   =   ',ncif
      print*,'         =>  FREE =     ',kfree
*. Here a scratch space of length KFREE-1 could be allocated
*. Check length of scratch
*
       RETURN
       END subroutine biotrn_mem
