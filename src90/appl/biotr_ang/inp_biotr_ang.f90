 
!###
! check and open .lsj file
!###
 
      SUBROUTINE OPEN_LSJF(TTFILE, ICI, JFILE, TRFILE, IULSJ, NAME) 
      use inout_C
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:11:08  11/17/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: ICI 
      INTEGER , INTENT(IN) :: IULSJ 
      CHARACTER , INTENT(IN) :: TTFILE*64 
      CHARACTER , INTENT(OUT) :: TRFILE*64 
      CHARACTER , INTENT(IN) :: JFILE(2)*24 
      CHARACTER , INTENT(IN) :: NAME(2)*24 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J 
!-----------------------------------------------
 
    9 FORMAT('  Transition between files: '/,2X,A24,/,2X,A24) 
 
 
      J = INDEX(TTFILE,' ') 
!      if (j .eq. 1) then
!        WRITE(ISCW,*) ' Names may not start with blanks'
!        go to 3
!      end if
!
      ICI = 1 
      OPEN(UNIT=IUJ(1), FILE=JFILE(1), STATUS='OLD', POSITION='asis') 
      OPEN(UNIT=IUJ(2), FILE=JFILE(2), STATUS='OLD', POSITION='asis') 
      TRFILE = TTFILE(1:J-1)//'.lsj' 
      OPEN(UNIT=IULSJ, FILE=TRFILE, STATUS='UNKNOWN', POSITION='asis') 
      WRITE (IULSJ, 9) NAME(1), NAME(2) 
 
      RETURN  
      END SUBROUTINE OPEN_LSJF 


 
 
!###
! prepare for calclation of elements
!###
 
      SUBROUTINE INI_CALC(IEM, LAM, NPAIR, NTERMS, IFL) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:11:08  11/17/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IEM 
      INTEGER , INTENT(IN) :: LAM 
      INTEGER , INTENT(IN) :: NPAIR 
      INTEGER , INTENT(OUT) :: NTERMS 
      INTEGER , INTENT(INOUT) :: IFL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IWRITE 
!-----------------------------------------------
 
   10 FORMAT(' ',41X,'_ '/,' < INITIAL STATE ||',A2,'(',I1,&
      ')|| FINAL STATE > = >  COEFF * WEIGHT(INITIAL,I) * WEIGHT(FINAL,J) * < N&
      &L||',A2,'(',I1,')||N''L''>') 
 
    4 FORMAT('+',41X,'_ '/,41X,'I,J'/,/,/,5X,'COEFF      I    J     < NL||',A2,&
         '(',I1,')||N''L''>'/,/) 
 
      WRITE (IWRITE, 10) IEM, LAM, IEM, LAM 
      WRITE (IWRITE, 4) IEM, LAM 
      WRITE (6, *) ' npair = ', NPAIR 
      NTERMS = 0 
      IF (IFL == 2) IFL = 3 
 
      RETURN  
      END SUBROUTINE INI_CALC 


 
!###
! initialize some more variables
!###
 
      SUBROUTINE INI_VAR(LBUF, MXL, LMAX, NTESTB) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:11:08  11/17/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: LBUF 
      INTEGER , INTENT(OUT) :: MXL 
      INTEGER , INTENT(OUT) :: NTESTB 
      INTEGER , INTENT(IN) :: LMAX(2) 
!-----------------------------------------------
 
!mrg  LBUF = 96*1027
      LBUF = 1 
      MXL = MAX0(LMAX(1),LMAX(2)) - 1 
!. Define TEST level for BIOTRN : 0 => silence, 1=> just a bit
!mrg  NTESTB = 1
      NTESTB = 15 
 
      RETURN  
      END SUBROUTINE INI_VAR 


 
!###
! initialize some other variables
!###
 
      SUBROUTINE INI_TYPE(IWRITE, IM, IFL, LAM, LAM2, VOK) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:11:08  11/17/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IWRITE 
      INTEGER , INTENT(OUT) :: IFL 
      INTEGER , INTENT(IN) :: LAM 
      INTEGER , INTENT(OUT) :: LAM2 
      LOGICAL , INTENT(OUT) :: VOK 
      CHARACTER , INTENT(INOUT) :: IM 
!-----------------------------------------------
!
 
    7 FORMAT('1',/,' ELECTRIC TRANSITION OF ORDER ',I3,/) 
    8 FORMAT('1',/,' MAGNETIC TRANSITION OF ORDER ',I3,/) 
      IF (IM=='*' .OR. IM==' ') STOP ' END OF CASE' 
      IF (IM == 'e') IM = 'E' 
      IF (IM == 'm') IM = 'M' 
      IF (IM == 'E') THEN 
         IFL = 1 
         WRITE (IWRITE, 7) LAM 
      ELSE 
         IFL = 2 
         WRITE (IWRITE, 8) LAM 
      ENDIF 
      LAM2 = LAM + LAM 
      VOK = .FALSE. 
 
!     .. let us compute velocity also for rel = .true.
!     if(.not.rel.and.ifl.eq.1.and.lam.le.2) vok = .true.
      IF (IFL==1 .AND. LAM<=2) VOK = .TRUE. 
!
 
      RETURN  
      END SUBROUTINE INI_TYPE 


 
 
!###
! allocate memory for orthogonality
!###
 
      subroutine mem_orth(norbi, norbf, north, ibugm, qiorth) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:11:08  11/17/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use orth_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: norbi 
      integer , intent(in) :: norbf 
      integer , intent(out) :: north 
      integer , intent(in) :: ibugm 
      INTEGER, POINTER :: QIORTH(:)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
 
      if (norbi/=0 .or. norbf/=0) then 
         write (6, *) ' norbi or norbf is different from 0...' 
         write (6, *) ' something wrong...' 
      endif 
      north = norbi*norbf 
      write (6, *) ' north = ', north 
      if (north > 0) then 
         if (ibugm /= 0) write (6, *) ' qiorth  allocation: north = ', north 
!         call alloc(qiorth,north,4)
         call orth 
      endif 
 
      return  
      end subroutine mem_orth 


 
 
!###
! set ncf
!###
 
      SUBROUTINE SET_NCF(LCFG, NCFG, MCFG, KCFG, NCF) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:11:08  11/17/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: LCFG 
      INTEGER , INTENT(IN) :: NCFG 
      INTEGER , INTENT(IN) :: MCFG 
      INTEGER , INTENT(IN) :: KCFG 
      INTEGER , INTENT(OUT) :: NCF(2) 
!-----------------------------------------------
!
 
      LCFG = NCFG 
      NCF(1) = MCFG 
      NCF(2) = KCFG 
 
      RETURN  
      END SUBROUTINE SET_NCF 


 
!###
! write header
!###
 
      SUBROUTINE WRITE_HDR(IWRITE) 
!
 
    5 FORMAT(/,/,/,20X,'========================'/,20X,'  T R A N S B I O  99 '&
        &/,20X,'========================'/,/) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:11:08  11/17/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IWRITE 
!-----------------------------------------------
 
      WRITE (IWRITE, 5) 
 
      RETURN  
      END SUBROUTINE WRITE_HDR 


 
!###
! set unit numbers
!###
      subroutine set_unit(iuls, iulsj) 
      use inout_C
      use inform_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:11:08  11/17/01  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(out) :: iuls 
      integer , intent(out) :: iulsj 
!-----------------------------------------------
 
!
!     Set machine dependent unit numbers:
!  IREAD  - The unit number of standard input
!  IWRITE - The unit number of standard output
!  ISCW   - The unit number of standard error (screen)
!  iuc(1) - The unit number for the <name>.c file for initial state
!  iuc(2) - The unit number for the <name>.c file for final   state
!  iuw(1) - The unit number for the <name>.w file for initial state
!  iuw(2) - The unit number for the <name>.w file for final   state
!  iul(1) - The unit number for the <name>.l file for initial state
!  iul(2) - The unit number for the <name>.l file for final   state
!  iuj(1) - The unit number for the <name>.j file for initial state
!  iuj(2) - The unit number for the <name>.j file for final   state
!  iuls   - The unit number for   'iul(1).iul(2).ls'    file (results)
!  iulsj  - The unit number for   'iuj(1).iuj(2).lsj'   file (results)
!
      iscw = 0 
      iscw2 = iscw 
      iread = 5 
      iwrite = 6 
      iwrite2 = iwrite 
      iuc(1) = 1 
      iuc(2) = 2 
      iuw(1) = 3 
      iuw(2) = 4 
      iul(1) = 8 
      iul(2) = 9 
      iuj(1) = 10 
      iuj(2) = 11 
      iuls = 13 
      iulsj = 14 
      return  
!     iut(1) = 15
!     iut(2) = 16
 
      end subroutine set_unit 


 
!###
! set debuging parameters
!###
 
      SUBROUTINE SET_DEBUG(IBUG1, IBUG2, IBUG3, IBUGM, NBUG6) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:11:08  11/17/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: IBUG1 
      INTEGER , INTENT(OUT) :: IBUG2 
      INTEGER , INTENT(OUT) :: IBUG3 
      INTEGER , INTENT(OUT) :: IBUGM 
      INTEGER , INTENT(OUT) :: NBUG6 
!-----------------------------------------------
 
      IBUG1 = 0 
      IBUG2 = 0 
      IBUG3 = 0 
      IBUGM = 0 
      NBUG6 = 0 
 
      RETURN  
      END SUBROUTINE SET_DEBUG 


 
!###
! input case
!##
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:11:08  11/17/01  
!...Switches:                     
      SUBROUTINE INP_ATOM(NAME, CFILE, WFILE, LFILE, JFILE, &
          LINTF, TFILE, TTFILE)
      use inout_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:11:08  11/17/01
!...Switches:
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER  :: NAME(2)*24
      CHARACTER , INTENT(OUT) :: CFILE(2)*24
      CHARACTER , INTENT(OUT) :: WFILE(2)*24
      CHARACTER , INTENT(OUT) :: LFILE(2)*24
      CHARACTER , INTENT(OUT) :: JFILE(2)*24
      CHARACTER , INTENT(OUT) :: LINTF(2)*24
      CHARACTER , INTENT(OUT) :: TFILE(2)*24
      CHARACTER , INTENT(OUT) :: TTFILE*64
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(8) :: Q
      INTEGER , DIMENSION(2) :: LSJ
      INTEGER  :: IDX(2)
      INTEGER :: J, I
      CHARACTER :: TRFILE*64
!-----------------------------------------------

      WRITE (ISCW, '(/A,A)') ' Name of Initial State'
      READ (IREAD, '(A)') NAME(1)
      WRITE (ISCW, '(/A,A)') ' Name of Final State'
      READ (IREAD, '(A)') NAME(2)

1     continue
      DO I = 1, 2
         J = INDEX(NAME(I),' ')
         IDX(I) = J
!         IF (J == 1) THEN
!           WRITE (ISCW, *) ' Names may not start with blanks'
!            GO TO 1
!         ENDIF
         CFILE(I) = NAME(I)(1:J-1)//'.c'
         call check_file(CFILE(I),ISCW);
!         WFILE(I) = NAME(I)(1:J-1)//'.w'
!         call check_file(WFILE(I),ISCW);
!         LFILE(I) = NAME(I)(1:J-1)//'.l'
!         JFILE(I) = NAME(I)(1:J-1)//'.j'
!         call check_file(JFILE(I),ISCW);
!         LINTF(I) = NAME(I)(1:J-1)//'.t'
!         TFILE(I)(1:J-1) = NAME(I)(1:J-1)
      END DO

!      TTFILE = TFILE(1)(1:IDX(1)-1)//'.'//TFILE(2)(1:IDX(2)-1)

!    3 j = index(ttfile,' ')
!      if (j .eq. 1) then
!        WRITE(ISCW,*) ' Names may not start with blanks'
!        go to 1
!      end if

      RETURN
      END SUBROUTINE INP_ATOM


 
!##
!    select print output
!##
      SUBROUTINE INP_PRINT(ISCW, IPRINT, TOL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:11:08  11/17/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: ISCW 
      REAL(DOUBLE)  :: TOL 
      LOGICAL  :: IPRINT 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NOD = 220 
      INTEGER, PARAMETER :: NWD = 60 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER :: PP 
!-----------------------------------------------
 
!      write(iscw,*) ' intermediate printing (y or n) ?  '
!      read(iread,'(a1)') pp
!      if ( pp .eq. 'y' .or. pp .eq. 'Y' ) then
!         Iprint = .true.
!         write(iscw,*) '  tolerance for printing ?  '
!         read(iread,*) tol
!      else
!         Iprint = .false.
!      endif
 
      RETURN  
      END SUBROUTINE INP_PRINT 


 
!##
! user input of the transition type
!###
      SUBROUTINE INP_TYPE(ISCW, IREAD, IM, LAM) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:11:08  11/17/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ISCW 
      INTEGER , INTENT(IN) :: IREAD 
      character  :: IM 
      INTEGER  :: LAM 
!-----------------------------------------------
!
 
      WRITE (ISCW, *) ' Type of transition ? (E1, E2, M1, M2, .. or *) ' 
 
      READ (IREAD, '(A1,I1)') IM, LAM 
 
      IF (LAM < 0) THEN 
         WRITE (ISCW, '(A1,I1)') ' Incorrect entry, type=', IM, LAM 
      ELSE IF (LAM == 0) THEN 
         LAM = 1 
      ENDIF 
      RETURN  

      END SUBROUTINE INP_TYPE 

!##
!   abort if memory allocation fails
!##

      subroutine mem_fail(iscw,sz,name,ierr)
      implicit none
!------------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: iscw, sz,ierr
      character(*) , INTENT(IN) :: name
!-----------------------------------------------
!
      write(iscw,'(A,A,A)') &
       & '... Error: allocate(', name , ') failed, exiting!... ';
      write(iscw,'(A,I2,A,I10)') &
         ' :: STAT = ', ierr, ', requested size = ', sz;
      stop

      end subroutine mem_fail


!##
!   check if the file exist
!##
 
      subroutine check_file(file_name,ISCW)
!------------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(*), INTENT(IN) :: file_name
      integer, intent(in) :: ISCW
 
!------------------------------------------------
!    LOCAL  variables
!-----------------------------------------------
      logical Y;
      integer idx,l;

      inquire(FILE=file_name, exist=Y)
      IDX = index(file_name,' ');
      L = len(file_name);

      if (.not.Y) then
            write (iscw,'(A,A)') file_name(1:L-IDX), &
                ' not found: the program is exiting!...'
         call exit(1)
      endif

      end subroutine check_file
 
 
