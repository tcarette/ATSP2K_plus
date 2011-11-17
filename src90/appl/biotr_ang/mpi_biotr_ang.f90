! Last update: September, 1997 by G. Gaigalas.
!
!  Things to do: - clean the old nonorthogonalities.
!                - check a Breit-Pauli case (I never did!).
!                - test inactive shells.
!                - add a threshold for contributions to the line strength,
!                  after counter-transformation.
!  BE AWARE: you should use fractional parentage modules to get consistency
!            with the Gedeminas-adapted angular codes !!!
!
!*************************************************************************
!
!  This program evaluates the gf values and transition probabilities
!  of electric and magnetic transitions in the LS coupling scheme
!  (length and velocity forms for E1/E2) or LSJ Breit-Pauli intermediate
!  coupling scheme (length form only).
!
!  M. R. Godefroid    Laboratoire de Chimie Physique Moleculaire
!                     Universite Libre de Bruxelles, Belgium
!  A. Hibbert         Department of Applied Mathematics
!                     Queen's University, Belfast, Northern Ireland
!  C. Froese Fischer  Department of Computer Science
!                     Vanderbilt University, Nashville, U.S.A.
!  P. Jonsson         Department of Physics,
!                     Lund Institute of Technology, Lund, Sweden
!  J. Olsen           Department of Theoretical Chemistry,
!                     Chemical Center, Lund, Sweden
!
!
!  References:
!  ----------
!  1. A. Hibbert et al, Comput. Phys. Commun. 51(1988)285
!  2. C. Froese Fischer et al, Comput. Phys. Commun. 64(1991)486-500
!  3. C. Froese Fischer and M. Godefroid, Comput. Phys. Commun. 64(1991)
!     501-519
!  4. P.A. Malmqvist, Int.J. of Quantum Chemistry, XXX, 479-94 (1986)
!  5. J. Olsen, M.R. Godefroid, P. Jonsson, P.A. Malmqvist and
!     C. Froese Fischer, Phys. Rev. A, submitted. (#AY5043)
!
!*************************************************************************
!
      PROGRAM BIOTR_ANG 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE DEBUG_C 
      USE DBG_C 
      USE INOUT_C 
      USE EMS_C 
      !USE NTRM_C 
      USE NOR_C 
      !USE RED_C 
      USE PARAM_C, ONLY: TOL 
      !USE FO_C 
      !USE RAS_C 
      !use inform_C
      use ndims_C, ONLY: NCFG
      !use nel_C
      use state_C
      use ovrlap_C
      !use medefn_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:25:13  11/17/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !USE initt_I 
      USE set_unit_I 
      USE initr_I 
      USE factrl_I 
      USE write_hdr_I 
      USE set_debug_I 
      USE inp_atom_I 
      USE inp_print_I 
      USE inp_type_I 
      USE nonh1_I 
      USE cfgin2_I 
      USE rasin_I 
      USE set_ncf_I 
      USE mem_orth_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NOD = 220 
!      INTEGER, PARAMETER :: LWORK = 6000000 
      INTEGER, PARAMETER :: MAXPROC = 9999 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(8) :: Q 
      INTEGER , DIMENSION(2) :: IDX, LSJ 
      INTEGER :: MYID, IULS, IULSJ, MCFG, KCFG, lcfg, U 
!      REAL(DOUBLE), DIMENSION(LWORK) :: WORK 
      LOGICAL :: IPRINT 
      CHARACTER :: IDSTRING*4, STARTDIR*128, PERMDIR*128, TMPDIR*128, FILE_CINT&
         *128, FILE_IH*128, FILE_ICO*128, FILE_CLST*128, FILE_ATOMC*128, IM, &
         YES, PP, M_OR_C 
      CHARACTER , DIMENSION(2) :: CFILE*24, WFILE*24, JFILE*24, NAME*24, & 
             LFILE*24, LINTF*24, TFILE*24 
      CHARACTER :: TTFILE*64, TRFILE*64, LIGNE*80, LINE*70, CONFIGI*64, &
         CONFIGF*64 
      CHARACTER , DIMENSION(8) :: ELC*3 
      CHARACTER, DIMENSION(15) :: COUPLE*3 
!-----------------------------------------------
!      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
! >>>>>>     <<<<<<<

    6 FORMAT(A1,I1) 
      MYID = 0 
 
! ... set unit numbers
      CALL SET_UNIT (IULS, IULSJ) 
 
      CALL INITR 
 
      CALL FACTRL (32) 
 
      IF (MYID == 0) CALL WRITE_HDR (IWRITE) 
 
      REL = .FALSE. 
 
      IF (MYID == 0) CALL SET_DEBUG (IBUG1, IBUG2, IBUG3, IBUGM, NBUG6) 
 
! ... get input files
      IF (MYID == 0)  &
        CALL INP_ATOM (NAME,CFILE,WFILE,LFILE,JFILE,LINTF,TFILE,TTFILE) 
 
! bcast file names !
 
!dbg  OPEN(UNIT=iut(1),STATUS='scratch',FORM='FORMATTED')
!dbg  OPEN(UNIT=iut(2),STATUS='scratch',FORM='FORMATTED')
!
! --- calculate the one-electron part for initial and final superpositions
 
      IF (MYID == 0) CALL INP_PRINT (ISCW, IPRINT, TOL) 
! bcast iprint, tol
      REL = .TRUE.                               ! compute only .lsj 
!bcast rel
 
!      IF (MYID == 0) CALL INP_TYPE (ISCW, IREAD, IM, LAM) 
!bcast IM,LAM
 
! open cfg lists
      OPEN(UNIT=IUC(1), FILE=CFILE(1), STATUS='OLD', POSITION='asis') 
      OPEN(UNIT=IUC(2), FILE=CFILE(2), STATUS='OLD', POSITION='asis') 
 
! angular data initial state
      WRITE (ISCW, *) &
         ' One-electron coupling coefficients for initial state...' 
      CALL NONH1 (1) 
      REWIND (IUC(1)) 
 
! angular data final state
      WRITE (ISCW, *) &
         ' One-electron coupling coefficients for final state...' 
      CALL NONH1 (2) 
      REWIND (IUC(2)) 
!
!
! check  cfg lists
      CALL CFGIN2 (MCFG, KCFG) 
 
!* --- specify the type of RAS calculation
      CALL RASIN 
!
! ... set some variables 
      CALL SET_NCF (LCFG, NCFG, MCFG, KCFG, NCF) 
!
!     allocate memory for orthogonality and call orth
 
      CALL MEM_ORTH (NORBI, NORBF, NORTH, IBUGM, QIORTH) 

! save program data for the next step biotr in a file _case.lst 
      U = 50;      
      open(unit=U,file='_case.lst',form='unformatted');
      write(ISCW,'(/A)') ' ...Saving program data in _case.lst'
      call save_data(U);
      close(U);

! save the integral and coefficient data in _intc.lst 
      open(unit=U,file='_intc.lst',form='unformatted');
      write(ISCW,'(A/)') ' ...Saving angular data in _intc.lst'
      call csf_W(U);
      close (U);

      STOP  
      END PROGRAM BIOTR_ANG 
