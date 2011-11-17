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
      PROGRAM BIOTR 
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE DEBUG_C
      USE DBG_C
      USE INOUT_C
      USE EMS_C
      USE NTRM_C
      USE NOR_C
      USE RED_C
      USE PARAM_C, ONLY: TOL
      USE FO_C
      USE RAS_C
      use inform_C
      use ndims_C, ONLY: NCFG
      use nel_C
      use state_C
      use ovrlap_C
      use medefn_C
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
      IMPLICIT NONE
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:47:57  11/18/01  
!...Switches:                     
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NOD = 220 
      INTEGER, PARAMETER :: LWORK = 1 
      INTEGER, PARAMETER :: MAXPROC = 9999 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LSJ(2) 
      INTEGER :: U
      INTEGER :: MYID, IULS, IULSJ, MCFG, KCFG, LCFG, ICI, LAM2,& 
                 NPAIR, LBUF, MXL, NTESTB,IDUMMY,NTEST 
      DOUBLE PRECISION :: WORK_(LWORK) 
      LOGICAL :: IPRINT 
      CHARACTER*1 :: IDSTRING*4, STARTDIR*128, PERMDIR*128, &
                     TMPDIR*128, FILE_CINT*128, FILE_IH*128, &
                     FILE_ICO*128, FILE_CLST*128, FILE_ATOMC*128, &
                     IM, YES, PP, M_OR_C 
      CHARACTER(len=24),dimension(2) :: CFILE,WFILE,JFILE,NAME,&
                     LFILE, LINTF, TFILE
      CHARACTER*1 :: TTFILE*64, TRFILE*64, LIGNE*80, LINE*70, CONFIGI*64, &
                     CONFIGF*64 
      !CHARACTER*1 :: ELC*3(8) 
      !CHARACTER*1 :: COUPLE*3(15) 
!-----------------------------------------------
!      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
! >>>>>>     <<<<<<<
!

      MYID = 0 

! ... set unit numbers

      CALL SET_UNIT (IULS, IULSJ)
      call check_file('_case.lst',ISCW);
      call check_file('_intc.lst',ISCW);

! ... write program header
      IF (MYID .EQ. 0) CALL WRITE_HDR (IWRITE)
      
! ... get input files <name>.j and <name>.w
      IF (MYID .EQ. 0) CALL INP_ATOM (NAME, CFILE, WFILE, LFILE, &
                        JFILE, LINTF, TFILE,TTFILE)

! compute only relativistic transitions
      REL = .TRUE.
!
! --- open .j and .lsj files according rel value
!
      CALL OPEN_LSJF (TTFILE, ICI, JFILE, TRFILE, IULSJ, NAME)

! ... get transition E1,..., M1,.., initilaize im, lam2, vok, ifl
      IF (MYID .EQ. 0) CALL INP_TYPE (ISCW,IWRITE,IREAD,IFL,LAM,LAM2,VOK,IM)
!bcast 

! change of the settings for debugging output in runtime (turned off)
!      IF (MYID .EQ. 0) CALL SET_DEBUG (IBUG1, IBUG2, IBUG3, IBUGM, NBUG6)
! print level (turned off)
!      IF (MYID == 0) CALL INP_PRINT (ISCW, IPRINT, TOL)

! bcast iprint, tol

! read program data from file _case.lst = U
      U = 50;
      open(unit=U,file='_case.lst',form='unformatted');
      write(ISCW,'(/A)') ' ...Reading program data saved in _case.lst'
      call read_data(U);
      close(U)

! read integral and coefficient data from _intc.lst = U
      open(unit=U,file='_intc.lst',form='unformatted');
      write(ISCW,'(A/)') ' ...Reading angular data saved in _intc.lst'
      call csf_R(U);
      close (U)

! bcast file names !
 
!dbg  OPEN(UNIT=iut(1),STATUS='scratch',FORM='FORMATTED')
!dbg  OPEN(UNIT=iut(2),STATUS='scratch',FORM='FORMATTED')
!

      call set_ncf(LCFG, NCFG, MCFG, KCFG, NCF)
 
! open for reading the .wfn files
      OPEN(UNIT=IUW(1), FILE=WFILE(1), STATUS='OLD', FORM='UNFORMATTED', &
         POSITION='asis') 
      OPEN(UNIT=IUW(2), FILE=WFILE(2), STATUS='OLD', FORM='UNFORMATTED', &
         POSITION='asis') 

      CALL READW2 

      CLOSE(IUW(1)) 
      CLOSE(IUW(2)) 
 
!
! --- calculate the overlap matrix between initial and final
!     orbitals
     CALL BRKT 
 
      CALL EIGVEC (ICI, CONFIGI, CONFIGF, LSJ) 
      IF (MYID .EQ. 0) WRITE (6, *) 'jv1,jv2', JV1(1), JV2(1) 
 
!
! --- allocate the /MULT/
      CALL ALMULT (NPAIR) 
!
! --- Here,we go : call routine to perform tranformations
!     of shells and CI coefficients so orbital bases become
!     biorthogonal
 
! set MXL, LMAX, NTESTB
      CALL INI_VAR (LBUF, MXL, LMAX, NTESTB) 
 
      IF (MYID .EQ. 0) WRITE (ISCW, *) ' Call BIOTRN...' 
      IF (IBUGM.NE.0 .AND. MYID.EQ.0) WRITE (6, *) ' wt1(1) = ', WT1(1), &
         ' wt2(1) = ', WT2(1) 
!

!. Here we are...
      CALL BIOTRN (P(1,1), NL(1,1), WT1(1), NCF(1), NVC(1), IDUMMY, &
         P(1,IWF(1)+1), NL(1,2), WT2(1), NCF(2), NVC(2), IDUMMY, NOD, &
         MXL, NINAC(1,1), LBUF, 15, 16, NTESTB) 
!
!. Here we were...
      IF (IBUGM.NE.0 .AND. MYID.EQ.0) WRITE (6, *) ' wt1(1) = ', WT1(1), &
         ' wt2(1) = ', WT2(1) 
!
! --- calculate the rotated slopes at the origin
      CALL SLOPE 

      IF (MYID .EQ. 0) WRITE (6, *) ' Test call of BRKT after BIOTRN' 
      CALL BRKT 
!
! --- calculate the overlap and radial one-electron transition integrals
!     velocity form only for E1/E2 using the non-relativistic option.
      WRITE (ISCW, *) ' Calculation of radial integrals...' 
      CALL RADINT (IM) 

!
! --- calculate the elements of vshell
 
      CALL INI_CALC (IEM(IFL), LAM, NPAIR, NTERMS, IFL) 
 
      CALL CALCUL (NPAIR); 

      CALL PROBAB (ICI, IULS, IULSJ, NPAIR, CONFIGI, CONFIGF, IM, IPRINT) 
 
!      call MPI_FINALIZE(ierr)
 
      STOP  
      END PROGRAM BIOTR 
 
