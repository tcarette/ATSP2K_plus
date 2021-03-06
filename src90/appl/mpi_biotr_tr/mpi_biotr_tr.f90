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
!-----------------------------------------------
!   MPI data
!-----------------------------------------------
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
      character*4 idstring
      character*128 startdir, permdir, tmpdir
      integer lclst,lentmp,lenstart,lenperm
! >>>>>>     <<<<<<<

!...Translated by Pacific-Sierra Research 77to90  4.3E  16:47:57  11/18/01  
!...Switches:                     
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NOD = 220 
      INTEGER, PARAMETER :: LWORK = 1 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LSJ(2) 
      INTEGER :: U,ii
      INTEGER :: IULS, IULSJ, MCFG, KCFG, LCFG, ICI, LAM2,& 
                 NPAIR, LBUF, MXL, NTESTB,IDUMMY,NTEST 
      DOUBLE PRECISION :: WORK_(LWORK) 
      LOGICAL :: IPRINT 
      CHARACTER*1 :: FILE_CINT*128, FILE_CLST*128, FILE_ATOMC*128, &
                     IM, YES, PP, M_OR_C 
      CHARACTER(len=128), dimension(2) :: FILE_W, FILE_J, FILE_C
      CHARACTER(len=24),dimension(2) :: CFILE,WFILE,JFILE,NAME,&
                     LFILE, LINTF, TFILE
      CHARACTER*1 :: TTFILE*64, TRFILE*64, LIGNE*80, LINE*70, CONFIGI*64, &
                     CONFIGF*64 
      integer :: w1, w2, f1, f2
!-----------------------------------------------
!      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
! >>>>>>     <<<<<<<
!

      MYID = 0 

! ... set unit numbers
!      CALL SET_UNIT (IULS, IULSJ)
!      call check_file('_case.lst',0);

      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      if (myid == 0) call check_file('_case.lst',0);
      write(idstring,'(I4.4)') myid

! ... write program header
      IF (MYID .EQ. 0) CALL WRITE_HDR (0)

      CALL SET_UNIT (IULS, IULSJ)
      
! ... get input files <name>.j and <name>.w
      IF (MYID .EQ. 0) CALL INP_ATOM (NAME, CFILE, WFILE, LFILE, &
                        JFILE, LINTF, TFILE,TTFILE)
      do ii = 1,2
         call MPI_BCAST(NAME(ii),24,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(CFILE(ii),24,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(WFILE(ii),24,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(LFILE(ii),24,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(JFILE(ii),24,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(LINTF(ii),24,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(TFILE(ii),24,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);
      end do
      call MPI_BCAST(TTFILE,64,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);

! compute only relativistic transitions
      REL = .TRUE.
!
! --- open .j and .lsj files according rel value
!
      if (myid == 0) &
            CALL OPEN_LSJF (TTFILE, ICI, JFILE, TRFILE, IULSJ, NAME)

! ... get transition E1,..., M1,.., initilaize im, lam2, vok, ifl
      IF (MYID .EQ. 0) &
          CALL INP_TYPE (ISCW,IWRITE,IREAD,IFL,LAM,LAM2,VOK,IM)
      call MPI_BCAST(LAM,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(IM,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(LAM2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(VOK,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(IFL,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);   
!bcast 

! change of the settings for debugging output in runtime (turned off)
!      IF (MYID .EQ. 0) CALL SET_DEBUG (IBUG1, IBUG2, IBUG3, IBUGM, NBUG6)
      call MPI_BCAST(IBUG1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(IBUG2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(IBUG3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(IBUGM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(NBUG6,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);

! print level (turned off)
!      IF (MYID == 0) CALL INP_PRINT (ISCW, IPRINT, TOL)
!      call MPI_BCAST(IPRINT,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr);
!      call MPI_BCAST(TOL,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);

! bcast iprint, tol

! read program data from file _case.lst = U
       U = 50;
       if (myid == 0) then
          open(unit=U,file='_case.lst',form='unformatted');
          write(ISCW,'(/A)') ' ...Reading program data from _case.lst'
       end if
         call read_data(U);
       if (myid == 0)  close(U)

! get directories
      startdir = ' '; permdir = ' '; tmpdir = ' ';
      call mpi_work_dir(startdir,permdir,tmpdir,lenstart,lenperm,lentmp)
!      lenstart = LEN_TRIM (startdir) !- 1
!      lenperm = LEN_TRIM (permdir) !- 1
!      lentmp = LEN_TRIM (tmpdir) !- 1

      file_clst  = tmpdir(1:lentmp)//'/_intc.lst.'//idstring
      lclst  = len_trim(file_clst);

      call check_file(file_clst(1:lclst),ISCW);

! read integral and coefficient data from _intc.lst = U
      open(unit=U,file=file_clst(1:lclst),form='unformatted');
      if (myid == 0) &
         write(ISCW,'(A/)') ' ...Reading angular data from _intc.lst'
      call csf_R(U);
      close (U)

! bcast file names !
 
      if (myid == 0 ) call set_ncf(LCFG, NCFG, MCFG, KCFG, NCF)
      call MPI_BCAST(LCFG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(MCFG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(KCFG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
 
! open for reading the .wfn files
      w1 = Len_trim(wfile(1))
      w2 = Len_trim(wfile(2));
      FILE_W(1) = startdir(1:lenstart)//'/'//wfile(1)(1:w1);
      FILE_W(2) = startdir(1:lenstart)//'/'//wfile(2)(1:w2);
      f1 = len_trim(FILE_W(1));
      f2 = len_trim(FILE_W(2));
      call check_file(FILE_W(1)(1:f1),ISCW);
      call check_file(FILE_W(2)(1:f2),ISCW);
      
      OPEN(UNIT=IUW(1), FILE=FILE_W(1)(1:f1), STATUS='OLD', &
           FORM='UNFORMATTED', POSITION='asis') 
      OPEN(UNIT=IUW(2), FILE=FILE_W(2)(1:f2), STATUS='OLD', &
           FORM='UNFORMATTED', POSITION='asis') 

      CALL READW2 

      CLOSE(IUW(1)) 
      CLOSE(IUW(2)) 
 
!
! --- calculate the overlap matrix between initial and final
!     orbitals
      CALL BRKT 

! open .c files 
      w1 = Len_trim(cfile(1))
      w2 = Len_trim(cfile(2));
      FILE_C(1) = startdir(1:lenstart)//'/'//cfile(1)(1:w1);
      FILE_C(2) = startdir(1:lenstart)//'/'//cfile(2)(1:w2);
      f1 = len_trim(FILE_C(1));
      f2 = len_trim(FILE_C(2));
      call check_file(FILE_C(1)(1:f1),ISCW);
      call check_file(FILE_C(2)(1:f2),ISCW);
      OPEN(UNIT=IUC(1), FILE=FILE_C(1)(1:f1), STATUS='OLD', &
           POSITION='asis')
      OPEN(UNIT=IUC(2), FILE=FILE_C(2)(1:f2), STATUS='OLD', &
           POSITION='asis')

! open .j files
      w1 = Len_trim(jfile(1))
      w2 = Len_trim(jfile(2));
      FILE_J(1) = startdir(1:lenstart)//'/'//jfile(1)(1:w1);
      FILE_J(2) = startdir(1:lenstart)//'/'//jfile(2)(1:w2);
      f1 = len_trim(FILE_J(1));
      f2 = len_trim(FILE_J(2));
      call check_file(FILE_J(1)(1:f1),ISCW);
      call check_file(FILE_J(2)(1:f2),ISCW);
      OPEN(UNIT=IUJ(1), FILE=FILE_J(1)(1:f1), STATUS='OLD', &
           POSITION='asis')
      OPEN(UNIT=IUJ(2), FILE=FILE_J(2)(1:f2), STATUS='OLD', &
           POSITION='asis')
      ici = 1
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

      if (myid == 0 ) &
         CALL PROBAB (ICI, IULS, IULSJ, NPAIR, CONFIGI, CONFIGF, IM, IPRINT) 
 
      call MPI_FINALIZE(ierr)
 
      STOP  
      END PROGRAM BIOTR 
 
