*
*     ------------------------------------------------------------------
*     3          P R O G R A M   L I S T I N G
*     ------------------------------------------------------------------
*
*
*     All comments in the program listing assume the radial function P
*     is the solution of an equation of the form
*
*      P" + ( 2Z/R - Y - L(L+1)/R**2 - E)P = X + T
*
*     where Y is called a potential function
*           X is called an exchange function, and
*           T includes contributions from off-diagonal energy parameter,
*             interactions between configurations, etc.
*
*     The program uses LOG(Z*R) as independent variable and
*                      P/DSQRT(R) as dependent variable.
*
*
*     ------------------------------------------------------------------
*                   M A I N    P R O G R AM
*     ------------------------------------------------------------------
*
*       The MAIN program controls the overall calculation and  allows
*   a series of cases to be processed as one run.  Each  case  itself
*   may  consist  of  a  series of atoms or ions in an iso-electronic
*   sequence.  In each case, all but the initial  estimates  for  the
*   first  are  obtained  by  scaling  the previous results using the
*   scaling of Sec.  (7-2).   Mixing coefficients are left unchanged.
*
*
*
        PROGRAM MCHF
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      include 'mpif.h'
      parameter (MAXPROC=100)
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)

 
        PARAMETER (NWD=94,NOD=220,NOFFD=800,MTERM=20,MEIG=20)
        DIMENSION IVAR(NWD)
*
      CHARACTER EL*3,ATOM*6,TERMs(20)*6, string*72
      LOGICAL PRINT,LD,STRONG, quest
      CHARACTER*8 ANS*1
      character whp*5, buffer*8

      integer	    iblock, nblock
      INTEGER       IN,OUT,PRI,OUC,OUF,OUH,OUD,IOU(7)
      COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,IUD,OUC,OUF,OUD,OUH,ISCW
      COMMON/LABEL/ EL(NWD),ATOM,TERMs
      LOGICAL       FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIEDKE
      COMMON/TEST/  FAIL,OMIT,EZERO,REL,ALL,TRACE
      COMMON/WAVE/  EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
      COMMON/PARAM/ H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      POINTER       (IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :     (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :     (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
      COMMON/NEL/   IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :      QMETH,QIEPTR

      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PVALUE,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
*
      POINTER (QWT,WT(1)),(QWP,WP(1)),  (IQWPTR,W(1))
      COMMON /CFGS/ETOTAL,QWT,QWP,IQWPTR

      pointer (ptmp,tmp(1))

      logical 			:: leigen(meig,mterm), lguess
      integer			:: unit_n
      character (len=128)	:: name_dotL
      double precision		:: eigst_weight(meig,mterm)   
      integer 			:: cf_tot(mterm) 
      !integer cf_tot
      Character*2 idstring
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm)
      character term_bl(mterm)*3
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess

      logical 			::clst_disk, clst_memory
      common/use_disk/clst_disk, clst_memory

      character*(128) file_c,file_y,file_cl,file_cfg,file_wo,
     :               file_wi,file_s,file_st,p_name
      integer llc,lly,llcl,llcf,llwi,llwo,lls,llst,lpn,lmpi_dir
      character*(128) mpi_dir,cwd,sh_command,tmpdir
      integer*2 serr


*
*
      DOUBLE PRECISION TIME(2)
      EQUIVALENCE (IUC,IOU(1)),(OUC,IOU(4))
*
*  ***** Define unit numbers and open files *********************
*                                                               *
*        UNIT NUMBERS AND FILE NAMES MAY BE MACHINE             *
*        DEPENDENT. CHECK THE FOLLOWING SECTION.                *
*                                                               *
*        IN - Standard input unit, normally the terminal        *
*        OUT- Standard output unit, normally the terminal       *
*        PRI- Printer output unit or file.                      *
*                                                               *


*      ... initialize mpi
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      write(idstring,'(I2.2)') myid

      IN = 5
      OUT = 6
      PRI = 3
      ISCW = 0

*>>>>>>>> specify files for input output >>>>>>>>>>>>>>>>>
!      call getwd(myid,mpi_dir,lmpi_dir,p_name,lpn)
!      mpi_dir = trim(mpi_dir)
!      serr = getcwd(cwd)
!      lcwd = len_trim(cwd)
!
!      if (serr.ne.0) then
!         print*, 'couldn''t get the current directory, exiting...'
!         call exit(23);
!      end if
    
      cwd = " "; mpi_dir = " "; tmpdir = " ";
      call mpi_work_dir(cwd, mpi_dir, tmpdir);

      lcwd = LEN_TRIM(cwd)  !- 1
      lenperm = LEN_TRIM(mpi_dir) !- 1
      lentmp = LEN_TRIM(tmpdir)  !- 1

      file_c   = tmpdir(1:lentmp)//'/c.lst.'//idstring
      file_y   = tmpdir(1:lentmp)//'/yint.lst.'//idstring
      file_cl  = cwd(1:lcwd)//'/cfg.inp'
      file_cfg = cwd(1:lcwd)//'/cfg.h'
      file_wo  = cwd(1:lcwd)//'/wfn.out'
      file_wi  = cwd(1:lcwd)//'/wfn.inp'
      file_s   = cwd(1:lcwd)//'/summry'
      file_st  = cwd(1:lcwd)//'/std.out'
      llc      = len_trim(file_c)
      lly      = len_trim(file_y)
      llcl     = len_trim(file_cl)
      llcf     = len_trim(file_cfg)
      llwi     = len_trim(file_wi)
      llwo     = len_trim(file_wo) 
      lls      = len_trim(file_s)
      llst     = len_trim(file_st)

      inquire(file = file_wi,exist=ld)
      if (ld) then
        iuf = 22
        open(unit=iuf,file=file_wi(1:llwi),status='old',
     :       form='unformatted')
      else
        iuf = 0
Ctc 19/03/2008 : crash of large calculations on several nodes of hydra@ulb
Ctc     print*, 'WARNING!!!! wfn.inp NOT FOUND!! mchf will continue..'
      WRITE(80+myid,*)
     :            'WARNING!!!! wfn.inp NOT FOUND!! mchf will continue..'
Ctc
      end if
      iuc = 21;iud = 23;ouc = 24;ouf = 25;
      open(unit=39,file=file_c(1:llc),status='unknown',
     :     form='unformatted');
      open(unit=23,file=file_y(1:lly),status='unknown',
     :     form='unformatted');
      open(unit=30,file=file_cl(1:llcl),status='unknown')
      open(unit=29,file=file_cfg(1:llcf),status='unknown')
      if (myid == 0) then
         open(unit=pri,file=file_s(1:lls),status='unknown')
      end if
      if (myid == 0)  then
         open (unit=ouf, file=file_wo(1:llwo), status='unknown',
     :                form = 'unformatted')
!         close (ouf)
      end if
      OPEN(UNIT=35,STATUS='SCRATCH',FORM='UNFORMATTED')

*>>>>>>>>>>>> end files <<<<<<<<<<<<<<<<<<<<<<<<<<
      if (myid.eq.0) then
         ISCW = 0
         IWRITE = 6
      else
         iscw = 40
         iwrite = 40
         open(unit=40,file=file_st(1:llst),status='unknown')
      end if

*<<< header inf >>>
 9    FORMAT(//20X,'==================================='/
     :         25X,' S M C H F_M P I  ... 2000   '/
     :         20X,'==================================='/)

      write(iscw,9)
      write(iscw,'(A,i4,A)')
     :  '                 ...snonh_mpi running on ',
     :                  nprocs,' processors...'
      write(iscw,*)

*
*  *****  WRITE OUT DIMENSION INFORMATION
*
      WRITE(iscw,99) 'NWD',NWD,'NO',NOD,'Lagrange Multipliers',(noffd)
99    FORMAT(//10X,'THE DIMENSIONS FOR THE CURRENT VERSION ARE:'/
     :       (10X,2(A6,'=',I3,4X),A,'=',I3/)/)
      write(iscw, '(/A/A/)') ' START OF CASE',' ============='

*<<<<<< end header inf >>>>>>>>


*  *****  INITIALIZE COMMON DATA ARRAYS
      CALL INITA
      CALL INITR

*     For testing purposes set to false
      lguess = .false.

*  ***** END OF INPUT/OUTPUT INTERFACE **************************
*     set logic variables controlling the usage of memory:
*     setting clst_memory = .true. attempts loading c.lst in memory
*     if spallcts returns NULL pointers and clst_disk = .true. then  
*     smchf continues with reading c.lst in diag.f and updatc.f 
*     
      clst_memory = .true.
      clst_disk = .true.
**

      FAIL = .FALSE.
      IJE(1:noffd) = 0
      if (myid == 0) then;
6       write(iscw,'(/A)') ' ATOM, Z in FORMAT(A, F) : '

        READ(IN,'(A)') STRING
        I = INDEX(STRING,',')
        IF ( I .EQ. 0) THEN
          write(iscw,*)' ATOM, and Z must be separated by commas '
          GO TO 6
        END IF
        ATOM = STRING(1:I-1)
        BUFFER = STRING(I+1:)
        J = INDEX(BUFFER,'.')
        IF (J .EQ. 0) THEN
          J = INDEX(BUFFER,' ')
          IF (J .NE. 0) BUFFER(J:J) = '.'
        END IF
        READ(BUFFER,'(F8.0)') Z
      end if
       
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_BCAST(Z,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ATOM,6,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
*
*  *****  DETERMINE DATA ABOUT THE PROBLEM
*
      call data(IVAR, eigst_weight, cf_tot,ptmp); 

      if (myid == 0 ) then;
        do i = 1, nblock;
          unit_n = 40 + i;
          name_dotL =  cwd(1:lcwd)//'/'//term_bl(i)(1:3)//'.l';
        open(unit_n,file=name_dotL,status='unknown',form='formatted');
        end do;
      end if; 
*
*  *****  SET PARAMETERS TO THEIR DEFAULT VALUE; NO,REL, and STRONG
*         were set in SUBROUTINE DATA
*
13    PRINT = .FALSE.
      CFGTOL = 1.D-8
      SCFTOL = 1.D-7
      NSCF = 200
      IC = 0
      ACFG = D0
      LD = .TRUE.
      TRACE = .FALSE.

      if (myid == 0) then
        write(iscw,'(/A)')'Default values for other parameters? (Y/N) '
        READ (IN,'(A)') ANS
        IF (ANS .EQ. 'N' .OR. ANS .EQ. 'n') THEN
*  *****  ADDITIONAL PARAMETERS
   50   write(iscw,'(/A)') ' Default values (NO,REL,STRONG) ? (Y/N) '
        READ(IN,'(A)') ANS
        IF (ANS .NE. 'Y' .AND. ANS .NE. 'y') THEN
          write(iscw,*) ' Enter values in FORMAT(I,L,L) '
          READ(IN,*) NO, REL, STRONG
          OMIT = .NOT. STRONG
          IF (NO .GT. (220)) THEN
            write(iscw,*)
     :            ' TOO MANY POINTS FOR EACH FUNCTION: MAX=(220)'
            STOP
          END IF
          ND = NO - 2
        END IF
        WRITE(iscw,61) ATOM,Z,NO,NWF,NWF-IB+1,NCFG,REL
61    FORMAT(/1X,A6,F6.0,I6,3I3,L3)
         write(iscw,'(/A,A)')' Default values for',
     :   ' PRINT, CFGTOL, SCFTOL ? (Y/N) '
          READ(IN,'(A)') ANS
          IF ( ANS .NE. 'Y' .AND. ANS .NE. 'y'  ) THEN
             write(iscw,'(A)')' Input free FORMAT(L, F, F) '
             READ(IN,*) PRINT, CFGTOL, SCFTOL
          END IF
          write(iscw,'(/A)') ' Default values for NSCF, IC ? (Y/N) '
          READ(IN,'(A)') ANS
          IF (ANS .NE. 'Y' .AND. ANS .NE. 'y' ) THEN
             write(iscw,'(A)') ' Input free FORMAT(I, I) '
             READ(IN,*) NSCF, IC
          END IF
          write(iscw,'(/A)') 'Default values for ACFG,LD,TRACE? (Y/N) '
          READ(IN,'(A)') ANS
          IF (ANS .NE. 'Y' .AND. ANS .NE. 'y') THEN
             write(iscw,'(A)') ' Input free FORMAT( F, L, L) '
             READ(IN,*) ACFG,LD,TRACE
          END IF
        END IF
       end if

      call MPI_BCAST(NO,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BARRIER(mpi_comm_world,ierr)
      call MPI_BCAST(REL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(STRONG,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(PRINT,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(CFGTOL,1,MPI_DOUBLE_PRECISION,0,
     :               MPI_COMM_WORLD,ierr)
      call MPI_BCAST(SCFTOL,1,MPI_DOUBLE_PRECISION,0,
     :               MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NSCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(IC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ACFG,1,MPI_DOUBLE_PRECISION,0,
     :               MPI_COMM_WORLD,ierr)
      call MPI_BCAST(LD,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(TRACE,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

      CALL SCF(IVAR,ACFG,SCFTOL,CFGTOL,LD,eigst_weight,cf_tot,tmp)

*     .. solve the problem

      if (myid == 0) then
        CALL OUTPUT(PRINT)
        CALL SUMMRY
*      call etime(time)
*6     write(iscw,'(//A/A//A,F8.3,A//)') ' END OF CASE',' ===========',
*     :      ' Total CPU time was ', TIME(1)/60, ' minutes'
        CLOSE(PRI)
      end if

      call MPI_FINALIZE(ierr)

      END
