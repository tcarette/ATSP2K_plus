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
*     Modified by:
*     G. Tachiev for simultaneous optimization of multiple blocks
*     August, 2000
*                  C O P Y R I G H T   2000
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
      logical 			:: leigen(meig,mterm), lguess
      integer			:: unit_n
      character (len=5)		:: name_dotL
      double precision 			:: eigst_weight(meig,mterm)   
      integer 			:: cf_tot(mterm) 
      Character*2 idstring
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm),
     :        nze_max(mterm)
      character term_bl(mterm)*3
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess,nze_max

      logical	:: hmx_memory, ico_memory, ih_memory, clst_memory
      common/memory_use/hmx_memory, ico_memory, ih_memory, clst_memory

*
*
      REAL TIME(2)
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

      IN = 5
      OUT = 6
      PRI = 3
      ISCW = 0

      OPEN(UNIT=PRI,FILE='summry',STATUS='UNKNOWN')

      ISCW = 0
      IWRITE = 6

*
*  *****  WRITE OUT HEADER
*
 9    FORMAT(//20X,'======================='/
     :         25X,'   M C H F  ... 2000   '/
     :         20X,'======================='/)
 
      write(iscw,9)
      write(iscw,*)

*
*  *****  WRITE OUT DIMENSION INFORMATION
*
      WRITE(iscw,99) 'NWD',NWD,'NO',NOD,'Lagrange Multipliers',(noffd)
99    FORMAT(//10X,'THE DIMENSIONS FOR THE CURRENT VERSION ARE:'/
     :       (10X,2(A6,'=',I3,4X),A,'=',I3/)/)
*
*  *****  INITIALIZE COMMON DATA ARRAYS
*

      CALL INITA
      CALL INITR
*                                                               *
*  ***** IN THE OPEN STATEMENTS CHECK FOR VALID FILE NAMES ******
*                                                               *

1     write(iscw, '(/A/A/)') ' START OF CASE',' ============='

      iuc = 21
      inquire(file='wfn.inp',exist=ld)
      if (ld) then
        iuf = 22
        open(unit=iuf,file='wfn.inp',status='old',form='unformatted')
      else
        iuf = 0
      end if
      iud = 23
*      open(unit=iud,file='yint.lst',status='old',form='unformatted')
*      open(unit=11,file='ih.lst',status='old',form='unformatted')
      open(unit=12,file='ico.lst',status='old',form='unformatted')
      open(unit=23,file='yint.lst',status='unknown',
     :     form='unformatted');
      ouc = 24
      ouf = 25
      open(unit=ouf,file='wfn.out',status='unknown',form='unformatted')
*     file cfg.h contains header information for the block
      open(unit=29,file ='cfg.h',status='old',form = 'formatted')
      open(unit=30,file ='cfg.inp',status='old',form = 'formatted')
      OPEN(unit=39,file='c.lst',status='old', form='unformatted')
      OPEN(UNIT=35,STATUS='SCRATCH',FORM='UNFORMATTED')

      lguess = .false.

*  ***** END OF INPUT/OUTPUT INTERFACE **************************
*

      FAIL = .FALSE.
      IJE(1:noffd) = 0
6     write(iscw,'(/A)') ' ATOM, Z in FORMAT(A, F) : '

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

*
*  *****  DETERMINE DATA ABOUT THE PROBLEM
*
*  initialize logical variables for memory use
      hmx_memory = .true.;
      ico_memory = .true.;
      ih_memory = .true.;
      clst_memory = .true.;

      call data(IVAR,eigst_weight,cf_tot); 

      do i = 1, nblock;
        unit_n = 40 + i;
        name_dotL = term_bl(i)//'.l';
        open(unit_n,file=name_dotL,status='unknown',form='formatted');
      end do;
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

*     .. solve the problem
      CALL SCF(IVAR,ACFG,SCFTOL,CFGTOL,LD,eigst_weight,cf_tot)
      close(11);
      close(12);
      close(13);
      CALL OUTPUT(PRINT);
      CALL SUMMRY;
*
*      call etime(time)
*6     write(iscw,'(//A/A//A,F8.3,A//)') ' END OF CASE',' ===========',
*     :      ' Total CPU time was ', TIME(1)/60, ' minutes'
      CLOSE(PRI)

      END
