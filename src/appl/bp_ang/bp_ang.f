*     ==================================================================
*
*      THIS PROGRAM EVALUATES THE BREIT-PAULI OPERATORS
*
*     ONE-ELECTRON OPERATOR
*     ELECTROSTATIC INTERACTION
*     SPIN-ORBIT INTERACTION
*     SPIN-OTHER-ORBIT INTERACTION
*     SPIN-SPIN INTERACTION
*
*       by A. Hibbert, R. Glass, and C. Froese Fischer
*       Comput. Phys. Commun. 64 (1991) 455
*
*       Modified for Gediminas g-angular structure
*       Modified for Unsorted version
*       Modified for restarting
*                               [Sep 1994]
*       Modified by G. Gaigalas for new libraries and orbit-orbit (1997)
*       Modified by G. Tachiev for sparse matrix data structures (2000)
*
*
*     ==================================================================
*
*       The PARAMETER values in this program define the following:
*               LSDIM - Number of coeffficients in a block
*               NTERMD- Maximum number of terms
*
*     ------------------------------------------------------------------
      PROGRAM bp_ang 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NTERMD=31)
*
*

      character*1 idstring
      character*128 startdir, permdir, tmpdir
      character*128 file_cint,file_ih,file_ico,file_clst,file_atomc

      CHARACTER*26 STRING,ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW
      CHARACTER QREL, QMASS
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(8),ISCW

      CHARACTER*128 input
      EXTERNAL INITT
      logical REL
*
*     Set machine dependent unit numbers
*               IREAD - The unit number of standard input
*               IWRITE- The unit number of standard output
*               ISCW  - The unit number of standard error (screen)
*               iuc   - The unit number for the <name>.c file
*
      iread = 4
      IOUT = 8
      in = 5

         ISCW = 0
         IWRITE = 6
      idstring = 's'
      
      startdir = '.'; permdir = '.'; tmpdir = '.';
      lenstart = 1; lenperm = 1; lentmp = 1;

      call inp_atom(ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW,
     :      REL,MASS,iscw,in)

      file_cint  = tmpdir(1:lentmp)//'/cint.lst.'//idstring
      file_ih    = tmpdir(1:lentmp)//'/ih.lst.'//idstring
      file_ico   = tmpdir(1:lentmp)//'/ico.lst.'//idstring
      file_clst  = tmpdir(1:lentmp)//'/c.lst.'//idstring
      file_atomc = permdir(1:lenperm)//'/'//atomc
      lcint  = len_trim(file_cint);
      lih    = len_trim(file_ih);
      lico   = len_trim(file_ico);
      lclst  = len_trim(file_clst);
      latomc = len_trim(file_atomc);

      OPEN(UNIT=IREAD,FILE=file_atomc(1:latomc),STATUS='UNKNOWN')
      OPEN(UNIT=IOUT, FILE=file_cint(1:lcint), STATUS='UNKNOWN',
     :     FORM='unformatted')
      OPEN(UNIT=11, FILE=file_ih(1:lih), STATUS='UNKNOWN',
     :     FORM='unformatted')
      OPEN(UNIT=12, FILE=file_ico(1:lico), STATUS='UNKNOWN',
     :     FORM='unformatted')
      OPEN (UNIT=50,FILE=file_clst(1:lclst), STATUS='UNKNOWN',
     :     FORM='UNFORMATTED')
*
*     Let's postpone the restart for the time being
*
      CALL INITA
      CALL BREVALA(MASS)

*      call etime(time)
*      write(iscw,'(//A/A//A,F8.3,A//)')
*     :       ' END OF CASE',' ===========',
*     :       ' Total time was ', time(1)/60, ' minutes' 
      END program bp_ang
