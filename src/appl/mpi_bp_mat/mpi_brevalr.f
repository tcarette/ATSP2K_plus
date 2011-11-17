***********************************************************************
*
*
* --- This SUBROUTINE computes the Breit-Pauli matrix elements
*     for a variety of operators:
*
*     ONE-ELECTRON OPERATOR (L-integrals)
*     ELECTROSTATIC INTERACTION (F, G, R) integrals)
*     SPIN-ORBIT INTERACTION  (Z-integrals)
*     SPIN-OTHER-ORBIT INTERACTION (N, V integrals)
*     SPIN-SPIN INTERACTION (S-integrals)
*
*     The matrix is decomposed into three parts:
*       Hnr -- non-relativistic integrals 
*       HZ  -- contributions from Z, N, and V integrals
*       HS  -- contributions from S integrals.
*
*
***********************************************************************
*
      Subroutine BREVALR(rel,skip,ec,nze,restart,onlydvd,iounew)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,NWCD=20,NWD=60,NTERMD=31,LINT=500)
*     MPI stuff ***********************************************
*
        INCLUDE 'mpif.h'
        parameter (MAXPROC=9999)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
****************************************************************

*
      CHARACTER ANS*2, SYMBOL*11
      DATA SYMBOL/'SPDFGHIKLMN'/
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON/IMAGNT/ IREL,ISTRICT,IELST
CG 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MSOO,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
CG 
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,iout,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk

      LOGICAL restart
      POINTER (qintptr,intptr(1)),(qpackn,ipackn(1)),(qvalue,value(1))
      COMMON /TABLE/qintptr,qpackn,qvalue,lmax
*     .. radial common
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON /NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
      POINTER (IQLSP,LSP(1))
      COMMON /NCFGS/ IQLSP,index(ntermd),IJK,flsj(ntermd,ntermd,2),
     :               termsh(ntermd),nterm
      POINTER (qjptr, jptr(1))
      DIMENSION NFLG(20)
      CHARACTER*3 EL(NWD)
      LOGICAL rel,skip
      logical onlydvd
      EXTERNAL INITT
  105 FORMAT (49H ISPORB=0 AND ISOORB=1 CAUSES THE PROGRAM TO FAIL,
     :  34H BECAUSE THE BLUME WATSON FORMULAE,/
     :  50H CANNOT BE USED FOR CLOSED SUBSHELLS.  TO OVERCOME,
     :  34H THIS, THE CODE HAS SET ISPORB = 1//)
   11 FORMAT(////' THE TYPE OF CALCULATION IS DEFINED BY ',
     :  'THE FOLLOWING PARAMETERS - '/
     : 5X,22H BREIT-PAULI OPERATORS,13X,8HIREL   =,I2/
     : 5X,27H PHASE CONVENTION PARAMETER,8X,8HICSTAS =,I2/)
   13 FORMAT(40H RELATIVISTIC OPERATORS TO BE INCLUDED -/5X,13H SPIN-ORB
     :IT (,I1,22H),  SPIN-OTHER-ORBIT (,I1,15H),  SPIN-SPIN (,I1,1H)/)
   24 FORMAT(36H0INITIAL DEBUG: IN 1-ELECTRON PART =,I2,2H ,,5X,
     : 20H IN 2-ELECTRON PART=,I2/16X,23HIN RECOUPLING PACKAGE =,I2/)
   42 FORMAT(//' MATRIX ELEMENTS CONSTRUCTED USING ',
     :       'THE SPHERICAL HARMONIC PHASE CONVENTION OF'/)
   43 FORMAT(/16X,47HCONDON AND SHORTLEY, THEORY OF ATOMIC STRUCTURE/
     :16X,47H-----------------------------------------------///)
   44 FORMAT(25X,42HFANO AND RACAH, IRREDUCIBLE TENSORIAL SETS/25X,42H--
     :----------------------------------------///)
   50 FORMAT(/20X,'============================'/
     :        20X,' B R E I T - P A U L I   C I'/
     :        20X,'============================'/)
   78 FORMAT(19H DEBUG PARAMETERS -/5X,16H NBUG6(TENSOR) =,I2/5X,
     : ' NBUG7(RELATIVISTIC OPERATORS - SO,SOO,SS) =',I2//)
*
*     
!      print *, 'iout,skip', iout,skip
!      if (myid == 0) then
!         call inp_type(irel,iread,iscw,ISPORB,ISOORB,ISPSPN,
!     :                   IORBORB,iwrite)
!      end if

      call MPI_BCAST(irel,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      IF (IREL.NE.1) ICOLOM=1
      skip = .true.
      if (irel .ne. 0) then
         skip = .false.
      end if
      IFULL = 0; ICSTAS = 1; ISPORB = 0; ISOORB = 0
      ISPSPN = 0; IORBORB=0

      call MPI_BCAST(ISPORB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(ISOORB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(ISPSPN,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(IORBORB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);

      read(iout) MSOO,nclosd, maxorb, lsdim, ncfg
      nwf = maxorb+nclosd
      call alctab(ncfg,nwf,skip)
      call alloc(qjptr,ncfg,4)
      call alloc(iqlsp, ncfg, 4)
      call alloc(qiajcmp,maxorb, 4)
      if(nclosd.ne.0) call alloc(qiajcld, nclosd, 4)
      read(iout) l(1:nwf)
!      print *, 'l(:)', l(1:nwf)
      read(iout) el(1:nwf)
!      print *, 'el(i)', el(1:nwf)
      read(iout) lmax,nnel, skip
!      print *, 'lmax,nnel,skip',lmax,nnel,skip
      read(iout) lsp(1:ncfg),jptr(1:ncfg)
!      print *,'lsp', lsp(1:ncfg),jptr(1:ncfg)
      read(iout) nterm, index(1:nterm)
!      print *, 'nterm',nterm, index(1:nterm)
*
*     .. read radial functions, 
      call readw(el,rel,ec,nnel,skip,.false.,iounew)
*
*     .. generate list of configurations
      call genintbr(nclosd,maxorb,lmax,iql,qintptr,qpackn,qvalue,
     :                 rel,skip,.true.)

      call genlst(ncfg,jptr)

      Return
      END
