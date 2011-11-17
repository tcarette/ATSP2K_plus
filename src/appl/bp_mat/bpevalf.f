***********************************************************************
*
*
* --- This SUBROUTINE computes the files for Breit-Pauli matrix elements
*     given angular data
*
*     The matrix is decomposed into three parts:
*       Hnr -- non-relativistic integrals 
*       HZ  -- contributions from Z, N, and V integrals
*       HS  -- contributions from S integrals.
*
*
***********************************************************************
*
      Subroutine bpevalf(rel,skip,ec,nze,restart,onlydvd,iounew)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,NWCD=20,NWD=60,NTERMD=31,LINT=500)
*
      CHARACTER SYMBOL*11
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
      CHARACTER*3 EL(NWD)
      LOGICAL rel,skip
      logical onlydvd
      EXTERNAL INITT

      IF (IREL.NE.1) ICOLOM=1
      skip = .true.
      if (irel .ne. 0) then
         skip = .false.
      end if
      IFULL = 0; ICSTAS = 1; ISPORB = 0; ISOORB = 0
      ISPSPN = 0; IORBORB=0

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
