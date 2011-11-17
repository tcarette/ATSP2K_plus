      subroutine inp_data(rel,skip,ec,nze,njv,minj,maxj,pflsj)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,NWCD=20,NWD=60,NTERMD=31,LINT=500)
      INCLUDE 'mpif.h'
      parameter (MAXPROC=9999)
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)

      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,iout,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk

      LOGICAL rel,skip
      logical onlydvd

      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MSOO,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS

      POINTER (qjptr, jptr(1))
      POINTER (IQLSP, LSP(1))
      CHARACTER*3 EL(NWD)
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON /NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)

      POINTER (qintptr,intptr(1)),(qpackn,ipackn(1)),(qvalue,value(1))
      COMMON /TABLE/qintptr,qpackn,qvalue,lmax
      COMMON /NCFGS/ IQLSP,index(ntermd),IJK,nterm,termsh(ntermd)
      pointer(pflsj,nze_flsj(1,1,1,1));

!      print *, 'iout,skip', iout,skip
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

      call READW(el,rel,ec,nnel,skip,onlydvd,iounew)

      call flsj_comp(ncfg,nterm,index,njv,maxj,minj,
     :               pflsj,lsp)

      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      return
      end
