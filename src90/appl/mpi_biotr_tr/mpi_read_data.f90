!
!     ------------------------------------------------------------------
!     read the program state of biotr_ang 
!     ------------------------------------------------------------------
      subroutine read_data(U)
      integer U;

      call closed_R(U);
      call CONSTS_R(U);
!      call CSF_R(U);
      call DBG_R(U);
      call DEBUG_R(U);
!      call diagnl_R(U);
!      call dimen_R(U);
      call ELT_R(U);
!      call EMS_R(U);
      call fact_R(U);
!      call FOUT_R(U);
      call fo_R(U);
      call INFORM_R(U);
      call INOUT_R(U);
      call MEDEFN_R(U);
      call ndims_R(U);
!      call nel_R(U);
      call nlorb_R(U);
      call non30_R(U);
      call NOR_R(U);
      call NTRM_R(U);
      call ovrlap_R(U);
      call param_R(U);
      call RAS_R(U);
      call RED_R(U);
      call SIGNF_R(U);
      call state_R(U);
!      call TRK_R(U);
      close(U);
      END subroutine read_data

! Define the routines saving data;
      subroutine closed_R(U);
         USE closed_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) NCLOSD, IBK;
         if (myid==0) read(U) B1ELC(1:4);
         call MPI_BCAST(NCLOSD,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IBK,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(B1ELC,4,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
      END subroutine closed_R

      subroutine CONSTS_R(U);
         ! the constants are in data statements;
         integer U;
      END subroutine consts_R;

! read integral data
      subroutine CSF_R(U);
         USE CSF_C

         integer U,IERR;

         read(U) size_c_1;
         read(U) size_int_1;
         nint1=size_int_1;
         ncoef1=size_c_1;
!         if (myid==0) write(0,'(I8,A)') size_c_1, &
!            ' L(i,j) coefficients for initial'
!         if (myid==0) write(0,'(I8,A)') size_int_1, &
!            ' L(i,j) integrals for initial'
         allocate(intptr1(size_int_1),STAT=IERR);
         if (ierr.ne.0) call mem_fail(6,size_int_1,'csf_R::intptr1',ierr);
         allocate(intgrl1(size_int_1),STAT=IERR);
         if (ierr.ne.0) call mem_fail(6,size_int_1,'csf_R::intgrl1',ierr);
         write(0,'(A)') ' ...Allocating data for cnn1, jann1, jbnn1 '
         allocate(cnn1(size_c_1),STAT=IERR);
         if (ierr.ne.0) call mem_fail(6,size_c_1,'csf_R::cnn1',ierr);
         allocate(jann1(size_c_1),STAT=IERR);
         if (ierr.ne.0) call mem_fail(6,size_c_1,'csf_R::jann1',ierr);
         allocate(jbnn1(size_c_1),STAT=IERR);
         if (ierr.ne.0) call mem_fail(6,size_c_1,'csf_R::jbnn1',ierr);
         write(0,'(A)') ' ...Reading cnn1, jann1, jbnn1 '
         read(U) intptr1(1:size_int_1),intgrl1(1:size_int_1);
         read(U) cnn1(1:size_c_1),jann1(1:size_c_1),jbnn1(1:size_c_1);
! associate pointers with arrays (not really needed in this code);
         qintgrl1=>intgrl1;
         qintptr1=>intptr1;
         qcnn1=>cnn1;
         qjann1=>jann1;
         qjbnn1=>jbnn1;

         read(U) size_c_2;
         read(U) size_int_2;
         ncoef2=size_c_2;
         nint2=size_int_2;
         write(0,'(I8,A)') size_c_2, &
            ' L(i,j) coefficients for final'
         write(0,'(I8,A)') size_int_2, &
            ' L(i,j) integrals  for final'
         allocate(intptr2(size_int_2),STAT=IERR);
         if (ierr.ne.0) call mem_fail(6,size_int_2,'csf_R::intptr2',ierr);
         allocate(intgrl2(size_int_2),STAT=IERR);
         if (ierr.ne.0) call mem_fail(6,size_int_2,'csf_R::intgrl2',ierr);
         write(0,'(A)') ' ...Allocating data for cnn2, jann2, jbnn2 '
         allocate(cnn2(size_c_2),STAT=IERR);
         if (ierr.ne.0) call mem_fail(6,size_c_2,'csf_R::cnn2',ierr);
         allocate(jann2(size_c_2),STAT=IERR);
         if (ierr.ne.0) call mem_fail(6,size_c_2,'csf_R::jann2',ierr);
         allocate(jbnn2(size_c_2),STAT=IERR);
         if (ierr.ne.0) call mem_fail(6,size_c_2,'csf_R::jbnn2',ierr);
         write(0,'(A)') ' ...Reading cnn2, jann2, jbnn2 '
         read(U) intptr2(1:size_int_2),intgrl2(1:size_int_2);
         read(U) cnn2(1:size_c_2),jann2(1:size_c_2),jbnn2(1:size_c_2);
! associate pointers with arrays (not really needed in this code);
         qintgrl2=>intgrl2;
         qintptr2=>intptr2;
         qcnn2=>cnn2;
         qjann2=>jann2;
         qjbnn2=>jbnn2;

      END subroutine csf_R
!
      subroutine DBG_R(U);
         USE DBG_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
      if (myid==0) read(U) IBUGM;
      call MPI_BCAST(IBUGM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine DBG_R;

      subroutine DEBUG_R(U);
         USE DEBUG_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) IBUG1, IBUG2, IBUG3, NBUG6, NBUG7, IFUL
      call MPI_BCAST(IBUG1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(IBUG2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(IBUG6,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(IBUG7,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(IFUL,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine diagnl_R(U);
      USE diagnl_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) IDIAG, JA, JB
      call MPI_BCAST(IDIAG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(JA,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(JB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine dimen_R(U);
         USE dimen_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) kfl1, kfl2, kfl3, kfl4, kfl5, kfl6, kfl7, mxihsh
      call MPI_BCAST(kfl1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(kfl2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(kfl3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(kfl4,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(kfl5,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(kfl6,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(kfl7,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      call MPI_BCAST(mxihsh,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine ELT_R(U);
         USE ELT_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U,I;
         if (myid==0) read(U) ELTENS(1:NWD);
         call MPI_BCAST(ELTENS,3*NWD,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine EMS_R(U);
         USE EMS_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) IEM(1:4);
         if (myid==0) read(U) IFL, JI, JF, LAM
         if (myid==0) read(U) REL, VOK
       call MPI_BCAST(IFL,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
       call MPI_BCAST(JI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
       call MPI_BCAST(JF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
       call MPI_BCAST(LAM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
       call MPI_BCAST(REL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr);
       call MPI_BCAST(VOK,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr);
       call MPI_BCAST(IEM,4,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine fact_R(U);
         USE fact_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) GAM(1:100);
         call MPI_BCAST(GAM,100,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine FOUT_R(U);
         USE fout_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) LCOUNT, NREC, IFLAG, LIJ, NIJ
         call MPI_BCAST(GAM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NREC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IFLAG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(LIJ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NIJ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine fo_R(U);
         USE fo_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) IWF(1:2),NCLOS(1:2),NOPEN(1:2);
         if (myid==0) read(U) MAXNFO
         call MPI_BCAST(IWF,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NCLOS,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NOPEN,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(MAXNFO,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine INFORM_R(U);
      USE INFORM_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) JSC(1:3);
         if (myid==0) read(U) IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,JSC0,ISCW; 
         call MPI_BCAST(JSC,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IREAD,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IWRITE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IOUT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(ISC1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(ISC2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(ISC3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(JSC0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(ISCW,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine INOUT_R(U);
         USE INOUT_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) IUC(1:2), IUW(1:2), IUL(1:2), IUJ(1:2), IUT(1:2)
         if (myid==0) read(U) IREAD, IWRITE, ISCW
         call MPI_BCAST(IUC,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IUW,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IUL,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IUJ,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IUT,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IREAD,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IWRITE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(ISCW,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine MEDEFN_R(U);
         USE MEDEFN_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) NJ(1:16),LJ(1:16),NOSH1(1:16), &
                              NOSH2(1:16),IJFUL(1:16)
         if (myid==0) read(U) J1QN1(1:31,1:3),J1QN2(1:31,1:3)
         if (myid==0) read(U) IHSH
         call MPI_BCAST(NJ,16,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(LJ,16,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NOSH1,16,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NOSH2,16,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IJFUL,16,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(J1QN1,31*3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(J1QN2,31*3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine ndims_R(U);
! used only in biotr_ang
         USE NDIMS_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) NCFG;
         call MPI_BCAST(NCFG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         allocate(NOCCSH(ncfg));
         allocate(NELCSH(8,ncfg));
         allocate(NOCORB(8,ncfg));
         allocate(J1QNRD(15,ncfg));
         if (myid==0) read(U) NOCCSH(1:NCFG);
         if (myid==0) read(U) NELCSH(1:8,1:NCFG),NOCORB(1:8,1:NCFG);
         if (myid==0) read(U) J1QNRD(1:15,1:NCFG);
         call MPI_BCAST(NOCCSH,NCFG,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NELCSH,8*NCFG,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NOCORB,8*NCFG,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(J1QNRD,15*NCFG,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine nel_R(U);
         !USE NEL_C, not comuted in biotr_ang;
         integer U;
      END subroutine

      subroutine nlorb_R(U);
         USE nlorb_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) nq1(1:NWD), lq1(1:NWD), nq2(1:NWD), lq2(1:NWD);
         call MPI_BCAST(nq1,NWD,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(lq1,NWD,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(nq2,NWD,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(lq1,NWD,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine non30_R(U);
         USE non30_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) MAXORB;
         call MPI_BCAST(MAXORB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         allocate(IAJCMP(MAXORB));
         allocate(LJCOMP(MAXORB));
         allocate(NJCOMP(MAXORB));
         allocate(IAJCLD(NWCD));
         allocate(LJCLSD(NWCD));
         if (myid==0) read(U) IAJCMP(1:MAXORB),LJCOMP(1:MAXORB), &
                              NJCOMP(1:MAXORB);
         if (myid==0) read(U) IAJCLD(1:NWCD),LJCLSD(1:NWCD);
         call MPI_BCAST(IAJCMP,MAXORB,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(LJCOMP,MAXORB,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NJCOMP,MAXORB,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IAJCLD,NWCD,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(LJCLSD,NWCD,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine NOR_R(U);
         USE NOR_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) NCOM, NCLOSI, NCLOSF, NORBI, NORBF, IWAR;
         call MPI_BCAST(NCOM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NCLOSI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NCLOSF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NORBI,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NORBF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IWAR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine NTRM_R(U);
         USE NTRM_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) NTERMS
         call MPI_BCAST(NTERMS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine ovrlap_R(U);
         USE OVRLAP_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) MU, NU, MUP, NUP, NONORT, NOVLPS, IROWMU, &
                  IROWNU, ICOLMU, ICOLNU, NORTH, IORDER, NCALLS, &
                  LMU, LNU, LMUP, LNUP, JMU, JNU, JMUP, JNUP

         call MPI_BCAST(MU,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NU,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(MUP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NUP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NONORT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NOVLPS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IROWMU,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IROWNU,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(ICOLMU,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(ICOLNU,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IORDER,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IORDER,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NCALLS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(LMU,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(LNU,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(LMUP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(LNUP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(JMU,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(JNU,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(JMUP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(JNUP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine param_R(U);
         USE PARAM_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) NO, ND, NWF, MASS, NCFG, IB, IC, ID, NSCF, NCLOSD
         if (myid==0) read(U) H, H1, H3, CH, EH, RHO, Z, TOL, D0, D1, &
                 D2, D3, D4, D5, D6, D8, D10, D12, D16, D30, FINE, RMASS 
         call MPI_BCAST(NO,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(ND,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NWF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(MASS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NCFG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(ID,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NSCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NCLOSD,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(H,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(H1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(H3,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(CH,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(EH,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(RHO,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(Z,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(TOL,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(D0,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(D1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(D2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(D3,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(D4,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(D5,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(D6,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(D7,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(D8,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(D10,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(D12,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(D16,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(D30,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(FINE,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(RMASS,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine RAS_R(U);
         USE RAS_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) NINAC(1:11,1:2);
         if (myid==0) read(U) NRAS1(1:11,1:10,1:2);
         if (myid==0) read(U) NRAS2(1:11,1:2);
         if (myid==0) read(U) NRAS3(1:11,1:10,1:2)
         if (myid==0) read(U) NL(1:11,1:2);
         if (myid==0) read(U) NGAS1(1:2);
         if (myid==0) read(U) NGAS3(1:2);
         if (myid==0) read(U) ITAB(1:NWD,1:2);
         if (myid==0) read(U) ELRAS(1:2*NWD);
         if (myid==0) read(U) ELRASI(1:nwd);
         if (myid==0) read(U) ELRASF(1:nwd);
         if (myid==0) read(U) LMAX(1:2);
         call MPI_BCAST(NINAC,2*11,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NRAS1,2*11*10,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NRAS2,2*11,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NRAS3,2*11*10,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NL,2*11,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NGAS1,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NGAS3,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(ITAB,2*NWD,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(ELRAS,3*2*NWD,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(ELRASI,3*NWD,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(ELRASF,3*NWD,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(LMAX,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine RED_R(U);
         USE RED_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) INDL(1:20), INDR(1:20), IORST(1:20);
         if (myid==0) read(U) LIS(1:16);
         if (myid==0) read(U) NCI1, KPL, KPR, NC, JIST, JFST, NCONTR
         call MPI_BCAST(INDL,20,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(INDR,20,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IORST,20,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(LIS,16,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NCI1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(KPL,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(KPR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NC,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(JIST,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(JFST,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(NCONTR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine SIGNF_R(U);
         USE SIGNF_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<        
         integer U;
         if (myid==0) read(U) SIGNFA;
         call MPI_BCAST(SIGNFA,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr);
      END subroutine

      subroutine state_R(U);
         use state_C;
         use inout_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U, n1, n2;
         character*39 RR;

         if (myid==0) read(U) NCF(1:2)
         call MPI_BCAST(NCF,2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         if (myid == 0) then
            read(iuj(1),'(A39,i7)') RR, n1
            rewind(iuj(1));
            if (n1.ne.ncf(1)) then
               write(0,'(A,A,i7,i7)') &
               ' Error! .j and .c initial state, have different sizes', & 
               ' exitting...', ncf(1), n1
               call exit(1);
            endif
   
            read(iuj(2),'(A39,i7)') RR, n2
            rewind(iuj(2));
            if (n2.ne.ncf(2)) then
               write(0,'(A,A,i7,i7)') &
               ' Error! .j and .c initial state, have different sizes', &
               ' exitting...', ncf(2), n2
               call exit(1);
            endif
          end if
          close(iuj(1));
          close(iuj(2));
      END subroutine

      subroutine TRK_R(U);
         USE TRK_C;
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
! >>>>>>     <<<<<<<
         integer U;
         if (myid==0) read(U) ID1(1:7), ID2(1:7), IK1(1:7), IK2(1:7);
         if (myid==0) read(U) BD1(1:3), BD2(1:3), BK1(1:3), BK2(1:3);
         call MPI_BCAST(IDL,7,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(ID2,7,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IK1,7,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(IK2,7,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(BD1,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(BD1,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(BD1,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
         call MPI_BCAST(BD1,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr);
      END subroutine


