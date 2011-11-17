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
         integer U;
         read(U) NCLOSD, IBK;
         read(U) B1ELC(1:4);
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
         write(0,'(I8,A)') size_c_1, &
            ' L(i,j) coefficients for initial'
         write(0,'(I8,A)') size_int_1, &
            ' L(i,j) integrals for initial'
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

      END subroutine
!
      subroutine DBG_R(U);
         USE DBG_C;
         integer U;
         read(U) IBUGM;
      END subroutine DBG_R;

      subroutine DEBUG_R(U);
         USE DEBUG_C;
         integer U;
         read(U) IBUG1, IBUG2, IBUG3, NBUG6, NBUG7, IFUL
      END subroutine

      subroutine diagnl_R(U);
      USE diagnl_C;
         integer U;
         read(U) IDIAG, JA, JB
      END subroutine

      subroutine dimen_R(U);
         USE dimen_C;
         integer U;
         read(U) kfl1, kfl2, kfl3, kfl4, kfl5, kfl6, kfl7, mxihsh
      END subroutine

      subroutine ELT_R(U);
         USE ELT_C;
         integer U,I;
         read(U) ELTENS(1:NWD);
      END subroutine

      subroutine EMS_R(U);
         USE EMS_C;
         integer U;
         read(U) IEM(1:4);
         read(U) IFL, JI, JF, LAM
         read(U) REL, VOK
      END subroutine

      subroutine fact_R(U);
         USE fact_C;
         integer U;
         read(U) GAM(1:100);
      END subroutine

      subroutine FOUT_R(U);
         USE fout_C;
         integer U;
         read(U) LCOUNT, NREC, IFLAG, LIJ, NIJ
      END subroutine

      subroutine fo_R(U);
         USE fo_C;
         integer U;
         read(U) IWF(1:2),NCLOS(1:2),NOPEN(1:2);
         read(U) MAXNFO
      END subroutine

      subroutine INFORM_R(U);
      USE INFORM_C;
         integer U;
         read(U) JSC(1:3);
         read(U) IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,JSC0,ISCW; 
      END subroutine

      subroutine INOUT_R(U);
         USE INOUT_C;
         integer U;
         read(U) IUC(1:2), IUW(1:2), IUL(1:2), IUJ(1:2), IUT(1:2)
         read(U) IREAD, IWRITE, ISCW
      END subroutine

      subroutine MEDEFN_R(U);
         USE MEDEFN_C;
         integer U;
         read(U) NJ(1:16),LJ(1:16),NOSH1(1:16),NOSH2(1:16),IJFUL(1:16)
         read(U) J1QN1(1:31,1:3),J1QN2(1:31,1:3)
         read(U) IHSH
      END subroutine

      subroutine ndims_R(U);
! used only in biotr_ang
         USE NDIMS_C;
         integer U;
         read(U) NCFG;
         allocate(NOCCSH(ncfg));
         allocate(NELCSH(8,ncfg));
         allocate(NOCORB(8,ncfg));
         allocate(J1QNRD(15,ncfg));
         read(U) NOCCSH(1:NCFG);
         read(U) NELCSH(1:8,1:NCFG),NOCORB(1:8,1:NCFG);
         read(U) J1QNRD(1:15,1:NCFG);
      END subroutine

      subroutine nel_R(U);
         !USE NEL_C, not comuted in biotr_ang;
         integer U;
      END subroutine

      subroutine nlorb_R(U);
         USE nlorb_C;
         integer U;
         read(U) nq1(1:NWD), lq1(1:NWD), nq2(1:NWD), lq2(1:NWD);
      END subroutine

      subroutine non30_R(U);
         USE non30_C;
         integer U;
         read(U) MAXORB;
         allocate(IAJCMP(MAXORB));
         allocate(LJCOMP(MAXORB));
         allocate(NJCOMP(MAXORB));
         allocate(IAJCLD(NWCD));
         allocate(LJCLSD(NWCD));
         read(U) IAJCMP(1:MAXORB),LJCOMP(1:MAXORB),NJCOMP(1:MAXORB);
         read(U) IAJCLD(1:NWCD),LJCLSD(1:NWCD);
      END subroutine

      subroutine NOR_R(U);
         USE NOR_C;
         integer U;
         read(U) NCOM, NCLOSI, NCLOSF, NORBI, NORBF, IWAR;
      END subroutine

      subroutine NTRM_R(U);
         USE NTRM_C;
         integer U;
         read(U) NTERMS
      END subroutine

      subroutine ovrlap_R(U);
         USE OVRLAP_C;
         integer U;
         read(U) MU, NU, MUP, NUP, NONORT, NOVLPS, IROWMU, &
                  IROWNU, ICOLMU, ICOLNU, NORTH, IORDER, NCALLS, &
                  LMU, LNU, LMUP, LNUP, JMU, JNU, JMUP, JNUP
      END subroutine

      subroutine param_R(U);
         USE PARAM_C;
         integer U;
         read(U) NO, ND, NWF, MASS, NCFG, IB, IC, ID, NSCF, NCLOSD
         read(U) H, H1, H3, CH, EH, RHO, Z, TOL, D0, D1, D2, D3, D4, D5, &
                 D6, D8, D10, D12, D16, D30, FINE, RMASS 
      END subroutine

      subroutine RAS_R(U);
         USE RAS_C;
         integer U;
         read(U) NINAC(1:11,1:2);
         read(U) NRAS1(1:11,1:10,1:2);
         read(U) NRAS2(1:11,1:2);
         read(U) NRAS3(1:11,1:10,1:2)
         read(U) NL(1:11,1:2);
         read(U) NGAS1(1:2);
         read(U) NGAS3(1:2);
         read(U) ITAB(1:NWD,1:2);
         read(U) ELRAS(1:2*NWD);
         read(U) ELRASI(1:nwd);
         read(U) ELRASF(1:nwd);
         read(U) LMAX(1:2);
      END subroutine

      subroutine RED_R(U);
         USE RED_C;
         integer U;
         read(U) INDL(1:20), INDR(1:20), IORST(1:20);
         read(U) LIS(1:16);
         read(U) NCI1, KPL, KPR, NC, JIST, JFST, NCONTR
      END subroutine

      subroutine SIGNF_R(U);
         USE SIGNF_C;
         integer U;
         read(U) SIGNFA;
      END subroutine

      subroutine state_R(U);
         use state_C;
         use inout_C;
         integer U, n1, n2;
         character*39 RR;

         read(U) NCF(1:2)

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


      END subroutine

      subroutine TRK_R(U);
         USE TRK_C;
         integer U;
         read(U) ID1(1:7), ID2(1:7), IK1(1:7), IK2(1:7);
         read(U) BD1(1:3), BD2(1:3), BK1(1:3), BK2(1:3);
      END subroutine


