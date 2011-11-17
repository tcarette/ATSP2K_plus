!
!     ------------------------------------------------------------------
!     save the program state of biotr_ang 
!     ------------------------------------------------------------------
      subroutine save_data(U)
      integer U;

      call closed_W(U);
      call CONSTS_W(U);
!      call CSF_W(U);
      call DBG_W(U);
      call DEBUG_W(U);
!      call diagnl_W(U);
!      call dimen_W(U);
      call ELT_W(U);
!      call EMS_W(U);
      call fact_W(U);
!      call FOUT_W(U);
      call fo_W(U);
      call INFORM_W(U);
      call INOUT_W(U);
      call MEDEFN_W(U);
      call ndims_W(U);
!      call nel_W(U);
      call nlorb_W(U);
      call non30_W(U);
      call NOR_W(U);
      call NTRM_W(U);
      call ovrlap_W(U);
      call param_W(U);
      call RAS_W(U);
      call RED_W(U);
      call SIGNF_W(U);
      call state_W(U);
!      call TRK_W(U);
      close(U);
      END subroutine save_data

! Define the routines saving data;
      subroutine closed_W(U);
         USE closed_C;
         integer U;
         write(U) NCLOSD, IBK;
         write(U) B1ELC(1:4);
      END subroutine closed_W

      subroutine CONSTS_W(U);
         ! the constants are in data statements;
         integer U;
      END subroutine consts_W;

! save integral data
      subroutine CSF_W(U);
         USE CSF_C
         integer U;
         write(U) size_c_1;
         write(U) size_int_1;
         write(U) intptr1(1:size_int_1),intgrl1(1:size_int_1);
         write(U) cnn1(1:size_c_1),jann1(1:size_c_1),jbnn1(1:size_c_1);
         write(U) size_c_2;
         write(U) size_int_2;
         write(U) intptr2(1:size_int_2),intgrl2(1:size_int_2);
         write(U) cnn2(1:size_c_2),jann2(1:size_c_2),jbnn2(1:size_c_2);
      END subroutine
!
      subroutine DBG_W(U);
         USE DBG_C;
         integer U;
         write(U) IBUGM;
      END subroutine DBG_W;

      subroutine DEBUG_W(U);
         USE DEBUG_C;
         integer U;
         write(U) IBUG1, IBUG2, IBUG3, NBUG6, NBUG7, IFUL
      END subroutine

      subroutine diagnl_W(U);
      USE diagnl_C;
         integer U;
         write(U) IDIAG, JA, JB
      END subroutine

      subroutine dimen_W(U);
         USE dimen_C;
         integer U;
         write(U) kfl1, kfl2, kfl3, kfl4, kfl5, kfl6, kfl7, mxihsh
      END subroutine

      subroutine ELT_W(U);
         use parameters_biotr_C
         USE ELT_C;
         integer U,I;
         write(U) ELTENS(1:NWD);
      END subroutine

      subroutine EMS_W(U);
         USE EMS_C;
         integer U;
         write(U) IEM(1:4);
         write(U) IFL, JI, JF, LAM
         write(U) REL, VOK
      END subroutine

      subroutine fact_W(U);
         USE fact_C;
         integer U;
         write(U) GAM(1:100);
      END subroutine

      subroutine FOUT_W(U);
         USE fout_C;
         integer U;
         write(U) LCOUNT, NREC, IFLAG, LIJ, NIJ
      END subroutine

      subroutine fo_W(U);
         USE fo_C;
         integer U;
         write(U) IWF(1:2),NCLOS(1:2),NOPEN(1:2);
         write(U) MAXNFO
      END subroutine

      subroutine INFORM_W(U);
      USE INFORM_C;
         integer U;
         write(U) JSC(1:3);
         write(U) IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,JSC0,ISCW; 
      END subroutine

      subroutine INOUT_W(U);
         USE INOUT_C;
         integer U;
         write(U) IUC(1:2), IUW(1:2), IUL(1:2), IUJ(1:2), IUT(1:2)
         write(U) IREAD, IWRITE, ISCW
      END subroutine

      subroutine MEDEFN_W(U);
         USE MEDEFN_C;
         integer U;
         write(U) NJ(1:16),LJ(1:16),NOSH1(1:16),NOSH2(1:16),IJFUL(1:16)
         write(U) J1QN1(1:31,1:3),J1QN2(1:31,1:3)
         write(U) IHSH
      END subroutine

      subroutine ndims_W(U);
! used only in biotr_ang
         USE NDIMS_C;
         integer U;
         write(U) NCFG;
         write(U) NOCCSH(1:NCFG);
         write(U) NELCSH(1:8,1:NCFG),NOCORB(1:8,1:NCFG);
         write(U) J1QNRD(1:15,1:NCFG);
      END subroutine

      subroutine nel_W(U);
         !USE NEL_C, not comuted in biotr_ang;
         integer U;
      END subroutine

      subroutine nlorb_W(U);
         USE nlorb_C;
         integer U;
         write(U) nq1(1:NWD), lq1(1:NWD), nq2(1:NWD), lq2(1:NWD);
      END subroutine

      subroutine non30_W(U);
         USE non30_C;
         integer U;
         write(U) MAXORB;
         write(U) IAJCMP(1:MAXORB),LJCOMP(1:MAXORB),NJCOMP(1:MAXORB);
         write(U) IAJCLD(1:NWCD),LJCLSD(1:NWCD);
      END subroutine

      subroutine NOR_W(U);
         USE NOR_C;
         integer U;
         write(U) NCOM, NCLOSI, NCLOSF, NORBI, NORBF, IWAR;
      END subroutine

      subroutine NTRM_W(U);
         USE NTRM_C;
         integer U;
         write(U) NTERMS
      END subroutine

      subroutine ovrlap_W(U);
         USE OVRLAP_C;
         integer U;
         write(U) MU, NU, MUP, NUP, NONORT, NOVLPS, IROWMU, &
                  IROWNU, ICOLMU, ICOLNU, NORTH, IORDER, NCALLS, &
                  LMU, LNU, LMUP, LNUP, JMU, JNU, JMUP, JNUP
      END subroutine

      subroutine param_W(U);
         USE PARAM_C;
         integer U;
         write(U) NO, ND, NWF, MASS, NCFG, IB, IC, ID, NSCF, NCLOSD
         write(U) H, H1, H3, CH, EH, RHO, Z, TOL, D0, D1, D2, D3, D4, D5, &
                 D6, D8, D10, D12, D16, D30, FINE, RMASS
      END subroutine

      subroutine RAS_W(U);
         use parameters_biotr_C
         USE RAS_C;
         integer U;
         write(U) NINAC(1:11,1:2)
         write(U) NRAS1(1:11,1:10,1:2)
         write(U) NRAS2(1:11,1:2);
         write(U) NRAS3(1:11,1:10,1:2)
         write(U) NL(1:11,1:2);
         write(U) NGAS1(1:2);
         write(U) NGAS3(1:2);
         write(U) ITAB(1:NWD,1:2);
         write(U) ELRAS(1:2*NWD);
         write(U) ELRASI(1:nwd);
         write(U) ELRASF(1:nwd);
         write(U) LMAX(1:2);
      END subroutine

      subroutine RED_W(U);
         USE RED_C;
         integer U;
         write(U) INDL(1:20), INDR(1:20), IORST(1:20);
         write(U) LIS(1:16);
         write(U) NCI1, KPL, KPR, NC, JIST, JFST, NCONTR
      END subroutine

      subroutine SIGNF_W(U);
         USE SIGNF_C;
         integer U;
         write(U) SIGNFA;
      END subroutine

      subroutine state_W(U);
         USE STATE_C;
         ! only ncf(I) computed in biotr_ang 
         integer U;
         write(U) NCF(1:2)
      END subroutine

      subroutine TRK_W(U);
         USE TRK_C;
         integer U;
         write(U) ID1(1:7), ID2(1:7), IK1(1:7), IK2(1:7);
         write(U) BD1(1:3), BD2(1:3), BK1(1:3), BK2(1:3);
      END subroutine


