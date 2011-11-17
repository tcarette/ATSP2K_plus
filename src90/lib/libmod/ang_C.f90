      MODULE caseop_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  12:15:26  11/16/01  
      INTEGER :: IOCASE 
      END MODULE caseop_C 

      MODULE consts_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  07:33:24  11/16/01  
      REAL(DOUBLE) :: ZERO, TENTH, HALF, ONE, TWO, THREE, FOUR, &
                      SEVEN, ELEVEN, EPS 
!
!
!     SET GLOBAL REAL CONSTANTS
!
      DATA ZERO, TENTH, HALF, ONE, TWO, THREE, FOUR, SEVEN, ELEVEN, EPS/ 0.0D00&
         , 0.1D00, 0.5D00, 1.0D00, 2.0D00, 3.0D00, 4.0D00, 7.0D00, 1.1D01, &
         1.0D-08/  

      END MODULE consts_C 

      MODULE debug_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  07:17:04  11/16/01  
      INTEGER :: IBUG1, IBUG2, IBUG3, NBUG6, NBUG7, IFULL 
      END MODULE debug_C 

      MODULE diagnl_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  07:17:04  11/16/01  
      INTEGER :: IDIAG, JA, JB 
      END MODULE diagnl_C 

      MODULE ems_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  11:44:13  11/16/01  
      INTEGER, DIMENSION(4) :: IEM 
      INTEGER :: IFL, JI, JF, LAM 
      LOGICAL :: REL, VOK 
      DATA IEM/2HE ,2HM ,2HMA, 2HMB/
      END MODULE ems_C 

      MODULE fact_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  10:07:18  11/16/01  
      REAL(DOUBLE), DIMENSION(100) :: GAM 
      END MODULE fact_C 

      MODULE inform_C 
!...Created by Pacific-Sierra Research 77to90  4.3E  07:17:04  11/16/01  
      INTEGER, DIMENSION(4) :: JSC,isc
      INTEGER :: IREAD, IWRITE, IOUT, ISC0, ISC1, ISC2, ISC3, JSC0 
      integer :: ISCW,iall
      integer :: IREAD2,IWRITE2,IREADF,ISCW2
      equivalence(ISCW,JSC(4))
      equivalence(ISC0,ISC(1))
      equivalence(ISC1,ISC(2))
      equivalence(ISC2,ISC(3))
      equivalence(ISC3,ISC(4))
      equivalence(IALL,JSC0)


      END MODULE inform_C 

      MODULE kampas_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  07:33:24  11/16/01  
      INTEGER, DIMENSION(2,20) :: IW1, IW2 
      INTEGER, DIMENSION(2,2,20,20) :: IWAA 
      REAL(DOUBLE), DIMENSION(2,20) :: RW1, RW2 
      REAL(DOUBLE), DIMENSION(2,2,20,20) :: RWAA 
      END MODULE kampas_C 
      MODULE kron_C 
!...Created by Pacific-Sierra Research 77to90  4.3E  09:30:21  11/16/01  
      INTEGER, DIMENSION(16,16) :: IDEL 
!
!
! --- READ IN OTHER INITIALIZATION DATA
!
      DATA IDEL/ 1, 16*0, 1, 16*0, 1, 16*0, 1, 16*0, 1, 16*0, 1, 16*0, 1, 16*0&
         , 1, 16*0, 1, 16*0, 1, 16*0, 1, 16*0, 1, 16*0, 1, 16*0, 1, 16*0, 1, 16&
         *0, 1/


      END MODULE kron_C 

      MODULE medefn_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  07:17:04  11/16/01  
      INTEGER, DIMENSION(16) :: NJ, LJ, NOSH1, NOSH2 
      INTEGER, DIMENSION(31,3) :: J1QN1, J1QN2 
      INTEGER, DIMENSION(16) :: IJFUL 
      INTEGER :: IHSH 
      integer, dimension(16,2) :: NOSH
      integer, dimension(31,3,2) :: J1QN
      equivalence(J1QN(1,1,1),J1QN1(1,1))
      equivalence(J1QN(1,1,2),J1QN2(1,1))
      equivalence(NOSH(1,1),NOSH1(1)),(NOSH(1,2),NOSH2(1))
      END MODULE medefn_C 

      MODULE medefn_a_C 
!...Created by Pacific-Sierra Research 77to90  4.3E  12:43:16  11/16/01  
      INTEGER, DIMENSION(16) :: NJ, LJ 
      INTEGER, DIMENSION(16,2) :: NOSH 
      INTEGER, DIMENSION(31,3,2) :: J1QN 
      INTEGER, DIMENSION(16) :: IJFUL 
      INTEGER :: IHSH 
      END MODULE medefn_a_C 

      MODULE mt15_C 
!...Created by Pacific-Sierra Research 77to90  4.3E  07:48:27  11/16/01  
      INTEGER, DIMENSION(1) :: M1 
      INTEGER, DIMENSION(7) :: M2 
      INTEGER, DIMENSION(17) :: M3 
      INTEGER, DIMENSION(47) :: M4 
      INTEGER, DIMENSION(73) :: M5 

      DATA M1/ 1060106/  
      DATA M2/ 1050206, 1050202, 1050210, 1050004, 1050008, 1070000, 1050012/  
      DATA M3/ 1040300, 1040304, 1040306, 1040308, 1040312, 1040102, 1040104, &
         2040104, 1060106, 2040106, 1040108, 2040108, 1040110, 2040110, 1040112&
         , 1040114, 1040116/  
      DATA M4/ 1030404, 1030406, 1030408, 30400, 1030412, 1050206, 2030206, &
         1030204, 2030204, 3030206, 1030208, 2030208, 4030206, 3030208, 1050202&
         , 2030202, 3030202, 1050210, 2030210, 3030210, 4030210, 1030212, &
         2030212, 1030214, 2030214, 1030216, 1030218, 1050004, 2030004, 3030004&
         , 1030006, 1050008, 2030008, 3030008, 4030004, 4030008, 1030010, &
         2030010, 1070000, 1050012, 2030000, 2030012, 3030012, 1030014, 1030016&
         , 2030016, 1030020/  
      DATA M5/ 20502, 20506, 20510, 1040300, 1020302, 2020302, 1040304, 2020304&
         , 3020304, 1040306, 2020306, 3020306, 4020306, 1040308, 2020308, &
         3020308, 4020308, 1020310, 2020310, 3020310, 1040312, 2020312, 3020312&
         , 1020314, 2020314, 1020316, 20318, 1040102, 2020102, 3020102, 4020102&
         , 1040104, 2040104, 3020104, 4020104, 5020104, 1060106, 2040106, &
         3020106, 4020106, 5020106, 6020106, 7020106, 1040108, 2040108, 3020108&
         , 4020108, 5020108, 6020108, 1040110, 2040110, 3020110, 4020110, &
         5020110, 6020110, 7020110, 1040112, 2020112, 3020112, 4020112, 5020112&
         , 1040114, 2020114, 3020114, 4020114, 5020114, 1040116, 2020116, &
         3020116, 1020118, 2020118, 1020120, 20122/  

      END MODULE mt15_C 

      MODULE mt67_C 
!...Created by Pacific-Sierra Research 77to90  4.3E  09:19:51  11/16/01  
      INTEGER, DIMENSION(119) :: M7, M6 
      INTEGER, DIMENSION(238) :: M76
      equivalence(M76(1),M7)
      equivalence(M76(120),M6)
      DATA M7/ 700, 20502, 504, 20506, 508, 20510, 512, 1040300, 2000300, &
         1020302, 2020302, 1040304, 2020304, 3020304, 4000304, 5000304, 6000304&
         , 1040306, 2020306, 3020306, 4020306, 5000306, 1040308, 2020308, &
         3020308, 4020308, 5000308, 6000308, 7000308, 1020310, 2020310, 3020310&
         , 4000310, 5000310, 1040312, 2020312, 3020312, 4000312, 5000312, &
         1020314, 2020314, 3000314, 1020316, 2000316, 3000316, 20318, 320, &
         1000100, 2000100, 1040102, 2020102, 3020102, 4020102, 5000102, 1040104&
         , 2040104, 3020104, 4020104, 5020104, 6000104, 7000104, 1060106, &
         2040106, 3020106, 4020106, 5020106, 6020106, 7020106, 8000106, 9000106&
         , 10000106, 1040108, 2040108, 3020108, 4020108, 5020108, 6020108, &
         7000108, 8000108, 9000108, 10000108, 1040110, 2040110, 3020110, &
         4020110, 5020110, 6020110, 7020110, 8000110, 9000110, 1040112, 2020112&
         , 3020112, 4020112, 5020112, 6000112, 7000112, 8000112, 9000112, &
         1040114, 2020114, 3020114, 4020114, 5020114, 6000114, 7000114, 1040116&
         , 2020116, 3020116, 4000116, 5000116, 1020118, 2020118, 3000118, &
         4000118, 1020120, 2000120, 20122, 124/
      DATA M6/ 10606, 1030404, 2010404, 3010404, 1030406, 2010406, 1030408, &
         2010408, 3010408, 10402, 1010410, 2010410, 30400, 1030412, 2010412, &
         10414, 10416, 1050206, 2030206, 6010206, 8010206, 1030204, 2030204, &
         3010204, 4010204, 3030206, 5010206, 1030208, 2030208, 4010208, 5010208&
         , 5010204, 4030206, 7010206, 9010206, 3030208, 6010208, 7010208, &
         1050202, 2030202, 3030202, 1050210, 2030210, 3030210, 4030210, 4010202&
         , 5010210, 6010210, 5010202, 6010202, 7010210, 8010210, 9010210, &
         1030212, 2030212, 3010212, 4010212, 5010212, 6010212, 1030214, 2030214&
         , 3010214, 4010214, 5010214, 6010214, 1030216, 2010216, 3010216, &
         1030218, 2010218, 3010218, 10220, 10222, 2010006, 3010006, 4010006, &
         1050004, 2030004, 3030004, 1030006, 1050008, 2030008, 3030008, 5010004&
         , 5010008, 6010004, 6010008, 7010008, 8010008, 4030004, 4030008, &
         1030010, 2030010, 10002, 3010010, 4010010, 1070000, 1050012, 2030000, &
         2030012, 3030012, 3010000, 4010012, 5010012, 4010000, 6010012, 7010012&
         , 1030014, 2010014, 3010014, 1030016, 2030016, 3010016, 4010016, &
         1010018, 2010018, 1030020, 2010020, 10024/

      END MODULE mt67_C 
      MODULE mt_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  10:32:45  11/16/01  
      INTEGER, DIMENSION(40) :: MT 

      DATA MT/ 100, 10000, 300, 20102, 104, 30000, 10202, 10004, 500, 100, &
         20302, 20102, 40104, 20104, 304, 104, 20306, 20106, 106, 20108, 308, &
         108, 20110, 112, 50000, 10000, 30202, 10202, 10404, 10204, 30004, &
         10004, 30206, 10206, 10006, 10208, 30008, 10008, 10210, 10012/

      END MODULE mt_C 

      MODULE ndims_C
        integer :: ncfg
        integer, dimension(:), pointer :: QNOC
        integer, dimension(:,:), pointer :: QNELCSH,QNOCORB,QJ1

        integer, dimension(:), allocatable, target :: NOCCSH
        integer, dimension(:,:), allocatable, target :: NELCSH
        integer, dimension(:,:), allocatable, target :: NOCORB,J1QNRD

      END MODULE ndims_C


      MODULE non30_C 
        integer :: maxorb
        integer, parameter :: nwcd = 20
        integer, parameter :: ncwd = 20
        integer, dimension(:), pointer :: qiajcmp,qljcomp,QNJCOMP
        integer, dimension(:), pointer :: QIAJCLD,QLJCLSD

        integer, dimension(:), allocatable, target :: iajcmp,ljcomp,NJCOMP
        integer, dimension(:), allocatable, target :: IAJCLD,LJCLSD

      END MODULE non30_C 

      MODULE ntrm_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  11:44:13  11/16/01  
      INTEGER :: NTERMS 
      END MODULE ntrm_C 

      MODULE occupation_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  12:38:29  11/16/01  
      INTEGER, DIMENSION(16,2) :: NCG, ICG 
      END MODULE occupation_C 
      MODULE operat_C 
!...Created by Pacific-Sierra Research 77to90  4.3E  07:48:27  11/16/01  
      INTEGER :: ICOLOM, ISOTOP, IORBORB 
      END MODULE operat_C 
      MODULE permat_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  12:15:26  11/16/01  
      INTEGER, DIMENSION(2,2) :: IRS 
      INTEGER, DIMENSION(20,20) :: IRL 
      REAL(DOUBLE), DIMENSION(2,2) :: RS 
      REAL(DOUBLE), DIMENSION(20,20) :: RL 
      END MODULE permat_C 

      MODULE rdint_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  11:44:13  11/16/01  
      REAL(DOUBLE), dimension(:,:), pointer :: IQRL, IQRV, IQOV
      REAL(DOUBLE), dimension(:,:), allocatable, target :: RLINT,RVINT,OVRLP
      END MODULE rdint_C 

      MODULE ribof_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  10:10:49  11/16/01  
      INTEGER, DIMENSION(238) :: IMPTF, IMGTF, IMPNF, IMGNF 
      DATA IMPTF/ 119*1, 119*120/
      DATA IMGTF/ 119*119, 119*238/
      DATA IMPNF/ 119*120, 119*1/
      DATA IMGNF/ 119*238, 119*119/

      END MODULE ribof_C 
      MODULE ribols3_C 
!...Created by Pacific-Sierra Research 77to90  4.3E  10:11:28  11/16/01  
      INTEGER, DIMENSION(90) :: IMPTLS3, IMGTLS3, IMPNLS3, IMGNLS3 

      DATA IMPTLS3/ 1, 9*2, 11, 11*12, 23, 13*24, 37, 15*38, 53, 17*54, 71, 19*&
         72/
      DATA IMGTLS3/ 1, 9*10, 11, 11*22, 23, 13*36, 37, 15*52, 53, 17*70, 71, 19&
         *90/
      DATA IMPNLS3/ 2, 9*1, 12, 11*11, 24, 13*23, 38, 15*37, 54, 17*53, 72, 19*&
         71/
      DATA IMGNLS3/ 10, 9*1, 22, 11*11, 36, 13*23, 52, 15*37, 70, 17*53, 90, 19&
         *71/

      END MODULE ribols3_C 
      MODULE ribols_C 
!...Created by Pacific-Sierra Research 77to90  4.3E  10:11:28  11/16/01  
      INTEGER, DIMENSION(40) :: IMPTLS, IMGTLS, IMPNLS, IMGNLS 

      DATA IMPTLS/ 1, 2, 3*3, 3*6, 16*9, 16*25/
      DATA IMGTLS/ 1, 2, 3*5, 3*8, 16*24, 16*40/
      DATA IMPNLS/ 2, 1, 3*6, 3*3, 16*25, 16*9/
      DATA IMGNLS/ 2, 1, 3*8, 3*5, 16*40, 16*24/

      END MODULE ribols_C 
      MODULE ribolsf_C 
!...Created by Pacific-Sierra Research 77to90  4.3E  10:11:28  11/16/01  
      INTEGER, DIMENSION(8) :: IMPTLSF, IMGTLSF, IMPNLSF, IMGNLSF 

      DATA IMPTLSF/ 301, 7*302/
      DATA IMGTLSF/ 301, 7*308/
      DATA IMPNLSF/ 302, 7*301/
      DATA IMGNLSF/ 308, 7*301/
      END MODULE ribolsf_C 

      MODULE savecom_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  11:44:13  11/16/01  
      INTEGER :: LRHO, LSIG 
      REAL(DOUBLE) :: CL2, CV2 
      END MODULE savecom_C 

      MODULE signf_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  11:44:13  11/16/01  
      REAL(DOUBLE) :: SIGNFA 
      END MODULE signf_C 

      MODULE skmt2_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  10:32:45  11/16/01  
      INTEGER, DIMENSION(8) :: MTF 
      INTEGER, DIMENSION(90) :: MT3 

      DATA MTF/ 60106, 70000, 50202, 50004, 50206, 50008, 50210, 50012/
      DATA MT3/ 80108, 90000, 70202, 70004, 70206, 70008, 70210, 70012, 70214, &
         70016, 100110, 110000, 90202, 90004, 90206, 90008, 90210, 90012, 90214&
         , 90016, 90218, 90020, 120112, 130000, 110202, 110004, 110206, 110008&
         , 110210, 110012, 110214, 110016, 110218, 110020, 110222, 110024, &
         140114, 150000, 130202, 130004, 130206, 130008, 130210, 130012, 130214&
         , 130016, 130218, 130020, 130222, 130024, 130226, 130028, 160116, &
         170000, 150202, 150004, 150206, 150008, 150210, 150012, 150214, 150016&
         , 150218, 150020, 150222, 150024, 150226, 150028, 150230, 150032, &
         180118, 190000, 170202, 170004, 170206, 170008, 170210, 170012, 170214&
         , 170016, 170218, 170020, 170222, 170024, 170226, 170028, 170230, &
         170032, 170234, 170036/

      END MODULE skmt2_C 
      MODULE terms_C 
!...Created by Pacific-Sierra Research 77to90  4.3E  07:48:27  11/16/01  
      INTEGER, DIMENSION(24) :: ITAB, JTAB 
      INTEGER, DIMENSION(333) :: NTAB 
      INTEGER :: NROWS 
!
! --- SETS QUANTUM NUMBERS OF TERMS WHICH CAN BE FORMED FROM
!     CONFIGURATIONS  L**Q . ONLY THE FIRST HALF OF THAT PART OF THE
!     TABLE, CORRESPONDING TO A GIVEN  L, IS INCLUDED, BECAUSE OF THE
!     SYMMETRY OF THE TABLE.  E.G. D**7 FORMS THE SAME TERMS AS D**3
!
!       The tables are set for a maximum value of L=9; the terms
!     for L>3 are assumed to be the same as those for L=3
!
!     S - SHELLS (ROWS 1 AND 2)
!
!     P - SHELLS (ROWS 3 TO 5)
!
!     D - SHELLS (ROWS 6 TO 10)
!
!     F - SHELLS (ROWS 11 AND 12)
!
!     G - SHELLS (ROWS 13 AND 14)
!
!     H - SHELLS (ROWS 15 AND 16)
!
!     I - SHELLS (ROWS 17 AND 18)
!
!     K - SHELLS (ROWS 19 AND 20)
!
!     L - SHELLS (ROWS 21 AND 22)
!
!     M - SHELLS (ROWS 23 AND 24)
!
      DATA NROWS/ 24/
!
!     THE ARRAYS I,J,N CORRESPOND TO THE ARRAYS ITAB,JTAB,NTAB
!
      DATA Itab/ 1, 1, 1, 3, 3, 1, 5, 8, 16, 16, 1, 7, 1, 7, 1, 7, 1, 7, 1, 7, 1, &
         7, 1, 7/
      DATA Jtab/ 0, 3, 6, 9, 18, 27, 30, 45, 69, 117, 165, 168, 189, 192, 213, 216&
         , 237, 240, 261, 264, 285, 288, 309, 312/
      DATA Ntab/ 1, 1, 2, 0, 1, 1, 1, 3, 2, 0, 1, 1, 2, 5, 1, 2, 3, 3, 1, 3, 2, 3&
         , 5, 2, 3, 1, 4, 1, 5, 2, 0, 1, 1, 2, 5, 1, 2, 9, 1, 2, 3, 3, 2, 7, 3&
         , 1, 5, 2, 3, 3, 2, 3, 5, 2, 3, 7, 2, 3, 9, 2, 3, 11, 2, 3, 3, 4, 3, 7&
         , 4, 0, 1, 1, 2, 5, 1, 2, 9, 1, 2, 3, 3, 2, 7, 3, 4, 1, 1, 4, 5, 1, 4&
         , 7, 1, 4, 9, 1, 4, 13, 1, 4, 3, 3, 4, 5, 3, 4, 7, 3, 4, 9, 3, 4, 11, &
         3, 4, 5, 5, 1, 5, 2, 3, 3, 2, 3, 5, 2, 3, 7, 2, 3, 9, 2, 3, 11, 2, 3, &
         3, 4, 3, 7, 4, 5, 1, 2, 5, 5, 2, 5, 7, 2, 5, 9, 2, 5, 13, 2, 5, 5, 4, &
         5, 9, 4, 5, 1, 6, 1, 7, 2, 2, 3, 3, 2, 7, 3, 2, 11, 3, 0, 1, 1, 2, 5, &
         1, 2, 9, 1, 2, 13, 1, 1, 9, 2, 2, 3, 3, 2, 7, 3, 2, 11, 3, 0, 1, 1, 2&
         , 5, 1, 2, 9, 1, 2, 13, 1, 1, 11, 2, 2, 3, 3, 2, 7, 3, 2, 11, 3, 0, 1&
         , 1, 2, 5, 1, 2, 9, 1, 2, 13, 1, 1, 13, 2, 2, 3, 3, 2, 7, 3, 2, 11, 3&
         , 0, 1, 1, 2, 5, 1, 2, 9, 1, 2, 13, 1, 1, 15, 2, 2, 3, 3, 2, 7, 3, 2, &
         11, 3, 0, 1, 1, 2, 5, 1, 2, 9, 1, 2, 13, 1, 1, 17, 2, 2, 3, 3, 2, 7, 3&
         , 2, 11, 3, 0, 1, 1, 2, 5, 1, 2, 9, 1, 2, 13, 1, 1, 19, 2, 2, 3, 3, 2&
         , 7, 3, 2, 11, 3, 0, 1, 1, 2, 5, 1, 2, 9, 1, 2, 13, 1/

      END MODULE terms_C 

      MODULE trk2_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  07:33:24  11/16/01  
      INTEGER, DIMENSION(7) :: ID3, ID4, IK3, IK4 
      REAL(DOUBLE), DIMENSION(3) :: BD3, BD4, BK3, BK4 
      END MODULE trk2_C 

      MODULE trk_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  07:33:24  11/16/01  
      INTEGER, DIMENSION(7) :: ID1, ID2, IK1, IK2 
      REAL(DOUBLE), DIMENSION(3) :: BD1, BD2, BK1, BK2 
      END MODULE trk_C 
