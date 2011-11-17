      MODULE buffer_C
!        use parameters_biotr_C
        integer :: ncfg
        integer, dimension(:), pointer :: qpackn, qjan, qjbn
        double precision, dimension(:), pointer :: qcn
        integer, dimension(:), allocatable, target ::ipackn,jbn,jan
        double precision, dimension(:), allocatable, target :: cn

      END MODULE buffer_C

      MODULE closed_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  00:07:05  11/17/01  
      INTEGER :: NCLOSD, IBK 
      REAL(DOUBLE), DIMENSION(4) :: B1ELC 
      END MODULE closed_C 

      MODULE csf_C
        use parameters_biotr_C
        integer :: nint1, nint2,ncoef1,ncoef2
        integer :: size_int_1, size_int_2,size_c_1, size_c_2
        integer, dimension(:), pointer :: qintgrl1,qintptr1
        integer, dimension(:), pointer :: qintgrl2,qintptr2
        integer, dimension(:), pointer :: qjann1,qjbnn1,qjann2,qjbnn2
        double precision, dimension(:), pointer :: qcnn1,qcnn2

        integer, dimension(:), allocatable, target :: intgrl1,intptr1
        integer, dimension(:), allocatable, target :: intgrl2,intptr2
        integer, dimension(:), allocatable, target :: jann1,jbnn1,jann2,jbnn2
        double precision, dimension(:), allocatable, target :: cnn1,cnn2
      END MODULE csf_C

      MODULE dbg_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  19:41:09  11/16/01  
      INTEGER :: IBUGM 
      END MODULE dbg_C 

      module dimen_C 
      USE vast_kind_param, ONLY:  double 
!...Created by Pacific-Sierra Research 77to90  4.3E  23:47:02  11/16/01  
      integer :: kfl1, kfl2, kfl3, kfl4, kfl5, kfl6, kfl7, mxihsh 
      end module dimen_C 

      MODULE elt_C 
      use parameters_biotr_C
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  19:41:09  11/16/01  
      CHARACTER, DIMENSION(NWD) :: ELTENS*3 
      END MODULE elt_C 

      MODULE fo_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  19:41:09  11/16/01  
      INTEGER, DIMENSION(2) :: IWF, NCLOS, NOPEN 
      INTEGER :: MAXNFO 
      END MODULE fo_C 

      MODULE fout_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  20:41:07  11/16/01  
      INTEGER :: LCOUNT, NREC, IFLAG, LIJ, NIJ 
      END MODULE fout_C 

      MODULE inout_C 
      use inform_C
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  19:41:09  11/16/01  
      INTEGER, DIMENSION(2) :: IUC, IUW, IUL, IUJ, IUT 
!      INTEGER :: IREAD, IWRITE, ISCW 
      END MODULE inout_C 

      module nlorb_C 
      USE vast_kind_param, ONLY:  double 
!...Created by Pacific-Sierra Research 77to90  4.3E  23:47:02  11/16/01  
      integer, parameter :: nwd = 128 
      integer, dimension(nwd) :: nq1, lq1, nq2, lq2 
      end module nlorb_C 

      MODULE nor_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  19:41:09  11/16/01  
      INTEGER :: NCOM, NCLOSI, NCLOSF, NORBI, NORBF, IWAR 
      END MODULE nor_C 

      MODULE ovrint_C 
      USE vast_kind_param, ONLY:  DOUBLE 
      integer, dimension(2) :: iovel,iover, iovep 
      END MODULE ovrint_C 

      MODULE ovrlap_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  20:41:07  11/16/01  
!      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
!     : ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
!     :     QIORTH

      integer, dimension(:), pointer :: QIORTH
      integer, dimension(:),allocatable, target :: iorth
      INTEGER :: MU, NU, MUP, NUP, NONORT, NOVLPS, IROWMU, IROWNU, ICOLMU, &
         ICOLNU, NORTH, IORDER, NCALLS, LMU, LNU, LMUP, LNUP, JMU, JNU, JMUP, &
         JNUP 
      END MODULE ovrlap_C 

      MODULE ras_C 
      use parameters_biotr_C
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  21:03:58  11/16/01  
      INTEGER, DIMENSION(11,2) :: NINAC 
      INTEGER, DIMENSION(11,10,2) :: NRAS1 
      INTEGER, DIMENSION(11,2) :: NRAS2 
      INTEGER, DIMENSION(11,10,2) :: NRAS3 
      INTEGER, DIMENSION(2) :: NGAS1, NGAS3 
      INTEGER, DIMENSION(11,2) :: NL 
      INTEGER, DIMENSION(NWD,2) :: ITAB 
      INTEGER, DIMENSION(2) :: LMAX 
      CHARACTER, DIMENSION(2*NWD) :: ELRAS*3 
      CHARACTER, DIMENSION(NWD) :: ELRASI*3, ELRASF*3 
      END MODULE ras_C 

      MODULE red_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  20:41:07  11/16/01  
      INTEGER, DIMENSION(20) :: INDL, INDR, IORST 
      INTEGER, DIMENSION(16) :: LIS 
      INTEGER :: NCI1, KPL, KPR, NC, JIST, JFST, NCONTR 
      END MODULE red_C 

      MODULE state_C 
      USE vast_kind_param, ONLY:  DOUBLE 
         !integer, dimension(:), pointer :: QIORTH
         !integer, dimension(:),allocatable, target :: iorth
         double precision, dimension(:), pointer :: QET1,QET2,QWT1,QWT2
         double precision, dimension(:,:), pointer :: QCFG1,QCFG2
         double precision, dimension(:),allocatable, target :: ET1,ET2,WT1,WT2
         double precision, dimension(:,:),allocatable, target :: CFG1,CFG2
         integer, dimension(:), pointer :: QLBL1,QJV1,QLBL2,QJV2
         integer, dimension(:),allocatable, target ::LBL1,JV1,LBL2,JV2
         INTEGER :: NJV(2),NVC(2),LGTH(2),NCF(2) 
      END MODULE state_C 
