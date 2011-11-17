      MODULE blume_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      REAL(DOUBLE), DIMENSION(4) :: COEFN2, COEFNK, COEFVK 
      END MODULE blume_C 

      MODULE nel_C 
      USE vast_kind_param, ONLY:  DOUBLE 
      use parameters_biotr_C
!...Created by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
        integer, dimension(:), pointer :: IQN,IQL,IQMAX_
        integer, dimension(:), allocatable, target :: N,L,MAX_
        REAL(DOUBLE), dimension(:), pointer :: IQAZ
        REAL(DOUBLE), dimension(:), allocatable, target :: AZ
        REAL(DOUBLE), dimension(:,:), pointer :: IQP
        REAL(DOUBLE), dimension(:,:),allocatable,target :: P
        integer, dimension(7) :: IQ
        character(LEN=3), dimension(NWD) :: EL 
        integer szX, szY
      END MODULE nel_C 

      MODULE param_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      INTEGER :: NO, ND, NWF, MASS, NCFG, IB, IC, ID, NSCF, NCLOSD 
      INTEGER :: MSOO
      REAL(DOUBLE) :: H, H1, H3, CH, EH, RHO, Z, TOL, D0, D1, D2, D3, D4, D5, &
         D6, D8, D10, D12, D16, D30, FINE, RMASS 
      EQUIVALENCE (MASS,MSOO)
      END MODULE param_C 

      MODULE radial_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      INTEGER, PARAMETER :: NOD = 220 
      REAL(DOUBLE), DIMENSION(NOD) :: R, RR, R2, YK, YR, X 
      END MODULE radial_C 
