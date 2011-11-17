      MODULE hyperstuff_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  13:56:18  11/16/01  
      INTEGER :: IAM, NODES 
      END MODULE hyperstuff_C 
      MODULE matrix_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  13:56:18  11/16/01  
      INTEGER, PARAMETER :: NZERMAX = 600000 
      INTEGER, PARAMETER :: NMAX = 9000 
      INTEGER, DIMENSION(NMAX) :: IND_COL 
      INTEGER, DIMENSION(NZERMAX) :: IND_ROW 
      REAL(DOUBLE), DIMENSION(NZERMAX) :: DMY_MATR 
      LOGICAL :: IUPPER 
      END MODULE matrix_C 
      MODULE temp_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  13:56:18  11/16/01  
      INTEGER, PARAMETER :: NMAX = 9000 
      REAL(DOUBLE), DIMENSION(NMAX) :: TM, TP 
      END MODULE temp_C 
