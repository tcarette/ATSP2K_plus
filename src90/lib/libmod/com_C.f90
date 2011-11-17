      MODULE blume_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  21:49:59  11/14/01  
      REAL(DOUBLE), DIMENSION(4) :: COEFN2, COEFNK, COEFVK 
      END MODULE blume_C 

      MODULE eav_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  21:49:59  11/14/01  
      REAL(DOUBLE), DIMENSION(10) :: CCA 
      REAL(DOUBLE), DIMENSION(35) :: CCB 
      END MODULE eav_C 

      MODULE enav_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  21:49:59  11/14/01  
      INTEGER, DIMENSION(15) :: KVALUE 
      INTEGER :: NINTS 
      REAL(DOUBLE), DIMENSION(15) :: COEFCT 
      END MODULE enav_C 

      MODULE fact_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  21:49:59  11/14/01  
      REAL(DOUBLE), DIMENSION(100) :: GAM 
      END MODULE fact_C 
