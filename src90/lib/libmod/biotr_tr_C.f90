      MODULE rydberg_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  15:58:56  11/18/01  
      REAL(DOUBLE) :: RYDBRG, ZMU, ENERGY 
      END MODULE rydberg_C 

      MODULE cor_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  15:58:56  11/18/01  
      REAL(DOUBLE) :: COREF, COREG, CORER, CORES, ZM1, ZM2, FZE, FZA, DR2 
      LOGICAL :: CORE, YM, YF 
      END MODULE cor_C 
      MODULE mult_C 

      USE vast_kind_param, ONLY:  DOUBLE 
        integer, dimension(:), pointer :: QIL,QIR,QJVL,QJVR 
        integer, dimension(:), allocatable, target :: IL,IR,JVL,JVR
        REAL(DOUBLE), dimension(:), pointer :: QSL,QSV 
        REAL(DOUBLE), dimension(:), allocatable, target :: SL,SV
        integer szX
      END MODULE mult_C 

      MODULE lsj_C 
      USE vast_kind_param, ONLY:  DOUBLE1=>DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  13:29:54  11/20/01  
      INTEGER :: LL1, LL2, IS1, IS2 
      END MODULE lsj_C 
      MODULE gtracloc_C 
      USE vast_kind_param, ONLY:  DOUBLE 
!...Created by Pacific-Sierra Research 77to90  4.3E  15:23:09  11/20/01  
      INTEGER :: LISTCUR, IELMCUR 
      END MODULE gtracloc_C 
