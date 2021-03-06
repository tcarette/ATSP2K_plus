      MODULE dyk_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  23:19:56  11/15/01  
      SUBROUTINE dyk (I, J, K) 
      INTEGER NOD 
      PARAMETER(NOD=220) 
      INTEGER NWD 
      PARAMETER(NWD=128) 
      INTEGER, INTENT(IN) :: I 
      INTEGER, INTENT(IN) :: J 
      INTEGER, INTENT(IN) :: K 
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE dzk_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      SUBROUTINE dzk (I, J, K) 
      INTEGER NOD 
      PARAMETER (NOD = 220) 
      INTEGER :: I 
      INTEGER :: J 
      INTEGER, INTENT(IN) :: K 
!VAST.../PARAM/ H(IN), H3(IN), EH(IN), Z(IN), ND(IN), D1(IN)
!VAST.../PARAM/ D2(IN), D4(IN), D5(IN), D6(IN), D8(IN)
!VAST.../RADIAL/ R(IN), R2(IN), YK(INOUT)
!VAST...Calls: L, P, AZ
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE ecore_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      SUBROUTINE ecore (EL, EC, REL) 
      USE vast_kind_param,ONLY: DOUBLE 
      INTEGER NOD 
      PARAMETER (NOD = 220) 
      REAL(DOUBLE), INTENT(OUT) :: EC 
      LOGICAL :: REL 
!VAST.../PARAM/ MSOO(IN), D2(IN), NCLOSD(IN)
!VAST...Calls: L, RK, ZCB, SN, HL
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE grad_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      REAL(KIND(0.0D0)) FUNCTION grad (I, J) 
      INTEGER NOD 
      PARAMETER (NOD = 220) 
      INTEGER, INTENT(IN) :: I 
      INTEGER, INTENT(IN) :: J 
!VAST.../PARAM/ H1(IN), Z(IN), ND(IN), D1(IN), D2(IN), D3(IN)
!VAST.../PARAM/ D4(IN), D5(IN), D6(IN)
!VAST.../RADIAL/ R(IN)
!VAST...Calls: L, P
!...This routine performs I/O.
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE hl_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      REAL(KIND(0.0D0)) FUNCTION hl (EL, I, J, REL) 
      INTEGER NOD 
      PARAMETER (NOD = 220) 
      CHARACTER (LEN = 3), DIMENSION(*), INTENT(IN) :: EL 
      INTEGER, INTENT(IN) :: I 
      INTEGER, INTENT(IN) :: J 
      LOGICAL, INTENT(IN) :: REL 
!VAST.../PARAM/ H(IN), H1(IN), Z(IN), ND(IN), D1(IN), D2(IN)
!VAST.../PARAM/ D3(IN), D4(IN), D5(IN), D6(IN)
!VAST.../RADIAL/ R(IN)
!VAST...Calls: L, P, RLSHFT
!...This routine performs I/O.
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE hlc_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      REAL(KIND(0.0D0)) FUNCTION hlc (EL, I, J, REL) 
      INTEGER NOD 
      PARAMETER (NOD = 220) 
      CHARACTER (LEN = 3), DIMENSION(*) :: EL 
      INTEGER :: I 
      INTEGER :: J 
      LOGICAL :: REL 
!VAST.../PARAM/ MSOO(IN), NCLOSD(IN)
!VAST...Calls: HL, L, RK, ZCB, SN
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE quadr_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      REAL(KIND(0.0D0)) FUNCTION quadr (I, J, KK) 
      INTEGER NOD 
      PARAMETER (NOD = 220) 
      INTEGER, INTENT(IN) :: I 
      INTEGER, INTENT(IN) :: J 
      INTEGER, INTENT(IN) :: KK 
!VAST.../PARAM/ H1(IN), Z(IN), D0(IN), D1(IN), D2(IN), D5(IN)
!VAST.../RADIAL/ R(IN), R2(IN)
!VAST...Calls: L, P, AZ
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE quads_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      REAL(KIND(0.0D0)) FUNCTION quads (I, J, KK) 
      INTEGER NOD 
      PARAMETER (NOD = 220) 
      INTEGER, INTENT(IN) :: I 
      INTEGER, INTENT(IN) :: J 
      INTEGER, INTENT(IN) :: KK 
!VAST.../PARAM/ H1(IN), Z(IN), D0(IN), D1(IN), D2(IN), D5(IN)
!VAST.../RADIAL/ R(IN), YK(IN)
!VAST...Calls: L, P
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE rk_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      REAL(KIND(0.0D0)) FUNCTION rk (I1, I2, I3, I4, K, REL) 
      INTEGER, INTENT(IN) :: I1 
      INTEGER, INTENT(IN) :: I2 
      INTEGER, INTENT(IN) :: I3 
      INTEGER, INTENT(IN) :: I4 
      INTEGER, INTENT(IN) :: K 
      LOGICAL :: REL 
!VAST.../PARAM/ Z(IN), MSOO(IN), D1(IN), D2(IN), RMASS(IN)
!VAST...Calls: YKF, QUADS, GRAD, QUADR, ZZ
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE rlshft_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      REAL(KIND(0.0D0)) FUNCTION rlshft (I1, I2) 
      INTEGER NOD 
      PARAMETER (NOD = 220) 
      INTEGER, INTENT(IN) :: I1 
      INTEGER, INTENT(IN) :: I2 
!VAST.../PARAM/ H(IN), H1(IN), Z(IN), D0(IN), D1(IN), D2(IN)
!VAST.../PARAM/ D4(IN), D5(IN), D12(IN), D16(IN), D30(IN)
!VAST.../PARAM/ FINE(IN)
!VAST.../RADIAL/ R(IN), RR(IN), R2(IN), YK(INOUT)
!VAST...Calls: L, P, AZ
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE tk_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      REAL(KIND(0.0D0)) FUNCTION tk (I, II, J, JJ, K) 
      INTEGER NOD 
      PARAMETER (NOD = 220) 
      INTEGER, INTENT(IN) :: I 
      INTEGER :: II 
      INTEGER, INTENT(IN) :: J 
      INTEGER :: JJ 
      INTEGER, INTENT(IN) :: K 
!VAST.../PARAM/ H(IN), H1(IN), D0(IN), D2(IN), D5(IN), D6(IN)
!VAST.../PARAM/ D8(IN), FINE(IN)
!VAST.../RADIAL/ YK(INOUT)
!VAST...Calls: DZK, YKK, L, P
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE ykf_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      SUBROUTINE ykf (I, J, K, REL) 
      INTEGER NOD 
      PARAMETER (NOD = 220) 
      INTEGER, INTENT(IN) :: I 
      INTEGER, INTENT(IN) :: J 
      INTEGER, INTENT(IN) :: K 
      LOGICAL, INTENT(IN) :: REL 
!VAST.../PARAM/ H3(IN), EH(IN), ND(IN), D4(IN), D30(IN), FINE(IN)
!VAST.../RADIAL/ YK(INOUT)
!VAST...Calls: ZK, P
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE ykk_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      SUBROUTINE ykk (I, J, K, KK) 
      INTEGER NOD 
      PARAMETER (NOD = 220) 
      INTEGER :: I 
!VAST...Dummy argument I is not referenced in this routine.
      INTEGER :: J 
!VAST...Dummy argument J is not referenced in this routine.
      INTEGER, INTENT(IN) :: K 
      INTEGER, INTENT(IN) :: KK 
!VAST.../PARAM/ H3(IN), EH(IN), NO(IN), ND(IN), D4(IN)
!VAST.../RADIAL/ YK(INOUT)
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE zeta_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      REAL(KIND(0.0D0)) FUNCTION zeta (I1, I2) 
      INTEGER NOD 
      PARAMETER (NOD = 220) 
      INTEGER :: I1 
      INTEGER :: I2 
!VAST.../BLUME/ COEFN2(IN), COEFNK(IN), COEFVK(IN)
!VAST.../PARAM/ Z(IN), FINE(IN), NCLOSD(IN)
!VAST...Calls: QUADR, L, SN, BWINT, VK
      END FUNCTION  
      END INTERFACE 
      END MODULE 
      MODULE zk_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      SUBROUTINE zk (I, J, K) 
      USE nel_C
      INTEGER NOD 
      PARAMETER (NOD = 220) 
      INTEGER, INTENT(IN) :: I 
      INTEGER, INTENT(IN) :: J 
      INTEGER, INTENT(IN) :: K 
!VAST.../PARAM/ H(IN), EH(IN), Z(IN), NO(IN), D1(IN)
!VAST.../RADIAL/ R(IN), RR(IN), YK(INOUT)
!VAST...Calls: L, P
      END SUBROUTINE  
      END INTERFACE 
      END MODULE 
      MODULE zz_I   
      INTERFACE
!...Generated by Pacific-Sierra Research 77to90  4.3E  22:22:18  11/14/01  
      REAL(KIND(0.0D0)) FUNCTION zz (I1, I2, I3, I4, K) 
      INTEGER :: I1 
      INTEGER :: I2 
      INTEGER :: I3 
      INTEGER :: I4 
      INTEGER, INTENT(IN) :: K 
!VAST...Calls: TK, L, UK, SN
      END FUNCTION  
      END INTERFACE 
      END MODULE 
