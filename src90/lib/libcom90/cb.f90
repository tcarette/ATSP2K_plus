!
!     -----------------------------------------------------------------
!          C B
!     -----------------------------------------------------------------
!
!
      REAL(KIND(0.0D0)) FUNCTION CB (L, LP, K) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE EAV_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:52:17  11/14/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rme_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: L 
      INTEGER  :: LP 
      INTEGER  :: K 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(0:4) :: ICBPTR 
      INTEGER :: L1, L2 
!-----------------------------------------------
      DATA ICBPTR/ 1, 6, 14, 23, 31/  
!
      IF (L <= LP) THEN 
         L1 = L 
         L2 = LP 
      ELSE 
         L1 = LP 
         L2 = L 
      ENDIF 
      IF (L2 <= 4) THEN 
         CB = CCB(ICBPTR(L1)+(K+L1-L2)/2+(L1+1)*(L2-L1)) 
      ELSE 
         CB = RME(L,LP,K)**2/(2*(2*L + 1)*(2*LP + 1)) 
      ENDIF 
      RETURN  
      END FUNCTION CB 
