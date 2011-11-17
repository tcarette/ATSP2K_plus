!=======================================================================
      SUBROUTINE MULTBC(N, K, M, C, TEMP, B) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!=======================================================================
!       called by: DVDRVR
!
!       Multiplies B(N,K)*C(K,M) and stores it in B(N,M)
!       Used for collapsing the expanding basis to current estimates,
!       when basis becomes too large, or for returning the results back
!
!       Subroutines called
!       DINIT, DGEMV, DCOPY
!-----------------------------------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:22  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dgemv_I 
      USE dcopy_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N 
      INTEGER  :: K 
      INTEGER  :: M 
      REAL(DOUBLE)  :: C(K*M) 
      REAL(DOUBLE)  :: TEMP(M) 
      REAL(DOUBLE)  :: B(N*K) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IROW 
!-----------------------------------------------
!-----------------------------------------------------------------------
      DO IROW = 1, N 
!              CALL DINIT(M,0.d0,TEMP,1)
         CALL DGEMV ('Transp', K, M, 1.D0, C, K, B(IROW), N, 0.D0, TEMP, 1) 
         CALL DCOPY (M, TEMP, 1, B(IROW), N) 
      END DO 
 
      RETURN  
      END SUBROUTINE MULTBC 
