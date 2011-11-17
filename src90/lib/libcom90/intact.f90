!
!     ------------------------------------------------------------------
!       I N T A C T
!     ------------------------------------------------------------------
!
      SUBROUTINE INTACT(L, LP, IEQUIV) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE ENAV_C 
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
      INTEGER , INTENT(IN) :: IEQUIV 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(0:4) :: ICBPTR 
      INTEGER :: I, K, L1, L2 
!-----------------------------------------------
      DATA ICBPTR/ 1, 6, 14, 23, 31/  
!
!     THIS SUBROUTINE GIVES THE INTERACTION ENERGY BETWEEN TWO SHELLS,
!     ONE WITH ORBITAL ANGULAR MOMENTUM  L , THE OTHER WITH ORBITAL
!     ANGULAR MOMENTUM  LP .   NOTICE THAT THE FIRST TERM OF THIS
!     INTERACTION ENERGY IS ALWAYS   F0(L,LP)   AND THIS IS NOT GIVEN
!     IN THIS SUBROUTINE.   THUS ONLY THE EXTRA TERMS ARE HERE PRODUCED.
!     FOR EQUIVALENT ELECTRONS (IEQUIV = 1) ,  THERE WILL BE  FK
!     INTEGRALS ONLY.   FOR NON-EQUIVALENT ELECTRONS (IEQUIV = 2) ,
!     THERE WILL BE  GK  INTEGRALS ONLY.
!
!     THE EXPRESSIONS FOR THE INTERACTION ENERGIES ARE GIVEN BY
!     R. D. COWAN, THE THEORY OF ATOMIC SPECTRA, EQUATIONS (6.38)
!      AND (6.39).
!
      I = 0 
      IF (IEQUIV == 1) THEN 
         DO K = 2, 2*L, 2 
            I = I + 1 
            KVALUE(I) = K 
            IF (L <= 4) THEN 
               COEFCT(I) = -CCA((L*(L-1)+K)/2) 
            ELSE 
               COEFCT(I) = -RME(L,L,K)**2/((2*L + 1)*(4*L + 1)) 
            ENDIF 
         END DO 
      ELSE 
         DO K = IABS(L - LP), L + LP, 2 
            I = I + 1 
            KVALUE(I) = K 
            IF (L <= LP) THEN 
               L1 = L 
               L2 = LP 
            ELSE 
               L1 = LP 
               L2 = L 
            ENDIF 
            IF (L2 <= 4) THEN 
               COEFCT(I) = -CCB(ICBPTR(L1)+(K+L1-L2)/2+(L1+1)*(L2-L1)) 
            ELSE 
               COEFCT(I) = -RME(L,LP,K)**2/(2*(2*L + 1)*(2*LP + 1)) 
            ENDIF 
         END DO 
      ENDIF 
      NINTS = I 
      RETURN  
      END SUBROUTINE INTACT 
