 
!
!     ------------------------------------------------------------------
!     V E C S U M
!     ------------------------------------------------------------------
!
      SUBROUTINE VECSUM(C, A, B, FACA, FACB, NDIM) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!
!     CACLULATE THE VECTOR C(I)=FACA*A(I)+FACB*B(I)
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:48:02  11/20/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NDIM 
      REAL(DOUBLE) , INTENT(IN) :: FACA 
      REAL(DOUBLE) , INTENT(IN) :: FACB 
      REAL(DOUBLE) , INTENT(OUT) :: C(1) 
      REAL(DOUBLE) , INTENT(IN) :: A(1) 
      REAL(DOUBLE) , INTENT(IN) :: B(1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
      REAL(DOUBLE) :: S 
!-----------------------------------------------
!
      IF (FACA/=0.0D0 .AND. FACB/=0.0D0) THEN 
         C(:NDIM) = FACA*A(:NDIM) + FACB*B(:NDIM) 
!
      ELSE IF (FACA==0.0D0 .AND. FACB/=0.0D0) THEN 
         C(:NDIM) = FACB*B(:NDIM) 
!
      ELSE IF (FACA/=0.0D0 .AND. FACB==0.0D0) THEN 
         C(:NDIM) = FACA*A(:NDIM) 
!
      ELSE IF (FACA==0.0D0 .AND. FACB==0.0D0) THEN 
         C(:NDIM) = 0.0D0 
      ENDIF 
!
      RETURN  
      END SUBROUTINE VECSUM 
