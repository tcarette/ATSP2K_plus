!
!     -------------------------------------------------------------
!      S L S
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                    (qls)         *
!     REDUCED MATRIX ELEMENT             (nl QLS::: a  :::nl QLS)  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      SUBROUTINE SLS(L, IT, LQ, LL, LS, ITS, LQS, LLS, LSS, S) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:56:36  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE subls_I 
      USE frma_I 
      USE rmeals_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: L 
      INTEGER  :: IT 
      INTEGER  :: LQ 
      INTEGER  :: LL 
      INTEGER  :: LS 
      INTEGER  :: ITS 
      INTEGER  :: LQS 
      INTEGER  :: LLS 
      INTEGER  :: LSS 
      REAL(DOUBLE)  :: S 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(5,2) :: IC 
      INTEGER :: IE1 
!-----------------------------------------------
      DATA IC/ 3, 5, 6, 8, 2, 9, 24, 25, 40, 5/  
      IF (L > 3) THEN 
         IF (L > 9) RETURN  
         IF (LSS == 1) THEN 
            CALL SUBLS (IT, ITS, L, S) 
         ELSE 
            CALL SUBLS (ITS, IT, L, S) 
            IE1 = LQ + LL + LS - LQS - LLS - LSS + 2*L 
            IF ((IE1/4)*4 /= IE1) S = -S 
         ENDIF 
      ELSE IF (L == 3) THEN 
         IF (IT > 300) THEN 
            IF (ITS < 300) THEN 
               WRITE (6, '(A)') ' ERROR IN    S L S ' 
               STOP  
            ENDIF 
            IF (LSS == 1) THEN 
               CALL SUBLS (IT, ITS, L, S) 
            ELSE 
               CALL SUBLS (ITS, IT, L, S) 
               IE1 = LQ + LL + LS - LQS - LLS - LSS + 2*L 
               IF ((IE1/4)*4 /= IE1) S = -S 
            ENDIF 
         ELSE 
            CALL FRMA (IT, LQ, LL, LS, ITS, LQS, LLS, LSS, S) 
         ENDIF 
      ELSE 
         CALL RMEALS (L, IT, LQ, LL, LS, ITS, LQS, LLS, LSS, S) 
      ENDIF 
      RETURN  
      END SUBROUTINE SLS 
