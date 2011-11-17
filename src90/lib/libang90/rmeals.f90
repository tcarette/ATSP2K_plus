!
!     -----------------------------------------------------------------
!      R M E A L S
!     -----------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                          December 1993   *
!
      SUBROUTINE RMEALS(L, J1, LQ, LL, LS, J2, LQS, LLS, LSS, COEF) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE RIBOLS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:13:23  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: L 
      INTEGER , INTENT(IN) :: J1 
      INTEGER , INTENT(IN) :: LQ 
      INTEGER , INTENT(IN) :: LL 
      INTEGER , INTENT(IN) :: LS 
      INTEGER , INTENT(IN) :: J2 
      INTEGER , INTENT(IN) :: LQS 
      INTEGER , INTENT(IN) :: LLS 
      INTEGER , INTENT(IN) :: LSS 
      REAL(DOUBLE) , INTENT(OUT) :: COEF 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(3,3) :: IP 
      INTEGER , DIMENSION(16,16) :: IDV, IDS 
      INTEGER :: JI1, JI2, JJ1, JJ2, I1, I2, IE1 
!-----------------------------------------------
      DATA IP/ 0, -24, 0, -24, -54, -30, 0, -30, 30/  
!
      DATA IDS/ 4*0, -60, 16*0, -12, 0, 8, 10*0, 384, -84, -90, 42, 2*0, 336, &
         -336, 8*0, 84, 30, 0, -12, -36, -18, -96, -30, -18, 5*0, -60, 0, -90, &
         3*0, 50, 0, -210, 3*0, 90, 4*0, -12, 42, 12, 0, 360, 270, -60, -18, 18&
         , -30, 270, -150, -66, 5*0, -36, -50, -270, 3*0, 24, 0, 360, 5*0, 8, 0&
         , 18, 0, -60, 0, -90, 0, 27, 5, -45, 0, -99, 4*0, 336, 96, -210, -18, &
         2*0, -1344, -84, 0, -108, 2*0, 132, 3*0, 336, -30, 0, -18, 24, -27, 84&
         , -105, 63, 243, 60, -165, 33, 4*0, -18, 0, 30, 0, -5, 0, 63, 105, -9&
         , 0, 2*99, 6*0, 270, -360, -45, -108, -243, 9, 8019, 396, 1521, -297, &
         -234, 4*0, -90, 150, 3*0, 60, 0, -396, 2*0, -132, 6*0, -66, 0, -99, 0&
         , 165, -99, 1521, 0, 3645, -3, 234, 8*0, -132, 33, 99, 297, -132, 3, &
         858, -546, 11*0, -234, 0, 234, 546, 910/  
!
      DATA IDV/ 34*1, 2*5, 4*1, 2*5, 8*1, 5, 5*1, 5, 28*1, 3*7, 3*1, 3*7, 7*1, &
         7, 5*1, 7, 9*1, 7, 1, 7, 3*1, 7, 1, 7, 4*1, 2*5, 4*1, 2*5, 8*1, 5, 5*1&
         , 5, 3*4, 1, 4, 11*1, 2*4, 20, 1, 4, 5, 6*1, 3*7, 1, 4, 20, 140, 7, 28&
         , 2*5, 5*1, 7, 5*1, 7, 9*1, 7, 1, 7, 1, 2*4, 28, 1, 308, 1, 11, 10*1, &
         2*5, 2*1, 2*5, 11*1, 5, 1, 11, 5, 11/  
!
      COEF = 0.00 
      JI1 = IMPTLS(J1) 
      JI2 = IMPNLS(J2) 
      IF (JI1 /= JI2) RETURN  
      IF (J1 < J2) THEN 
         JJ1 = J1 
         JJ2 = J2 
      ELSE 
         JJ1 = J2 
         JJ2 = J1 
      ENDIF 
      SELECT CASE (L)  
      CASE (0)  
         COEF = -2.00 
      CASE (1)  
         I1 = JJ1 - 2 
         I2 = JJ2 - 5 
         IF (IP(I1,I2) >= 0) THEN 
            COEF = DSQRT(DBLE(IP(I1,I2))) 
         ELSE 
            COEF = -DSQRT((-DBLE(IP(I1,I2)))) 
         ENDIF 
      CASE (2)  
         I1 = JJ1 - 8 
         I2 = JJ2 - 24 
         IF (IDS(I1,I2) >= 0) THEN 
            COEF = DSQRT(DBLE(IDS(I1,I2))/DBLE(IDV(I1,I2))) 
         ELSE 
            COEF = -DSQRT((-DBLE(IDS(I1,I2))/DBLE(IDV(I1,I2)))) 
         ENDIF 
      END SELECT 
      IF (J1 > J2) THEN 
         IE1 = LQ + LL + LS - LQS - LLS - LSS + 2*L 
         IF ((IE1/4)*4 /= IE1) COEF = -COEF 
      ENDIF 
      RETURN  
      END SUBROUTINE RMEALS 
