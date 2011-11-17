!
!     -----------------------------------------------------------------
!      R M E W P L S
!     -----------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                          December 1993   *
!
      SUBROUTINE RMEWPLS(J1, J2, K1, K2, K3, COEF) 
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
      INTEGER , INTENT(IN) :: J1 
      INTEGER , INTENT(IN) :: J2 
      INTEGER , INTENT(IN) :: K1 
      INTEGER , INTENT(IN) :: K2 
      INTEGER , INTENT(IN) :: K3 
      REAL(DOUBLE) , INTENT(OUT) :: COEF 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(3,3) :: IP120N, IP120L, IP111N, IP111L, IP021N, &
         IP021L 
      INTEGER :: JI1, JI2, I1, I2 
!-----------------------------------------------
      DATA IP120N/ 4*0, 60, 90, 0, 90, 0/  
      DATA IP120L/ 2*0, 40, 0, -90, 0, -40, 0, 70/  
      DATA IP111N/ 0, -72, 0, 72, 108, -90, 0, -90, 0/  
      DATA IP111L/ 0, 72, 0, -72, 108, -90, 0, -90, 0/  
      DATA IP021N/ 2*0, -40, 0, -90, 0, 40, 0, 70/  
      DATA IP021L/ 4*0, 60, 90, 0, 90, 0/  
      COEF = 0.00 
      JI1 = IMPTLS(J1) 
      JI2 = IMPTLS(J2) 
      IF (JI1 /= JI2) RETURN  
      IF (K1==1 .AND. K2==2 .AND. K3==0) THEN 
         IF (J1<6 .AND. J2<6) THEN 
            I1 = J1 - 2 
            I2 = J2 - 2 
            IF (IP120N(I1,I2) >= 0) THEN 
               COEF = DSQRT(DBLE(IP120N(I1,I2))) 
            ELSE 
               COEF = -DSQRT((-DBLE(IP120N(I1,I2)))) 
            ENDIF 
         ELSE IF (J1>5 .AND. J2>5) THEN 
            I1 = J1 - 5 
            I2 = J2 - 5 
            IF (IP120L(I1,I2) >= 0) THEN 
               COEF = DSQRT(DBLE(IP120L(I1,I2))) 
            ELSE 
               COEF = -DSQRT((-DBLE(IP120L(I1,I2)))) 
            ENDIF 
         ENDIF 
      ELSE IF (K1==1 .AND. K2==1 .AND. K3==1) THEN 
         IF (J1<6 .AND. J2<6) THEN 
            I1 = J1 - 2 
            I2 = J2 - 2 
            IF (IP111N(I1,I2) >= 0) THEN 
               COEF = DSQRT(DBLE(IP111N(I1,I2))) 
            ELSE 
               COEF = -DSQRT((-DBLE(IP111N(I1,I2)))) 
            ENDIF 
         ELSE IF (J1>5 .AND. J2>5) THEN 
            I1 = J1 - 5 
            I2 = J2 - 5 
            IF (IP111L(I1,I2) >= 0) THEN 
               COEF = DSQRT(DBLE(IP111L(I1,I2))) 
            ELSE 
               COEF = -DSQRT((-DBLE(IP111L(I1,I2)))) 
            ENDIF 
         ENDIF 
      ELSE IF (K1==0 .AND. K2==2 .AND. K3==1) THEN 
         IF (J1<6 .AND. J2<6) THEN 
            I1 = J1 - 2 
            I2 = J2 - 2 
            IF (IP021N(I1,I2) >= 0) THEN 
               COEF = DSQRT(DBLE(IP021N(I1,I2))) 
            ELSE 
               COEF = -DSQRT((-DBLE(IP021N(I1,I2)))) 
            ENDIF 
         ELSE IF (J1>5 .AND. J2>5) THEN 
            I1 = J1 - 5 
            I2 = J2 - 5 
            IF (IP021L(I1,I2) >= 0) THEN 
               COEF = DSQRT(DBLE(IP021L(I1,I2))) 
            ELSE 
               COEF = -DSQRT((-DBLE(IP021L(I1,I2)))) 
            ENDIF 
         ENDIF 
      ENDIF 
      RETURN  
      END SUBROUTINE RMEWPLS 
