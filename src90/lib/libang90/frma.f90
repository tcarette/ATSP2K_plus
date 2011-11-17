!     ..........................................................   :
!                                                                  :
!          Block                                                   :
!                  Standard Quantities  -  S Q L S F               :
!                         Part Two                                 :
!                                                                  :
!     Written by G. Gaigalas,                                      :
!                Institute of Theoretical Physics and Astronomy    :
!                Vilnius,  Lithuania                               :
!                                                       May 1995   :
!                                                                  :
!     ..........................................................   :
!
!
!     -----------------------------------------------------------------
!      F R M A
!     -----------------------------------------------------------------
!
      SUBROUTINE FRMA(J1, LQ, LL, LS, J2, LQS, LLS, LSS, COEF) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE RIBOF_C 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:10:49  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE frma01_I 
      USE frma02_I 
      USE frma03_I 
      USE frma04_I 
      USE frma05_I 
      USE frma06_I 
      USE frma07_I 
      USE frma08_I 
      USE frma09_I 
      USE frma10_I 
      USE frma11_I 
      USE frma12_I 
      USE frma13_I 
      USE frma14_I 
      USE frma15_I 
      USE frma16_I 
      USE frma17_I 
      USE frma18_I 
      USE frma19_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
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
      INTEGER , DIMENSION(78) :: IV 
      INTEGER :: ISKA, JI1, JI2, JJ1, JJ2, IQ, IL, IS, IQS, ILS, ISS, IVAR, IE1 
!-----------------------------------------------
      DATA IV/ 1, 49, 147, 1, 539, 49, 33, 1, 11, 1176, 12936, 2646, 4116, &
         181104, 296352, 543312, 7392, 2, 162, 22176, 20790, 66528, 3234, 0, 0&
         , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
         980628, 1, 2, 196, 1176, 51744, 7392, 25872, 37044, 5488, 4116, &
         1992144, 121968, 2646, 1008, 1, 0, 0, 0, 0, 0, 0, 0, 12603360, 491040&
         , 11319, 142688, 7546, 0, 0, 0, 0/  
      ISKA = 0 
      JI1 = IMPTF(J1) 
      JI2 = IMPNF(J2) 
      IF (JI1 /= JI2) RETURN  
      IF (J1 < J2) THEN 
         JJ1 = J1 
         JJ2 = J2 
         IQ = LQ 
         IL = LL 
         IS = LS 
         IQS = LQS 
         ILS = LLS 
         ISS = LSS 
      ELSE 
         JJ1 = J2 
         JJ2 = J1 
         IQ = LQS 
         IL = LLS 
         IS = LSS 
         IQS = LQ 
         ILS = LL 
         ISS = LS 
      ENDIF 
      IF (JJ1 < 12) THEN 
         CALL FRMA01 (JJ1, JJ2, ISKA) 
         IVAR = IV(JJ1) 
      ELSE IF (JJ1 < 18) THEN 
         CALL FRMA02 (JJ1, JJ2, ISKA) 
         IVAR = IV(JJ1) 
      ELSE IF (JJ1 < 24) THEN 
         CALL FRMA03 (JJ1, JJ2, ISKA) 
         IVAR = IV(JJ1) 
      ELSE IF (JJ1 < 29) THEN 
         CALL FRMA04 (JJ1, JJ2, ISKA, IV(JJ1)) 
         IVAR = IV(JJ1) 
      ELSE IF (JJ1 < 34) THEN 
         CALL FRMA05 (JJ1, JJ2, ISKA, IV(JJ1)) 
         IVAR = IV(JJ1) 
      ELSE IF (JJ1 < 40) THEN 
         CALL FRMA06 (JJ1, JJ2, ISKA, IV(JJ1)) 
         IVAR = IV(JJ1) 
      ELSE IF (JJ1 < 47) THEN 
         CALL FRMA07 (JJ1, JJ2, ISKA, IV(JJ1)) 
         IVAR = IV(JJ1) 
      ELSE IF (JJ1 < 57) THEN 
         CALL FRMA08 (JJ1, JJ2, ISKA) 
         IVAR = IV(JJ1) 
      ELSE IF (JJ1 < 63) THEN 
         CALL FRMA09 (JJ1, JJ2, ISKA) 
         IVAR = IV(JJ1) 
      ELSE IF (JJ1 < 70) THEN 
         CALL FRMA10 (JJ1, JJ2, ISKA, IV(JJ1)) 
         IVAR = IV(JJ1) 
      ELSE IF (JJ1 < 75) THEN 
         CALL FRMA11 (JJ1, JJ2, ISKA) 
         IVAR = IV(JJ1) 
      ELSE IF (JJ1 < 79) THEN 
         CALL FRMA12 (JJ1, JJ2, ISKA, IV(JJ1)) 
         IVAR = IV(JJ1) 
      ELSE IF (JJ1 < 86) THEN 
         CALL FRMA13 (JJ1, JJ2, ISKA, IVAR) 
      ELSE IF (JJ1 < 91) THEN 
         CALL FRMA14 (JJ1, JJ2, ISKA, IVAR) 
      ELSE IF (JJ1 < 96) THEN 
         CALL FRMA15 (JJ1, JJ2, ISKA, IVAR) 
      ELSE IF (JJ1 < 101) THEN 
         CALL FRMA16 (JJ1, JJ2, ISKA, IVAR) 
      ELSE IF (JJ1 < 107) THEN 
         CALL FRMA17 (JJ1, JJ2, ISKA, IVAR) 
      ELSE IF (JJ1 < 113) THEN 
         CALL FRMA18 (JJ1, JJ2, ISKA, IVAR) 
      ELSE IF (JJ1 < 120) THEN 
         CALL FRMA19 (JJ1, JJ2, ISKA, IVAR) 
      ENDIF 
      IF (ISKA > 0) THEN 
         COEF = DSQRT(DBLE(ISKA)/DBLE(IVAR)) 
      ELSE IF (ISKA == 0) THEN 
         COEF = ZERO 
         RETURN  
      ELSE 
         COEF = -DSQRT((-DBLE(ISKA)/DBLE(IVAR))) 
      ENDIF 
      IF (J1 > J2) THEN 
         IE1 = LQ + LL + LS - LQS - LLS - LSS + 2*3 
         IF ((IE1/4)*4 /= IE1) COEF = -COEF 
      ENDIF 
      RETURN  
      END SUBROUTINE FRMA 
