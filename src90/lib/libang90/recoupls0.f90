!
!     --------------------------------------------------------------
!     R E C O U P L S 0
!     --------------------------------------------------------------
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      LOGICAL FUNCTION RECOUPLS0 (K, JA1, JA2, JA3, JA4, KA) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE MEDEFN_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:00:59  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: K 
      INTEGER , INTENT(IN) :: JA1 
      INTEGER , INTENT(IN) :: JA2 
      INTEGER , INTENT(IN) :: JA3 
      INTEGER , INTENT(IN) :: JA4 
      INTEGER , INTENT(IN) :: KA 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IJ, IA1, IA2 
!-----------------------------------------------
      RECOUPLS0 = .TRUE. 
      IF (IHSH == 1) RETURN  
      IF (JA1==1 .AND. JA2==2) GO TO 1 
      IF (KA /= 0) GO TO 5 
!
!  CASES WHEN :          KA = 0
!                  OR    JA1 = JA2
!                  OR    JA1 = 1    JA2 = 2
!
    1 CONTINUE 
      DO I = 1, IHSH 
         IJ = IHSH + I - 1 
         IF (I /= 1) THEN 
            IF (J1QN1(IJ,K) /= J1QN2(IJ,K)) RECOUPLS0 = .FALSE. 
         ENDIF 
         IF (KA == 0) GO TO 9 
         IF (I == JA1) CYCLE  
         IF (I == JA2) CYCLE  
    9    CONTINUE 
         IF (I==JA1 .AND. K==1) CYCLE  
         IF (I==JA2 .AND. K==1) CYCLE  
         IF (I==JA3 .AND. K==1) CYCLE  
         IF (I==JA4 .AND. K==1) CYCLE  
         IF (J1QN1(I,K) == J1QN2(I,K)) CYCLE  
         RECOUPLS0 = .FALSE. 
      END DO 
      RETURN  
!
!  OTHER CASES
!
    5 CONTINUE 
      IA1 = JA1 - 1 
      IA2 = JA2 - 1 
      IF (JA1 == 1) IA1 = JA1 
      DO I = 1, IHSH 
         IJ = IHSH + I 
         IF (I /= IHSH) THEN 
            IF (I<IA1 .OR. I>=IA2) THEN 
               IF (J1QN1(IJ,K) /= J1QN2(IJ,K)) RECOUPLS0 = .FALSE. 
            ENDIF 
         ENDIF 
         IF (I == JA1) CYCLE  
         IF (I == JA2) CYCLE  
         IF (KA==2 .AND. I==JA3) CYCLE  
         IF (KA==3 .AND. I==JA3) CYCLE  
         IF (KA==3 .AND. I==JA4) CYCLE  
         IF (J1QN1(I,K) == J1QN2(I,K)) CYCLE  
         RECOUPLS0 = .FALSE. 
      END DO 
      RETURN  
      END FUNCTION RECOUPLS0 
