!
!     -------------------------------------------------------------
!      C 1 E 0 S M
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                 ---         ---  *
!                                                 I  Q   1  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
!                                                 I  QM  0  CM  I  *
!                                                 ---         ---  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                          December 1993   *
!
      SUBROUTINE C1E0SM(Q, QM, C, CM, A) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:48:27  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ittk_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(IN) :: Q 
      REAL(DOUBLE) , INTENT(IN) :: QM 
      REAL(DOUBLE) , INTENT(IN) :: C 
      REAL(DOUBLE) , INTENT(IN) :: CM 
      REAL(DOUBLE) , INTENT(OUT) :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IIQ, IIC, IS, IG 
!-----------------------------------------------
      A = ZERO 
      IIQ = TWO*Q + TENTH 
      IIC = TWO*C + TENTH 
      IF (ITTK(IIQ,IIC,2) == 0) RETURN  
      IF (ABS(QM - CM) <= EPS) THEN 
         IF (Q + TENTH >= ABS(QM)) THEN 
            IF (C + TENTH >= ABS(CM)) THEN 
               IF (ABS(QM) > EPS) GO TO 5 
               IS = Q + C + ONE + TENTH 
               IF ((IS/2)*2 /= IS) GO TO 4 
    5          CONTINUE 
               IG = Q - C + TWO + TENTH 
               IF (IG > 0) THEN 
                  IF (IG <= 3) THEN 
                     SELECT CASE (IG)  
                     CASE DEFAULT 
                        A = -SQRT(((C + CM + ONE)*(C - CM + ONE))/((C + ONE)*(&
                           TWO*C + THREE))) 
                        RETURN  
                     CASE (2)  
                        A = CM/SQRT(C*(C + ONE)) 
                        RETURN  
                     CASE (1)  
                        A = SQRT(((C + CM)*(C - CM))/((TWO*C - ONE)*C)) 
                        RETURN  
                     END SELECT 
                  ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
    4 CONTINUE 
      A = ZERO 
      RETURN  
      END SUBROUTINE C1E0SM 
