!
!     -------------------------------------------------------------
!      C 1 E 1 S M
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                 ---         ---  *
!                                                 I  Q   1  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
!                                                 I  QM  1  CM  I  *
!                                                 ---         ---  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                          December 1993   *
!
      SUBROUTINE C1E1SM(Q, QM, SM, C, CM, A) 
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
      REAL(DOUBLE) , INTENT(IN) :: SM 
      REAL(DOUBLE) , INTENT(IN) :: C 
      REAL(DOUBLE) , INTENT(IN) :: CM 
      REAL(DOUBLE) , INTENT(OUT) :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IIQ, IIC, IE 
      REAL(DOUBLE), DIMENSION(2) :: GC 
!-----------------------------------------------
      GC(1) = ONE 
      GC(2) = -ONE 
      A = ZERO 
      IIQ = TWO*Q + TENTH 
      IIC = TWO*C + TENTH 
      IF (ITTK(IIQ,IIC,2) == 0) RETURN  
      IF (ABS(QM + SM - CM) > EPS) RETURN  
      IF (Q + TENTH < ABS(QM)) RETURN  
      IF (C + TENTH < ABS(CM)) RETURN  
      IE = 0 
      IF (ABS(SM - ONE) < EPS) IE = 1 
      IF (ABS(SM + ONE) < EPS) IE = 2 
      IF (IE == 0) RETURN  
      IF (ABS(Q + ONE - C) >= EPS) THEN 
         IF (ABS(Q - C) < EPS) GO TO 2 
         IF (ABS(Q - ONE - C) > EPS) RETURN  
         A = SQRT((C - GC(IE)*CM+ONE)*(C-GC(IE)*CM+TWO)/((TWO*C+TWO)*(TWO*C+&
            THREE))) 
         RETURN  
      ENDIF 
      A = SQRT((C + GC(IE)*CM-ONE)*(C+GC(IE)*CM)/((TWO*C-ONE)*TWO*C)) 
      RETURN  
    2 CONTINUE 
      A = -GC(IE)*SQRT((C + GC(IE)*CM)*(C-GC(IE)*CM+ONE)/((C+ONE)*TWO*C)) 
      RETURN  
      END SUBROUTINE C1E1SM 
