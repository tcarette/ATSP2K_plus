!
!     -------------------------------------------------------------
!      C L E 0 S M
!     -------------------------------------------------------------
!
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING              *
!                                                 ---         ---  *
!                                                 I  Q   S  C   I  *
!     CLEBSCH - GORDAN COEFFICIENT:               I             I  *
!                                                 I  QM  0  CM  I  *
!                                                 ---         ---  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                          December 1993   *
!
      SUBROUTINE CLE0SM(Q, QM, S, C, CM, A) 
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
      USE c1e0sm_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE)  :: Q 
      REAL(DOUBLE)  :: QM 
      REAL(DOUBLE) , INTENT(IN) :: S 
      REAL(DOUBLE)  :: C 
      REAL(DOUBLE)  :: CM 
      REAL(DOUBLE)  :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IIQ, IIC, IIS 
!-----------------------------------------------
      A = ZERO 
      IIQ = TWO*Q + TENTH 
      IIC = TWO*C + TENTH 
      IIS = TWO*S + TENTH 
      IF (ITTK(IIQ,IIC,IIS) == 0) RETURN  
      IF (S >= EPS) THEN 
         CALL C1E0SM (Q, QM, C, CM, A) 
         RETURN  
      ENDIF 
      IF (Q + TENTH < ABS(QM)) RETURN  
      IF (C + TENTH < ABS(CM)) RETURN  
      IF (ABS(Q - C) > EPS) RETURN  
      IF (ABS(QM - CM) > EPS) RETURN  
      A = ONE 
      RETURN  
      END SUBROUTINE CLE0SM 
