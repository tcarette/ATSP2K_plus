!
!     -------------------------------------------------------------
!      T R A N S I T I O N
!     -------------------------------------------------------------
!
!     CALCULATE TRANSITION OPERATORS                               *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville           September 1997   *
!
      SUBROUTINE TRANSITION(L1, L2, I, A) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE1=>DOUBLE 
      use consts_C
      use medefn_C
      use ems_C
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:28:09  11/20/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: L1 
      INTEGER  :: L2 
      INTEGER , INTENT(OUT) :: I 
      REAL(DOUBLE) , INTENT(OUT) :: A 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!      INTEGER :: IFL 
!-----------------------------------------------
!      LOGICAL REL,VOK
!      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
!      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
!     :     J1QN2(31,3),IJFUL(16)
!      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      I = 200 
      IF (IFL /= 4) THEN 
         A = -DSQRT(DBLE(TWO*J1QN1(2*IHSH-1,2))) 
      ELSE 
         A = -DSQRT(DBLE(J1QN1(2*IHSH-1,2)*J1QN1(2*IHSH-1,3))) 
      ENDIF 
      RETURN  
      END SUBROUTINE TRANSITION 
