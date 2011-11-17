!
!     -------------------------------------------------------------
!      A C O N S T
!     -------------------------------------------------------------
!                                                                  *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vilnius, LITHUANIA                              January 1997 *
!
      SUBROUTINE ACONST(IQ, IP, IR, IS, IT, A) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:09:47  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IQ 
      INTEGER , INTENT(IN) :: IP 
      INTEGER , INTENT(IN) :: IR 
      INTEGER , INTENT(IN) :: IS 
      INTEGER , INTENT(IN) :: IT 
      REAL(DOUBLE) , INTENT(OUT) :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: AI 
!-----------------------------------------------
      AI = DBLE(IR + IQ - IP)*DBLE(IP - IR + IQ)*DBLE(IP + IR - IQ + 2)*DBLE(IP&
          + IR + IQ + 2)*DBLE(IT + IQ - IS)*DBLE(IS - IT + IQ)*DBLE(IS + IT - &
         IQ + 2)*DBLE(IS + IT + IQ + 2) 
      A = DSQRT(AI/DBLE(256)) 
      RETURN  
      END SUBROUTINE ACONST 
