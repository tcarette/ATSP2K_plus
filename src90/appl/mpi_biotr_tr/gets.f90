!
!     ------------------------------------------------------------------
!     G E T S
!     ------------------------------------------------------------------
!
      SUBROUTINE GETS(S, NSHLI, NSHLF) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!
! Obtain overlap matrix between all I shells and
! all F shells
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:40:05  11/20/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE quadr_I 
      USE wrtmat_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NSHLI 
      INTEGER  :: NSHLF 
      REAL(DOUBLE)  :: S(NSHLI,NSHLF) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: II, IF, NTEST 
!-----------------------------------------------
!
!
      DO II = 1, NSHLI 
         DO IF = 1, NSHLF 
            S(II,IF) = QUADR(II,IF + NSHLI,0) 
         END DO 
      END DO 
!
      NTEST = 0 
      IF (NTEST /= 0) THEN 
         WRITE (6, *) ' S matrix from GETS' 
         CALL WRTMAT (S, NSHLI, NSHLF, NSHLI, NSHLF) 
      ENDIF 
!
      RETURN  
      END SUBROUTINE GETS 
