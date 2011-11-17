!
!     --------------------------------------------------------------
!      D L S A 2
!     --------------------------------------------------------------
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE DLSA2(K, JA1, JA2, KA, IRE, IAT, REC) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:15:35  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ixjtik_I 
      USE sixj_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: K 
      INTEGER , INTENT(IN) :: JA1 
      INTEGER , INTENT(IN) :: JA2 
      INTEGER  :: KA 
      INTEGER , INTENT(IN) :: IRE 
      INTEGER , INTENT(OUT) :: IAT 
      REAL(DOUBLE) , INTENT(OUT) :: REC 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA1, IB1, IA2, IB2, IT2, IT2S, N1, J2, N2, IFAZ 
      REAL(DOUBLE) :: A2 
!-----------------------------------------------
      REC = ZERO 
      IA1 = J1QN1(JA1,K) - 1 
      IB1 = J1QN2(JA1,K) - 1 
      IA2 = J1QN1(JA2,K) - 1 
      IB2 = J1QN2(JA2,K) - 1 
      IF (JA1==1 .AND. JA2==2) THEN 
         IT2 = IA1 
         IT2S = IB1 
         N1 = IHSH + 1 
         J2 = J1QN1(N1,K) - 1 
      ELSE 
         N1 = IHSH + JA2 - 1 
         J2 = J1QN1(N1,K) - 1 
         N2 = IHSH + JA2 - 2 
         IT2 = J1QN1(N2,K) - 1 
         IT2S = J1QN2(N2,K) - 1 
      ENDIF 
      IF (IRE == 0) THEN 
         IF (IXJTIK(KA,IB2,IA2,J2,IT2,IT2S) == 0) RETURN  
         IAT = 1 
      ELSE 
         CALL SIXJ (KA, IB2, IA2, J2, IT2, IT2S, 0, A2) 
         REC = A2*DSQRT(DBLE((IB2 + 1)*(IT2 + 1))) 
         IFAZ = J2 + IT2S + IA2 + KA 
         IF ((IFAZ/4)*4 /= IFAZ) REC = -REC 
         IAT = 1 
      ENDIF 
      RETURN  
      END SUBROUTINE DLSA2 
