!     ..........................................................   :
!                                                                  :
!          Block                                                   :
!                       R E C O U P L S                            :
!                                                                  :
!     For Calculation Recoupling Coeficients                       :
!                                                                  :
!     Written by G. Gaigalas,                                      :
!                  Department  of  Computer Science,               :
!                  Vanderbilt University,  Nashville               :
!                                                  February 1994   :
!                                                                  :
!     ..........................................................   :
!
!     --------------------------------------------------------------
!     D L S A 1
!     --------------------------------------------------------------
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE DLSA1(K, JA1, KA, IRE, IAT, REC) 
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
      INTEGER  :: KA 
      INTEGER , INTENT(IN) :: IRE 
      INTEGER , INTENT(OUT) :: IAT 
      REAL(DOUBLE) , INTENT(OUT) :: REC 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA1, IB1, K1, J1, IT1, IT1S, JJ, IFAZ 
      REAL(DOUBLE) :: A1 
!-----------------------------------------------
      REC = ZERO 
      IA1 = J1QN1(JA1,K) - 1 
      IB1 = J1QN2(JA1,K) - 1 
      IF (JA1 > 2) THEN 
         K1 = IHSH + JA1 - 2 
         J1 = J1QN1(K1,K) - 1 
         K1 = K1 + 1 
         IT1 = J1QN1(K1,K) - 1 
         IT1S = J1QN2(K1,K) - 1 
      ELSE 
         JJ = 1 
         IF (JA1 == 1) JJ = 2 
         J1 = J1QN1(JJ,K) - 1 
         K1 = IHSH + 1 
         IT1 = J1QN1(K1,K) - 1 
         IT1S = J1QN2(K1,K) - 1 
      ENDIF 
      IF (IRE == 0) THEN 
         IF (IXJTIK(KA,IB1,IA1,J1,IT1,IT1S) == 0) RETURN  
         IAT = 1 
      ELSE 
         CALL SIXJ (KA, IB1, IA1, J1, IT1, IT1S, 0, A1) 
         A1 = A1*DSQRT(DBLE((IA1 + 1)*(IT1S + 1))) 
         IFAZ = J1 + IT1 + IB1 + KA 
         IF ((IFAZ/4)*4 /= IFAZ) A1 = -A1 
         REC = A1 
         IAT = 1 
         IF (JA1 == 1) THEN 
            IFAZ = IA1 + IB1 + 2*J1 - IT1 - IT1S 
            IF ((IFAZ/4)*4 /= IFAZ) REC = -REC 
         ENDIF 
      ENDIF 
      RETURN  
      END SUBROUTINE DLSA1 
