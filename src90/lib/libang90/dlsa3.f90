!
!     --------------------------------------------------------------
!      D L S A 3
!     --------------------------------------------------------------
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE DLSA3(K, JA1, JA2, KA, IRE, IAT, REC) 
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
      INTEGER :: IMAX, IMIN, I, JI, K1, ITI, ITIS, ITI1, ITI1S, IFAZ 
      REAL(DOUBLE) :: AA, A3 
!-----------------------------------------------
      REC = ZERO 
      AA = ONE 
      IMAX = JA2 - 1 
      IMIN = JA1 + 1 
      IF (JA1 == 1) IMIN = IMIN + 1 
      IF (IMIN < JA2) THEN 
         DO I = IMIN, IMAX 
            JI = J1QN1(I,K) - 1 
            K1 = IHSH + I - 2 
            ITI = J1QN1(K1,K) - 1 
            ITIS = J1QN2(K1,K) - 1 
            K1 = K1 + 1 
            ITI1 = J1QN1(K1,K) - 1 
            ITI1S = J1QN2(K1,K) - 1 
            IF (IRE == 0) THEN 
               IF (IXJTIK(KA,ITIS,ITI,JI,ITI1,ITI1S) == 0) RETURN  
            ELSE 
               CALL SIXJ (KA, ITIS, ITI, JI, ITI1, ITI1S, 0, A3) 
               A3 = A3*DSQRT(DBLE((ITI + 1)*(ITI1S + 1))) 
               IFAZ = KA + JI + ITI + ITI1S 
               IF ((IFAZ/4)*4 /= IFAZ) A3 = -A3 
               AA = AA*A3 
            ENDIF 
         END DO 
      ENDIF 
      REC = AA 
      IAT = 1 
      RETURN  
      END SUBROUTINE DLSA3 
