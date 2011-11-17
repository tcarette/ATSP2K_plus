!
!     --------------------------------------------------------------
!     R E C O U P L S 2
!     --------------------------------------------------------------
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE RECOUPLS2(K, JA1, JA2, KA, IRE, IAT, REC) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE MEDEFN_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:00:59  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dlsa2_I 
      USE dlsa1_I 
      USE dlsa3_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: K 
      INTEGER  :: JA1 
      INTEGER  :: JA2 
      INTEGER  :: KA 
      INTEGER  :: IRE 
      INTEGER  :: IAT 
      REAL(DOUBLE) , INTENT(OUT) :: REC 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA1, IB1, IA2, IB2, ISKR 
      REAL(DOUBLE) :: S, SS, RE 
!-----------------------------------------------
      IAT = 0 
      S = DBLE(J1QN1(JA1,K)) 
      SS = DBLE(J1QN1(JA2,K)) 
      REC = ONE/DSQRT(S*SS) 
      IF (IRE/=0 .AND. KA==0) THEN 
         IAT = 1 
      ELSE 
         IA1 = J1QN1(JA1,K) - 1 
         IB1 = J1QN2(JA1,K) - 1 
         IA2 = J1QN1(JA2,K) - 1 
         IB2 = J1QN2(JA2,K) - 1 
         IAT = 0 
         CALL DLSA2 (K, JA1, JA2, KA, IRE, IAT, RE) 
         IF (IAT == 0) RETURN  
         REC = RE*REC*DSQRT(DBLE(IA2 + 1))/DSQRT(DBLE((KA + 1)*(IB2 + 1))) 
         IF (JA1==1 .AND. JA2==2) RETURN  
         IAT = 0 
         CALL DLSA1 (K, JA1, KA, IRE, IAT, RE) 
         IF (IAT == 0) RETURN  
         REC = RE*REC 
         ISKR = JA2 - JA1 
         IF (JA1 == 1) ISKR = JA2 - 1 - JA1 
         IF (ISKR <= 1) RETURN  
         IAT = 0 
         CALL DLSA3 (K, JA1, JA2, KA, IRE, IAT, RE) 
         REC = RE*REC 
      ENDIF 
      RETURN  
      END SUBROUTINE RECOUPLS2 
