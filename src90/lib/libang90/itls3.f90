!
!     -------------------------------------------------------------
!      I T L S 3
!     -------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      INTEGER FUNCTION ITLS3 (IK, ID, KG1, KG2, BK, BD, IBT, BT, ITP, ITG, IQ) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
      USE RIBOLS_C 
      USE RIBOLSF_C 
      USE RIBOF_C 
      USE RIBOLS3_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:51:45  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ittk_I 
      USE mes_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: KG1 
      INTEGER  :: KG2 
      INTEGER , INTENT(OUT) :: ITP 
      INTEGER , INTENT(OUT) :: ITG 
      INTEGER , INTENT(IN) :: IQ 
      INTEGER  :: IK(7) 
      INTEGER  :: ID(7) 
      INTEGER , INTENT(OUT) :: IBT(7) 
      REAL(DOUBLE)  :: BK(3) 
      REAL(DOUBLE) , INTENT(IN) :: BD(3) 
      REAL(DOUBLE) , INTENT(OUT) :: BT(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ITK, ITD, ITP1, ITG1 
!-----------------------------------------------
      ITLS3 = 0 
      IF (ID(3) > 9) RETURN  
      IF (ITTK(ID(5),IK(5),KG1) == 0) RETURN  
      IF (ITTK(ID(6),IK(6),KG2) == 0) RETURN  
      ITK = IK(1) 
      ITD = ID(1) 
      IF (ID(3) < 3) THEN 
         ITP1 = IMPTLS(ITK) 
         ITP = IMPNLS(ITD) 
         IF (ITP1 /= ITP) RETURN  
         ITG1 = IMGTLS(ITK) 
         ITG = IMGNLS(ITD) 
      ELSE IF (ID(3) == 3) THEN 
         IF (ITK > 300) THEN 
            IF (ITD < 300) CALL MES (53) 
            IF (ID(4) > 2) CALL MES (13) 
            IF (IK(4) > 2) CALL MES (13) 
            ITK = ITK - 300 
            ITD = ITD - 300 
            ITP1 = IMPTLSF(ITK) 
            ITP = IMPNLSF(ITD) 
            IF (ITP1 /= ITP) RETURN  
            ITG1 = IMGTLSF(ITK) 
            ITG = IMGNLSF(ITD) 
         ELSE 
            IF (ITD > 300) CALL MES (53) 
            ITP1 = IMPTF(ITK) 
            ITP = IMPNF(ITD) 
            IF (ITP1 /= ITP) RETURN  
            ITG1 = IMGTF(ITK) 
            ITG = IMGNF(ITD) 
         ENDIF 
      ELSE 
         IF (ID(4) > 2) CALL MES (13) 
         IF (IK(4) > 2) CALL MES (13) 
         ITP1 = IMPTLS3(ITK) 
         ITP = IMPNLS3(ITD) 
         IF (ITP1 /= ITP) RETURN  
         ITG1 = IMGTLS3(ITK) 
         ITG = IMGNLS3(ITD) 
      ENDIF 
      IF (ITG1 /= ITG) RETURN  
      ITLS3 = 1 
      IBT(2) = ID(2) 
      IBT(3) = ID(3) 
      IBT(4) = ID(4) + IQ 
      BT(3) = BD(3) + HALF*DBLE(IQ) 
      RETURN  
      END FUNCTION ITLS3 
