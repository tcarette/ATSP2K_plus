!
!     -------------------------------------------------------------
!      I T L S
!     -------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      INTEGER FUNCTION ITLS (IK, ID, KG, KGG, BD, IBT, BT, KG1, KGG1, ITP, ITG&
         , IQ) 
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
      INTEGER , INTENT(IN) :: KG 
      INTEGER , INTENT(IN) :: KGG 
      INTEGER  :: KG1 
      INTEGER  :: KGG1 
      INTEGER , INTENT(OUT) :: ITP 
      INTEGER , INTENT(OUT) :: ITG 
      INTEGER , INTENT(IN) :: IQ 
      INTEGER  :: IK(7) 
      INTEGER  :: ID(7) 
      INTEGER , INTENT(OUT) :: IBT(7) 
      REAL(DOUBLE) , INTENT(IN) :: BD(3) 
      REAL(DOUBLE) , INTENT(OUT) :: BT(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ITK, ITD, ITP1, ITG1 
!-----------------------------------------------
      ITLS = 0 
      IF (ID(3) > 9) RETURN  
      KG1 = 2*KG 
      KGG1 = 2*KGG 
      IF (ITTK(ID(5),IK(5),KG1) == 0) RETURN  
      IF (ITTK(ID(6),IK(6),KGG1) == 0) RETURN  
      ITK = IK(1) 
      ITD = ID(1) 
      IF (ID(3) < 3) THEN 
         ITP1 = IMPTLS(ITK) 
         ITP = IMPTLS(ITD) 
         IF (ITP1 /= ITP) RETURN  
         ITG1 = IMGTLS(ITK) 
         ITG = IMGTLS(ITD) 
      ELSE IF (ID(3) == 3) THEN 
         IF (ITK > 300) THEN 
            IF (ITD < 300) CALL MES (51) 
            IF (ID(4) > 2) CALL MES (11) 
            IF (IK(4) > 2) CALL MES (11) 
            ITK = ITK - 300 
            ITD = ITD - 300 
            ITP1 = IMPTLSF(ITK) 
            ITP = IMPTLSF(ITD) 
            IF (ITP1 /= ITP) RETURN  
            ITG1 = IMGTLSF(ITK) 
            ITG = IMGTLSF(ITD) 
         ELSE 
            IF (ITD > 300) CALL MES (51) 
            ITP1 = IMPTF(ITK) 
            ITP = IMPTF(ITD) 
            IF (ITP1 /= ITP) RETURN  
            ITG1 = IMGTF(ITK) 
            ITG = IMGTF(ITD) 
         ENDIF 
      ELSE 
         IF (ID(4) > 2) CALL MES (11) 
         IF (IK(4) > 2) CALL MES (11) 
         ITP1 = IMPTLS3(ITK) 
         ITP = IMPTLS3(ITD) 
         IF (ITP1 /= ITP) RETURN  
         ITG1 = IMGTLS3(ITK) 
         ITG = IMGTLS3(ITD) 
      ENDIF 
      IF (ITG1 /= ITG) RETURN  
      ITLS = 1 
      IBT(2) = ID(2) 
      IBT(3) = ID(3) 
      IBT(4) = ID(4) + IQ 
      BT(3) = BD(3) + HALF*DBLE(IQ) 
      RETURN  
      END FUNCTION ITLS 
