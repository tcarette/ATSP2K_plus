!
!     -------------------------------------------------------------
!      S A V E N O N
!     -------------------------------------------------------------
!                                                                  *
!     THIS SUBROUTINE FOR            B I O R T O G                 *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville           September 1997   *
!
      SUBROUTINE SAVENON(I, A, KL, LA, LB, LC, LD, JA, JB, IPTR) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      use non30_C
      use rdint_C
      USE CONSTS_C 
      USE EMS_C 
      USE NTRM_C 
      USE SAVECOM_C 
      USE SIGNF_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:05:53  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rmetr_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I 
      INTEGER  :: KL 
      INTEGER  :: LA 
      INTEGER  :: LB 
      INTEGER  :: LC 
      INTEGER  :: LD 
      INTEGER  :: JA 
      INTEGER  :: JB 
      INTEGER  :: IPTR 
      REAL(DOUBLE)  :: A 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: CL, CV, RMET 
!-----------------------------------------------
!      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
!     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
!     :       (QLJCLSD,LJCLSD(1))
!      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
!      POINTER(IQRL,RLINT(MAXORB,1)),(IQRV,RVINT(MAXORB,1)),
!     :       (IQOV,OVLP(MAXORB,1))
!      COMMON/RDINT/IQRL,IQRV,IQOV
      IF (I == 200) THEN 
         CL = ZERO 
         CV = ZERO 
         NTERMS = NTERMS + 1 
         RMET = RMETR(LRHO,LSIG) 
         CL = A*RMET 
         IF (VOK) CV = CL*RVINT(LB,LD) 
         CL = CL*RLINT(LB,LD) 
         CL2 = CL2 + CL 
         CV2 = CV2 + CV 
      ELSE 
         A = A*SIGNFA 
         CALL SAVELS (I, A, KL, LA, LB, LC, LD, JA, JB, IPTR) 
      ENDIF 
      RETURN  
      END SUBROUTINE SAVENON 
