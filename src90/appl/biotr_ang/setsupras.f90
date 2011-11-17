!
!     --------------------------------------------------------------
!     S E T S U P R A S
!     --------------------------------------------------------------
!
! --- Sets up the RAS information of shells.
!
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville           September 1997   *
!
      SUBROUTINE SETSUPRAS(IK,NCLO) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE RAS_C 
      use non30_C
      use inout_C
      use closed_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  00:14:13  11/17/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IK,NCLO 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LM, I, J 
!-----------------------------------------------
!      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2),iut(2)
!      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
!     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
!     :       (QLJCLSD,LJCLSD(1))
!      cOMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      LM = 0 
      NINAC(1:11,IK) = 0 
      NRAS2(1:11,IK) = 0 
      IF (NCLOSD /= 0) THEN 
         DO I = 1, NCLOSD 
            J = LJCLSD(I) + 1 
            LM = MAX0(J,LM) 
            NINAC(J,IK) = NINAC(J,IK) + 1 
         END DO 
      ENDIF 
      IF (MAXORB /= 0) THEN 
         DO I = 1, MAXORB 
            J = LJCOMP(I) + 1 
            LM = MAX0(J,LM) 
            NRAS2(J,IK) = NRAS2(J,IK) + 1 
         END DO 
      ENDIF 
      IF (LM > 11) THEN 
         WRITE (ISCW, *) ' l-value too large in RASIN (SETSUPRAS)  ' 
         STOP  
      ELSE 
         LMAX(IK) = LM 
      ENDIF 
      RETURN  
      END SUBROUTINE SETSUPRAS 
