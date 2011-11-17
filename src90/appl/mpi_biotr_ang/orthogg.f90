!
!     ------------------------------------------------------------------
!       O R T H O G G
!     ------------------------------------------------------------------
!
      SUBROUTINE ORTHOGG(LET) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE SIGNF_C 
 
      use inform_C
      use debug_C
      use medefn_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  00:01:38  11/17/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: LET 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: N5, N6, N7, I, L1, L2, L3, N72 
!-----------------------------------------------
!      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(4),IALL,JSC(3),ISCW
!      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
!      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
!     :     J1QN2(31,3),IJFUL(16)
!
!     THIS SUBROUTINE CHECKS FOR POSSIBLE ORTHOGONALITY DUE TO
!     COUPLING DIFFERENCES OR UNEVEN PARITY
!
  102 FORMAT(' ORTHOGONALITY IN COUPLING SCHEMES OF CONFIGURATIONS') 
  103 FORMAT(' THE TWO CONFIGURATIONS HAVE DIFFERING NUMBERS OF ELECTRONS') 
  104 FORMAT(' THE TWO CONFIGURATIONS HAVE DIFFERING TOTAL PARITY') 
      N5 = 0 
      N6 = 0 
      N7 = 0 
      N5 = SUM(NOSH1(:IHSH)) 
      N6 = SUM(NOSH2(:IHSH)) 
      N7 = DOT_PRODUCT(LJ(:IHSH),NOSH1(:IHSH)-NOSH2(:IHSH)) 
!
!     CHECK ON NUMBER OF ELECTRONS
!
      IF (N5 - N6 /= 0) THEN 
         WRITE (IWRITE, 103) 
         LET = 0 
!
!     CHECK ON PARITY
!
      ELSE 
         IF (N7 - (N7/2)*2 /= 0) THEN 
            WRITE (IWRITE, 104) 
            LET = 0 
         ELSE 
            N72 = N7/2 
            SIGNFA = 1.D0 
            IF (N72 - (N72/2)*2 /= 0) SIGNFA = -SIGNFA 
!
! --- COUPLING ORTHOGONALITY TEST FOR FIRST TWO SHELLS
!
            LET = 1 
         ENDIF 
      ENDIF 
      RETURN  
      END SUBROUTINE ORTHOGG 
