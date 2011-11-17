!
!     ------------------------------------------------------------------
!     G S T A T E
!     ------------------------------------------------------------------
!
      SUBROUTINE GSTATE(NFIRST, NLAST) 
      use inout_C
      use ndims_C
      use non30_C
      use nor_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  21:56:52  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE lval_I 
      USE numval_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NFIRST 
      INTEGER , INTENT(IN) :: NLAST 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NWD = 128 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(15) :: J1QN, J2QN, J3QN 
      INTEGER :: IRD, I, N, J, K, NCOM1, NOR11, JJ, M, N1 
      CHARACTER , DIMENSION(15) :: JCQN, JQNG, J3QNG 
      CHARACTER , DIMENSION(2) :: LABEL*8 
      CHARACTER :: LINE*72 
!-----------------------------------------------
!
!      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2),iut(2)
!      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
!     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
!      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
!     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
!     :       (QLJCLSD,LJCLSD(1))
!      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
!      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
!      COMMON /NOR/NCOM,NCLOSI,NCLOSF,NORBI,NORBF,IWAR
      DATA LABEL/ 'INITIAL', 'FINAL'/  
!
!      DATA DEFINING THE STATE IS READ IN AND PRINTED OUT.
!
    5 FORMAT(8(1X,A3,'(',I2,')')) 
   55 FORMAT(8(1X,A3,1X,I2,1X)) 
   36 FORMAT(15(1X,A1,A1,A1)) 
   24 FORMAT(/,/,' INITIAL STATE CONFIGURATIONS:-') 
   25 FORMAT(/,'     ',I5,'.  ',10(1X,A3,'(',I2,')')) 
   26 FORMAT(11X,10(1X,4X,I1,A1,I1)) 
   27 FORMAT(22X,9(1X,4X,A1,A1,I1)) 
   28 FORMAT(' ----------------------------  '/) 
   29 FORMAT(/,/,' FINAL STATE CONFIGURATIONS:-') 
   30 FORMAT(2X,'ELECTRON ',A3,' NOT FOUND IN THE LIST OF ELECTRONS',&
         ' FOR THE ',A8,' STATE') 
      IF (NFIRST == 1) THEN 
!        WRITE(IWRITE,24)
         IRD = IUC(1) 
      ELSE 
!        WRITE(IWRITE,29)
         IRD = IUC(2) 
      ENDIF 
      WRITE (IWRITE, 28) 
      DO I = NFIRST, NLAST 
         READ (IRD, '(A72)') LINE 
         N = 0 
         J = 2 
   60    CONTINUE 
         IF (LINE(J:J+2)/='   ' .AND. N<8) THEN 
            N = N + 1 
            J = J + 8 
            GO TO 60 
         ENDIF 
         NOCCSH(I) = N 
         READ (LINE, 55) (NOCORB(J,I),NELCSH(J,I),J=1,N) 
         K = I 
         IF (NFIRST /= 1) K = I - NFIRST + 1 
!dbg  WRITE(IWRITE,25) K,(NOCORB(J,I),NELCSH(J,I),J=1,N)
         NCOM1 = NCOM + 1 
         NOR11 = NCOM1 + NORBI 
         L61: DO J = 1, N 
            DO JJ = 1, MAXORB 
               IF (NFIRST==1 .AND. JJ>=NOR11) EXIT  
               IF (NFIRST/=1 .AND. JJ>=NCOM1 .AND. JJ<NOR11) CYCLE  
               IF (NOCORB(J,I) /= IAJCMP(JJ)) CYCLE  
               NOCORB(J,I) = JJ 
               CYCLE  L61 
            END DO 
!
!        ELECTRON NOT FOUND IN THE LIST
!
            WRITE (IWRITE, 30) NOCORB(J,I), LABEL(NFIRST) 
            STOP  
         END DO L61 
         M = 2*N - 1 
         N1 = N + 1 
!GG      READ(IRD,6)      (J3QN(J),JCQN(J),J1QN(J),J=1,M)
         READ (IRD, 36) (J3QNG(J),JCQN(J),JQNG(J),J=1,M) 
!dbg  WRITE(IWRITE,26) (J3QNG(J),JCQN(J),J1QN(J),J=1,N)
!dbg      IF(N.EQ.1) GO TO 64
!dbg  WRITE(IWRITE,27) (J3QN(J),JCQN(J),J1QN(J),J=N1,M)
!dbg   64 CONTINUE
         DO J = 1, M 
            J2QN(J) = 2*LVAL(JCQN(J)) + 1 
!GG
            J1QN(J) = NUMVAL(JQNG(J)) 
            J3QN(J) = ICHAR(J3QNG(J)) - ICHAR('1') + 1 
!GG
            J1QNRD(J,I) = (J3QN(J)*64+J2QN(J))*64 + J1QN(J) 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE GSTATE 
