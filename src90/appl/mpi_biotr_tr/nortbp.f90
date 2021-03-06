!
!     ------------------------------------------------------------------
!     N O R T B P
!     ------------------------------------------------------------------
!
      SUBROUTINE NORTBP(JA, JB) 
      use inout_C
      use ovrlap_C
      use ndims_C
      use non30_C
      use nor_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:25:16  11/20/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: JA 
      INTEGER , INTENT(IN) :: JB 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(2) :: ILNO, IRNO 
      INTEGER :: N1, N2, IL, IR, I, NI, J, NJ, NA, NB, NC, K, ISTO, &
                 NMU, NMUP, NNU, NNUP 
!-----------------------------------------------
!
!      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2),iut(2)
!      POINTER (QIORTH,IORTH(1))
!      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
!     : ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
!     : QIORTH
!      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
!     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
!      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
!     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
!     :       (QLJCLSD,LJCLSD(1))
!      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
!      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
!      COMMON /NOR/NCOM,NCLOSI,NCLOSF,NORBI,NORBF,IWAR
!
  101 FORMAT(/,&
         ' INCORRECT NON-ORTHOGONALITY SET UP IN THE MATRIX ELEMENT  -  (',I2,&
         '/V/',I2,')') 
!
      N1 = NOCCSH(JA) 
      N2 = NOCCSH(JB) 
      IWAR = 0 
      JMU = 0 
      JNU = 0 
      JMUP = 0 
      JNUP = 0 
      IL = 0 
      IR = 0 
!
! --- BEGIN SEARCH FOR NON-ORTHOGONAL SUBSHELLS IN THIS MATRIX ELEMENT
!
      DO I = 1, N1 
         NI = NOCORB(I,JA) 
         L2: DO J = 1, N2 
            NJ = NOCORB(J,JB) 
            IF (NI == NJ) CYCLE  L2 
            NA = MIN0(NI,NJ) 
            IF (NA <= NCOM) CYCLE  L2 
            NB = MAX0(NI,NJ) 
!
!   THE FOLLOWING CARD IS NOT THE SAME THAN IN NORTBP OF NON
!
            NC = NORBF*(NA - NCOM - 1) + NB - NCOM - NORBI 
            IF (IORTH(NC) /= 1) CYCLE  L2 
            IF (IL == 0) GO TO 4 
            IF (ILNO(IL) == I) GO TO 14 
            IF (IL == 2) GO TO 100 
    4       CONTINUE 
            IL = IL + 1 
            ILNO(IL) = I 
   14       CONTINUE 
            IF (IR == 0) GO TO 7 
            DO K = 1, IR 
               IF (IRNO(K) /= J) CYCLE  
               CYCLE  L2 
            END DO 
            IF (IR == 2) GO TO 100 
    7       CONTINUE 
            IR = IR + 1 
            IRNO(IR) = J 
         END DO L2 
      END DO 
      IF (IL /= 0) THEN 
         IF (IR /= 1) THEN 
            IF (IRNO(1) > IRNO(2)) THEN 
               ISTO = IRNO(1) 
               IRNO(1) = IRNO(2) 
               IRNO(2) = ISTO 
            ENDIF 
         ENDIF 
         JMU = ILNO(1) 
         IF (IL /= 1) JNU = ILNO(2) 
         JMUP = IRNO(1) 
         IF (IR == 1) GO TO 10 
         JNUP = IRNO(2) 
         GO TO 10 
      ENDIF 
      NOVLPS = 0 
      RETURN  
  100 CONTINUE 
      WRITE (IWRITE, 101) JA, JB 
      IWAR = 1 
      RETURN  
   10 CONTINUE 
      NMU = NOCORB(JMU,JA) 
      NMUP = NOCORB(JMUP,JB) 
      LMU = LJCOMP(NMU) 
      LMUP = LJCOMP(NMUP) 
      IF (JNU/=0 .OR. JNUP/=0) THEN 
         IF (JNU == 0) JNU = JMU 
         IF (JNUP == 0) JNUP = JMUP 
         NNU = NOCORB(JNU,JA) 
         NNUP = NOCORB(JNUP,JB) 
         LNU = LJCOMP(NNU) 
         LNUP = LJCOMP(NNUP) 
         NOVLPS = 2 
         RETURN  
      ENDIF 
      NOVLPS = 1 
      RETURN  
      END SUBROUTINE NORTBP 
