!
!     ------------------------------------------------------------------
!       C F G N 1
!     ------------------------------------------------------------------
!
      SUBROUTINE CFGN1() 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      use inform_C
      use non30_C
      use ndims_C
      use ovrlap_C
      use closed_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:47:32  11/17/01  
!...Switches:                     
!      PARAMETER (NWD=128,NWCD=20)
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE cfgo1_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NWD = 128 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(2) :: IEL, IBUFF 
      INTEGER :: IASTER, IBLNK, NWF, I, J, M1, N, JJ, I1, J1, M, IJ 
      CHARACTER , DIMENSION(3,NWD) :: JAJCMP 
      CHARACTER :: BUFFER*8 
!-----------------------------------------------
!
!       Read the configurations for a state and determine the
!       non-orthogonal orbitals
!
      DATA IASTER, IBLNK/ 4H*   , 4H    /  
!
!      CALL CFGO1 (NCFG,MAXORB,NCLOSD,QIAJCMP,QLJCOMP,QNJCOMP,
!     :                 QNOC,QNELCSH,QNOCORB,QJ1,QIAJCLD,QLJCLSD)
      CALL CFGO1(NCLOSD); 
      NWF = MAXORB 
 
      print*,nclosd
!
!     Now deal with the non-orthogonality
!
      IF (NWF > 1) THEN 
         WRITE (6, *) ' alloc, qiorth: nwf*(nwf-1)/2 = ', NWF*(NWF - 1)/2 
!        call alloc(qiorth,NWF*(NWF-1)/2,8)
         ALLOCATE (IORTH(NWF*(NWF-1)/2)) 
      ENDIF 
!
!  ---  SEPARATE THE ELECTRON LABEL CHARACTERS AND LEFT JUSTIFY
!
      DO I = 1, MAXORB 
         WRITE (BUFFER, '(A3)') IAJCMP(I) 
         READ (BUFFER, '(3A1)') (JAJCMP(J,I),J=1,3) 
         IF (JAJCMP(1,I) /= ' ') CYCLE  
         JAJCMP(1,I) = JAJCMP(2,I) 
         JAJCMP(2,I) = JAJCMP(3,I) 
         JAJCMP(3,I) = ' ' 
      END DO 
!
!  ---  INITIALIZE THE ORTHOGONALITY ARRAY
!
      M1 = (MAXORB*(MAXORB - 1))/2 
      IORTH(:M1) = -1 
!
!  ---  SET ORBITALS IN THE SAME CONFIGURATION TO BE ORTHOGONAL
!
      WRITE (6, *) ' NCFG = ', NCFG 
      DO I = 1, NCFG 
         N = NOCCSH(I) 
         DO J = 1, N - 1 
            DO JJ = J + 1, N 
               I1 = NOCORB(J,I) 
               J1 = NOCORB(JJ,I) 
               IF (J1 > I1) THEN 
                  M = I1 
                  I1 = J1 
                  J1 = M 
               ENDIF 
               IORTH(J1+((I1-1)*(I1-2))/2) = 0 
            END DO 
         END DO 
      END DO 
!
! --- DETERMINE THE NON-ORTHOGONAL ORBITALS
!
      NORTH = 0 
      DO J = 1, MAXORB - 1 
         DO I = J + 1, MAXORB 
            IJ = J + ((I - 1)*(I - 2))/2 
            IF (.NOT.(JAJCMP(2,I)==JAJCMP(2,J) .AND. JAJCMP(3,I)/=' ' .AND. &
               JAJCMP(3,J)/=' ' .AND. JAJCMP(3,I)/=JAJCMP(3,J) .AND. IORTH(IJ)&
               /=0)) CYCLE  
            NORTH = NORTH + 1 
            IORTH(IJ) = 1 
         END DO 
      END DO 
!
      READ (IREAD, *, END=90) 
   79 CONTINUE 
      READ (IREAD, '(2(1X,A3))', END=90) IBUFF(1), IBUFF(2) 
      IF (IBUFF(1)/=IASTER .AND. IBUFF(1)/=IBLNK) THEN 
         L80: DO I = 1, 2 
            DO J = 1, MAXORB 
               IF (IBUFF(I) /= IAJCMP(J)) CYCLE  
               IEL(I) = J 
               CYCLE  L80 
            END DO 
            WRITE (*, '(A,A3,A)') ' ELECTRON ', IBUFF(I), ' NOT FOUND' 
            STOP  
         END DO L80 
         IF (IEL(1) > IEL(2)) THEN 
            I = IEL(1) 
            IEL(1) = IEL(2) 
            IEL(2) = I 
         ENDIF 
         IJ = IEL(1) + ((IEL(2)-1)*(IEL(2)-2))/2 
         IF (IORTH(IJ) == 1) NORTH = NORTH - 1 
         IORTH(IJ) = 0 
         WRITE (IWRITE, '(1X,A3,A,A3)') IBUFF(1), ' is orthogonal to ', IBUFF(2&
            ) 
         GO TO 79 
      ENDIF 
   90 CONTINUE 
      RETURN  
      END SUBROUTINE CFGN1 
