!
!     ------------------------------------------------------------------
!     A N A L Y S 2
!     ------------------------------------------------------------------
!
      SUBROUTINE ANALY2(MCFG, KCFG, LIST) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      use parameters_biotr_C
      USE ELT_C 
      !USE DEBUG_C 
      USE DBG_C 
      USE INOUT_C 
      USE FO_C 
      USE NOR_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  19:43:33  11/16/01  
!...Switches:                     
!
!        This routine analyzes the format of the configuration input
!        data, for two sets, not necessarily orthogonal.
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE find_I 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: MCFG 
      INTEGER , INTENT(OUT) :: KCFG 
      CHARACTER , INTENT(INOUT) :: LIST(*)*3 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(8) :: IEL 
      INTEGER , DIMENSION(2) :: NORB, ICFG 
      INTEGER , POINTER, DIMENSION(:,:) :: QAFTER 
      INTEGER , ALLOCATABLE, TARGET, DIMENSION(:,:) :: AFTER 
      INTEGER :: IERR, I, J, NCFG, ISTATE, IRD, ICLOSD, NCLO, IORB, K, I1, I2, &
         J1, J2, IORIG, LAST, II, IORD, IL, LASTEL, NOR11 
      REAL(DOUBLE) :: STAT 
      CHARACTER :: LINE*72 
      CHARACTER , DIMENSION(NWD,2) :: OF*3 
      CHARACTER , DIMENSION(NWD) :: ELC*3 
      CHARACTER , DIMENSION(8) :: EL*3 
      CHARACTER , DIMENSION(NWD,2) :: ELBLC*3 
      CHARACTER , DIMENSION(2) :: LABEL*7 
      CHARACTER :: ANS*6 
!-----------------------------------------------
!      POINTER (QAFTER,AFTER(3*NWD,3*NWD))
 
      DATA LABEL/ 'Initial', 'Final  '/  
!
    1 FORMAT(A72) 
    6 FORMAT(A7,' State : '/,'-------------') 
    4 FORMAT(/,'   closed shells: '/(1X,18(1X,A3))) 
    5 FORMAT('   open shells: '/(1X,18(1X,A3))) 
    7 FORMAT(18(1X,A3)) 
    8 FORMAT(/,' The following biorthonormal set is built for TENSOR:'/(1X,18(&
         1X,A3))) 
      IF (IBUGM /= 0) WRITE (6, *) ' qafter  allocation: 9*nwd*nwd= ', 9*NWD*&
         NWD 
      ALLOCATE (AFTER(3*NWD,3*NWD), STAT=IERR) 
      if (ierr.ne.0) call mem_fail(6,3*NWD*3*NWD,'analy2::after',ierr);
      qafter=>after;
      AFTER(:3*NWD,:3*NWD) = 0 
!
      NCFG = 0 
      DO ISTATE = 1, 2 
         IRD = IUC(ISTATE) 
!
!  ---  Determine input format and number of closed shells
!
         READ (IRD, '(A30,I3,I4)') LINE, ICLOSD, IWF(ISTATE) 
         READ (IRD, '(A72)') LINE 
         NCLO = 0 
         J = 2 
   10    CONTINUE 
         IF (LINE(J:J+2) /= '   ') THEN 
            NCLO = NCLO + 1 
            ELBLC(NCLO,ISTATE) = LINE(J:J+2) 
            J = J + 4 
            IF (J < 72) GO TO 10 
         ENDIF 
         NCLOS(ISTATE) = NCLO 
!
!  --- We have the processed configuration list format: skip
!
         IF (IWF(ISTATE) > ICLOSD) READ (IRD, '(20(1X,A3))') (ELC(J),J=ICLOSD&
             + 1,IWF(ISTATE)) 
!
!  ---  Determine the number of configurations and electrons
!
         IORB = 0 
   20    CONTINUE 
         READ (IRD, 1, END=55) LINE 
         IF (LINE(1:1)/='*' .AND. LINE(2:2)/='*') THEN 
!
!  ------  A new configuration has been read; find the electrons
!
            NCFG = NCFG + 1 
            J = 2 
            I = 0 
!mrg 30      IF (LINE(J:J+2) .NE. '   ' .AND. I.LT.(5)) THEN
   30       CONTINUE 
            IF (LINE(J:J+2)/='   ' .AND. I<8) THEN 
!
!  --------- An electron has been found; is it a new one?
!
               I = I + 1 
               EL(I) = LINE(J:J+2) 
               K = 1 
   40          CONTINUE 
               IF (K <= IORB) THEN 
                  IF (OF(K,ISTATE) /= EL(I)) THEN 
                     K = K + 1 
                     IF (K > NWD) THEN 
                        WRITE (IWRITE, *) ' TOO MANY ELECTRONS: MAX=', NWD 
                        STOP  
                     ENDIF 
                     GO TO 40 
                  ELSE 
                     IEL(I) = K 
                  ENDIF 
               ELSE 
!
!  ------------  A new electron has been found; add it to the list
!
                  IORB = K 
                  OF(IORB,ISTATE) = EL(I) 
                  IEL(I) = K 
               ENDIF 
               J = J + 8 
               GO TO 30 
            ENDIF 
!
!  ------  Add data to the AFTER matrix
!
            DO I1 = 2, I 
               J1 = NWD*ISTATE + IEL(I1) 
               AFTER(J1,NWD*ISTATE+IEL(:I1-1)) = 1 
            END DO 
            READ (IRD, *) 
            GO TO 20 
         ENDIF 
   55    CONTINUE 
         NORB(ISTATE) = IORB 
         ICFG(ISTATE) = NCFG 
         NOPEN(ISTATE) = IORB 
         WRITE (IWRITE, 6) LABEL(ISTATE) 
         WRITE (IWRITE, 4) (ELBLC(I,ISTATE),I=1,NCLOS(ISTATE)) 
         WRITE (IWRITE, 5) (OF(I,ISTATE),I=1,NOPEN(ISTATE)) 
         WRITE (6, *) 
      END DO 
!
!  ---   set parameters
!
      NORBI = NORB(1) 
      NORBF = NORB(2) 
!mrg
      IWF(1) = NCLOS(1) + NOPEN(1) 
      IWF(2) = NCLOS(2) + NOPEN(2) 
!mrg
 
      MAXNFO = IWF(1) + IWF(2) 
      WRITE (6, *) ' total number of radial distributions = ', MAXNFO 
      WRITE (6, *) ' total number of initial orbitals     = ', IWF(1) 
      WRITE (6, *) ' total number of final   orbitals     = ', IWF(2) 
!
!  ---  determine the common initial/final state orbitals
!
!
      ELC(:NORBI) = OF(:NORBI,1) 
      NCOM = NORBI 
!
!  ---  add others from final state
!
      L54: DO I = 1, NORBF 
         DO J = 1, NORBI 
            IF (OF(I,2) /= OF(J,1)) CYCLE  
            CYCLE  L54 
         END DO 
         NCOM = NCOM + 1 
         IF (NCOM > NWD) STOP ' Too many common electrons: MAX=(60)' 
         ELC(NCOM) = OF(I,2) 
      END DO L54 
!
      WRITE (ISCW, '(2/)') 
!
! --- Transfer electrons to common orthogonal set
!
      DO ISTATE = 1, 2 
         IORIG = NWD*ISTATE 
         LAST = NORB(ISTATE) 
         DO I = 1, NCOM 
!
! --- Find electron and transfer AFTER information
!
            J = 1 
  202       CONTINUE 
            IF (J <= LAST) THEN 
               IF (ELC(I) /= OF(J,ISTATE)) THEN 
                  J = J + 1 
                  GO TO 202 
               ELSE 
                  II = IORIG + J 
                  DO K = 1, IORIG + LAST 
                     IF (AFTER(I,K) == 0) AFTER(I,K) = AFTER(II,K) 
                     IF (AFTER(K,I) == 0) AFTER(K,I) = AFTER(K,II) 
                     AFTER(II,K) = 2 
                     AFTER(K,II) = 2 
                  END DO 
                  NORB(ISTATE) = NORB(ISTATE) - 1 
               ENDIF 
            ELSE 
               WRITE (ISCW, *) ' Common electron ', ELC(I), ' not found in ', &
                  LABEL(ISTATE), ' state' 
            ENDIF 
         END DO 
      END DO 
!
!  ---  Check if the ordering of the electrons is inconsistent
!
      DO I = 1, NWD*3 
         EL(1) = FIND(I,OF,ELC) 
         DO J = 1, NWD*3 
            EL(2) = FIND(J,OF,ELC) 
            IF (AFTER(I,J)/=1 .OR. AFTER(J,I)/=1) CYCLE  
            WRITE (ISCW, *) ' The order of ', EL(1), ' and ', EL(2), &
               ' is inconsistent' 
            STOP  
         END DO 
      END DO 
!
!  ---  Reorder the electrons to satisfy the after relations found
!         in the different configurations
!
      IORD = 1 
   70 CONTINUE 
      IF (IORD <= NCOM) THEN 
!
!  ------  Search for a row with no 1's in the NCOM rows
!
         L71: DO I = 1, NCOM 
            DO J = 1, NWD*2 + NORBF 
               IF (AFTER(I,J) /= 1) CYCLE  
               CYCLE  L71 
            END DO 
!
!  ---------  The current row contains all 0's or 2's
!
            IF (AFTER(I,I) == 2) CYCLE  L71 
!
!  ------------  We have the next electron; delete the corresponding
!                  rows and columns from the AFTER matrix
!
            LIST(IORD) = ELC(I) 
            IORD = IORD + 1 
            AFTER(I,:NWD*2+NORBF) = 2 
            AFTER(:NWD*2+NORBF,I) = 2 
            GO TO 70 
         END DO L71 
      ENDIF 
      IF (IORD /= NCOM + 1) THEN 
!
!        SEARCH FOR THE ELECTRON NOT INCLUDED
!
         DO I = 1, NCOM 
            IF (AFTER(I,I) == 2) CYCLE  
            DO J = NWD + 1, NWD*2 + NORBF 
               IF (AFTER(I,J) /= 1) CYCLE  
               WRITE (ISCW, *) ELC(I), ' cannot be included in common set' 
               IL = 1 
               IF (J > NWD*2) IL = 2 
               WRITE (ISCW, *) ' Occurs AFTER ', FIND(J,OF,ELC), ' in ', LABEL(&
                  IL), ' state' 
               STOP  
            END DO 
         END DO 
      ENDIF 
!
!  ---  ORDER THE REMAINING ELECTRONS FOR THE INITIAL AND FINAL STATE
!
      LAST = NCOM 
 
      LASTEL = NORBI 
      DO ISTATE = 1, 2 
         LAST = LAST + NORB(ISTATE) 
  304    CONTINUE 
         IF (IORD <= LAST) THEN 
            IORIG = NWD*ISTATE 
            L301: DO I = IORIG + 1, IORIG + LASTEL 
               DO J = 1, IORIG + LASTEL 
                  IF (AFTER(I,J) /= 1) CYCLE  
                  CYCLE  L301 
               END DO 
!
!           The current row contains no 1's
!
               IF (AFTER(I,I) == 2) CYCLE  L301 
!
!               We have the next electron
!
               IF (IORD > NWD) THEN 
                  WRITE (IWRITE, *) ' Too many electrons: MAX=', NWD 
                  STOP  
               ENDIF 
               LIST(IORD) = OF(I-IORIG,ISTATE) 
               IORD = IORD + 1 
               AFTER(I,:IORIG+LASTEL) = 2 
               AFTER(:IORIG+LASTEL,I) = 2 
               GO TO 304 
            END DO L301 
         ENDIF 
         LASTEL = NORBF 
      END DO 
!
      NORBI = NORB(1) 
      NORBF = NORB(2) 
      NCLOSI = NCLOS(1) 
      NCLOSF = NCLOS(2) 
      MCFG = ICFG(1) 
      KCFG = ICFG(2) - MCFG 
      IF (NCOM > 0) WRITE (IWRITE, 8) (LIST(I),I=1,NCOM) 
      ELTENS(:NCOM) = LIST(:NCOM) 
      NOR11 = NCOM + NORBI 
      DEALLOCATE (AFTER) 
      RETURN  
      END SUBROUTINE ANALY2 
