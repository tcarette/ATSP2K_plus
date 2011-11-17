!
!     ------------------------------------------------------------------
!       A N A L Y S E 1
!     ------------------------------------------------------------------
!
      SUBROUTINE ANALY1(IREAD, IWRITE, NCLOSD, MAXORB, N, NCFG, LIST) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:38:59  11/16/01  
!...Switches:                     
!
!        This routine analyzes the format of the configuration input
!        data and determines a consistent ordering of the electrons
!
!     IMPLICIT INTEGER (Q)
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NWD = 128 
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IREAD 
      INTEGER , INTENT(IN) :: IWRITE 
      INTEGER , INTENT(OUT) :: NCLOSD 
      INTEGER , INTENT(OUT) :: MAXORB 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(OUT) :: NCFG 
      CHARACTER , INTENT(OUT) :: LIST(NWD)*3 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , POINTER, DIMENSION(:,:) :: QAFTER 
      INTEGER , ALLOCATABLE, TARGET, DIMENSION(:,:) :: AFTER 
      INTEGER , DIMENSION(8) :: IEL 
      INTEGER :: IERR, I, J, NCI, NWI, K, I1, I2, IORD 
      REAL(DOUBLE) :: STAT 
      CHARACTER :: LINE*72 
      CHARACTER , DIMENSION(NWD) :: OF*3 
      CHARACTER , DIMENSION(8) :: EL*3 
      CHARACTER :: HEAD*30 
!-----------------------------------------------
!      POINTER(QAFTER,AFTER((NWD),1))
!
    1 FORMAT(A72) 
!
      !call alloc(qafter,nwd*nwd,4)
      ALLOCATE (AFTER(NWD,NWD), STAT=IERR) 
 
!
      AFTER = 0 
 
!
! --- Check format of configuration file
!
      REWIND IREAD 
      READ (IREAD, '(A30,I3,I4)') HEAD, NCI, NWI 
!
!  ---  Determine the number of common closed subshells
!
      READ (IREAD, '(A72)') LINE 
      NCLOSD = 0 
      J = 2 
   10 CONTINUE 
      IF (LINE(J:J+2) /= '   ') THEN 
         NCLOSD = NCLOSD + 1 
         J = J + 4 
         IF (J < 72) GO TO 10 
      ENDIF 
!
! --- if not clist format
!
   12 CONTINUE 
      IF (NWI > NCI) THEN 
         READ (IREAD, '(A72)') LINE 
         NWI = NWI - 20 
         GO TO 12 
      ENDIF 
!
!  ---  Determine the number or configurations and electrons
!
!
      MAXORB = 0 
      NCFG = N 
   20 CONTINUE 
      READ (IREAD, 1, END=55) LINE 
      IF (LINE(1:1)/='*' .AND. LINE(2:2)/='*') THEN 
!
!  ------  A new configuration has been read; find the electrons
!
         NCFG = NCFG + 1 
         J = 2 
         I = 0 
   30    CONTINUE 
         IF (LINE(J:J+2)/='   ' .AND. I<8) THEN 
!
!  --------- An electron has been found; is it a new one?
!
            I = I + 1 
            EL(I) = LINE(J:J+2) 
            K = 1 
   40       CONTINUE 
            IF (K <= MAXORB) THEN 
               IF (OF(K) /= EL(I)) THEN 
                  K = K + 1 
                  IF (K > NWD) THEN 
                     WRITE (IWRITE, *) 'SET NWD Larger: NWD=', NWD 
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
               MAXORB = K 
               OF(MAXORB) = EL(I) 
               IEL(I) = K 
            ENDIF 
            J = J + 8 
            GO TO 30 
         ENDIF 
!
!  ------  Add data to the AFTER matrix
!
         DO I1 = 2, I 
            AFTER(IEL(I1),IEL(:I1-1)) = 1 
         END DO 
         READ (IREAD, *) 
         IF (I > 8) READ (IREAD, *) 
         GO TO 20 
      ENDIF 
!
!  ---  Check if the ordering of the electrons is inconsistent
!
   55 CONTINUE 
      DO I = 1, MAXORB 
         DO J = 1, MAXORB 
            IF (AFTER(I,J)/=1 .OR. AFTER(J,I)/=1) CYCLE  
            WRITE (IWRITE, *) ' The order of ', OF(I), ' and ', OF(J), &
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
      IF (IORD <= MAXORB) THEN 
!
!  ------  Search for a row with no 1's
!
         L71: DO I = 1, MAXORB 
            DO J = 1, MAXORB 
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
            LIST(IORD) = OF(I) 
            IORD = IORD + 1 
            AFTER(I,:MAXORB) = 2 
            AFTER(:MAXORB,I) = 2 
            GO TO 70 
         END DO L71 
      ENDIF 
      IF (ALLOCATED(AFTER)) THEN 
!
!      call dalloc(qafter,nwd*nwd)
         DEALLOCATE (AFTER) 
      ENDIF 
      RETURN  
      END SUBROUTINE ANALY1 
