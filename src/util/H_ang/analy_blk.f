*
*     ------------------------------------------------------------------
*       A N A L Y S E_B L K
*     ------------------------------------------------------------------
*
      SUBROUTINE ANALY_BLK(NCLOSD,MAXORB,NB,NBsize,list)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
*        This routine analyzes the format of the configuration input
*        data for the different blocks, determines the number of blocks,
*        the size of each block, and a consistent ordering of the electrons
*
*        NB             - number of blocks
*        NBsize(1:NB)   - size of each block
*        Nclosd         - number of closed shells
*        Maxorb         - number of orther orbitals
*        List(1:maxorb) - list of orbitals 
*
      PARAMETER (NWD=70, NBD=20)
      CHARACTER LIST(NWD)*3, LINE*72, OF(NWD)*3, EL(8)*3, LINEP*72
      INTEGER NBSIZE(NBD)
      POINTER(QAFTER,AFTER((NWD),1))
      INTEGER AFTER,IEL(8)
*
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,JSC0,JSC(4)
*
      call alloc(qafter,nwd*nwd,4)
*
      DO 2 I = 1,(NWD)
         DO 3 J = 1,(NWD)
            AFTER(I,J) = 0
  3      CONTINUE
  2   CONTINUE
*
* --- Skip header
*
      REWIND IREAD
      READ(IREAD,*)
      READ(IREAD,'(A72)') Line
      LINEP = LINE
*
*  ---  Determine the number of common closed subshells
*
      NCLOSD = 0
      J = 2
 10   IF (LINE(J:J+2) .NE. '   ' ) THEN
         NCLOSD = NCLOSD + 1
         J = J+4
         IF (J .LT. 72) GO TO 10
      END IF
*
*  ---  Determine the number of blocks, configurations in each block
*       and number of electrons
*
*
      NB = 0
      MAXORB = 0
 15   NCFG=0
 20   READ(IREAD,'(A72)',END=55) LINE
      IF (LINE(1:1) .NE. '*'  ) THEN
*
*  ------  A new configuration has been read; find the electrons
*
         NCFG = NCFG + 1
         J = 2
         I = 0
 30      IF (LINE(J:J+2) .NE. '   ' .AND. I.LT.(8)) THEN
*
*  --------- An electron has been found; is it a new one?
*
            I = I+1
            EL(I) = LINE(J:J+2)
            K = 1
 40         IF (K .LE. MAXORB) THEN
               IF ( OF(K) .NE. EL(I) ) THEN
                  K = K+1
                  IF (K .GT. (NWD)) THEN
                     WRITE(IWRITE,*) 'SET NWD Larger: NWD=',NWD
                     STOP
                  END IF
                  GO TO 40
                 ELSE
                  IEL(I) = K
               END IF
              ELSE
*
*  ------------  A new electron has been found; add it to the list
*
               MAXORB = K
               OF(MAXORB) = EL(I)
               IEL(I) = K
            END IF
            J = J+8
            GO TO 30
         END IF
*
*  ------  Add data to the AFTER matrix
*
         DO 50 I1 = 2,I
            DO 51 I2 = 1,I1-1
               AFTER(IEL(I1),IEL(I2)) = 1
 51         CONTINUE
 50      CONTINUE
         READ(IREAD,*)
         GO TO 20
      ELSE
*       .. we have reached an * indicating the end of the current block
        NB = NB + 1
        NBsize(nb) = NCFG
*       .. skip next two lines (might be end-of-file)
        READ(IREAD,*,END=55)
        READ(IREAD, '(A72)',END=55) LINE
        IF (LINE .ne. LINEP) then
          WRITE(IWRITE,'(A,I8,A)') 'Closed shells for Block',nb+1,
     :          'not the same as for first block'
          STOP
        END IF
	GO TO 15
      END IF
*
*  ---  Check if the ordering of the electrons is inconsistent
*
 55   DO 60 I = 1,MAXORB
         DO 61 J = 1,MAXORB
            IF (AFTER(I,J) .EQ. 1 .AND. AFTER(J,I) .EQ. 1) THEN
                WRITE(IWRITE,*) ' The order of ',OF(I),' and ',
     :                OF(J),' is inconsistent'
                STOP
            END IF
 61      CONTINUE
 60   CONTINUE
*
*  ---  Reorder the electrons to satisfy the after relations found
*         in the different configurations
*
      IORD = 1
 70   IF (IORD .LE. MAXORB ) THEN
*
*  ------  Search for a row with no 1's
*
         DO 71 I = 1,MAXORB
            DO 72 J = 1,MAXORB
               IF (AFTER(I,J) .EQ. 1 ) GO TO 71
 72         CONTINUE
*
*  ---------  The current row contains all 0's or 2's
*
            IF (AFTER(I,I) .NE. 2 ) THEN
*
*  ------------  We have the next electron; delete the corresponding
*                  rows and columns from the AFTER matrix
*
               LIST(IORD) = OF(I)
               IORD = IORD+1
               DO 73 J = 1,MAXORB
                  AFTER(I,J) = 2
                  AFTER(J,I) = 2
 73            CONTINUE
               GO TO 70
            END IF
 71      CONTINUE
      END IF
*
      call dalloc(qafter,nwd*nwd)
      RETURN
      END
