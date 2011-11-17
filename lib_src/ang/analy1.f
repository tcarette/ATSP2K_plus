*      
*     Routines for MCHF_LIB_ANG
*
*                   C O P Y R I G H T -- 1994
*
*     Computer Physics Communication, Vol. 64, 399-405 (1991)
*
*     NOTE:  This file contains routines for programs using one atomic 
*            state and two atomic states as in transitions. The COMMONS
*            are different, thus there will be warning messages if this
*	     file compiled without first using fsplit.
*     ------------------------------------------------------------------
*	A N A L Y S E 1
*     ------------------------------------------------------------------
*
      SUBROUTINE ANALY1(IREAD,IWRITE,NCLOSD,MAXORB,N,NCFG,NOCCSH,LIST,
     :                   NCD)
*
*        This routine analyzes the format of the configuration input
*        data and determines a consistent ordering of the electrons
*
      PARAMETER (NWD=30)
      INTEGER NOCCSH(NCD),AFTER(NWD,NWD),IEL(5)
      CHARACTER LIST(NWD)*3, LINE*72, OF(NWD)*3, EL(5)*3
*
  1   FORMAT(A72)
*
      DO 2 I = 1,(NWD)
         DO 3 J = 1,(NWD)
            AFTER(I,J) = 0
  3      CONTINUE
  2   CONTINUE
*
*  ---  Determine the number of common closed subshells
*
      READ(IREAD,'(/A72)' ) LINE
      NCLOSD = 0
      J = 2
 10   IF (LINE(J:J+2) .NE. '   ' ) THEN
         NCLOSD = NCLOSD + 1
         J = J+4
         IF (J .LT. 72) GO TO 10
      END IF
*
*  ---  Determine the number or configurations and electrons
*
      MAXORB = 0
      NCFG = N
 20   READ(IREAD,1,END=55) LINE
      IF (LINE(1:1) .NE. '*'  .AND. LINE(2:2) .NE. '*' ) THEN
*
*  ------  A new configuration has been read; find the electrons
*
         NCFG = NCFG + 1
         IF (NCFG .GT. NCD )
     :	    WRITE(IWRITE,'(A,I5)') ' TOO MANY CONFIGURATIONS: MAX=',NCD
         J = 2
         I = 0
 30      IF (LINE(J:J+2) .NE. '   ' .AND. I.LT.(5)) THEN
*
*  --------- An electron has been found; is it a new one?
*
            I = I+1
            EL(I) = LINE(J:J+2)
            K = 1
 40         IF (K .LE. MAXORB) THEN
               IF ( OF(K) .NE. EL(I) ) THEN
                  K = K+1
                  IF (K .GT. NWD) THEN
		     WRITE(IWRITE,*) ' TOO MANY ELECTRONS: MAX=',NWD
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
         NOCCSH(NCFG) = I
*
*  ------  Add data to the AFTER matrix
*
         DO 50 I1 = 2,I
            DO 51 I2 = 1,I1-1
               AFTER(IEL(I1),IEL(I2)) = 1
 51         CONTINUE
 50      CONTINUE
         READ(IREAD,*)
         IF (I .GT. 5) READ(IREAD,*)
         GO TO 20
      END IF
*
*  ---  Check if the ordering of the electrons is inconsistent
*
 55   DO 60 I = 1,MAXORB
         DO 61 J = 1,MAXORB
            IF (AFTER(I,J) .EQ. 1 .AND. AFTER(J,I) .EQ. 1) THEN
                WRITE(IWRITE,*) ' The order of ',OF(I),' and ',
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
      RETURN
      END
