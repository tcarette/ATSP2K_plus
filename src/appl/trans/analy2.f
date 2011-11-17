*
*     ------------------------------------------------------------------
*	A N A L Y S 2
*     ------------------------------------------------------------------
*
      SUBROUTINE analy2(MCFG,KCFG,LIST,LORTH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
*        This routine analyzes the format of the configuration input
*        data, for two sets, not necessarily orthogonal.
*
      PARAMETER (NWD=80)
      INTEGER AFTER,IEL(8),NORB(2),NCLOS(2),ICFG(2)
      CHARACTER*3 LIST(*), LINE*72, OF(NWD,2), ELC(NWD), EL(8), FIND
      CHARACTER*7 LABEL(2)
      CHARACTER*6 ANS
      LOGICAL LORTH
      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2)
      COMMON /NOR/NCOM,NCLOSI,NCLOSF,NORBI,NORBF,IWAR
Cmrg
      COMMON /FO/iclsd(2),iwfn(2)
Cmrg
      POINTER (QAFTER,AFTER(3*NWD,3*NWD))
      DATA LABEL/'Initial','Final  '/
*
  1   FORMAT(A72)
    4 FORMAT(/10H THERE ARE,I3,' INITIAL STATE ORBITALS AS FOLLOWS: '/
     :      (1X,18(1X,A3)))
    5 FORMAT(/10H THERE ARE,I3,' FINAL STATE ORBITALS AS FOLLOWS: '/
     :      (1X,18(1X,A3)))
    6 FORMAT(' List common orbitals, terminating with a blank orbital.'/
     :       ' Upper and lower case characters must match.'/
     :       ' Fixed format (18(1X,A3)) as inicated below:'/
     :       ' AAA AAA AAA AAA AAA AAA AAA .... etc (up to 18/line)')
    7 FORMAT(18(1X,A3))
    8 FORMAT(/10H THERE ARE,I3,' COMMON ORBITALS AS FOLLOWS: '/
     :      (1X,18(1X,A3)))
*
      if(ibug1.ne.0) print*,' qafter  allocation: 9*nwd*nwd= ',9*nwd*nwd
      call alloc(qafter,9*NWD*NWD,4)
      DO 2 I = 1,(3*NWD)
         DO 3 J = 1,(3*NWD)
            AFTER(I,J) = 0
  3      CONTINUE
  2   CONTINUE
*
      NCFG = 0
      DO 100 ISTATE = 1,2
        IRD = iuc(istate)
*
*  ---  Determine input format and number of closed shells
*
      print*,' ird = ',ird
      READ(IRD,'(A30,I3,I4)') line,iclosd,iwf
Cmrg
      iclsd(istate) = iclosd
      iwfn(istate) = iwf
Cmrg
      READ(IRD,'(A72)' ) LINE
      NCLO = 0
      J = 2
 10   IF (LINE(J:J+2) .NE. '   ' ) THEN
         NCLO = NCLO + 1
         J = J+4
         IF (J .LT. 72) GO TO 10
      END IF
      NCLOS(ISTATE) = NCLO
*
*  --- We have the processed configuration list format: skip
*
      if (iwf .gt. iclosd)
     :    READ(IRD,'(20(1X,A3))') (ELC(j),j=iclosd+1,iwf)
*
*  ---  Determine the number of configurations and electrons
*
      IORB = 0
 20   READ(IRD,1,END=55) LINE
      IF (LINE(1:1) .NE. '*'  .AND. LINE(2:2) .NE. '*' ) THEN
*
*  ------  A new configuration has been read; find the electrons
*
         NCFG = NCFG + 1
         J = 2
         I = 0
Cmrg 30      IF (LINE(J:J+2) .NE. '   ' .AND. I.LT.(5)) THEN
 30      IF (LINE(J:J+2) .NE. '   ' .AND. I.LT.(8)) THEN
*
*  --------- An electron has been found; is it a new one?
*
            I = I+1
            EL(I) = LINE(J:J+2)
            K = 1
 40         IF (K .LE. IORB) THEN
               IF ( OF(K,ISTATE) .NE. EL(I) ) THEN
                  K = K+1
                  IF (K .GT. (NWD)) THEN
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
               IORB = K
               OF(IORB,ISTATE) = EL(I)
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
               J1 = (NWD)*ISTATE + IEL(I1)
               J2 = (NWD)*ISTATE + IEL(I2)
               AFTER(J1,J2) = 1
 51         CONTINUE
 50      CONTINUE
         READ(IRD,*)
         GO TO 20
      END IF
 55   NORB(ISTATE) = IORB
      ICFG(ISTATE) = NCFG
  100 CONTINUE
*
*  ---   set parameters
*
      NORBI = NORB(1)
      NORBF = NORB(2)
*
*  ---  determine the common initial/final state orbitals
*
      WRITE(ISCW,4) NORBI,(OF(I,1),I=1,NORBI)
      WRITE(ISCW,5) NORBF,(OF(I,2),I=1,NORBF)
      ANS = 'Y'
      IF (.NOT. LORTH) THEN
         WRITE(ISCW,'(/A)')
     :     ' Initial & final state orbitals an orthonormal set ? (Y/N) '
         READ(IREAD,'(A1)') ANS
      END IF
      IF ( ANS .EQ. 'Y' .OR. ANS .EQ. 'y') THEN
         DO 53 I = 1,NORBI
            ELC(I) = OF(I,1)
   53    CONTINUE
         NCOM = NORBI
*
*  ---  add others from final state
*
         DO 54 I = 1,NORBF
            DO 56 J = 1,NORBI
               IF (OF(I,2) .EQ. OF(J,1)) GO TO 54
   56       CONTINUE
            NCOM = NCOM + 1
            IF (NCOM .GT. (NWD))
     :         STOP ' Too many common electrons: MAX=(60)'
            ELC(NCOM) = OF(I,2)
   54    CONTINUE
      ELSE
         WRITE(ISCW,6)
         READ(IREAD,7) (ELC(I),I=1,18)
         IF (ELC(18) .NE. '   ') READ(IREAD,7) (ELC(I),I=19,(NWD))
         NCOM = 0
   52    IF (ELC(NCOM+1) .NE. '   ') THEN
            NCOM = NCOM + 1
            IF (NCOM .LT. (NWD)) GO TO 52
         END IF
      END IF
      WRITE(ISCW,'(//)')
*
* --- Transfer electrons to common orthogonal set
*
      DO 200 ISTATE = 1,2
         IORIG = (NWD)*ISTATE
         LAST = NORB(ISTATE)
         DO 201 I=1,NCOM
*
* --- Find electron and transfer AFTER information
*
         J = 1
  202    IF (J .LE. LAST) THEN
            IF (ELC(I) .NE. OF(J,ISTATE) ) THEN
               J = J+1
               GO TO 202
              ELSE
               II = IORIG + J
               DO 210 K = 1, IORIG+LAST
                  IF (AFTER(I,K) .EQ. 0) AFTER(I,K) = AFTER(II,K)
                  IF (AFTER(K,I) .EQ. 0) AFTER(K,I) = AFTER(K,II)
                  AFTER(II,K) = 2
                  AFTER(K,II) = 2
  210          CONTINUE
               NORB(ISTATE) = NORB(ISTATE) - 1
            END IF
           ELSE
            WRITE(ISCW,*) ' Common electron ',ELC(I),' not found in ',
     :             LABEL(ISTATE),' state'
         END IF
  201    CONTINUE
  200 CONTINUE
*
*  ---  Check if the ordering of the electrons is inconsistent
*
      DO 60 I = 1,(NWD)*3
         EL(1) = FIND(I,OF,ELC)
         DO 61 J = 1,(NWD)*3
            EL(2) = FIND(J,OF,ELC)
            IF (AFTER(I,J) .EQ. 1 .AND. AFTER(J,I) .EQ. 1) THEN
                WRITE(ISCW,*) ' The order of ',EL(1),' and ',
     :                EL(2),' is inconsistent'
                STOP
            END IF
 61      CONTINUE
 60   CONTINUE
*
*  ---  Reorder the electrons to satisfy the after relations found
*         in the different configurations
*
      IORD = 1
 70   IF (IORD .LE. NCOM ) THEN
*
*  ------  Search for a row with no 1's in the NCOM rows
*
         DO 71 I = 1,NCOM
            DO 72 J = 1,(NWD)*2+NORBF
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
               LIST(IORD) = ELC(I)
               IORD = IORD+1
               DO 74 J = 1,(NWD)*2+NORBF
                  AFTER(I,J) = 2
                  AFTER(J,I) = 2
   74          CONTINUE
               GO TO 70
            END IF
 71      CONTINUE
      END IF
      IF (IORD .NE. NCOM+1) THEN
*
*        SEARCH FOR THE ELECTRON NOT INCLUDED
*
      DO 73 I = 1,NCOM
         IF (AFTER(I,I) .NE. 2) THEN
         DO 75 J = (NWD)+1,(NWD)*2+NORBF
            IF (AFTER(I,J) .EQ. 1) THEN
               WRITE(ISCW,*)
     :           ELC(I),' cannot be included in common set'
               IL = 1
               IF ( J .GT. (NWD)*2 ) IL = 2
               WRITE(ISCW,*) ' Occurs AFTER ',FIND(J,OF,ELC),' in ',
     :                     LABEL(IL),' state'
               STOP
            END IF
   75    CONTINUE
         END IF
   73 CONTINUE
      END IF
*
*  ---  ORDER THE REMAINING ELECTRONS FOR THE INITIAL AND FINAL STATE
*
      LAST = NCOM

      LASTEL = NORBI
      DO 300 ISTATE = 1,2
         LAST = LAST + NORB(ISTATE)
  304    IF (IORD .LE. LAST) THEN
         IORIG = (NWD)*ISTATE
         DO 301 I = IORIG+1, IORIG+LASTEL
            DO 302 J = 1,IORIG+LASTEL
               IF (AFTER(I,J) .EQ. 1) GO TO 301
  302       CONTINUE
*
*           The current row contains no 1's
*
            IF (AFTER(I,I) .NE. 2) THEN
*
*               We have the next electron
*
                IF (IORD.GT.(NWD)) THEN
		  WRITE(IWRITE,*) ' Too many electrons: MAX=',NWD
		  STOP
		END IF
                LIST(IORD) = OF(I-IORIG,ISTATE)
                IORD = IORD+1
                DO 303 J = 1,IORIG+LASTEL
                   AFTER(I,J) = 2
                   AFTER(J,I) = 2
  303           CONTINUE
                GO TO 304
             END IF
  301    CONTINUE
         END IF
         LASTEL = NORBF
  300 CONTINUE
*
      NORBI = NORB(1)
      NORBF = NORB(2)
      NCLOSI = NCLOS(1)
      NCLOSF = NCLOS(2)
      MCFG = ICFG(1)
      KCFG = ICFG(2) - MCFG
      IF (NCOM .GT. 0) WRITE(IWRITE,8) NCOM,(LIST(I),I=1,NCOM)
      WRITE(IWRITE,4) NORBI,(LIST(I),I=NCOM+1,NCOM+NORBI)
      NOR11 = NCOM + NORBI
      WRITE(IWRITE,5) NORBF,(LIST(I),I=NOR11+1,NOR11+NORBF)
      call dalloc(qafter,9*nwd*nwd)
      RETURN
      END
