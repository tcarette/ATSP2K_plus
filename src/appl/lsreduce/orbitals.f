
*     ------------------------------------------------------------------
*       O R B I T A L S
*     ------------------------------------------------------------------
*
      SUBROUTINE ORBITALS(maxorb,el,qiajcmp,qljcomp,
     :                    qnjcomp,qiajcld,qljclsd,nb)
*
*       Process the lists of closed shells, orbitals and set parameters
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
ctc      PARAMETER (NWD=60,NWCD=20,LSDIM=30000)
      PARAMETER (NWD=70,NWCD=20,LSDIM=1000000)
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      CHARACTER EL(NWD)*3, string*72, buffer*8
      CHARACTER*1 JAJCLD(3,NWCD),JAJCMP(3,NWD)
*
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,JSC0,JSC(4)
      COMMON /CLOSED/B1ELC(4),NCLOSD,IBK
*
    3 FORMAT(18(1X,A3))
    4 FORMAT(3A1)
    8 FORMAT(A3)
   23 FORMAT(/10H THERE ARE,I3,21H ORBITALS AS FOLLOWS://
     : 5X,21(1X,A3):/5X,21(1X,A3))
   28 FORMAT(/10H THERE ARE ,I3,31H CLOSED SUBSHELLS COMMON TO ALL ,
     :  27H CONFIGURATIONS AS FOLLOWS: //
     :  5X, 21(1X,A3))
*
* --- ALLOCATE MEMORY: NWFD = MAXORB
*
      NWFD = MAXORB
      call alloc(qiajcmp,nwfd,4)
      call alloc(qljcomp,nwfd,4)
      call alloc(qnjcomp,nwfd,4)
      call alloc(qiajcld,nwcd,4)
      call alloc(qljclsd,nwcd,4)
*
* ---  We have the EL list from analyz_blk
*
      DO 30 I = 1,MAXORB
         READ(EL(I),8) IAJCMP(I)
         READ(EL(I),4) (JAJCMP(J,I),J=1,3)
30    CONTINUE
      WRITE(IWRITE,23) MAXORB,(IAJCMP(I),I=1,MAXORB)
      DO 60 I=1,MAXORB
      IF (JAJCMP(1,I) .EQ. ' ') THEN
         JAJCMP(1,I) = JAJCMP(2,I)
         JAJCMP(2,I) = JAJCMP(3,I)
         JAJCMP(3,I) = ' '
      ENDIF
      LJCOMP(I) = LVAL(JAJCMP(2,I))
      NJCOMP(I) = ICHAR(JAJCMP(1,I)) - ICHAR('1') + 1
   60 CONTINUE
*
* --- We know the number of closed shells but not their properties
*
      REWIND(IREAD)
      READ(IREAD,'(A72)')
*
* --- READ IN THE COMMON SET OF CLOSED SUBSHELLS
*
      READ(IREAD,3) (EL(I),I=1,NCLOSD)
      DO 70 I=1,NCLOSD
         READ(EL(I),8) IAJCLD(I)
         READ(EL(I),4) (JAJCLD(J,I),J=1,3)
         J = 3
         IF (JAJCLD(1,I) .NE. ' ') J = 2
         LJCLSD(I) = LVAL(JAJCLD(J,I))
 70   CONTINUE
      WRITE(IWRITE,28) NCLOSD,(IAJCLD(I),I=1,NCLOSD)
*
*  ---  SEPARATE THE ELECTRON LABEL CHARACTERS AND LEFT JUSTIFY
*
       DO 10 I = 1,MAXORB
          WRITE(BUFFER,'(A3)') IAJCMP(I)
          READ(BUFFER,'(3A1)') (JAJCMP(J,I),J=1,3)
          IF (JAJCMP(1,I) .EQ. ' ') THEN
             JAJCMP(1,I) = JAJCMP(2,I)
             JAJCMP(2,I) = JAJCMP(3,I)
             JAJCMP(3,I) = ' '
          END IF
 10    CONTINUE
*     .. write out intial data about the problem
      write(iout) nclosd, maxorb, nb, lsdim
      if (nclosd .gt. 0) then
         write(string,'(24A3)') (iajcld(i),i=1,nclosd)
      end if
      write(iout) string
      do i = 1, maxorb, 24
         m = min(maxorb,i+23)
         write(string,'(24A3)') (iajcmp(j),j=i,m)
         write(iout) string
       end do

*     ..write out intial data about the problem
2     FORMAT(20(1x,A3))
      NWF = MAXORB + NCLOSD
        write(25,'(i4,a15)') nclosd, 'Closed Shells:'
        WRITE(25,3) (iajcld(j),j=1,nclosd)
        write(25,'(i4,a16)') nwf-nclosd, 'Other Orbitals:'
        WRITE(25,3) (iajcmp(j),j=1,maxorb)
*
      RETURN
      END
