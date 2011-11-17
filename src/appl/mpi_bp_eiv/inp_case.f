      subroutine inp_case(QREL,QMASS,ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW,
     :      REL,MASS,iscw,iread,ATOMNAME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220, NWD=60)
      LOGICAL REL, SKIP

      CHARACTER*26 STRING,ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW,ATOMNAME
      CHARACTER QREL, QMASS

1     WRITE(ISCW,'(/A,A)') ' Enter ATOM, relativistic (Y/N)',
     :               ' with mass correction (Y/N)'
      READ(IREAD,'(A)') STRING
      I = INDEX(STRING,',')
      IF (I .eq. 0)     THEN
         WRITE(ISCW,*) ' Separate with commas'
         GO TO 1
      ELSE IF (I .GT. 22) THEN
         WRITE(ISCW,*) ' ATOM name can have at most 22 characters'
         GO TO 1
      END IF

!      WRITE(ISCW,'(/A,A)') ' Enter file with initial estimates'

      QREL  = STRING(I+1:I+1)
      QMASS = STRING(I+3:I+3)
      ATOMC = STRING(1:I-1)//'.c'
      ATOML = STRING(1:I-1)//'.l'
      ATOMJ = STRING(1:I-1)//'.j'
      ATOMW = STRING(1:I-1)//'.w'
      ATOMNEW = STRING(1:I-1)//'.new'
      ATOMNAME = string(1:I-1)

      IF (QREL.EQ.'N' .OR. QREL.EQ.'n') REL = .FALSE.
      MASS = 0
      IF (QMASS.EQ.'Y' .OR. QMASS.EQ.'y') THEN
         WRITE(ISCW,'(A)') ' Gradient or Slater form? (G/S):'
         READ(IREAD,'(A1)') QMASS
         MASS = 1
         IF (QMASS.EQ.'S' .OR. QMASS.EQ.'s') MASS = 2
      END IF

      return
      end



