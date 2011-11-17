*
*     ------------------------------------------------------------------
*	I F N M N X
*     ------------------------------------------------------------------
*
      FUNCTION IFNMNX(IVEC,NEL,IMXMN)
*
* Smallest or largest value of integer array
*
* IMXMN = 1 => Largest value
* IMXMN = 2 => Smallest
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IVEC(NEL)
*
      IVAL = IVEC(1)
      IF ( IMXMN .EQ. 1 ) THEN
        DO 100 IEL = 2, NEL
           IF(IVEC(IEL).GT. IVAL ) IVAL = IVEC(IEL)
  100   CONTINUE
      ELSE IF ( IMXMN .EQ. 2 ) THEN
        DO 200 IEL = 2, NEL
           IF(IVEC(IEL).LT. IVAL ) IVAL = IVEC(IEL)
  200   CONTINUE
      ELSE
        WRITE(6,*) ' Stop in IFNMNX '
        WRITE(6,*) ' Improper calue of IMXMN ', IMXMN
        STOP 'IFNMNX'
      END IF
*
      IFNMNX = IVAL
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' Value returned from IFNMNX', IFNMNX
      END IF
*
      RETURN
      END
