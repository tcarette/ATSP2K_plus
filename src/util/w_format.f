*       Program to FORMAT wavefunctions
*
*       Created by C. Froese Fischer June 16, 1987
*       Vanderbilt University

      PROGRAM  W_FORMAT
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER AT*6,TT*6,EL1*3,NEW*3,NAME*24
      DIMENSION PT(220)
*
      OPEN(UNIT=2,FILE='wfn.fmt',STATUS='UNKNOWN')
      OPEN(UNIT=3,FILE='wfn.inp',STATUS='OLD',
     :   FORM='UNFORMATTED')
      IIN = 3
      IOUT=2
3     FORMAT(A6,A6,3X,A3,I6,F6.2,3(1PE18.10))
2     READ(IIN,END=5) AT,TT,EL1,M,ZT,ETI,EKI,AZI,(PT(J),J=1,M)
*      READ(IIN,END=5)(PT(J),J=1,M)
         WRITE(IOUT,3) AT,TT,EL1,M,ZT,ETI,EKI,AZI
         WRITE(IOUT,'(4(1PE18.10))')(PT(J),J=1,M)
      GO TO 2
5     ENDFILE 2
      CLOSE(UNIT=3)
      CLOSE(UNIT=2)
      END

