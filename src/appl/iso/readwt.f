*-----------------------------------------------------------------------
*     R E A D W T
*-----------------------------------------------------------------------
*
*     If the wavefunction expansion has been obtained from an
*     MCHF calculation, the weights are J independent and are read from 
*     the file <name>.c . If the wavefunction expansion has been 
*     obtained from a CI calculation, the weights  
*     are read from the file <name>.l .NR is the ordernumber of the
*     main cfg component in the CI expansion
*     and determines where in the file the weights are read.
*     The weights are saved in an array WT(K) where K is the number of
*     the configurations in the expansion  

      SUBROUTINE  READWT(MAXORB,IMCHF,NCFG,NR)
C      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER CFGNR
      CHARACTER*80 LINE
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC0,ISC1,ISC2,ISC3,
     : IALL,JSC(3),ISCW
      COMMON/RYDBERG/RYDBRG,ZUM,ENERGY

      POINTER(QWT,WT(1))

      COMMON /NAN/ QWT

*     Read from <name>.c

*     allocate memory

      call alloc(qwt,ncfg,8)

      IF (IMCHF.EQ.1) THEN
         REWIND(UNIT=IREAD)
         DO 5 K=1,3+((MAXORB-1)/20)
            READ(IREAD,*)
5        CONTINUE
         DO 10 K=1,NCFG
            READ(IREAD,100) WT(K)
            READ(IREAD,*)
10       CONTINUE
100      FORMAT(T65,F15.12)
Cww100      FORMAT(T65,F10.7)
         RETURN
      ENDIF

*     Weights from CI calc. If IMCHF=2 read from <name>.l 
*     else if IMCHF=3 read from <name>.j

      J=0
      READ(63,200) NNCFG
200   FORMAT(T40,I7)
Cww old ci format 200   FORMAT(T38,I4)
      IF (NNCFG.NE.NCFG) THEN
         WRITE(0,*) ' Unequal numer of cfg in <name>.c and <name>.j(l)'
         STOP
      ENDIF

*     Determine how many lines the array of weights occupy
*     and the number of weights on the last line.

      IF (MOD(NNCFG,7).EQ.0) THEN
         NRLINES=NNCFG/7
         NRLL=0
      ELSE
         NRLINES=NNCFG/7+1
         NRLL=MOD(NNCFG,7)
      ENDIF

*     Find the right array of weights in the file.
*     Firstly the right J quantum number should be found.

20    READ(63,*)
      READ(63,*)
      READ(63,300) JJ,NUMBER
300   FORMAT(T8,I6,T23,I4)

      IF (IMCHF.EQ.2) THEN
         JJ=0
         MAX=0
         MIN=0 
      ENDIF
      IF (JJ.GT.MAX) THEN
         DO 30 K=1,NUMBER*(NRLINES+2)
            READ(63,*)
30       CONTINUE
         GO TO 20
      ENDIF    

      J=J+1
      N=0
      ENERGY=0.0

*     Secondly the right configuration should be found.

      DO 40 I=1,NUMBER
      READ(63,*)
      READ(63,400) CFGNR,ENERG
      IF (CFGNR.EQ.NR) THEN 
	 ENERGY=ENERG
         N=N+1 
         DO 50 K=1,NNCFG/7
                  READ(63,500) (WT(L),L=7*(K-1)+1,7*K)
                  WRITE(*,*) (WT(L),L=7*(K-1)+1,7*K)
50       CONTINUE
         IF (NRLL.NE.0) THEN
                 READ(63,500) (WT(L),L=NNCFG-NRLL+1,NNCFG)
                 write(*,*) (WT(L),L=NNCFG-NRLL+1,NNCFG)
         ENDIF
      ELSE  
         DO 60 K=1,NRLINES
            READ(63,*)
60       CONTINUE 
      ENDIF 
40    CONTINUE
400   FORMAT(I6,F16.9)
CGG400   FORMAT(I6)
500   FORMAT(7F11.8)
CGG500   FORMAT(7F15.12)
Cww500   FORMAT(7F10.7)
6     FORMAT(A,I4,A)
      IF (N.EQ.0) THEN
         IF (IMCHF.EQ.2) THEN
         WRITE(0,*) 'There is no eigenvector in the file <name>.l that' 
         WRITE(0,6) 'has cfg',NR,' as the dominating component.'
         WRITE(0,*) 'Change the number of the dominating component in'
         WRITE(0,*) 'the wanted eigenvector to the right one according'
         WRITE(0,*) 'to your identification and rerun the program.'
         STOP
         ELSE 
         WRITE(0,6) 'There is no eigenvector with 2*J=',JJ,' in the'
        WRITE(0,6) 'file <name>.j that has cfg',NR,' as the dominating'
         WRITE(0,*) 'component. Change the number of the dominating'
         WRITE(0,*) 'component in the wanted eigenvector to the right'
         WRITE(0,*) 'one according to your identification and rerun'
         WRITE(0,*) 'the program.'
         STOP
         ENDIF
      ENDIF
      IF (N.GT.1) THEN
         IF (IMCHF.EQ.2) THEN
         WRITE(0,*) 'There is more than one eigenvector in the file'
         WRITE(0,6) '<name>.l that has cfg',NR,' as the dominating'
         WRITE(0,*) 'component. Discard the unwanted eigenvectors'
         WRITE(0,*) 'and change the digit giving the number of eigen-'
         WRITE(0,*) 'vectors to the right one and rerun the program.'
         STOP
         ELSE
         WRITE(0,6) 'There is more than one eigenvector with 2*J=',JJ
         WRITE(0,6) 'in the file <name>.j that has cfg',NR,' as'
         WRITE(0,*) 'the dominating component. Discard the unwanted'
         WRITE(0,*) 'eigenvectors and change the digit giving the'
         WRITE(0,*) 'number of eigenvectors to the right one and' 
         WRITE(0,*) 'rerun the program.'
         STOP
         ENDIF
      ENDIF
      IF (JJ.EQ.MIN) THEN
         GO TO 999
      ENDIF
      GO TO 20
999   RETURN
      END
