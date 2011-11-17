*-----------------------------------------------------------------------
*     R E A D W T
*-----------------------------------------------------------------------
*
*     This subroutine reads the J dependent weights of the configura-
*     tions. If the wavefunction expansion has been obtained from an
*     MCHF calculation, the weights are J independent and are read from 
*     the file <name>.c . If the wavefunction expansion has been 
*     obtained from a CI calculation, the weights are J dependent and 
*     are read from the file <name>.j .NR is the ordernumber of the
*     main cfg component in the CI expansion
*     and determines where in the file the weights are read.
*     The weights are saved in an array WT(K,J) where K is the number of
*     the configurations in the expansion and J defines the J 
*     quantumnumber.

      SUBROUTINE  READWT(MAXORB,IMCHF,NCFG,NR,JJMAX,JJMIN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER CFGNR
      CHARACTER*80 LINE

      POINTER(QWT,WT(20,1))

      COMMON /NAN/ QWT
      COMMON/INFORM/IREADC,IWRITE,IOUT,IREADJ,IREADW,ISC(6),ISCW

*     J independent weights. Read from <name>.c

*     allocate memory

      call alloc(qwt,20*ncfg,8)

      MAX=JJMAX
      MIN=JJMIN

      IF (IMCHF.EQ.1) THEN
         REWIND(UNIT=IREADC)
         DO 5 K=1,3+((MAXORB-1)/20)
            READ(IREADC,*)
5        CONTINUE
         DO 10 K=1,NCFG
            READ(IREADC,100) WT(1,K)
            READ(IREADC,*)
10       CONTINUE
100      FORMAT(T65,F11.8)
cww100      FORMAT(T65,F10.7)
         RETURN
      ENDIF

*     Weights from CI calc. If IMCHF=2 read from <name>.l 
*     else if IMCHF=3 read from <name>.j

      J=0
      READ(IREADJ,200) NNCFG
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

20    READ(IREADJ,*,end=888) 
      READ(IREADJ,*,end=888) 
      READ(IREADJ,300,end=888) JJ,NUMBER
300   FORMAT(T8,I6,T23,I4)

      IF (IMCHF.EQ.2) THEN
         MAX=0
         MIN=0 
      ENDIF
      IF (JJ.GT.MAX) THEN
         DO 30 K=1,NUMBER*(NRLINES+2)
            READ(IREADJ,*)
30       CONTINUE
         GO TO 20
      ENDIF    

      J=J+1
      N=0

*     Secondly the right configuration should be found.

      DO 40 I=1,NUMBER
      READ(IREADJ,*)
      READ(IREADJ,400) CFGNR
      IF (CFGNR.EQ.NR) THEN 
         N=N+1 
         DO 50 K=1,NNCFG/7
                  READ(IREADJ,500) (WT(J,L),L=7*(K-1)+1,7*K)
                  WRITE(*,*) (WT(J,L),L=7*(K-1)+1,7*K)
50       CONTINUE
         IF (NRLL.NE.0) THEN
                 READ(IREADJ,500) (WT(J,L),L=NNCFG-NRLL+1,NNCFG)
                 write(*,*) (WT(J,L),L=NNCFG-NRLL+1,NNCFG)
         ENDIF
      ELSE  
         DO 60 K=1,NRLINES
            READ(IREADJ,*)
60       CONTINUE 
      ENDIF 
40    CONTINUE
400   FORMAT(I6)
500   FORMAT(7F11.8)
cww500   FORMAT(7F10.7)
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
888   if(n.eq.0) then
        write(0,*) 'I have searched the .j file, and did not find'
        write(0,*) 'any Eigenvalue with ',nr,' as major component.'
        stop
      end if
999   RETURN
      END
