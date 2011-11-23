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

      SUBROUTINE  READWT(NCFG,NR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER CFGNR

*      POINTER(QWT,WT(20,1))
CAB      POINTER(QWT,WT(1))
      POINTER(QWT,WT(*))

      COMMON /NAN/ QWT
      COMMON/INFORM/IREADC,IWRITE,IOUT,IREADJ,IREADW,ISC(6),ISCW

*     allocate memory

      call alloc(qwt,20*ncfg,8)
      N=0
*     Weights from <name>.l 
      READ(IREADJ,200) NNCFG
200   FORMAT(T40,I7)
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

      READ(IREADJ,*,end=888) 
      READ(IREADJ,*,end=888) 
      READ(IREADJ,300,end=888) JJ,NUMBER
300   FORMAT(T8,I6,T23,I4)


       READ(IREADJ,*)
       READ(IREADJ,400) CFGNR
       N=N+1 
       DO 50 K=1,NNCFG/7
           READ(IREADJ,500) (WT(L),L=7*(K-1)+1,7*K)
*           WRITE(*,'3x,7(F8.5,1x)') (WT(L),L=7*(K-1)+1,7*K)
50     CONTINUE
           IF (NRLL.NE.0) THEN
               READ(IREADJ,500) (WT(L),L=NNCFG-NRLL+1,NNCFG)
*               write(*,'3x,7(F8.5,1x)') (WT(L),L=NNCFG-NRLL+1,NNCFG)
           ENDIF
400    FORMAT(I6)
500    FORMAT(7F11.8)
6      FORMAT(A,I4,A)
          
       IF (N.EQ.0) THEN
          WRITE(0,*) 'There is no eigenvector in the file <name>.l that' 
          WRITE(0,*) 'has cfg',NR,' as the dominating component.'
          WRITE(0,*) 'Change the number of the dominating component in'
          WRITE(0,*) 'the wanted eigenvector to the right one according'
          WRITE(0,*) 'to your identification and rerun the program.'
          STOP
       ENDIF
       IF (N.GT.1) THEN
           WRITE(0,*) 'There is more than one eigenvector in the file'
           WRITE(0,*) '<name>.l that has cfg',NR,' as the dominating'
           WRITE(0,*) 'component. Discard the unwanted eigenvectors'
           WRITE(0,*) 'and change the digit giving the number of eigen-'
           WRITE(0,*) 'vectors to the right one and rerun the program.'
           STOP
       ENDIF
       RETURN
          
888   if(n.eq.0) then
        write(0,*) 'I have searched the .j file, and did not find'
        write(0,*) 'any Eigenvalue with ',nr,' as major component.'
        stop
      end if
      
      END
