*     ------------------------------------------------------------- 
* 
*     A PROGRAM FOR ROTATING A SET OF RADIAL ORBITALS
* 
*                By P. Jonsson and M. Godefroid
* 
*                1993 
* 
*     ------------------------------------------------------------ 
* 
      PROGRAM ROTATE 
      IMPLICIT REAL*8(A-H,O-Z) 
      PARAMETER (NOD=220,NWD=60) 
      CHARACTER*1 END,PP,ORBSYM
      CHARACTER*3 EL(75),ELC(3)
      CHARACTER*6 ATOM,TERM,ALABEL(2)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
      DIMENSION IPT(3),ET(75),EK(75),PROT(NOD)
      LOGICAL PRINT
      PARAMETER (IREAD=5,IWRITE=6)
*
*
C
C     Open wfn files
C
      OPEN(UNIT=1,FILE='wfn.inp',STATUS='OLD',
     :     FORM='UNFORMATTED')
      OPEN(UNIT=2,FILE='wfn.rot',STATUS='UNKNOWN',
     :     FORM='UNFORMATTED')

      CALL INITR
C      write(0,*) ' Enter tolerance for keeping an orbital'
C      read(5,*) tol
C
      WRITE(*,*)
     : ' Enter orbitals symmetry (s,p,d ..) to be rotated'
      READ(*,'(A1)') ORBSYM
      LROT= LVAL(ORBSYM)
      write(*,*)
      WRITE(*,*) ' Original orbitals : '
      write(*,*) ' ------------------- '
      write(*,*)
C
C     Read the radial functions.
C
      NUMBERW = 0 
      I = 1
111   READ(1,END=112) ATOM,TERM,EL(I),MR,Z,ET(I),EK(I),AZ(I),
     :     (P(J,I),J=1,MR)
      IF (EL(I)(1:1) .NE. ' ') THEN
         N(I) = ICHAR(EL(I)(1:1)) - ICHAR('1') + 1
         L(I) = LVAL(EL(I)(2:2))
      ELSE
         N(I) = ICHAR(EL(I)(2:2)) - ICHAR('1') + 1
         L(I) = LVAL(EL(I)(3:3))
      ENDIF
      NUMBERW = NUMBERW + 1
      MM = MR+1
      DO 113 J = MM,NO
         P(J,I) = D0
113   CONTINUE
      MAX(I) = MR
C
C     Initialize radial grid
C
      IF (I.EQ.1) THEN
         DO 2 J=1,NO
            R(J)=EXP(RHO)/Z
            RR(J)=R(J)*R(J)
            R2(J)=SQRT(R(J))
            RHO=RHO+H
2        CONTINUE
      ENDIF
      rad = qquadr(i,i,1)
      write(*,1001) el(i),rad,qquadr(i,i,0)
      I = I + 1
      GOTO 111  
112   CONTINUE
C
C     Find first orbital of right symmetry
C
      write(*,*) ' How many do you want to skip ? '
      read(*,*) iskip
      IC = 0
      DO 114 J = 1,NUMBERW
         IF (L(J).EQ.LROT) THEN
            IFIRST=J
            IC = IC + 1
            IF (IC.EQ.iskip+1) GOTO 126
         ENDIF
114   CONTINUE
126   write(*,1002) ifirst
C
C     Rotate the orbital
C
      nfo = 1
      DO 115 J = IFIRST+1,NUMBERW
         IF (L(J).EQ.LROT) THEN
            nfo = nfo + 1
            DO 116 K = 1,NO
               P(K,IFIRST) = P(K,IFIRST) + P(K,J)
116         CONTINUE
            AZ(IFIRST) = AZ(IFIRST) +  AZ(J)
            maxrot = max0(max(ifirst),max(j))
         ENDIF
115   CONTINUE
      print*,' coucou '
      if (nfo.eq.1) then
        write(*,*) ' Only one orbital of this l-symmetry... stop...'
        call exit
      end if
      print*,' nfo = ',nfo,' el = ',el,' ifirst = ',ifirst
      WRITE(*,1004) el(ifirst),nfo-1
      write(*,*)
      WRITE(*,*) ' Rotated  orbitals : '
      write(*,*) ' ------------------- '
      write(*,*)
C
C     Normalize the rotated orbital and save it
C
      QNORM = DSQRT(QQUADR(IFIRST,IFIRST,0))
      DO 117 J = 1,NO
         P(J,IFIRST) = P(J,IFIRST)/QNORM
         PROT(J) = P(J,IFIRST)
117   CONTINUE
      AZ(IFIRST) = AZ(IFIRST)/DSQRT(DFLOAT(nfo))
      AZROT = AZ(IFIRST) 
      QNORM = DSQRT(QQUADR(IFIRST,IFIRST,0))
C
C     Read the radial functions again and replace the 
C     ifirst orbital with the rotated one.
C
      REWIND(UNIT=1)
      I = 1
      M = 1
   12 READ(1,END=13) ATOM,TERM,EL(I),MR,Z,ET(I),EK(I),AZ(I),
     :     (P(J,I),J=1,MR)
      IF (EL(I)(1:1) .NE. ' ') THEN
         N(I) = ICHAR(EL(I)(1:1)) - ICHAR('1') + 1
         L(I) = LVAL(EL(I)(2:2))
      ELSE
         N(I) = ICHAR(EL(I)(2:2)) - ICHAR('1') + 1
         L(I) = LVAL(EL(I)(3:3))
      ENDIF
      MM = MR+1
      DO 24 J = MM,NO
24    P(J,I) = D0
      MAX(I)= MR

      IF (I.EQ.IFIRST) THEN
         DO 118 J = 1,NO
            P(J,I) = PROT(J)
118      CONTINUE
         AZ(I) = AZROT
         MAX(I) = maxrot
      ENDIF
C
C  *****  ORTHOGONALIZE TO INNER FUNCTIONS
C
      IM = I - 1
      DO 6 II =1,IM
      IF (L(i) .ne. L(ii)) GO TO 6
      If( (EL(i)(1:1) .NE. ' ' .and. EL(ii)(1:1) .NE. ' ')
     :  .and. (EL(i)(3:3) .NE. EL(ii)(3:3))) Go to 6
      PN = QQUADR(I,II,0)
      IF (ABS(ABS(PN) - 1.d0) .lt. 1.D-08) GO TO 12 
      IF ( DABS(PN) .GT. 1.D-08 ) THEN
         PNN = DSQRT(1.d0-PN*PN)
         IF (P(50,I) - PN*P(50,II) .LT. D0) PNN = -PNN 
         write(*,*) i,ii,pn
	 mr = max0(max(i),max(ii)) 
         DO 25 J = 1,MR
            P(J,I) =(P(J,I) - PN*P(J,II))/PNN
25       CONTINUE
         AZ(I) = (AZ(I) - PN*AZ(II))/PNN
      END IF
6     CONTINUE
       pnn = QQUADR(I,I,0)
       pnn = dsqrt(pnn)
C       if (pnn .lt. tol) go to 9
       do 30 j = 1, MR
	  P(J,I) = P(J,I)/pnn
30     continue
       pnn = QQUADR(I,I,0)
       rad = qquadr(i,i,1)
       write(*,*) ' Orbital ',el(i), ' with radius',rad, ' <|> = ',pnn
       WRITE(2) ATOM,TERM,EL(I),MR,Z,ET(I),EK(I),AZ(I),
     :     (P(J,I),J=1,MR)
       I = I + 1
9     IF ( I .LT. NWD) THEN
	GO TO 12
      else
	write(*,*) ' TOO many orbitals'
      end if
13    continue
      
      REWIND(UNIT=2)
      I = 1
      M = 1
  912 READ(2,END=913) ATOM,TERM,EL(I),MR,Z,ET(I),EK(I),AZ(I),
     :     (P(J,I),J=1,MR)
      IF (EL(I)(1:1) .NE. ' ') THEN
         N(I) = ICHAR(EL(I)(1:1)) - ICHAR('1') + 1
         L(I) = LVAL(EL(I)(2:2))
      ELSE
         N(I) = ICHAR(EL(I)(2:2)) - ICHAR('1') + 1
         L(I) = LVAL(EL(I)(3:3))
      ENDIF
      MM = MR+1
      DO 924 J = MM,NO
924    P(J,I) = D0
      MAX(I)=MR
      i = i + 1
      GOTO 912
913   CONTINUE
      write(*,*)
      WRITE(*,*) ' Overlap matrix after rotation' 
      write(*,*)
      DO 812 K = 1,NUMBERW
         DO 813 I = 1,K
           if (l(k).ne.l(i)) go to 813
            write(*,1000) el(k),el(i),qquadr(k,i,0)
813      CONTINUE
812   CONTINUE
      write(*,*)
      WRITE(*,*) ' The rotated orbitals are written on wfn.rot !'
1000  FORMAT(1H ,'<',A3,'|',A3,'> = ',F18.10)
1001  FORMAT(1H ,' <r> of original ',A3,' electron = ',F18.10
     :,' <|> = ',F18.10)
1002  FORMAT(//,1H ,'The first orbital of this symmetry has # ',I3)
1004  FORMAT(//1H ,'The orbital ',A3,' are  mixed (identical weights)
     : with the ',I4,//,' orbitals of the same symmetry. The latter
     : will be orthogonalized to the inner functions.',//)
      end
*
*     ------------------------------------------------------------------
*              Q U A D R
*     ------------------------------------------------------------------
*
*                                   kk
*       Evaluates the integral of  r   P (r) P (r) with respect to r
*                                       i     j
*
      DOUBLE PRECISION  FUNCTION QQUADR(I,J,KK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NOD=220,NWD=60)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),P(NOD,NWD),YK(NOD),
     :   YR(NOD),X(NOD),AZ(NWD),L(NWD),MAX(NWD),N(NWD)
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
      DIMENSION G(NOD)
*
      K = KK + 2
      LI = L(I)
      LJ = L(J)
      DEN = LI + LJ + 1 + K
      ZR = Z*R(4)
      BI = (P(4,I)/(AZ(I)*R2(4)*R(4)**LI) - D1+ZR/(LI+1) )/ZR**2
      BJ = (P(4,J)/(AZ(J)*R2(4)*R(4)**LJ) - D1+ZR/(LJ+1) )/ZR**2
      ALPHA= (D1/(LI + 1) + D1/(LJ + 1))/(DEN + D1)
      ZR = Z*R(1)
      BETA = (DEN+D1)*ALPHA**2 - D2*(BI+BJ+D1/((LI+1)*(LJ+1)))/(DEN+D2)
      D = P(1,I)*P(1,J)*R(1)**K*(((BETA*ZR+ALPHA)*ZR+D1)/(DEN*H1)+D5)
      DD = D0
      M = MIN0(MAX(I),MAX(J))
      DO 1 JJ = 1,M
        G(JJ) = P(JJ,I)*P(JJ,J)*R(JJ)**K
    1 CONTINUE
      DO 2 JJ= 3,M,2
        D = D + G(JJ)
    2 CONTINUE
      DO 3 JJ = 2,M,2
        DD = DD + G(JJ)
    3 CONTINUE
      QQUADR = H1*(D + D2*DD)
      RETURN
      END

