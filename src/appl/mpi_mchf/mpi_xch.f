        SUBROUTINE XCH(I,IOPT)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
        PARAMETER (NWD=70,NOD=220,NOFFD=1600)
*
*     MPI stuff ***********************************************
*
	INCLUDE 'mpif.h'
	parameter (MAXPROC=100)
	common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)	
	common /PVM/ istart,ifinish
****************************************************************
        dimension yhelp(2*NOD)
****************************************************************

        CHARACTER EL*3,ATOM*6,TERM*6
        INTEGER       IN,OUT,PRI,OUC,OUF,OUH,OUD
        COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,IUD,OUC,OUF,OUD,OUH,ISCW
        COMMON/LABEL/ EL(NWD),ATOM,TERM
        LOGICAL       FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON/TEST/  FAIL,OMIT,EZERO,REL,ALL,TRACE
       COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
        COMMON/PARAM/  H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,
     :                NCFG,IB,IC,ID,
     :                D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,
     :                NSCF,NCLOSD,RMASS
        COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      POINTER (QWT,WT(1)),(QWP,WP(1)),  (IQWPTR,W(1))
      COMMON /CFGS/ETOTAL,QWT,QWP,IQWPTR
*
        POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :     (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :     (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
        COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :             QMETH,QIEPTR
*
      POINTER (pkval,kval(1)),(pcoef,coef(1)),
     :        (pih,ih(1)),(pjh,jh(1)),
     :        (pcoeff,coeff(1)),(pnijptr,nijptr(1)),
     :        (pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PCOEF,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
*
        LOGICAL SAME,EXIT
*

        DO 1 J=1,NO
  1     X(J) = D0
        DO 2 J = 1+myid,NWF,nprocs
           IF ((I.LE.NCLOSD .AND. I.NE.J) .OR.
     :         (I.GT.NCLOSD .AND. J.LE.NCLOSD))  THEN
              DO 4 K = IABS(L(I)-L(J)),L(I)+L(J),2
                 C = - D2*CB(L(I),L(J),K)*SUM(J)
                 CALL YKF(J,I,K,REL)
                 DO 6 JJ = 1,NO
                    X(JJ) = X(JJ) + C*YK(JJ)*P(JJ,J)
  6              CONTINUE
  4           CONTINUE
           END IF
  2     CONTINUE
        SUMI = SUM(I)
        IF (I .LE. NCLOSD) GO TO 51
*
        IBEGIN = INTPTR(1)+1
        IEND = INTPTR(2)
        DO 7 INT = IBEGIN+myid,IEND,nprocs
           call unpacki(2,int,kv,iel1,iel2,iel3,iel4)
           IE1 = 0
           IF (IEL1 .EQ. I) THEN
              IE1 = IEL1
              IE2 = IEL2
           ELSE IF (IEL2 .EQ. I) THEN
              IE1 = IEL2
              IE2 = IEL1
           END IF
           IF (IE1 .NE. 0) THEN
              C = D2*COEF(INT)/SUMI
              CALL YKF(IE1,IE2,KV,REL)
              DO 8 JJ = 1,NO
                 X(JJ) = X(JJ) + C*YK(JJ)*P(JJ,IE2)
 8            CONTINUE
           END IF
 7      CONTINUE
*
*       gnore coef with abs(c) smaller than zlimit
        zlimit = 1.D-10
        inonzero = 0
        izero = 0

        IBEGIN = INTPTR(4) + 1
        IEND = INTPTR(5)
        DO 50 INT = IBEGIN+myid,IEND,nprocs
           call unpacki(5,int,kk,i1,i2,j1,j2)
           IF ((I1-I)*(I2-I) .EQ. 0 .OR. (J1-I)*(J2-I) .EQ. 0) THEN
              C = COEF(INT)/SUMI
              CC = C
*   skip all zero coefficients
         if (abs(c).gt.zlimit) then
              CC = C
*
*
*  ***** COUNT THE NUMBER OF OCCURRENCES OF I
*
              IK = 0
              IF (I1 .EQ. I) IK = IK + 1
              IF (I2 .EQ. I) IK = IK + 1
              IF (J1 .EQ. I) IK = IK + 1
              IF (J2 .EQ. I) IK = IK + 1
              EXIT = .FALSE.
              DO 11 II2=1,2
              DO 12 II1=1,2
              GO TO (10, 20, 30, 40) IK
10          CONTINUE
*
*  ***** I OCCURS JUST ONCE IN RK
*
              IF (I1 .NE. I) GO TO 13
              GO TO 16
20            CONTINUE
*
*  ***** I OCCURS TWICE IN THE RK INTEGRAL
*
              IF (I1 .NE. I) GO TO 13
              IF (J1 .EQ. I) GO TO 17
*
*  ***** TEST IF THE PAIR (I1,J1) = PAIR (I2,J2)
*
              ICODE1 = 100*I1 + J1
              ICODE2 = 100*I2 + J2
              ICODE3 = 100*J2 + I2
              SAME = ICODE1 .EQ. ICODE2 .OR. ICODE1 .EQ. ICODE3
              IF ( .NOT. SAME ) GO TO 15
              GO TO 17
30            CONTINUE
*
*  ***** I OCCURS THREE TIMES IN THE RK INTEGRAL
*
*
              IF (I1 .EQ. I) GO TO 13
              CALL YKF(I2, J2, KK, REL)
              DO 33 J = 1,NO
33              X(J) = X(J) + CC*P(J,I1)*YK(J)
              CALL YKF(I1, J1, KK, REL)
              CC = D2*CC
              DO 34 J = 1,NO
34               X(J) = X(J) + CC*P(J,I2)*YK(J)
              GO TO 50
*
*  ***** I OCCURS FOUR TIMES IN RK INTEGRAL
*
40            CC = D4*CC
              GO TO 16
17            CC = D2*CC
16            EXIT = .TRUE.
!              call MPI_BCAST(EXIT,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
15            CALL YKF(I2,J2,KK,REL)
              DO 14 J=1,NO
14               X(J) = X(J) +CC*P(J,J1)*YK(J)
              IF (EXIT) GO TO 50
13            III = I1
              I1= I2
              I2= III
              III = J1
              J1 = J2
12            J2 = III
              III = I1
              I1 = J1
              J1 = III
              III = I2
              I2= J2
11            J2= III
           END IF
        end if
*

50      CONTINUE
*
51      IBEGIN = INTPTR(5) + 1
        IEND = INTPTR(6)
        DO 60 INT = IBEGIN+myid,IEND,nprocs
          call unpacki(6,int,kv,iel1,iel2,iel3,iel4)
*         ... Include only if off-diagonal ...
          IF (IEL1 .NE. IEL2) THEN
           I1 = IEL1
           I2 = IEL2
           IF (I1 .NE. I) THEN
              ITEMP = I1
              I1 = I2
              I2 = ITEMP
           END IF
           IF (I1 .EQ. I) THEN
              C = COEF(INT)/SUMI
*  skip all coef = 0
        if (abs(c).gt.zlimit) then
              CALL DIFF(I2)
              DO 62 J = 1,NO
                 X(J) = X(J) + C*YK(J)/R(J)
 62           CONTINUE
              DO 64 II = 1,NCLOSD
                 CC = -D2*(4*L(II)+2)*C
                 CALL YKF(II,II,0,REL)
                 DO 65 J = 1,NO
                    X(J) = X(J) + CC*YK(J)*P(J,I2)
 65              CONTINUE
                 DO 66 K = IABS(L(I)-L(II)),L(I)+L(II),2
                    CCC = CC*CB(L(I),L(II),K)
                    CALL YKF(I2,II,K,REL)
                    DO 67 J = 1,NO
                       X(J) = X(J) - CCC*YK(J)*P(J,II)
 67                 CONTINUE
 66              CONTINUE
 64           CONTINUE
           END IF
        end if

           IF (I .LE. NCLOSD) THEN
              C = -D2*COEF(INT)
* skip for all coef = 0
        if (abs(c).gt.zlimit) then
              CALL YKF(I1,I2,0,REL)
              CC = D2*C
              DO 61 J = 1,NO
                X(J) = X(J) + CC*YK(J)*P(J,I)
 61           CONTINUE
              DO 63 K = IABS(L(I)-L(I1)),L(I)+L(I1),2
                CC = C*CB(L(I),L(I1),K)
                CALL YKF(I2,I,K,REL)
                DO 68 J = 1,NO
                   X(J) = X(J) - CC*YK(J)*P(J,I1)
 68             CONTINUE
                CALL YKF(I1,I,K,REL)
                DO 69 J = 1,NO
                   X(J) = X(J) - CC*YK(J)*P(J,I2)
 69             CONTINUE
 63           CONTINUE
           END IF
          END IF
        end if
 60     CONTINUE
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
* Global summation of both x array and YR from routine POTL.  X is 
* temporarily appended at the end of YR. Then call gdsum for 2*N0 vars 
* and with work array yhelp. After this the result is put back to x. 
* WARN: in grange and solve POTL should be called before XCH

        call dcopy(N0,x,1,YR(N0+1),1)
        
        call gdsummpi(YR, 2*no, yhelp)

        call dcopy(N0,YR(N0+1),1,x,1)
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*
        GO TO (75,76,77),IOPT
 76     DO 78 J = 1,NO
 78        X(J) = X(J)/R(J)
        GO TO 75
 77     DO 79 J =1,NO
 79        X(J) = R(J)*X(J)
        DO 74 J = 1,NWF
           IF (J .NE. I) THEN
           C = E(I,J)
           IF (DABS(C) .LE. 1.D-20 ) GO TO 74
           DO 73 JJ = 1,NO
 73        X(JJ) = X(JJ) + C*P(JJ,J)*RR(JJ)
           END IF
 74     CONTINUE
*
*  *****  CHECK IF EXCHANGE IS ZERO: IF SO, METHOD 2 SHOULD BE USED.
*
 75     IF (METH(I) .EQ. 2 .OR. METH(I) .GT. 3) RETURN
        IF ( DABS(X(1)) + DABS(X(2)) + DABS(X(3)) .EQ. D0 ) METH(I) = 2
        END
