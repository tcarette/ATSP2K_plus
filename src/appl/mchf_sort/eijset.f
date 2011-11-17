*
*
        SUBROUTINE EIJSET(I,J,E)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
        PARAMETER (NWD=60,NOD=220,NOFFD=800)
*
*
       COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
        COMMON/PARAM/ H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,
     :                NCFG,IB,IC,ID,
     :                D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,
     :                NSCF,NCLOSD,RMASS
        POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :     (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :     (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
        COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :             QMETH,QIEPTR
*
	iscw=0
        IBEGIN = 1
        IF (I .GT. 1) IBEGIN = IEPTR(I-1)+1
        IEND = IEPTR(I)
        DO 10 II = IBEGIN,IEND
           IF (IJE(II) .EQ. J) THEN
              EIJ(II) = E
              RETURN
           END IF
 10     CONTINUE
*
* ***** J-value not found - enter into list
*
        IF (IJE(noffd) .NE. 0)
     :   write(iscw,*)'Too many off-diagonal energy parameters'
*
*  ***** Find point at which the insertion should be made
*
        IEND = IEPTR(I)
        IF (IEND .NE. 0) THEN
           IP = 1
           IF (I .GT. 1) IP = IEPTR(I-1)+1
 30        IF (IJE(IP) .LT. J .AND. IP .LE. IEND) THEN
              IP = IP + 1
              GO TO 30
           END IF
        ELSE
           IP = 1
        END IF
*
* *****  IP is the location in which EIJ should be stored
*        Move other data
*
 
        DO 40 JJ = (noffd)-1,IP,-1
           IJE(JJ+1) = IJE(JJ)
           EIJ(JJ+1) = EIJ(JJ)
 40     CONTINUE
*
* ***** Space has been made - insert data
*
        IJE(IP) = J
        EIJ(IP) = E
*
* ***** Update pointers
*
        DO 50 II = I,NWF
           IEPTR(II) = IEPTR(II) +1
 50     CONTINUE
        END
