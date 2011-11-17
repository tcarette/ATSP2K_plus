*
*     ------------------------------------------------------------------
*    3-12      G R A N G E
*     ------------------------------------------------------------------
*
*       Controls the calculation of off-diagonal energy parameters.
*   It searches for all pairs (i,j) which are constrained through  an
*   orthogonality requirement.   When one of the pair , say P
*                                                            i
*   must be orthogonal not only to P  but also to P  where n = n ,
*                                   j              k        j   k
*   a system of equations must be solved, the exact form depending on
*   whether or not any of the functions are part of the frozen  core.
*   When  only one pair with a given set of principal quantum numbers
*   is present, ELAGR(I,J) is used to  determine  the  off-  diagonal
*   energy  parameters  as  long  as  |q  -q | > 0.05.  Otherwise Eq.
*                                       i   j
*   (7-10) modified for configuration interaction is used.
*
*
      SUBROUTINE GRANGE(ivar)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
        PARAMETER (NWD=70,NOD=220,NOFFD=1600,NODD=100)
        INTEGER ivar(NWD)
*
*     MPI stuff ***********************************************
*
	INCLUDE 'mpif.h'
	parameter (MAXPROC=100)
	common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)	
	common /PVM/ istart,ifinish
****************************************************************

        CHARACTER EL*3,ATOM*6,TERM*6
        INTEGER       IN,OUT,PRI,OUC,OUF,OUH,OUD
        COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,IUD,OUC,OUF,OUD,OUH,ISCW
        LOGICAL       FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON/LABEL/ EL(NWD),ATOM,TERM
        COMMON/TEST/  FAIL,OMIT,EZERO,REL,ALL,TRACE
       COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
        COMMON/PARAM/ H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,
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
      POINTER (pkval,kval(1)),(pvalue,value(1)),
     :        (pih,ih(1)),(pjh,jh(1)),
     :        (pcoeff,coeff(1)),(pnijptr,nijptr(1)),
     :        (pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PVALUE,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
*
        COMMON AC(NODD,NODD),BC(NODD),JV(NODD),IV(NODD)
        LOGICAL  FIRST,both(NODD)
*
*       CLEAR THE ARRAY FOR CHANGING THE WT ARRAY
 
*
       nit = nwf-ib+1
*
*  *****  FOR EACH l COMPUTE OFF-DIAGONAL ENERGY PARAMETERS
*
        DO 10 IL = 0,7
           IJ = 0
           DO 11 iiv = 1,nit
             i = ivar(iiv)
             IF ( L(I) .NE. IL ) GO TO 11
             DO 12 J = 1,nwf
               IF (l(j).eq. il .and. i.ne.j) THEN
*                .. we have an orthogonality condition
*                   should it be included in this order?
                 jjv = 0
                 do 15 jj = 1,nit
                   if (ivar(jj) .eq. j) jjv = jj
   15            continue
                 if (i .lt. j) then
                    iin = i
                    jin = j
                 else
                    iin = j
                    jin = i
                 end if
                 if (jjv .lt. iiv .and.
     :              abs(e(iin,jin)) .gt. 1.d-10) then
*                  .. include
                   IJ = IJ + 1
                   IF ( IJ.GT.(NODD)) then
                     write(iscw,*) 'TOO MANY LAGRANGE MULTIPLIERS'
                     stop
                   end if
                   IV(IJ) = I
                   JV(IJ) = J
                   if (jjv .eq.0) then
                      both(ij) = .false.
                   else
                      both(ij) = .true.
                   end if
                 END IF
               END IF
12           CONTINUE
11         CONTINUE
*
*  ***** IJ IS THE NUMBER OF LAGRANGE MULTIPLIERS FOR l = IL
*
           IF (IJ .EQ. 0) GO TO 10
* * * * * * *
           if (myid.eq.0) then
              call dinit(ij,D0,bc,1)
              DO 14 III = 1,IJ
                 call dinit(ij,D0,AC(1,iii),1)
14            CONTINUE
           endif
* * * * * * *
           DO 16 iiV = 1,nit
              i = ivar(iiv)
              IF ( L(I) .NE. IL ) GO TO 16
                 FIRST = .TRUE.
                 DO 18 II = 1,IJ
                    J = 0
                    IF ( IV(II) .EQ. I) THEN
                       J = JV(II)
                     ELSE IF ( JV(II) .EQ. I) THEN
                       J = IV(II)
                    END IF
                    IF ( J .NE. 0) THEN
                       IF (FIRST) THEN
                          CALL POTL(I)
                          CALL XCH(I,2)
* * * * * * *>
                          if (myid.eq.0) call dcopy(no,yr,1,yk,1)
                          FIRST = .FALSE.
                       END IF
* * * * * * *
                       if (myid.eq.0) then
                         DO 22 JJ = 1,NO
                            YR(JJ) = P(JJ,J)
22                       CONTINUE
                         BC(II) = BC(II) +
     :                    HL(EL,I,J,REL)-D2*QUADS(I,J,1)-QUAD(J,NO,YR,X)
                       endif
* * * * * * *
                    END IF
18              CONTINUE
16         CONTINUE
* * * * * * *
           if (myid.eq.0) then
           DO 24 II = 1,IJ
              DO 26 III = 1,II
                 IF ( II .EQ. III) THEN
                    AC(II,II) = D1/SUM(IV(II))
                    IF (both(ii)) THEN
                       AC(II,II) = AC(II,II) + D1/SUM(JV(II))
                    END IF
                  ELSE IF (IV(II) .EQ. IV(III) .AND.
     :                    E(JV(II),JV(III)) .EQ. D0 ) THEN
                     AC(II,III) = QUADR(JV(II),JV(III),0)/SUM(IV(II))
                     AC(III,II) = AC(II,III)
                  ELSE IF (JV(II) .EQ. JV(III) .AND. both (ii)
     :              .AND. E(IV(II),IV(III)) .EQ. D0) THEN
                     AC(II,III) = QUADR(IV(II),IV(III),0)/SUM(JV(II))
                     AC(III,II) = AC(II,III)
                  END IF
26            CONTINUE
24         CONTINUE
           CALL LINEQN(NODD,IJ,AC,BC)
          endif

	 call MPI_BCAST(bc,ij,MPI_DOUBLE_PRECISION,0,
     $                      MPI_COMM_WORLD,ierr)

         DO 28 II=1,IJ
              CALL EIJSET(IV(II),JV(II),BC(II)/SUM(IV(II)))
              IF ( both(ii) )
     :            CALL EIJSET(JV(II),IV(II),BC(II)/SUM(JV(II)))
28         CONTINUE
10    CONTINUE
*
*  *****  PRINT THE OFF-DIAGONAL ENERGY PARAMETERS
*
* * * * * * *
        if (myid.eq.0) then
        DO 30 iiv =1,nit
           i = ivar(iiv)
           DO 32 J = 1,nwf
              jjv = 0
              do 33 jj = 1,nit
                if (ivar(jj) .eq. j) jjv = jj
33            continue
              if (jjv .lt. iiv .and. abs(e(i,j)) .gt. 1.d-10) then
                 WRITE(OUT,35) EL(I),EL(J),E(I,J),EL(J),EL(I),E(J,I)
35               FORMAT(7X,2(3X,'E(',2A3,') =',F12.5))
              END IF
32         CONTINUE
30      CONTINUE
        endif
* * * * * * *

        RETURN
        END
