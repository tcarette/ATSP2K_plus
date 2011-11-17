*
*     ------------------------------------------------------------------
*    3-32      S C F
*     -----------------------------------------------------------------
*
*       This routine controls the SCF procedure described in Chapter
*   7.  If certain input parameters are zero (or blank) they will  be
*   set to their default value.
*
*          Parameter       Default Value
*          --------        -------------
*          CFGTOL          1.D-10
*          SCFTOL          1.D-7
*          IC              (NWF + 1 - IB)/4 + 3
*          NSCF            12
*
*   The self-consistency convergence criterion is
*
*          Z2 = DSQRT( SCFTOL*(Z*NWF/2) )
*
*   It is increased by a factor two at the end of each iteration whereas
*   CFGTOL is increased by DSQRT(2).
*
*
      SUBROUTINE SCF(IVAR,ACFG,SCFTOL,CFGTOL,LD,eigst_weight,cf_tot)

        IMPLICIT DOUBLE PRECISION(A-H,O-Z)

*
        PARAMETER (NWD=60,NOD=220,NOFFD=800,meig=20,mterm=20)
        DIMENSION IVAR(NWD)
*
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
      POINTER (pkval,kval(1)),(pvalue,value(1)),
     :        (pih,ih(1)),(pjh,jh(1)),
     :        (pcoeff,coeff(1)),(pnijptr,nijptr(1)),
     :        (pinptr,inptr(1))
        COMMON/STATE/INTPTR(6),IDIM,PKVAL,PVALUE,NZE,PIH,PJH,
     :               NCODIM,PCOEFF,PNIJPTR,PINPTR

       POINTER (pico,ico(1))
       COMMON/ICOFF/pico

      POINTER (qjptr,jptr(1))
      COMMON/COLS/qjptr


*      POINTER (qhmx,hmx(1)),(qtm,tm(1)),(qtp,tp(1)),
*     :        (qdiag,hii(1)),(qiwork,iwork(1))
*      common/spd/qhmx,qtm,qtp,qdiag,qiwork 

      Logical leigen(meig,mterm), lguess
      double precision			:: eigst_weight(meig,mterm)
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm),
     :        nze_max(mterm)
      character term_bl(mterm)*3
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess,nze_max
      LOGICAL LAST,LD,CONV,ECONV

      integer cf_tot(nblock), NCFG_(1)

*
*  *****  SET THE SCF CONVERGENCE PARAMETER TO AN OPTIMISTIC VALUE
*
      TOL = 1.d-8
      Z2 = SCFTOL*DSQRT(Z*NWF)
      write(OUT,15)
15    format(//)
      write(OUT,16) OMIT,ACFG,SCFTOL,NO,REL
   16 FORMAT(10X,44HWEAK ORTHOGONALIZATION DURING THE SCF CYCLE=,L4/
     :       10X,44HACCELERATING PARAMETER FOR MCHF ITERATION  =,F5.2/
     :       10X,44HSCF CONVERGENCE TOLERANCE (FUNCTIONS)      =,1PD9.2
     :      /10X,44HNUMBER OF POINTS IN THE MAXIMUM RANGE      =,I4/
     :       10X,44HRELATIVISTIC DIAGONAL  ENERGY CORRECTIONS  =,L4)
      nit = nwf-ib+1
*
*  *****  SET ITERATION PARAMETERS
*
      IPR = 0
      ECONV = .FALSE.
      LAST = .FALSE.
      if (nit == 0) LAST = .true.
      DP1 = D0
      ETOTAL = D0
      EC = D0
      ICYCLE = 0
      call update
      
      IF (IB.GT.NWF)GO TO 17

      CALL DIAG(ECONV,ACFG,CFGTOL,LAST,icycle,eigst_weight,cf_tot)
*
*  *****  PERFORM NSCF SELF-CONSISTENT FIELD ITERATIONS
*
9     DO 100 I = 1,NSCF    
      ICYCLE = ICYCLE + 1

      if ((mod(icycle,10)).eq.0) then
      	write (iscw,*) ICYCLE,' iterations, output written to wfn.out '	
*       	FAIL = .TRUE.
        call output(print)
	rewind(unit=ouf)
      end if
 
      Z2 = SCFTOL*DSQRT(Z*NWF)
 
      WRITE(OUT,7) ICYCLE,CFGTOL,Z2
7     FORMAT(//10X,17HITERATION NUMBER ,I3/10X,16H----------------//
     1 10X,50HCONVERGENCE CRITERIA:ENERGY  (CFGTOL)            =,1PD9.1/
     2 11X,49H                   :FUNCTION(SCFTOL*SQRT(Z*NWF))=,1PD9.1/)

      DP1 = D0
      IF (IB .GT. NWF) GO TO 17
      CALL GRANGE(ivar)
*
*  *****  SOLVE EACH DIFFERENTIAL EQUATION IN TURN
*
         WRITE(OUT,14)
14       FORMAT(/20X,' EL',9X,'ED',13X,'AZ',11X,'NORM',7X,'DPM')

      DO 2 JV = 1,nit
         j = ivar(jv)
         CALL DE(J,ivar)
         IF ( FAIL ) RETURN
         DP = DPM(J)*DSQRT(SUM(J))
         IF ( DP1 .GE. DP ) GO TO 2
         DP1 = DP
         JJ = J
2     CONTINUE
      IF (( ID .EQ. 1) .AND. DP1 .LT. Z2) GO TO 6
      IF ( IC .LE. 0) GO TO 6
*
*  *****  SOLVE IC DIFFERENTIAL EQUATIONS EACH TIME SELECTING THE
*  *****  ONE WITH THE LARGEST DPM
*
      DO 4 II =1,IC
        CALL DE(JJ,ivar)
        IF ( FAIL ) RETURN
        DP1 = D0
        DO 5 JV = 1,nit
          J = ivar(jv)
          DP = DSQRT(SUM(J))*DPM(J)
          IF ( DP1 .GT. DP ) GO TO 5
          JJ = J
          DP1 = DP
5       CONTINUE
        IF (DP1 .LT. Z2) GO TO 6
4     CONTINUE
6     CALL ORTHOG(ivar)
      IF (DP1 .LE. Z2 )  LAST = .TRUE.
      CALL UPDATE
      IF ( LAST ) GO TO 17
      IF ( I .EQ. NSCF ) GO TO 1

12    CALL DIAG(ECONV,ACFG,CFGTOL,LAST,icycle,eigst_weight,cf_tot)
*
*  *****  IF FUNCTIONS APPEAR TO HAVE CONVERGED,SOLVE EACH AGAIN, IN
*  *****  TURN, AND TEST AGAIN
*
      CONV = ECONV .or. DP1 .LE. Z2
      IF (CONV) LAST =.TRUE.
*
*  *****  INCREASE THE CONVERGENCE CRITERION FOR SELF-CONSISTENCY
*

1     CONTINUE

      WRITE(OUT,8) EL(JJ),DP1
8     FORMAT(/ 6X,34HLEAST SELF-CONSISTENT FUNCTION IS ,A3,
     1           27H :WEIGHTED MAXIMUM CHANGE =,1PD10.2)

100   continue

      
18    write(iscw, 13)
13    FORMAT(10X/' SCF ITERATIONS HAVE CONVERGED TO THE ABOVE ACCURACY')
      WRITE(PRI,13)

      FAIL = .TRUE.
*
*  *****  PERFORM FINAL CALCULATIONS
*
17    ACFG = D0

      CALL DIAG(ECONV,ACFG,CFGTOL,.true.,icycle,eigst_weight,cf_tot)
      NIT = NWF - IB + 1
      
*     WRITE(PRI, 105) NIT, DP1, CFGTOL
105   FORMAT(//10X,'NUMBER OF FUNCTIONS ITERATED          =',I6/
     1         10X,'MAXIMUM WEIGHTED CHANGE IN FUNCTIONS  =',D10.2/
     2         10X,'TOLERANCE FOR THE MCHF ITERATION      =',D10.2)

      END
