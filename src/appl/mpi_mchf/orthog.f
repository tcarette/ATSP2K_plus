*
*     ------------------------------------------------------------------
*    3-22      O R T H O G
*     ------------------------------------------------------------------
*
*       This routine orthogonalizes the set of radial functions when an
*   orthogonality constraint applies.  A Gram-Schmidt type of  process
*   is used.  When more than one radial function with a given (nl) is
*   present, it may be necessary to solve a 2x2 system of equations.
*
*
      SUBROUTINE ORTHOG(ivar)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
Ctc 19/03/2008 : crash of large calculations on several nodes of hydra@ulb
        include 'mpif.h'
        parameter (MAXPROC=100)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
Ctc
*
        PARAMETER (NWD=94,NOD=220,NOFFD=800)
        INTEGER ivar(NWD)
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
*     .. orthogonalize only orbitals that were varied
*        in the order of variation
      nit = nwf-ib+1
      DO 10 iV= 1,nit
        i = ivar(iv)
        DO 20 ii = 1,nwf
          if (e(i,ii) .ne. 0.D0 .and. ii .ne. i) then
*           .. we have an orthogonality condition
            jv = 0
            do 22 iiv = 1,nit
              if (ivar(iiv) .eq. ii) jv = iiv
22          continue
            if (jv .lt. iv) then
*           .. orbital i is orthogonalized to ii
            PN = QUADR(I,II,0)
            IF ( DABS(PN) .GT. 1.D-15 ) THEN
*             .. only write if sufficiently large
Ctc 19/03/2008 : crash of large calculations on several nodes of hydra@ulb
Ctc           if (DABS(PN) .gt. 1.d-10) WRITE(out,24) EL(ii),EL(I),PN
             if (DABS(PN) .gt. 1.d-10) WRITE(80+myid,24) EL(ii),EL(I),PN
Ctc
24            FORMAT(6X,'<',A3,'|',A3,'>=',1PD8.1)
              VARIED(I) = .TRUE.
              M = MAX0(MAX(I),MAX(II))
              DO 25 J = 1,M
25              P(J,I) =P(J,I) - PN*P(J,II)
*             .. let's modify the AZ also
	      AZ(i) = (AZ(i) - PN*AZ(ii))
	      Max(i) = m
            END IF
            END IF
          end if
20      CONTINUE
	m = max(i)
	PN = 1.d0/sqrt(quadr(i,i,0))
	if (P(50,i) .lt. d0) PN=-PN
	do 27 j = 1,m
	    P(j,i) = P(j,i)*PN
27      continue
	Az(i) = Az(i)*PN
30      IF (abs(p(m,i)) .lt. 1.d-12) then
	   p(m,i) = 0.d0
	   m = m-1
	   if (m .gt. 100) go to 30
	end if
	max(i) = m
10    CONTINUE
      END
