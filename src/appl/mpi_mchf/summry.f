*
*     ------------------------------------------------------------------
*    3-35      S U M M R Y
*     ------------------------------------------------------------------
*
*       The results of a calculation are summarized.   These include
*   the following for each electron:
*
*          E(NL)   - diagonal energy parameter
*          AZ(NL)  - starting parameter, P(r)/r**(l+1) as r -> 0.
*          SIGMA   - screening parameter as defined by Eq. (6-  ).
*          1/R**3  - expected value of <1/r**3>
*          1/R     - expected value of <1/r>
*          R       - expected mean radius
*          R**2    - expected value of <r**2>
*          I(NL)   - -(1/2)<nl|L|nl>
*          KE      - I(NL) + Z <r>
*          REL     - Relativistic shift (mass-velocity, Darwin term,
*                    spin-spin contact term)
*
*   These results are followed by:
*
*          TOTAL ENERGY--RELATIVISTIC OR NON-RELATIVISTIC (ET)
*          KINETIC ENERGY-- NON-RELATIVISTIC (EN)
*          POTENTIAL ENERGY (EP) = ET - EN
*          RATIO                 = - EP/EN
*                      k   k   k
*   The values of all F , G , R  and <nl|L|n'l> integrals which enter
*   into the calculation are printed, but only if OUD > 0.
*
*
      SUBROUTINE SUMMRY
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
Ctc 19/03/2008 : crash of large calculations on several nodes of hydra@ulb
        include 'mpif.h'
        parameter (MAXPROC=100)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
Ctc
*
        PARAMETER (NWD=70,NOD=220,NOFFD=800)
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
      POINTER (pkval,kval(1)),(pcoef,coef(1)),
     :        (pih,ih(1)),(pjh,jh(1)),
     :        (pcoeff,coeff(1)),(pnijptr,nijptr(1)),
     :        (pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PCOEF,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
*
      COMMON R1(NWD),RM1(NWD),EK(NWD),SS(3)
*
      WRITE(PRI,9) ATOM,TERM
9     FORMAT(/// 24X,5HATOM ,A6,3X,5HTERM ,A6//63X,13HMEAN VALUE OF,21X,
     1   22HONE ELECTRON INTEGRALS /5X,2HNL,7X,5HE(NL),9X,'I(nl)',
     2   5X,'-L(nl)/2', 5X, 'RelS', 5X, 'S(nl)', 5X, 'Az(nl)')
 
      PI = ACOS(-D1)
      EN = D0
      REL = .FALSE.
*
*  *****  COMPUTE AND PRINT ONE-ELECTRON PARAMETERS
*
      DO 10 I = 1,NWF
      R1(I) = QUADR(I,I,1)
      EK(I) = -D5*HL(EL,I,I,REL)
      RM1(I) = QUADR(I,I,-1)
      EKINP = EK(I) + Z*RM1(I)
      EN = EN+ SUM(I)*EKINP
      RH = 3*N(I)*N(I) - L(I)*(L(I) + 1)
      SC = Z - D5*RH/R1(I)
      S(I) = SC
      RELS = RLSHFT(I,I)
      WRITE (3,15)EL(I),E(I,I),EK(I),EKINP,RELS,S(I),AZ(I)
15    FORMAT(1X,A3,F14.7,3F13.6,F7.2,F13.5)
10    CONTINUE
*
*  *****  Compute Moments
*
      WRITE(PRI,8) 'Delta(R)'
 8    FORMAT(//2X,'nl',6X,A8,5X,'1/R**3',7X,'1/R',9X,'R',8X,'R**2')
      DO 11 I = 1,NWF
      RM3 = 0
      IF (L(I) .NE. 0) RM3 = QUADR(I,I,-3)
      RP2 = QUADR(I,I,2)
      RZ = 0.
      IF ( L(I) .EQ. 0) RZ = AZ(I)**2/(4.*PI)
      WRITE(PRI,16) EL(I),RZ,RM3,RM1(I),R1(I),RP2,SUM(I)
16    FORMAT(1X,A3,F14.3,F13.4,F11.5,F10.5,F11.5,F13.10)
11    CONTINUE
*
*  *****  ADD CONTRIBUTION FROM THE 'L' INTEGRALS
*
      IBEGIN = INTPTR(5) + 1
      IEND = INTPTR(6)
      DO 32 I = IBEGIN,IEND
        call unpacki(6,i,kv,iel1,iel2,iel3,iel4)
        IF (IEL1 .NE. IEL2) THEN
          CONT = HL(EL,IEL1,IEL2,REL)-D2*Z*QUADR(IEL1,IEL2,-1)
          EN = EN + CONT*COEF(I)
      END IF
32    CONTINUE
      EPOTL = ETOTAL - EN
      RATIO =-EPOTL/EN
Ctc 19/03/2008 : crash of large calculations on several nodes of hydra@ulb
Ctc   WRITE(OUT,26) ETOTAL,EPOTL,EN,RATIO
      WRITE(80+myid,26) ETOTAL,EPOTL,EN,RATIO
Ctc
      WRITE(PRI,26) ETOTAL,EPOTL,EN,RATIO
26    FORMAT(//5X,'ENERGY (a.u.)'/5X,'------'/
     : 10X,' Total              ',F16.9/
     : 10X,' Potential          ',F16.9/
     : 10X,' Kinetic            ',F16.9/
     : 10X,' Ratio              ',F16.9)

      if (qwt.ne.0) call dalloc(qwt,ncfg)
      if (iqp.ne.0)  call dalloc(iqp,nod*nwf)
      if (iqn.ne.0) call dalloc(iqn,nwf)
      if (iql.ne.0) call dalloc(iql,nwf)
      if (iqaz.ne.0) call dalloc(iqaz,nwf)
      if (iqmax.ne.0) call dalloc(iqmax,nwf)
      if (qvard.ne.0) call dalloc(qvard,nwf)
      if (qsum.ne.0) call dalloc(qsum,nwf)
      if (qs.ne.0) call dalloc(qs,nwf)
      if (qdpm.ne.0) call dalloc(qdpm,nwf)
      if (qmeth.ne.0) call dalloc(qmeth,nwf)
      if (qieptr.ne.0) call dalloc(qieptr,nwf)
      if (pkval.ne.0) call dalloc(pkval,idim)
      if (pcoef.ne.0) call dalloc(pcoef,idim)
      if (pih.ne.0) call dalloc(pih,nze)
      if (pcoeff.ne.0) call dalloc(pcoeff,1.5*ncodim)
      if (pinptr.ne.0) call dalloc(pinptr,1.5*ncodim)
*
13    RETURN
      END
