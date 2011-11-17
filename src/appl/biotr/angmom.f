*
*     ------------------------------------------------------------------
*       A N G M O M
*     ------------------------------------------------------------------
*
      SUBROUTINE ANGMOM(NEW,NZERO,IFIRST)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=128, NWCD=20)
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(4),IALL,JSC(3),ISCW
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON /fout/lcount,nrec,iflag,lij,nij
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      POINTER(QIORTH,IORTH(1))
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     : ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     :     QIORTH
      COMMON/RED/INDL(20),INDR(20),IORST(20),NCI1,KPL,KPR,NC,LIS(16),
     : JIST,JFST,NCONTR
      CHARACTER*30 FORMAT(3)
      CHARACTER*1 NCHAR(8)
*
      DATA NCHAR/'1','2','3','4','5','6','7','8'/
      DATA FORMAT/'(2H <,I3,5H |H| ,I2,6H > = <, ',
     :              '  (A3,1H(,I2,1H)) ,5H |H| ,   ',
     :                '(A3,1H(,I2,1H)),2H >,/)     '/
*
* --- THIS PROGRAMME CONSIDERS EITHER SUPERPOSITION OF CONFIGURATIONS OR
*     MULTI-CONFIGURATIONAL HARTREE-FOCK WAVE FUNCTIONS.  USING THE
*     RESULT THAT THE TWO-ELECTRON HAMILTONIAN MATRIX ELEMENT
*     (PSI/V/PSIP)  CAN BE WRITTEN AS A SUM OF SLATER INTEGRALS, THE
*     PRESENT CODE  -  WEIGHTS  -  CALCULATES THE COEFFICIENTS OF THESE
*     INTEGRALS.  PSI AND PSIP ARE ALLOWED TO RUN OVER NCFG CONFIGURATNS
*
*
* --- CONSIDER (PSI/V/PSIP) AS PSI AND PSIP RUN OVER ALL CONFIGURATIONS
*
      NC = 0
      NFIRST = NCFG - NEW + 1
      DO 1 JA=NFIRST,NCFG
        if(mod(ja,100).eq.0) write(ISCW,*) '   ja = ',ja
        DO 2 JB=1,JA
          IFLAG=0
          IDIAG=0
          IF (NORTH .NE. 0) CALL NORTBPn(JA,JB)
*          N1=NOCCSH(JA)
*          N2=NOCCSH(JB)
*          IF (IFULL .NE. 0) THEN
*          FORMAT(2)(2:2) = NCHAR(N1)
*          FORMAT(2)(30:30) = NCHAR(N2)
*          WRITE(IWRITE,'(///)')
*          WRITE(IWRITE,FORMAT) JA,JB,
*     :        (IAJCMP(NOCORB(J,JA)),NELCSH(J,JA),J=1,N1),
*     :        (IAJCMP(NOCORB(J,JB)),NELCSH(J,JB),J=1,N2)
*          END IF
*
* --- SET UP DEFINING QUANTUM NUMBERS FOR EACH MATRIX ELEMENT
*
          CALL SETUP(JA,JB,LET)
          IF (LET .NE. 0) THEN
*          IF(IBUG1.GT.0.OR.IBUG2.GT.0) CALL VIJOUT(JA,JB)
*
* --- TEST ON number of electrons, parity, coupling
*
            CALL ORTHOGG(LET)
            IF(LET.NE.0) THEN
*
* --- IF NO SUCH ORTHOGONALITY IS EXHIBITED, CALCULATE WEIGHTS OF SLATER
*     INTEGRALS
*
              CALL LMATRIX
              IF (IFLAG .NE. 0) NIJ = NIJ + 1
            ENDIF
          ENDIF
    2   CONTINUE
    1 CONTINUE
*
      nwf = maxorb
      if (nwf .gt. 1) then
        print*,' dalloc, qiorth: nwf*(nwf-1)/2 = ',(nwf*(nwf-1))/2
        call dalloc(qiorth,(nwf*(nwf-1))/2)
      end if
      print*,' dalloc, qnoc: ncfg = ',ncfg
      call dalloc(qnoc,ncfg)
      print*,' dalloc, qnelcsh: 8*ncfg = ',8*ncfg
      call dalloc(qnelcsh,8*ncfg)
      print*,' dalloc, qnocorb: 8*ncfg = ',8*ncfg
      call dalloc(qnocorb,8*ncfg)
      print*,' dalloc, qj1    :15*ncfg = ',15*ncfg
      call dalloc(qj1,15*ncfg)
      print*,' dalloc, qljclcd: nwcd = ',nwcd
      call dalloc(qljclsd,nwcd)
      qnelcsh=0;
      qnocorb=0;
      qj1=0;
      qljclsd=0
*
*
      RETURN
      END
