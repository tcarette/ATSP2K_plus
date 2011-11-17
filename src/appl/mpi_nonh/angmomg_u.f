*
*     ------------------------------------------------------------------
*       A N G M O M G
*     ------------------------------------------------------------------
*
      SUBROUTINE ANGMOMG(NEW,NZERO,IFIRST)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=70, NWCD=20)
      PARAMETER (LSDIM=30000)
*
      POINTER (qcn,cn(1)),(qinptr,inptr(lsdim)),(qpackn,ipackn(1)),
     :        (qlused,lused(1)),(qnijptr,nijptr(lsdim)),
     :        (qjan,jan(lsdim)),(qjbn,jbn(lsdim)),(qico,ico(1)),
     :        (qintptr,idummy(1)) 
      COMMON /buffer/qcn,qinptr,qpackn,qlused,qintptr,lmax,qnijptr,
     :               qjan,qjbn,qico
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(4),IALL,JSC(3),ISCW
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      POINTER  (qjptr, jptr(1))
      COMMON /fout/n,ntot,idum(6),nrec(8),iflag,lij,nij,qjptr,cf_tot
      COMMON/DIAGNL/IDIAG,JA,JB
      CHARACTER*30 FORMAT(3)
      CHARACTER*1 NCHAR(8)
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

      LOGICAL lused
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
      if (mod(jb,100) .eq. 0) write(ISCW,*) '   jb =',jb
      NFIRST = NCFG - NEW + 1
      irow = max(jb,nfirst)
      if (jb .gt. nzero) then
        last = jb
      else
        last = ncfg
      end if
      Do 6 JA = irow, last
*     DO 6 JA=irow,NCFG
*     IF (JB .GT. NZERO .AND. IFIRST .EQ. 1 .AND. JA .NE. JB ) GO TO 6
*     IF (JB.GT.NZERO .AND. JB.LT.NFIRST .AND. IFIRST.EQ.0) GO TO 6
*     write(ISCW,*) '     ja =',ja
      IFLAG=0
      IDIAG=0
      IF(JA.EQ.JB) IDIAG=1
      IF (NORTH .NE. 0) THEN
        WRITE(ISCW,'(A)') ' this prog. com. with orthogonal orbitals'
        STOP
      ENDIF
      N1=NOCCSH(JA)
      N2=NOCCSH(JB)
*      IF (IFULL .NE. 0) THEN
*         FORMAT(2)(2:2) = NCHAR(N1)
*         FORMAT(2)(30:30) = NCHAR(N2)
*         WRITE(IWRITE,'(///)')
*         WRITE(IWRITE,FORMAT) JA,JB,
*     :        (IAJCMP(NOCORB(J,JA)),NELCSH(J,JA),J=1,N1),
*     :        (IAJCMP(NOCORB(J,JB)),NELCSH(J,JB),J=1,N2)
*      END IF
*
* --- SET UP DEFINING QUANTUM NUMBERS FOR EACH MATRIX ELEMENT
*
C      CALL SETUP(JA,JB)
      CALL SHELLS(JA,JB,LET)
      IF(LET.EQ.0) GO TO 6
      IF(IBUG1.GT.0.OR.IBUG2.GT.0) CALL VIJOUT(JA,JB)
*
* --- TEST ON POSSIBLE RECOUPLING ORTHOGONALITY
*
      CALL ORTHOGG(LET)
      IF(LET.EQ.0) GO TO 6
*
* --- IF NO SUCH ORTHOGONALITY IS EXHIBITED, CALCULATE WEIGHTS OF SLATER
*     INTEGRALS
*
       CALL NONRELAT
      IF (IFLAG .NE. 0) then
        if (lij/lsdim .eq. 1) then
*          .. a record needs to be written to iout
           write(iout) lsdim,(jan(i),i=1,lsdim), (ico(i),i=1,lsdim)
           lij = 0
        end if
        NIJ = NIJ + 1
        LIJ = LIJ + 1
        jan(lij) = ja
        jbn(lij) = jb
        ico(lij) = ntot
      endif
    6 CONTINUE
*
      END
