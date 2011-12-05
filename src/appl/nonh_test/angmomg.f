*     ------------------------------------------------------------------
*       A N G M O M G
*     ------------------------------------------------------------------
*
      SUBROUTINE ANGMOMG(NEW,NZERO,IFIRST,nih)
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
      COMMON /fout/n,ntot,idum(6),nrec(8),iflag,lij,nij,qjptr
      COMMON/DIAGNL/IDIAG,JA,JB
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
        NIJ = NIJ + 1
        LIJ = LIJ + 1
        if (nih .eq. LSDIM) then
           write(11) nih, (jan(i),i=1,nih)
           write(12) nih, (ico(i),i=1,nih)
           nih = 0
        end if
        nih = nih + 1
        jan(nih) = ja
        jbn(nih) = jb
        ico(nih) = ntot
      endif
    6 CONTINUE
*
      END
