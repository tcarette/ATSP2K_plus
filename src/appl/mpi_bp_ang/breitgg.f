*
*     ------------------------------------------------------------------
*       B R E I T G G
*     ------------------------------------------------------------------
*
      SUBROUTINE BREITGG(NEW,NZERO,IFIRST,idg,skip,nze)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (MXIHSH=16)
*
      PARAMETER (LSDIM=30000)
      POINTER (qcn,cn(lsdim)),(qinptr,inptr(lsdim)),
     :        (qnijptr,nijptr(lsdim)),(qjan,jan(lsdim)),
     :        (qjbn,jbn(lsdim)),(qintptr,intptr(0:2*lmax+1,7)),
     :        (qpackn,ipackn(1)),(qlused,lused(1)),(qico,ico(1))
      COMMON /buffer/qcn,qinptr,qpackn,qlused,qintptr,lmax,qnijptr,
     :               qjan,qjbn,qico
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON /fout/n,ntot,iflag,nih,nij,qjptr
      COMMON/IMAGNT/ IREL,ISTRICT,IELST
      COMMON /INFORM/IREAD,IWRITE,IOUT,ISC(8),ISCW
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),
     :      J1QN1(31,3),J1QN2(31,3),IJFUL(16)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      POINTER(QACMULT,ACMULT(1))
      COMMON /SPORB/ QACMULT
CG 
      COMMON/STEGG/IX,IGGG,IRHO,ISIG,IRHOP,ISIGP

!!    Add so that these operators can be manipulated (CFF:March_14, 2001)
      COMMON /OPERAT/ ICOLOM,ISOTOP,IORBORB
      COMMON /BREIT/ ISPORB,ISOORB,ISPSPN

CG 
      LOGICAL INCL,skip
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
*      NFIRST = NCFG - NEW + 1
*      irow = max(jb,nfirst)
*      if (jb .gt. nzero) then
*	last = jb
*      else
*	last = ncfg
*      end if

      NFIRST = NCFG - NEW + 1
      irow = max(jb,nfirst)
      if (jb .gt. nzero) then
!       turn off the two-body relativistic operators (CFF: March_14,2001)
        isoorb = 0
        ispspn = 0
        iorborb = 0
      end if
      last = ncfg
      m = last-irow+1
   
      IFLAG = 0
      nih = 0
      Do 6 JA = irow, last
*       clear acmult  (do we know why?)
	do 60 iel = 1,maxorb
	  acmult(iel) = 0.d0
   60   continue
      INCL = .TRUE.
      IDIAG = 0
      IF(JA.EQ.JB) THEN
        IDIAG=1
        IF(IDG.EQ.1) INCL = .TRUE.
      ENDIF
      IF(INCL) THEN
*
* ... Set up defining quantum numbers for each matrix element.
*
        CALL SHELLS(JA,JB,LET)
        IF(LET.NE.0) THEN
          IF(IHSH.GT.MXIHSH) STOP
*
* ... TEST ON POSSIBLE RECOUPLING ORTHOGONALITY.
*
          CALL SETUPGG
          IF(IX.LE.4) THEN
            CALL ORTHOG(LET,INCL)
            IF(LET.NE.0) THEN
              IF (IFULL.NE.0) WRITE(IWRITE,77)
              IF(IBUG1.NE.0 .OR. IBUG2.NE.0) CALL VIJOUT(JA,JB)
              CALL NONBP
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      if (iflag .eq. 1) then
*       .. this matrix element contributed
       nih = nih+1
       jan(nih) = ja
       ico(nih) = nij
       iflag = 0
      end if
    6 CONTINUE
      nze = nze+nih
*
   77 FORMAT(///30X,'MULTIPLYING FACTOR',11X,'TYPE OF INTEGRAL')
      END
