*
*     ------------------------------------------------------------------
*       B R E I T G G
*     ------------------------------------------------------------------
*
      SUBROUTINE BREITGG(NEW,NZERO,IFIRST,idg,irfst,skip,nze)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER irfst(*)
      PARAMETER (MXIHSH=16)
*
      POINTER (qh,h(ncfg,3)),(qjan,jan(ncfg,3))
      COMMON /buffer/qh, qjan, nrow(3),iflag(3)
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON/IMAGNT/ IREL,ISTRICT,IELST
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,
     ;             iouhn,iouhz,iouhs,iouhm,iLS,idisk
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
      integer big
      data    big/268435456/
CG 
      COMMON/STEGG/IX,IGGG,IRHO,ISIG,IRHOP,ISIGP
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
      NFIRST = NCFG - NEW + 1
      irow = max(jb,nfirst)
      if (jb .gt. nzero) then
	last = jb
      else
	last = ncfg
      end if
      m = last-irow+1
      do 10 i = 1,m
	h(i,1) = 0.d0
  10  continue
      if (irel .ne. 0) then
	do 20 i = 1,m
 	  h(i,2) = 0.d0
	  h(i,3) = 0.d0
  20    continue
      end if
      nrow(1) = 0
      nrow(2) = 0
      nrow(3) = 0
      Do 6 JA = irow, last
*     write(ISCW,*) '   ja =',ja
*       clear acmult
	do 60 iel = 1,maxorb
	  acmult(iel) = 0.d0
   60   continue
      INCL = .TRUE.
*     IF (JB.LE.NZERO.AND.IRFST(JB).EQ.1) INCL = .TRUE.
      IFLAG(1) = 0
      IFLAG(2) = 0
      IFLAG(3) = 0
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
    6 CONTINUE
*       .. write out the contents for the current column
      m = nrow(1)
      if (idisk .eq. 0) then
         if (nze+m .gt. big) then
            nze=ncfg
            idisk = 1
         else
            nze = nze + m
         end if
      end if
      write(iouhn) jb,m,(h(i,1),i=1,m),(jan(i,1),i=1,m)
      if (.not. skip) then
        m = nrow(2)
        write(iouhz) jb,m,(h(i,2),i=1,m),(jan(i,2),i=1,m)
        m = nrow(3)
        write(iouhs) jb,m,(h(i,3),i=1,m),(jan(i,3),i=1,m)
      end if
*
   77 FORMAT(///30X,'MULTIPLYING FACTOR',11X,'TYPE OF INTEGRAL')
      END
