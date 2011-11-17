!
!     ------------------------------------------------------------------
!       A N G M O M
!     ------------------------------------------------------------------
!
      SUBROUTINE ANGMOM(NEW, NZERO, IFIRST) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE FOUT_C 
      USE RED_C 
      use inform_C
      use debug_C
      use diagnl_C
      use fout_C
      use ndims_C
      use non30_C
      use ovrlap_C
      !use medefn_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:23:02  11/17/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE nortbpn_I 
      USE setupb_I 
      USE orthogg_I 
      USE lmatrix_I 
      IMPLICIT NONE

!-----------------------------------------------
!   MPI data
!-----------------------------------------------
! >>>>>>     <<<<<<<
      include 'mpif.h'
      integer, parameter :: MAXPROC=999
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      integer :: myid, nprocs, ierr, istat
      character*4 idstring
      character*128 startdir, permdir, tmpdir
      integer lclst,lentmp,lenstart,lenperm
! >>>>>>     <<<<<<<

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NEW 
      INTEGER  :: NZERO 
      INTEGER  :: IFIRST 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NFIRST, LET, NWF 
      CHARACTER , DIMENSION(3) :: FORMAT*30 
      CHARACTER, DIMENSION(8) :: NCHAR 
!-----------------------------------------------
!     !PARAMETER (NWD=128, NWCD=20)
!      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(4),IALL,JSC(3),ISCW
!      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
!      COMMON/DIAGNL/IDIAG,JA,JB
!      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
!     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
!      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
!     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
!     :       (QLJCLSD,LJCLSD(1))
!      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
!      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
!      POINTER(QIORTH,IORTH(1))
!      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
!     : ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
!     :     QIORTH
!
      DATA NCHAR/ '1', '2', '3', '4', '5', '6', '7', '8'/  
      DATA FORMAT/ '(2H <,I3,5H |H| ,I2,6H > = <, ', &
         '  (A3,1H(,I2,1H)) ,5H |H| ,   ', '(A3,1H(,I2,1H)),2H >,/)     '/  
!
! --- THIS PROGRAMME CONSIDERS EITHER SUPERPOSITION OF CONFIGURATIONS OR
!     MULTI-CONFIGURATIONAL HARTREE-FOCK WAVE FUNCTIONS.  USING THE
!     RESULT THAT THE TWO-ELECTRON HAMILTONIAN MATRIX ELEMENT
!     (PSI/V/PSIP)  CAN BE WRITTEN AS A SUM OF SLATER INTEGRALS, THE
!     PRESENT CODE  -  WEIGHTS  -  CALCULATES THE COEFFICIENTS OF THESE
!     INTEGRALS.  PSI AND PSIP ARE ALLOWED TO RUN OVER NCFG CONFIGURATNS
!
!
! --- CONSIDER (PSI/V/PSIP) AS PSI AND PSIP RUN OVER ALL CONFIGURATIONS
!
      NC = 0 
      NFIRST = NCFG - NEW + 1 
      DO JA = myid + NFIRST, NCFG, nprocs 
         IF (MOD(JA,1000) == 0) WRITE (ISCW, '(A,I6)') '   ja = ', JA 
         IF (JA == ncfg) write(ISCW,'(A,I5)') '   JA = ',JA
         DO JB = 1, JA 
            IFLAG = 0 
            IDIAG = 0 
            IF (NORTH /= 0) CALL NORTBPN (JA, JB) 
!          N1=NOCCSH(JA)
!          N2=NOCCSH(JB)
!          IF (IFULL .NE. 0) THEN
!          FORMAT(2)(2:2) = NCHAR(N1)
!          FORMAT(2)(30:30) = NCHAR(N2)
!          WRITE(IWRITE,'(///)')
!          WRITE(IWRITE,FORMAT) JA,JB,
!     :        (IAJCMP(NOCORB(J,JA)),NELCSH(J,JA),J=1,N1),
!     :        (IAJCMP(NOCORB(J,JB)),NELCSH(J,JB),J=1,N2)
!          END IF
!
! --- SET UP DEFINING QUANTUM NUMBERS FOR EACH MATRIX ELEMENT
!
            CALL SETUPB (JA, JB, LET) 
            IF (LET == 0) CYCLE  
!          IF(IBUG1.GT.0.OR.IBUG2.GT.0) CALL VIJOUT(JA,JB)
!
! --- TEST ON number of electrons, parity, coupling
!
            CALL ORTHOGG (LET) 
            IF (LET == 0) CYCLE  
!
! --- IF NO SUCH ORTHOGONALITY IS EXHIBITED, CALCULATE WEIGHTS OF SLATER
!     INTEGRALS
!
            CALL LMATRIX 
            IF (IFLAG == 0) CYCLE  
            NIJ = NIJ + 1 
         END DO 
      END DO 
!
      NWF = MAXORB 
      IF (NWF > 1) THEN 
         if (myid==0) &
           WRITE (6, *) ' dalloc, qiorth: nwf*(nwf-1)/2 = ', (NWF*(NWF - 1))/2 
          DEALLOCATE (IORTH) 
      ENDIF 

      if (myid==0) then
         WRITE (6, *) ' dalloc, qj1    :15*ncfg = ', 15*NCFG
         WRITE (6, *) ' dalloc, qnocorb: 8*ncfg = ', 8*NCFG
         WRITE (6, *) ' dalloc, qnelcsh: 8*ncfg = ', 8*NCFG
         WRITE (6, *) ' dalloc, qnoc: ncfg = ', NCFG
      end if

         DEALLOCATE (NOCCSH) 
         DEALLOCATE (NELCSH) 
         DEALLOCATE (NOCORB) 
         DEALLOCATE (J1QNRD) 
!      WRITE (6, *) ' dalloc, qljclcd: nwcd = ', NWCD 
!         DEALLOCATE (LJCLSD) 
!
      RETURN  
      END SUBROUTINE ANGMOM 
