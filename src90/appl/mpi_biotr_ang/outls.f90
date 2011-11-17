!
!     ------------------------------------------------------------------
!       O U T L S
!     ------------------------------------------------------------------
!
      SUBROUTINE OUTLS(INTGRL, NINT, INTPTR, CNN, NCOEF, JANN, JBNN) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CLOSED_C 
      USE DEBUG_C 
      USE FOUT_C 
      USE INFORM_C 
!
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:57:27  11/17/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE qsort_I 
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
      INTEGER , INTENT(OUT) :: NINT 
      INTEGER , INTENT(OUT) :: NCOEF 
      INTEGER , INTENT(OUT) :: INTGRL(*) 
      INTEGER , INTENT(OUT) :: INTPTR(*) 
      INTEGER , INTENT(OUT) :: JANN(*) 
      INTEGER , INTENT(OUT) :: JBNN(*) 
      REAL(DOUBLE) , INTENT(OUT) :: CNN(*) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: LSDIM = 30000 
      INTEGER, PARAMETER :: LSTACK = 128 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , POINTER, DIMENSION(:) :: QPACKN, QJAN, QJBN, QIPT 
      INTEGER , ALLOCATABLE, TARGET, DIMENSION(:) :: IPACKN, JBN, JAN, IPT 
      INTEGER , DIMENSION(LSTACK) :: ISTACK 
      INTEGER :: NCC, IOU, II, IST, JJ, J, M, IERR, ICASE, NC, LAST, K
      REAL(DOUBLE), POINTER, DIMENSION(:) :: QCN 
      REAL(DOUBLE), ALLOCATABLE, TARGET, DIMENSION(:) :: CN 
      CHARACTER, DIMENSION(8) :: INT 
!-----------------------------------------------
!
!     The list of unique integrals L(j,i) is formed, in the order of
!     increasing symmetry, with j .le. i.  With each integral (INTGRL)
!     there is a pointer (INTPTR) that points to the last element in the
!     array of coefficients for that integral.
!
!
!
!
!
!     ..local arrays
!      POINTER (qcn,cn(1)),(qpack,ipackn(1)),(qjan,jan(1)),
!     :        (qjbn,jbn(1)),(qipt,ipt(1))
!
      DATA INT/ 'F', 'G', 'R', 'R', 'R', 'R', 'R', 'L'/  
!
!===  Begin processing data
!
      NCC = LCOUNT + NREC*LSDIM 

      IOU = ISC(4) 
 
      if (myid == 0) WRITE (ISCW, '(A/(4I15))') &
              'Number of L-integrals', NCC 
 
      ALLOCATE (CN(NCC),STAT=ierr) 
      if (ierr.ne.0) call mem_fail(6,ncc,'outls::cn',ierr);
      ALLOCATE (IPACKN(NCC),STAT=ierr) 
      if (ierr.ne.0) call mem_fail(6,ncc,'outls::IPACKN',ierr);
      ALLOCATE (JAN(NCC),STAT=ierr) 
      if (ierr.ne.0) call mem_fail(6,ncc,'outls::JAN',ierr);
      ALLOCATE (JBN(NCC),STAT=ierr) 
      if (ierr.ne.0) call mem_fail(6,ncc,'outls::JBN',ierr);
      ALLOCATE (IPT(NCC),STAT=ierr) 
      if (ierr.ne.0) call mem_fail(6,ncc,'outls::IPT',ierr);
      qcn=>cn;
      qpackn=>ipackn;
      qjan=>jan;
      qjbn=>jbn;
      qipt=>ipt;
 
!
      II = 0 
!
      NINT = 0 
!        .. read data
      IST = 0 
      DO JJ = 1, NREC 
         READ (IOU) (CN(J+IST),J=1,LSDIM), (IPACKN(J+IST),J=1,LSDIM), (JAN(J+&
            IST),J=1,LSDIM), (JBN(J+IST),J=1,LSDIM) 
         IST = IST + LSDIM 
      END DO 
      M = LCOUNT 
      IF (M /= 0) READ (IOU) (CN(J+IST),J=1,M), (IPACKN(J+IST),J=1,M), (JAN(J+&
         IST),J=1,M), (JBN(J+IST),J=1,M) 
      CLOSE(IOU) 
!
      CALL QSORT (NCC, IPACKN, IPT, ISTACK, LSTACK, IERR) 
      IF (IERR == 1) THEN 
         WRITE (ISCW, *) ' Stack dimension not large enough for sort', &
            'CASE = ', ICASE, 'N = ', NC 
         CALL EXIT (1) 
      ENDIF 
!
!        Form the list of integrals with pointers to the data
!
      LAST = 0 
  110 CONTINUE 
      J = LAST + 1 
      LAST = J 
      IF (J <= NCC) THEN 
!
!           Find  last item in the list with this integral
!
  120    CONTINUE 
         LAST = LAST + 1 
         IF (LAST <= NCC) THEN 
            IF (IPACKN(IPT(J)) == IPACKN(IPT(LAST))) GO TO 120 
         ENDIF 
         LAST = LAST - 1 
         NINT = NINT + 1 
         INTGRL(NINT) = IPACKN(IPT(LAST)) 
         INTPTR(NINT) = LAST 
         GO TO 110 
      ENDIF 
      NCOEF = NCC 
!
!       .. transfer the sorted data
!
      CNN(:NCOEF) = CN(IPT(:NCOEF)) 
      JANN(:NCOEF) = JAN(IPT(:NCOEF)) 
      JBNN(:NCOEF) = JBN(IPT(:NCOEF)) 


!
!      .. deallocate the local arrays.
 
      DEALLOCATE (CN) 
      DEALLOCATE (IPACKN) 
      DEALLOCATE (JAN) 
      DEALLOCATE (JBN) 
      DEALLOCATE (IPT) 
      RETURN  
 
      END SUBROUTINE OUTLS 
