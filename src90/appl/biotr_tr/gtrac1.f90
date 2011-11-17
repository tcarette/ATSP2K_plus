!
!     ------------------------------------------------------------------
!     G T R A C 1
!     ------------------------------------------------------------------
!
      SUBROUTINE GTRAC1(I, L, IFIRST, NFOUND, BUF, IBUF, LBUF, JLEI, LU, NTESTG&
         ) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:21:18  11/20/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!      USE gtracxvn_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I 
      INTEGER  :: L 
      INTEGER  :: IFIRST 
      INTEGER  :: NFOUND 
      INTEGER  :: IBUF 
      INTEGER  :: LBUF 
      INTEGER  :: JLEI 
      INTEGER , INTENT(IN) :: LU 
      INTEGER  :: NTESTG 
      REAL(DOUBLE)  :: BUF 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ILIST 
!-----------------------------------------------
!
! Obtain racah coefficient for given I and l starting
! from element IFIRST
!
! JLEI = 0 => include all terms with given I
! JLEI = 1 => include only terms with J<I
!
!
!     write(6,*) ' Old GTRACX '
!     CALL  GTRACX(I,L,IFIRST,NFOUND,BUF,
!    &                IBUF,LBUF,JLEI,LU,NTESTG)
!
!     write(6,*) ' new gtracx '
!     CALL  GTRACXN(I,L,IFIRST,NFOUND,BUF,
!    &                IBUF,LBUF,JLEI,LU,NTESTG)
! For new version where the coupling coefficients reside in core,
! the program must know which of the two lists that should be read.
! This is done by assuming that LU is as ' in the bad old days'
! i.e
!         LU = 15 => Left (initial state)
!         LU = 16 => right (final state)
      IF (LU == 15) THEN 
         ILIST = 1 
      ELSE IF (LU == 16) THEN 
         ILIST = 2 
      ELSE 
         WRITE (6, *) ' PROBLEM in GTRAC1, Unrecognized LU = ', LU 
         STOP 'GTRAC1: Undefined LU' 
      ENDIF 
!
      CALL GTRACXVN (I, L, IFIRST, NFOUND, BUF, IBUF, LBUF, JLEI, &
                     ILIST, NTESTG) 
!
      RETURN  
      END SUBROUTINE GTRAC1 
