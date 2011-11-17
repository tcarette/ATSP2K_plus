!
!     ------------------------------------------------------------------
!     C F G I N 2
!     ------------------------------------------------------------------
!
      SUBROUTINE CFGIN2(MCFG, KCFG) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      use parameters_biotr_C
      !USE DEBUG_C 
      USE DBG_C 
      USE INOUT_C 
      USE NOR_C 
      !USE FO_C 
      USE RAS_C 
      use ndims_C
      use non30_C
!
! --- Read two sets of configurations and determine the orthogonality
! --- conditions between them
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:37:01  11/17/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE analy2_I 
      USE lval_I 
      USE gstate_I 
      USE cfgtst_I 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: MCFG 
      INTEGER  :: KCFG 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, ICLOSDI, IWFI, ICLOSDF, IWFF, ierr 
      CHARACTER , DIMENSION(NWD) :: EL*3, ELC*3 
      CHARACTER , DIMENSION(NWD,3) :: JAJCMP 
      CHARACTER , DIMENSION(2) :: INPUT*24 
      CHARACTER :: HEADI*72, HEADF*72, HEADER*72 
      CHARACTER, DIMENSION(2) :: LABEL*7 
!-----------------------------------------------
      DATA LABEL/ 'Initial', 'Final  '/  
!      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
!     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
!      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
!     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
!     :       (QLJCLSD,LJCLSD(1))
!      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
!      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
!
    3 FORMAT(20(1X,A3)) 
    7 FORMAT(A30,I3,I4) 
   22 FORMAT(/,/,' STATE ',' (WITH',I5,' CONFIGURATIONS):'/,' ',31('-'),/) 
   23 FORMAT(/,' THERE ARE',I3,' ORBITALS AS FOLLOWS:'/,/,5X,21(1X,A3),:,/,5X,&
         21(1X,A3)) 
      CALL ANALY2 (MCFG, KCFG, EL) 
!
      REWIND (UNIT=IUC(1)) 
      REWIND (UNIT=IUC(2)) 
!
      MAXORB = NCOM + NORBI + NORBF 
      NCFG = MCFG + KCFG 
      IF (IBUGM /= 0) THEN 
!
! --- allocate the memory
!
         WRITE (6, *) ' qiajcmp allocation: maxorb   = ', MAXORB 
         WRITE (6, *) ' qljcomp allocation: maxorb   = ', MAXORB 
         WRITE (6, *) ' qnjcomp allocation: maxorb   = ', MAXORB 
         WRITE (6, *) ' qnoc    allocation: ncfg     = ', NCFG 
         WRITE (6, *) ' qnelcsh allocation: 8*ncfg   = ', 8*NCFG 
         WRITE (6, *) ' qnocorb allocation: 8*ncfg   = ', 8*NCFG 
         WRITE (6, *) ' qj1     allocation: 15*ncfg   = ', 15*NCFG 
      ENDIF 
 
    !  if (.not.allocated(IAJCMP)) ALLOCATE (IAJCMP(MAXORB),stat=ierr)
      allocate (IAJCMP(MAXORB),stat=ierr)
!      if (ierr.ne.0) call mem_fail(6,MAXORB,'cfgin2::IAJCMP',ierr);
      if (.not.allocated(LJCOMP)) ALLOCATE (LJCOMP(MAXORB),stat=ierr) 
!      if (ierr.ne.0) call mem_fail(6,MAXORB,'cfgin2::LJCOMP',ierr);
      if (.not.allocated(NJCOMP)) ALLOCATE (NJCOMP(MAXORB),stat=ierr) 
!      if (ierr.ne.0) call mem_fail(6,MAXORB,'cfgin2::NJCOMP',ierr);
      if (.not.allocated(NOCCSH)) ALLOCATE (NOCCSH(NCFG),stat=ierr) 
!      if (ierr.ne.0) call mem_fail(6,NCFG,'cfgin2::NOCCSH',ierr);
      if (.not.allocated(NELCSH)) ALLOCATE (NELCSH(8,NCFG),stat=ierr) 
!      if (ierr.ne.0) call mem_fail(6,8*NCFG,'cfgin2::NELCSH',ierr);
      if (.not.allocated(NOCORB)) ALLOCATE (NOCORB(8,NCFG),stat=ierr) 
!      if (ierr.ne.0) call mem_fail(6,8*NCFG,'cfgin2::NOCORB',ierr);
      if (.not.allocated(J1QNRD)) ALLOCATE (J1QNRD(15,NCFG),stat=ierr) 
!      if (ierr.ne.0) call mem_fail(6,15*NCFG,'cfgin2::J1QNRD',ierr);
      qljcomp=>LJCOMP;
      qiajcmp=>IAJCMP;
      qnjcomp=>NJCOMP;
      qnoc=>NOCCSH;
      qnelcsh=>NELCSH;
      qnocorb=>NOCORB;
      qj1=>J1QNRD;
  
 
!
! --- set up the electrons
!
      READ (EL, '(A3)') (IAJCMP(I),I=1,MAXORB) 
      READ (EL, '(3A1)') ((JAJCMP(I,J),J=1,3),I=1,MAXORB) 
!
! --- set up of ljcomp
!
      DO I = 1, MAXORB 
         IF (JAJCMP(I,1) == ' ') THEN 
            JAJCMP(I,1) = JAJCMP(I,2) 
            JAJCMP(I,2) = JAJCMP(I,3) 
            JAJCMP(I,3) = ' ' 
         ENDIF 
         LJCOMP(I) = LVAL(JAJCMP(I,2)) 
         NJCOMP(I) = ICHAR(JAJCMP(I,1)) - ICHAR('1') + 1 
      END DO 
!
! ---- check common closed shells
!
      IF (NCLOSI /= NCLOSF) STOP &
         ' Common closed shells not the same in the two states' 
!
      READ (IUC(1), 7) HEADI, ICLOSDI, IWFI 
      READ (IUC(2), 7) HEADF, ICLOSDF, IWFF 
      HEADER = HEADI(1:30)//'=>'//HEADF(1:30) 
!
! --- check closed shells further
!
      READ (IUC(1), 3) (ELC(I),I=1,NCLOSI) 
      READ (IUC(2), 3) (EL(I),I=1,NCLOSF) 
      DO I = 1, NCLOSF 
         J = 1 
    2    CONTINUE 
         IF (EL(I) == ELC(J)) CYCLE  
         J = J + 1 
         IF (J <= NCLOSI) THEN 
            GO TO 2 
         ELSE 
            STOP ' Common closed sub-shells not the same' 
         ENDIF 
      END DO 
!
! --- Check if electron lists are present
!
      IF (IWFI > ICLOSDI) READ (IUC(1), 3) (ELC(I),I=ICLOSDI + 1,IWFI) 
      IF (IWFF > ICLOSDF) READ (IUC(2), 3) (EL(I),I=ICLOSDF + 1,IWFF) 
!
! --- get initial state configurations
!
      CALL GSTATE (1, MCFG) 
      CALL GSTATE (MCFG + 1, NCFG) 
!
!  --- check the data
!
      !CALL CFGTST(NCFG,QLJCOMP,QNOC,QNELCSH,QNOCORB,QJ1)
      CALL CFGTST 
      REWIND (UNIT=IUC(1)) 
      REWIND (UNIT=IUC(2)) 
      RETURN  
      END SUBROUTINE CFGIN2 
