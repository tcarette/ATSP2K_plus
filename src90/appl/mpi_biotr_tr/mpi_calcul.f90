!
!     -------------------------------------------------------------
!       C A L C U L
!     -------------------------------------------------------------
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville           September 1997   *
!
      SUBROUTINE CALCUL(NPAIR) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      use debug_C
      use inout_C
      use medefn_C
      use ems_C
      use ovrlap_C
      use nor_C
      use state_C
      use mult_C
      use consts_C
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:27:03  11/20/01
!...Switches:
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE nortbp_I
      USE setup_oc_I
      USE ittk_I
      USE nontrans_I
      USE fline_I
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

    5 FORMAT(/,/,' JI =',I4,'  JF =',I4) 
    6 FORMAT(/,' ','NORTH = ',I3,/,' JMU = ',I3,2X,'JNU = ',I3,/,' JMUP= ',I3,&
         2X,'JNUP= ',I3,/,' NOVLPS = ',I3,/) 
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NPAIR 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KA, KB, IREZ, IIFLL, JII, JFI, NCONTR, JFF, NC, &
                 LET, N1, ITIK, K, KL, KR 
      REAL(DOUBLE) :: CL2, CV2, FWW 
!-----------------------------------------------
!      IBUG1 = 1
!      NBUG6 = 1
      KA = LAM 
      KB = 0 
      IREZ = 1 
      IF (IFL == 3) IREZ = 2 
      WRITE (ISCW, *) ' Double sum over CSF ...' 
      DO IIFLL = 1, IREZ 
         IF (IIFLL == 2) THEN 
            IFL = 4 
            KA = LAM - 1 
            KB = 1 
         ENDIF 
         DO JII = myid + 1, NCF(1), nprocs 
            JI = JII 
            IF (MOD(JI,1000) == 0) WRITE (ISCW, *) '     JA = ', JI 
            IF (JI==NCF(1)) write(ISCW,'(A,I5)') '   JA = ',JI
            DO JFI = 1, NCF(2) 
               JF = JFI 
               IF (NPAIR == 0) THEN 
                  WRITE (6, *) ' ji = ', JI, ' jf = ', JF 
                  STOP  
               ENDIF 
               NCONTR = 0 
               IF (NBUG6 /= 0) WRITE (IWRITE, 5) JI, JF 
               JFF = JF + NCF(1) 
               NOVLPS = 0 
               JMUP = 0 
               JNUP = 0 
               JMU = 0 
               JNU = 0 
!
               NC = 0 
!
               IF (NORTH /= 0) CALL NORTBP (JI, JFF) 
               IF (IWAR == 1) CYCLE  
               IF (NBUG6 /= 0) WRITE (IWRITE, 6) NORTH, JMU, JNU, JMUP, JNUP, &
                  NOVLPS 
!
! --- set up the occupation and coupling arrays
!
               CALL SETUP_OC (JI, JFF, LET) 
               IF (LET == 0) CYCLE  
!
! --- test selection rules delta(S)=0 for Ek and M1
!                          delta(L)=0 for M1
               N1 = 2*IHSH - 1 
               ITIK = 1 
               IF (ITTK(J1QN1(N1,2)-1,J1QN2(N1,2)-1,2*KA)==0.OR. &
                  ITTK(J1QN1(N1,3)-1,J1QN2(N1,3)-1,2*KB)==0) ITIK = 0 
               IF (ITIK == 0) CYCLE  
               CALL NONTRANS (KA, KB, CL2, CV2) 
              !print*,KA,KB,CL2,CV2, ' nontrans'
!
! --- calculate the contribution of <JI/ O /JF> to the line
!     strengths for the npair (J,J') found
               IF (DABS(CL2)<=ZERO .AND. DABS(CV2)<=ZERO) CYCLE  
               DO K = 1, NPAIR 
                  KL = (IL(K)-1)*NCF(1) + JI 
                  KR = (IR(K)-1)*NCF(2) + JF 
                  FWW = FLINE(K)*WT1(KL)*WT2(KR) 
                  SL(K) = SL(K) + CL2*FWW 
                  IF (VOK) SV(K) = SV(K) + CV2*FWW 
                  IF (IBUG1 == 0) CYCLE  
                  WRITE (6, *) ' pair = ', K, ' wt1 = ', WT1(KL), ' wt2 = ', &
                     WT2(KR) 
                  WRITE (6, *) '              sl(pair) = ', SL(K) 
               END DO 
            END DO 
         END DO 
      END DO 
      call mpi_allr_dp(SL,nvc(1)*nvc(2));
      if (VOK) call mpi_allr_dp(SV,nvc(1)*nvc(2));
      if (myid == 0) then
               write(6,'(5F12.8)') SL
               write(6,'(5F12.8)') sv
      end if

      RETURN  
      END SUBROUTINE CALCUL 
