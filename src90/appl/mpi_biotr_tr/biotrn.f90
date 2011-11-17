!
!     ------------------------------------------------------------------
!     B I O T R N
!     ------------------------------------------------------------------
!
      SUBROUTINE BIOTRN(PI, NLI, CII, NCSFI, NCII, IREOI, PF, NLF, CIF, &
         NCSFF, NCIF, IREOF, NGRID, MXL, NINSHL, LBUF, LUI, LUF, NTESTG) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      use inout_C
! Main routine for  employing biorthogonal rotations for RAS
! type wave functions, allowing the calculation of
! transition moments between two RAS states.
!
! Written in Bruxelles, some days in May 1993
!
! The task of this part of the code is to change
! two sets of orbitals into biorthogonal orbitals
! and counterrotate  the CI coefficients.
!
! The routine is interfaced to the MCHF program
! through the above parameter list + the subroutines
!
!               GETS  : Obtain overlap between 2 P functions
!               GTRAC1 : Read a buffer of coupling coefficients
!
! =====
! Input
! =====
!NOTE: SCR_ is not used, instead SCR is allocated further below.
! PI : Input radial shells for initial state
! NLI : Number of shells per L for initial state
! CII : CI expansion coeffcients of Initial state
!       Organized as CII(ICSF,IVEC)
! NCSFI : Number of CSF's for initial state
! NCII : Number of CI vectors for Initial vector
! IIREO : NOT ACTIVE !!!
!        Reorder array of initial  shells, from MCHF order to
!        and ordering like
!            Loop over L
!              Inactive shells
!              RAS1 Shells
!              RAS2 Shells
!              RAS3 Shells
!
! PF : Input radial shells for final   state
! NLF : Number of shells per L for final   state
! CIF : CI expansion coeffcients of final   state
!       Organized as CII(ICSF,IVEC)
! NCSFF : Number of CSF's for final   state
! NCIF : Number of CI vectors for final   vector
! IFREO : Reorder array of final shells, from MCHF order to
!        and ordering like the one given for IIREO
!
! NGRID : Number of gridpoints
! MXL   : Largest L value
! NINSHL: Number of inactive shells  per L
! LBUF  : Length of buffer to be used for accessing Racah coefs
! LUI   : Unit number for Initial state racah coefficients
! LUF   : Unit number for Final state Racah Coeficients
! LSCR : Total length of scratch space.
! NTESTG : Global print flag : = 0 => complete silence
!                              = 1 => test and print overlap matrices ,
!                              .gt. 1 => We hope you now what you are doing
!
! ======
! Output
! ======
! PI : Initial shells in biorthogonal basis
! CII : CI vectors for initials states vorresponding to
!       biorthogonal basis
! PF : Final   shells in biorthogonal basis
! CIF : CI vectors for final states corresponding to
!      biorthoginal ewxpansion
!
! NOTE : In the CURRENT version no REORDERING takes
!        place inside BIOTRN. What this means is
!        something wise men still discuss
!
! Inactive Shells : The pfunctions of the inactive shells must be
!                    be supplied to the program, and NLI,NLF
!                   must refer to the total number of
!                   occupied shells ( inactive+active)
!                   The Racah coefficients over the
!                   inactive shells should however
!                   not be supplied
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  16:54:58  11/18/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ifnmnx_I 
      USE ielsum_I 
      USE gets_I 
      USE copvec_I 
      USE wrtmat_I 
      USE invmat_I 
      USE ulla_I 
      USE trpmat_I 
      USE matml4_I 
      USE scalve_I 
      USE setvec_I 
      USE pamtmt_I 
!      USE citra_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NCSFI 
      INTEGER  :: NCII 
      INTEGER  :: IREOI 
      INTEGER  :: NCSFF 
      INTEGER  :: NCIF 
      INTEGER  :: IREOF 
      INTEGER  :: NGRID 
      INTEGER , INTENT(IN) :: MXL 
      INTEGER  :: LBUF 
      INTEGER  :: LUI 
      INTEGER  :: LUF 
      INTEGER , INTENT(IN) :: NTESTG 
      INTEGER  :: NLI(MXL) 
      INTEGER  :: NLF(MXL) 
      INTEGER , INTENT(IN) :: NINSHL(MXL) 
      REAL(DOUBLE)  :: PI(NGRID,*) 
      REAL(DOUBLE)  :: CII(NCSFI,NCII) 
      REAL(DOUBLE)  :: PF(NGRID,*) 
      REAL(DOUBLE)  :: CIF(NCSFF,NCIF) 
      REAL(DOUBLE), dimension(:), allocatable, target :: SCR
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NTESTL, NTEST, ILI, ILF, NLIMX, NLFMX, NLIFMX, &
         NLTI, NLTF, KFREE, KLPI, KLPF, KLPIF, KLSTOT, KLSIF, KLSIFI, &
         KLCI, KLCF, KLSCR, KBUF, KIBUF, KLCISC, KLCIF, L, NI, NF, NIFMN, &
         III, JJJ, KLPMX, KLPMN, NMX, NMN, KLCMX, KLCMN, NDIFF, J, KLTI, &
         I, NINL, KLTF, ierr, LSCR 
      REAL(DOUBLE) :: TII, TIII 
!-----------------------------------------------
!
!
      NTEST = 0;
      NTESTL = 0;
      call biotrn_mem(NLI, NCSFI, NCII, NLF, NCSFF, NCIF, NGRID, MXL, &
         LBUF, NTESTG,KFREE,KLPI,KLPF,KLPIF,KLSTOT,KLSIF,KLSIFI, &
         KLCI,KLCF,KLSCR,KBUF,KIBUF,KLCISC,KLCIF,NLIMX,NLFMX,NLIFMX, &
         NLTI,NLTF)

! allocate memory for the scratch arrray
      write(iscw,'(/A,I8/)') ' Allocating scratch array, n = ', &
        kfree
      allocate(scr(kfree),STAT=ierr);
      if (ierr.ne.0) call mem_fail(iscw,kfree,'biotrn::SCR',ierr);

      LSCR = KFREE
!. Obtain overlap matrix
      CALL GETS (SCR(KLSTOT), NLTI, NLTF) 
 
      DO L = 0, MXL 
         WRITE (ISCW, *) '   L = ', L 
         WRITE (ISCW, *) '       Orbital rotation...' 
         WRITE (6, *) '   L = ', L 
         WRITE (6, *) '       Orbital rotation...' 
         IF (NTEST >= 1) THEN 
            WRITE (6, *) 
            WRITE (6, *) &
               ' BIOTRN : Information on transformations of shells with L =', L 
            WRITE (6, *) 
         ENDIF 
!
! =========================================================
! 1 : Obtain Biorthogonal forms of initial and final shells
! =========================================================
!
!. The overlap matrix can be written
!
!     * * * * * * *
!     *       *   *
!     *   S   * X *
!     *       *   *
!     * * * * * * *
!
! where S is the quadratic subblock.
! The basis functions corresponding to the quadratic
! subblock are made biorthogonal with an UL decomposition.
! The remaining basis
! functions are made biorthogonal by choosing  the
! transformation
!            * * *
!            *   *
!            * Y *
!            *   *
!            * * *
!            * 1 *
!            * * *
!
!
! With Y = -S-1*X
!
         NI = NLI(L+1) 
         NF = NLF(L+1) 
!
!.1.1 : Obtain shells of given L in proper order for
!       biorthogonal treatment
!
         IF (L == 0) THEN 
            ILI = 1 
            ILF = 1 
         ELSE 
            ILI = ILI + NLI(L+1-1) 
            ILF = ILF + NLF(L+1-1) 
         ENDIF 
!mrg Remember: that is the good place: nowhere else!
!mrg    if(l.ne.X) go to 1000
!mrg
! I shells
!. Keep : Should maybe be used someday
!       DO 101 I = 1, NLTI
!         IF(IREOI(I).GE. ILI. AND.
!    &    IREOI(I).LE. ILI+NI-1) THEN
!           IEFF = IREOI(I) - ILI + 1
!           CALL COPVEC(PI(1,I),
!    &           SCR(KLPI-1+(IEFF-1)*NGRID+1), NGRID)
!          IF(ntest.ge.10)
!    &     write(6,*) ' Loop 101 I IEFF ',i,IEFF
!         END IF
! 101   CONTINUE
         CALL COPVEC (PI(1,ILI), SCR(KLPI), NI*NGRID) 
! F shells
!       DO 102 I = 1, NLTF
!         IF(IREOF(I).GE. ILF. AND.
!    &    IREOF(I).LE. ILF+NF-1) THEN
!           IEFF = IREOF(I) - ILF + 1
!
!         IF(ntest.ge.10)
!    &      write(6,*) ' Loop 102 I IEFF ',i,IEFF
!           CALL COPVEC(PF(1,I),
!    &           SCR(KLPF-1+(IEFF-1)*NGRID+1), NGRID)
!         END IF
! 102   CONTINUE
         CALL COPVEC (PF(1,ILF), SCR(KLPF), NF*NGRID) 
!
! 1.2 obtain biorthogonal of the first min(ni,nf) shells
!
         NIFMN = MIN(NI,NF) 
!
!. Overlap matrix SIF = Integral (PI(I)*PF(J))
!
         DO III = 1, NIFMN 
            SCR(KLSIF+III-1:NIFMN*(NIFMN-1)+KLSIF+III-1:NIFMN) = SCR(NLTI*(ILF-&
               1)+KLSTOT-2+III+ILI:(NIFMN-2+ILF)*NLTI+KLSTOT-2+III+ILI:NLTI) 
         END DO 

         IF (NTEST >= 15) THEN 
            WRITE (6, *) ' Overlap matrix ' 
            CALL WRTMAT (SCR(KLSIF), NIFMN, NIFMN, NIFMN, NIFMN) 
         ENDIF 
!
! Obtain upper triangular CI and CF so CI(T) S CF = 1
! or CF CI(T) = S-1, which corresponds to an UL decomposition
!
!. Invert S
!
!             INVMAT(A,B,MATDIM,NDIM)
         CALL COPVEC (SCR(KLSIF), SCR(KLSIFI), NIFMN**2) 
         CALL INVMAT (SCR(KLSIFI), SCR(KLCI), NIFMN, NIFMN) 
 
!. UL decompose
         CALL COPVEC (SCR(KLSIFI), SCR(KLSIF), NIFMN**2) 
         CALL ULLA (SCR(KLSIF), SCR(KLCF), SCR(KLCI), NIFMN, SCR(KLSCR)) 
         CALL TRPMAT (SCR(KLCI), NIFMN, NIFMN, SCR(KLSCR)) 
         CALL COPVEC (SCR(KLSCR), SCR(KLCI), NIFMN**2) 
!
!. The transformation matrix between the first NIFMX
!. shells is now known, biorthogonalize remaining orbitals
         IF (NI/=NF .AND. NI/=0 .AND. NF/=0) THEN 
            IF (NI > NF) THEN 
               KLPMX = KLPI 
               KLPMN = KLPF 
               NMX = NI 
               NMN = NF 
               KLCMX = KLCI 
               KLCMN = KLCF 
            ELSE 
               KLPMX = KLPF 
               KLPMN = KLPI 
               NMX = NF 
               NMN = NI 
               KLCMX = KLCF 
               KLCMN = KLCI 
            ENDIF 
            NDIFF = NMX - NMN 
! Y = -S-1 * X
!. overlap X between remaining orbitals and the other set
            IF (NI > NF) THEN 
! I columns F rows
               DO III = NMN + 1, NMX 
                  SCR(KLSIF+(III-NMN-1)*NF:KLSIF+(III-NMN)*NF-1) = SCR(NLTI*(&
                     ILF-1)+KLSTOT-2+III+ILI:(NF-2+ILF)*NLTI+KLSTOT-2+III+ILI:&
                     NLTI) 
               END DO 
            ELSE IF (NF > NI) THEN 
! F columns I rows
               DO JJJ = NMN + 1, NMX 
                  SCR(KLSIF+(JJJ-NMN-1)*NI:KLSIF+(JJJ-NMN)*NI-1) = SCR(KLSTOT-1&
                     +(JJJ+ILF-2)*NLTI+ILI:NI+KLSTOT+(JJJ+ILF)*NLTI-2*(NLTI+1)+&
                     ILI) 
               END DO 
            ENDIF 
!
            IF (NI > NF) THEN 
               CALL TRPMAT (SCR(KLSIFI), NMN, NMN, SCR(KLSCR)) 
               CALL COPVEC (SCR(KLSCR), SCR(KLSIFI), NMN**2) 
            ENDIF 
            CALL MATML4 (SCR(KLSCR), SCR(KLSIFI), SCR(KLSIF), NMN, NDIFF, NMN, &
               NMN, NMN, NDIFF, 0) 
            CALL SCALVE (SCR(KLSCR), -1.0D0, NMN*NDIFF) 
            CALL COPVEC (SCR(KLSCR), SCR(KLSIF), NMN*NDIFF) 
! Construct complete CMX
            CALL SETVEC (SCR(KLSCR), 0.0D0, NMX**2) 
            DO J = 1, NMX 
               IF (J <= NIFMN) THEN 
                  CALL COPVEC (SCR(KLCMX+(J-1)*NIFMN), SCR(KLSCR+(J-1)*NMX), &
                     NMN) 
               ELSE 
                  CALL COPVEC (SCR(KLSIF+(J-NMN-1)*NMN), SCR(KLSCR+(J-1)*NMX), &
                     NMN) 
                  SCR(KLSCR-1+(J-1)*NMX+J) = 1.0D0 
               ENDIF 
            END DO 
!
            CALL COPVEC (SCR(KLSCR), SCR(KLCMX), NMX**2) 
         ENDIF 
!. The two upper triangular matrices CI and CF are now known
!. Rotate the shells
         CALL MATML4 (SCR(KLPIF), SCR(KLPI), SCR(KLCI), NGRID, NI, NGRID, NI, &
            NI, NI, 0) 
         CALL COPVEC (SCR(KLPIF), SCR(KLPI), NI*NGRID) 
 
         CALL COPVEC (SCR(KLPI), PI(1,ILI), NI*NGRID) 
         CALL MATML4 (SCR(KLPIF), SCR(KLPF), SCR(KLCF), NGRID, NF, NGRID, NF, &
            NF, NF, 0) 
         CALL COPVEC (SCR(KLPIF), SCR(KLPF), NF*NGRID) 
         CALL COPVEC (SCR(KLPF), PF(1,ILF), NF*NGRID) 
!
         IF (NTEST >= 1) THEN 
            WRITE (6, *) ' Test of overlap of biorthonormal functions' 
! F columns I rows
            DO JJJ = 1, NF 
               SCR(KLSIF+(JJJ-1)*NI:KLSIF+JJJ*NI-1) = SCR(KLSTOT-1+(JJJ+ILF-2)*&
                  NLTI+ILI:NI+KLSTOT+(JJJ+ILF)*NLTI-2*(NLTI+1)+ILI) 
            END DO 
            CALL MATML4 (SCR(KLSCR), SCR(KLCI), SCR(KLSIF), NI, NF, NI, NI, NI&
               , NF, 1) 
            CALL MATML4 (SCR(KLSIF), SCR(KLSCR), SCR(KLCF), NI, NF, NI, NF, NF&
               , NF, 0) 
            WRITE (6, *) &
               ' new overlap matrix ( should be 1 on diagonal, 0 elsewhere )' 
            CALL WRTMAT (SCR(KLSIF), NI, NF, NI, NF) 
         ENDIF 
 
         IF (NTEST >= 1) THEN 
            WRITE (6, *) 
            WRITE (6, *) ' Orbital Rotation matrix for I state' 
            CALL WRTMAT (SCR(KLCI), NI, NI, NI, NI) 
            WRITE (6, *) ' Orbital Rotation matrix for F state' 
            CALL WRTMAT (SCR(KLCF), NF, NF, NF, NF) 
            WRITE (6, *) 
         ENDIF 
!
! ======================================
! 2. Counter rotate the CI coefficients
! ======================================
!
!. Initial state
!. CI rotation parameters
         WRITE (ISCW, *) '       CI counter-rotation...' 
         KLTI = KLSIF 
         CALL PAMTMT (SCR(KLCI), SCR(KLTI), SCR(KLSCR), NI) 
         DO I = 1, NI 
            TII = SCR(KLTI-1+(I-1)*NI+I) 
            TIII = 1.0D0/TII 
            CALL SCALVE (SCR(KLTI+(I-1)*NI), TIII, I - 1) 
         END DO 
         NINL = NINSHL(L+1) 
         CALL CITRA (CII, NCSFI, NCII, L, NI, SCR(KLTI), LBUF, LUI,  &
            NINL, SCR(KLCIF), SCR(KLCISC), SCR(KBUF), SCR(KIBUF), NTESTL) 
         CALL COPVEC (SCR(KLCIF), CII, NCII*NCSFI) 
!. Final state
!. CI rotation parameters
         KLTF = KLSIF 
         CALL PAMTMT (SCR(KLCF), SCR(KLTF), SCR(KLSCR), NF) 
         DO I = 1, NF 
            TII = SCR(KLTF-1+(I-1)*NF+I) 
            TIII = 1.0D0/TII 
            CALL SCALVE (SCR(KLTF+(I-1)*NF), TIII, I - 1) 
         END DO 
!
         CALL CITRA (CIF, NCSFF, NCIF, L, NF, SCR(KLTF), LBUF, LUF, &
            NINL, SCR(KLCIF), SCR(KLCISC), SCR(KBUF), SCR(KIBUF), NTESTL) 
         CALL COPVEC (SCR(KLCIF), CIF, NCIF*NCSFF) 
!. End of loop over L
      END DO 
!
      RETURN  
      END SUBROUTINE BIOTRN 
