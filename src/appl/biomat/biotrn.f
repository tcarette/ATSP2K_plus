*
*     ------------------------------------------------------------------
*	B I O T R N
*     ------------------------------------------------------------------
*
      SUBROUTINE BIOTRN(PI,NLI,CII,NCSFI,NCII,IREOI,
     &                  PF,NLF,CIF,NCSFF,NCIF,IREOF,
     &                  NGRID,MXL,NINSHL,
     &                  LBUF,LUI,LUF,
     &                  NTESTG)
*
* Main routine for  employing biorthogonal rotations for RAS
* type wave functions, allowing the calculation of
* transition moments between two RAS states.
*
* Written in Bruxelles, some days in May 1993
*
* The task of this part of the code is to change
* two sets of orbitals into biorthogonal orbitals
* and counterrotate  the CI coefficients.
*
* The routine is interfaced to the MCHF program
* through the above parameter list + the subroutines
*
*               GETS  : Obtain overlap between 2 P functions
*               GTRAC1 : Read a buffer of coupling coefficients
*
* =====
* Input
* =====
*
* PI : Input radial shells for initial state
* NLI : Number of shells per L for initial state
* CII : CI expansion coeffcients of Initial state
*       Organized as CII(ICSF,IVEC)
* NCSFI : Number of CSF's for initial state
* NCII : Number of CI vectors for Initial vector
* IIREO : NOT ACTIVE !!!
*        Reorder array of initial  shells, from MCHF order to
*        and ordering like
*            Loop over L
*              Inactive shells
*              RAS1 Shells
*              RAS2 Shells
*              RAS3 Shells
*
* PF : Input radial shells for final   state
* NLF : Number of shells per L for final   state
* CIF : CI expansion coeffcients of final   state
*       Organized as CII(ICSF,IVEC)
* NCSFF : Number of CSF's for final   state
* NCIF : Number of CI vectors for final   vector
* IFREO : Reorder array of final shells, from MCHF order to
*        and ordering like the one given for IIREO
*
* NGRID : Number of gridpoints
* MXL   : Largest L value
* NINSHL: Number of inactive shells  per L
* LBUF  : Length of buffer to be used for accessing Racah coefs
* LUI   : Unit number for Initial state racah coefficients
* LUF   : Unit number for Final state Racah Coeficients
* LSCR : Total length of scratch space.
* NTESTG : Global print flag : = 0 => complete silence
*                              = 1 => test and print overlap matrices , print header and  t matrix
*                              .gt. 1 => We hope you now what you are doing
*
* ======
* Output
* ======
* PI : Initial shells in biorthogonal basis
* CII : CI vectors for initials states vorresponding to
*       biorthogonal basis
* PF : Final   shells in biorthogonal basis
* CIF : CI vectors for final states corresponding to
*      biorthoginal ewxpansion
*
* NOTE : In the CURRENT version no REORDERING takes
*        place inside BIOTRN. What this means is 
*        something wise men still discuss
*
* Inactive Shells : The pfunctions of the inactive shells must be
*                    be supplied to the program, and NLI,NLF
*                   must refer to the total number of
*                   occupied shells ( inactive+active)
*                   The Racah coefficients over the
*                   inactive shells should however
*                   not be supplied
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INTEGER LSCR
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2),iut(2)
      DIMENSION PI(NGRID,*),NLI(MXL),CII(NCSFI,NCII)
      DIMENSION PF(NGRID,*),NLF(MXL),CIF(NCSFF,NCIF)
C     DIMENSION IREOI(*),IREOF(*)
      DIMENSION NINSHL(MXL)
*. Scratch

      double precision, dimension(:), allocatable,target :: scr
*

      call biotrn_mem(NLI, NCSFI, NCII, NLF, NCSFF, NCIF, NGRID, MXL, 
     :    LBUF, NTESTG,KFREE,KLPI,KLPF,KLPIF,KLSTOT,KLSIF,KLSIFI,
     :    KLCI,KLCF,KLSCR,KBUF,KIBUF,KLCISC,KLCIF,NLIMX,NLFMX,NLIFMX,
     : NLTI,NLTF)

      LSCR = KFREE

! allocate memory for the scratch arrray
      write(iscw,'(/A,I8/)') ' Allocating scratch array, n = ',
     :   kfree
      allocate(scr(kfree));

*. Obtain overlap matrix
      CALL GETS(SCR(KLSTOT),NLTI,NLTF)

      DO 1000 L = 0, MXL
        write(iscw,*) '   L = ',L  
        write(iscw,*) '       Orbital rotation...'
        write(6,*) '   L = ',L  
        write(6,*) '       Orbital rotation...'
        IF(NTEST.GE.1) WRITE(6,*) 
        IF(NTEST.GE.1) WRITE(6,*) 
     &  ' BIOTRN : Information on transformations of shells with L =',L
        IF(NTEST.GE.1) WRITE(6,*) 
*
* =========================================================
* 1 : Obtain Biorthogonal forms of initial and final shells
* =========================================================
*
*. The overlap matrix can be written
*
*     * * * * * * *
*     *       *   *
*     *   S   * X *
*     *       *   *
*     * * * * * * *
*
* where S is the quadratic subblock.
* The basis functions corresponding to the quadratic
* subblock are made biorthogonal with an UL decomposition.
* The remaining basis
* functions are made biorthogonal by choosing  the
* transformation
*            * * *
*            *   *
*            * Y *
*            *   *
*            * * *
*            * 1 *
*            * * *
*
*
* With Y = -S-1*X
*
        NI = NLI(L+1)
        NF = NLF(L+1)
*
*.1.1 : Obtain shells of given L in proper order for
*       biorthogonal treatment
*
        IF(L .EQ. 0 ) THEN
          ILI = 1
          ILF = 1
         ELSE
          ILI = ILI + NLI(L+1-1)
          ILF = ILF + NLF(L+1-1)
        END IF
Cmrg Remember: that is the good place: nowhere else!
Cmrg    if(l.ne.X) go to 1000
Cmrg
* I shells
C. Keep : Should maybe be used someday
C       DO 101 I = 1, NLTI
C         IF(IREOI(I).GE. ILI. AND.
C    &    IREOI(I).LE. ILI+NI-1) THEN
C           IEFF = IREOI(I) - ILI + 1
C           CALL COPVEC(PI(1,I),
C    &           SCR(KLPI-1+(IEFF-1)*NGRID+1), NGRID)
C          IF(ntest.ge.10)
C    &     write(6,*) ' Loop 101 I IEFF ',i,IEFF
C         END IF
C 101   CONTINUE
        CALL COPVEC(PI(1,ILI),SCR(KLPI),NI*NGRID)
* F shells
C       DO 102 I = 1, NLTF
C         IF(IREOF(I).GE. ILF. AND.
C    &    IREOF(I).LE. ILF+NF-1) THEN
C           IEFF = IREOF(I) - ILF + 1
C
C         IF(ntest.ge.10)
C    &      write(6,*) ' Loop 102 I IEFF ',i,IEFF
C           CALL COPVEC(PF(1,I),
C    &           SCR(KLPF-1+(IEFF-1)*NGRID+1), NGRID)
C         END IF
C 102   CONTINUE
        CALL COPVEC(PF(1,ILF),SCR(KLPF),NF*NGRID)
*
* 1.2 obtain biorthogonal of the first min(ni,nf) shells
*
        NIFMN = MIN(NI,NF)
*
*. Overlap matrix SIF = Integral (PI(I)*PF(J))
*
       DO 51 III = 1, NIFMN
       DO 51 JJJ = 1, NIFMN
        SCR(KLSIF+(JJJ-1)*NIFMN+III-1) = 
     &  SCR(KLSTOT-1+(JJJ+ILF-1-1)*NLTI+III+ILI-1)
   51  CONTINUE
       IF(NTEST.GE.15) THEN
         WRITE(6,*) ' Overlap matrix '
         CALL wrtmat(SCR(KLSIF),NIFMN,NIFMN,NIFMN,NIFMN)
       END IF
*
* Obtain upper triangular CI and CF so CI(T) S CF = 1
* or CF CI(T) = S-1, which corresponds to an UL decomposition
*
*. Invert S
*
C             INVMAT(A,B,MATDIM,NDIM)
         CALL COPVEC(SCR(KLSIF),SCR(KLSIFI),NIFMN**2)
         Call INVMAT(SCR(KLSIFI),SCR(KLCI),NIFMN,NIFMN)

*. UL decompose
         CALL COPVEC(SCR(KLSIFI),SCR(KLSIF),NIFMN**2)
         CALL ULLA(SCR(KLSIF),SCR(KLCF),SCR(KLCI),
     &             NIFMN,SCR(KLSCR))
         CALL TRPMAT(SCR(KLCI),NIFMN,NIFMN,SCR(KLSCR))
         CALL COPVEC(SCR(KLSCR),SCR(KLCI),NIFMN**2)
*
*. The transformation matrix between the first NIFMX
*. shells is now known, biorthogonalize remaining orbitals
         IF(NI.NE.NF.AND.NI.NE.0.AND.NF.NE.0) THEN
           IF(NI.GT.NF) THEN
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
           END IF
           NDIFF = NMX - NMN
* Y = -S-1 * X
*. overlap X between remaining orbitals and the other set
       IF(NI.GT.NF) THEN
* I columns F rows
       DO 52 III = NMN+1, NMX
       DO 52 JJJ = 1, NF
        SCR(KLSIF+(III-NMN-1)*NF+JJJ-1) = 
     &  SCR(KLSTOT-1+(JJJ+ILF-1-1)*NLTI+III+ILI-1)
   52  CONTINUE
      ELSE IF (NF.GT.NI) THEN
* F columns I rows
       DO 53 JJJ = NMN+1, NMX
       DO 53 III = 1, NI
        SCR(KLSIF+(JJJ-NMN-1)*NI+III-1) = 
     &  SCR(KLSTOT-1+(JJJ+ILF-1-1)*NLTI+III+ILI-1)
   53  CONTINUE
      END IF
*
           IF(NI.GT.NF) THEN
             CALL TRPMAT(SCR(KLSIFI),NMN,NMN,SCR(KLSCR))
             CALL COPVEC(SCR(KLSCR),SCR(KLSIFI),NMN ** 2 )
           END IF
           CALL MATML4(SCR(KLSCR),SCR(KLSIFI),SCR(KLSIF),
     &                 NMN,NDIFF,NMN,NMN,NMN,NDIFF,0)
           CALL SCALVE(SCR(KLSCR),-1.0D0,NMN*NDIFF)
           CALL COPVEC(SCR(KLSCR),SCR(KLSIF),NMN*NDIFF)
* Construct complete CMX
           CALL SETVEC(SCR(KLSCR),0.0D0,NMX**2)
           DO 300 J = 1, NMX
            IF(J.LE.NIFMN) THEN
             CALL COPVEC(SCR(KLCMX+(J-1)*NIFMN),
     &                   SCR(KLSCR+(J-1)*NMX),NMN)
            ELSE
             CALL COPVEC(SCR(KLSIF+(J-NMN-1)*NMN),
     &                   SCR(KLSCR+(J-1)*NMX),NMN)
             SCR(KLSCR-1+(J-1)*NMX+J) = 1.0D0
            END IF
  300      CONTINUE
*
           CALL COPVEC(SCR(KLSCR),SCR(KLCMX),NMX**2)
         END IF
*. The two upper triangular matrices CI and CF are now known
*. Rotate the shells
         CALL MATML4(SCR(KLPIF),SCR(KLPI),SCR(KLCI),
     &               NGRID,NI,NGRID,NI,NI,NI,0)
         CALL COPVEC(SCR(KLPIF),SCR(KLPI),NI*NGRID)
      
         CALL COPVEC(SCR(KLPI),PI(1,ILI),NI*NGRID)
         CALL MATML4(SCR(KLPIF),SCR(KLPF),SCR(KLCF),
     &               NGRID,NF,NGRID,NF,NF,NF,0)
         CALL COPVEC(SCR(KLPIF),SCR(KLPF),NF*NGRID)
         CALL COPVEC(SCR(KLPF),PF(1,ILF),NF*NGRID)
*
      if(ntest .GE. 1 ) then
        write(6,*) ' Test of overlap of biorthonormal functions'
* F columns I rows
       DO 54 JJJ =1, NF       
       DO 54 III = 1, NI
        SCR(KLSIF+(JJJ-1)*NI+III-1) = 
     &  SCR(KLSTOT-1+(JJJ+ILF-1-1)*NLTI+III+ILI-1)
   54  CONTINUE
       CALL MATML4(SCR(KLSCR),SCR(KLCI),SCR(KLSIF),
     &             NI,NF,NI,NI,NI,NF,1)
       CALL MATML4(SCR(KLSIF),SCR(KLSCR),SCR(KLCF),
     &             NI,NF,NI,NF,NF,NF,0)
        write(6,*)  
     &  ' new overlap matrix ( should be 1 on diagonal, 0 elsewhere )'
        call wrtmat(scr(klsif),ni,nf,ni,nf)
      end if

        IF(NTEST.GE.1) THEN
         WRITE(6,*)
         WRITE(6,*) ' Orbital Rotation matrix for I state'
         CALL WRTMAT(SCR(KLCI),NI,NI,NI,NI)
         WRITE(6,*) ' Orbital Rotation matrix for F state'
         CALL WRTMAT(SCR(KLCF),NF,NF,NF,NF)
         WRITE(6,*)
        END IF
*
* ======================================
* 2. Counter rotate the CI coefficients
* ======================================
*
*. Initial state
*. CI rotation parameters
        write(iscw,*) '       CI counter-rotation...'
       KLTI = KLSIF
       CALL PAMTMT(SCR(KLCI),SCR(KLTI),SCR(KLSCR),NI)
       DO 301 I = 1, NI
         TII = SCR(KLTI-1+(I-1)*NI+I)
         TIII = 1.0D0 / TII
         CALL SCALVE(SCR(KLTI+(I-1)*NI),TIII,I-1)
  301  CONTINUE
        NINL = NINSHL(L+1)
        CALL CITRA(CII,NCSFI,NCII,L,NI,SCR(KLTI),LBUF,LUI,NINL,
     &             SCR(KLCIF),SCR(KLCISC),SCR(KBUF),
     &             SCR(KIBUF),NTESTL )
       CALL COPVEC(SCR(KLCIF),CII,NCII*NCSFI)
*. Final state
*. CI rotation parameters
       KLTF = KLSIF
       CALL PAMTMT(SCR(KLCF),SCR(KLTF),SCR(KLSCR),NF)
       DO 302 I = 1, NF
         TII = SCR(KLTF-1+(I-1)*NF+I)
         TIII = 1.0D0 / TII
         CALL SCALVE(SCR(KLTF+(I-1)*NF),TIII,I-1)
  302  CONTINUE
*
        CALL CITRA(CIF,NCSFF,NCIF,L,NF,SCR(KLTF),LBUF,LUF,NINL,
     &             SCR(KLCIF),SCR(KLCISC),SCR(KBUF),
     &             SCR(KIBUF),NTESTL )
       CALL COPVEC(SCR(KLCIF),CIF,NCIF*NCSFF)
*. End of loop over L
 1000 CONTINUE
*
       RETURN
       END
