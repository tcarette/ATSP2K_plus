*
*      A program to compute the density function of an MCHF atomic
*      state function, the natural orbitals and the density at the
*      nucleus.
*
*
*      A. Borgoo,       VUB Bruxelles
*      G. Gaigalas      ITPA Vilnius
*      O. Scharf,       ULB Bruxelles
*      M.R. Godefroid,  ULB Bruxelles
*
*=======================================================================    
*
*	List of Changes:
*
*=======================================================================    
*
*	Restrictions:
*
*     max nwd=60 different orbitals
*     max nod=220 size of the radial points
*     possible orbitals orbitals="spdfghiklmn"
*
*     precision epsilo=1.D-16
*
*=======================================================================    

********************************************
*                                          *
      PROGRAM DENSITY 
*                                          *
********************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
    
*     RESTRICTIONS:     
      PARAMETER (EPSILO=1.D-16,NWD=60,NOD=220)
      CHARACTER*11 ORBITALS
      PARAMETER (ORBITALS='spdfghiklmn')
      
*     EXTERNAL FUNCTIONS:
      EXTERNAL INITT
Cmrg  NAG Routine for the Eigenvalues, Eigenvectors      
Cmrg  EXTERNAL f02faf
*     LAPACK Routine for the Eigenvalues, Eigenvectors      
      EXTERNAL dsyev


      CHARACTER*24 NAME(50)
*      CHARACTER IEM(4)*2
      CHARACTER*3 NOPNAME(NWD)
      CHARACTER*3 ELC(NWD), ATOM*6, TERM*6
      CHARACTER*3 ELNAME
      
      INTEGER DOMINANT
      INTEGER NOORBITALPRINT
      INTEGER NOMATRIXPRINT
      INTEGER NONATURALORBITAL
      INTEGER DENSITYDIM
      INTEGER MM
      INTEGER CO
  
      REAL *8 NAZ
      REAL *8 NOP

      DIMENSION EGVL(NWD),EGVC(NWD,NWD),WORK(3*NWD) 
      DIMENSION DENSITYPLOT(NOD+1)
      DIMENSION PC(NOD,NWD),NC(NWD),LC(NWD),ORB(NOD)
      DIMENSION MM(NWD)
      DIMENSION AZC(NWD)
      DIMENSION NAZ(NWD)
      DIMENSION DENSITYMATRIX( NWD, NWD )
      DIMENSION FACTORMATRIX( NWD, NWD )
      DIMENSION CFACTORMATRIX( NWD )
      DIMENSION NOP(NOD,NWD)
*      dimension orba1(10,20),dipa1(10,20),conta1(10,20),quadb1(10,20)
**      dimension vol1(10,20),nbound(2,10),radint1(nwd,nwd)
**      dimension azsqr1(nwd,nwd)
*      dimension cfgcontri(10000,20,4)
**      DIMENSION ORBA(20),DIPA(20),CONTA(20),QUADB(20),VOL,VOLF1(20)
**      DIMENSION ORBF1(20),DIPF1(20),CONTF1(20),QUADF1(20),WTJIJF(20)


**      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
***      ,ISC(6),ISCW
***      COMMON/DIAGNL/IJDIAG,JA,JB
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :J1QN2(31,3),IJFUL(16)
CAB      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
CAB     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
CAB      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
CAB     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
CAB     :       (QLJCLSD,LJCLSD(1))
CAB      POINTER(QWT,WT(1))
      POINTER(QNOC,NOCCSH(*)),(QNELCSH,NELCSH(8,*)),
     :       (QNOCORB,NOCORB(8,*)),(QJ1,J1QNRD(15,*))
      POINTER(QIAJCMP,IAJCMP(*)),(QLJCOMP,LJCOMP(*)),
     :       (QNJCOMP,NJCOMP(*)),(QIAJCLD,IAJCLD(*)),
     :       (QLJCLSD,LJCLSD(*))
      POINTER(QWT,WT(*))
      COMMON /NAN/ QWT
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     1 ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP
      COMMON/OVRINT/IOVEL(2),IOVER(2),IOVEP(2)
      COMMON/FACT/GAM(100)
      COMMON/NTRM/NTERMS
      
*     THIS IS IMPORTANT:
      COMMON /PARAT/D0,D1,D2,D3,D4,D5,D6,D10,H,H1,NO,ND
      COMMON/INFORM/IREADC,IWRITE,IOUT,IREADJ,IREADW
      COMMON/HYPER/VHY(20),NHY,IRHY(20),ISHY(20),JRHY(20),JSHY(20),
     1 IMUHY(20),IMUPHY(20),JMUHY(20),JMUPHY(20),
     2 INUHY(20),INUPHY(20),JNUHY(20),JNUPHY(20)
*     The radial data for nonclosed shells
      COMMON/ADATA/P(NOD,NWD),R(NOD),RR(NOD),R2(NOD),
     :       ZED,AZ(NWD),L(NWD),MAX(NWD)
      COMMON/ADATA2/ATOM,TERM,ELNAME(NWD)

      
*                       Data statements

*      DATA AMHZ,BMHZ,GS/47.70534D0,234.9647D0,2.002319D0/
*      DATA IEM/2HL1,2HCS,2HS1,2HC2/
*      data orba1,dipa1,conta1,quadb1/800*0.D0/
*      data cfgcontri/800000*0.D0/

*                       Output formats


*     The following section concerns input/output and may be 
*     system dependent. Check allowed unit numbers and file name
*     conventions - modify, if necessary

      IWRITE=6
      IREADC=4
      IREADJ=2
      IREADW=8
      IMATRIXOUT=11
      IWRITEW=9
      IWRITENW=15
      IOUT=3
      IORBITALOUT=12

**********************************************************
*
*  START OF THE PROGRAM ...
*
**********************************************************

*     HEADLINE
      WRITE(0,*) 'Density calculation, Summer 2009 '

*     GET USER INPUT
      CALL USERINPUT( NAME(1),NOORBITALPRINT,
     :           NOMATRIXPRINT,NONATURALORBITAL,
     :           iprintall )
      
*     GET THE FILE NAMES
      K=INDEX(NAME(1),' ')
      NAME(1)(K:K+1)='.c'
      NAME(3)=NAME(1)(1:K-1)//'.w'
      NAME(4)=NAME(1)(1:K-1)//'.d'
      NAME(6)=NAME(1)(1:K-1)//'.l'
      NAME(7)=NAME(1)(1:K-1)//'.plt'
      NAME(8)=NAME(1)(1:K-1)//'.matrix'
      NAME(9)=NAME(1)(1:K-1)//'.n'
      NAME(10)=NAME(1)(1:K-1)//'.nw'
 
*     OPEN REQUIRED FILES
      OPEN(UNIT=IREADC, FILE=NAME(1), STATUS='OLD')
      OPEN(UNIT=IREADJ, FILE=NAME(6), STATUS='OLD')
      OPEN(UNIT=IREADW, FILE=NAME(3), STATUS='OLD',FORM='UNFORMATTED')
      
*     START WRITING TO THE DENSITY FILE .d 
      OPEN(UNIT=IOUT,   FILE=NAME(4), STATUS='UNKNOWN')
      WRITE(IOUT,*) "# Electron density"
      WRITE(IOUT,*) "# Input file: ",NAME(1)

********************************************************
*     
*     SETTING UP THE CALCULATION
*
********************************************************
      WRITE(*,*) 
      WRITE(*,*) "ANALYSING THE CALCULATION"
      WRITE(*,*) "========================="
      WRITE(*,*) 
      WRITE(*,*) "ACCURACY IS SET TO ", EPSILO
      
*     The integral over the density (Number of electrons)	
      DENSITYN=0.0 
   
*     Set factorials
      NFACTS=32
      CALL FACTRL(NFACTS)

*     Read the configuration list and print it.
      CALL CFGN1(NCLOSD)
      WRITE(*,*) ""

*     Read the weights of the configurations and store them
*     in a vector WT(K1). 
      CALL READWT(NCFG,1)

*     Check if the read weights are correct
      CNORM = D0
      DO 933 K = 1,NCFG
         CNORM = CNORM + WT(K)*WT(K)
933   CONTINUE
      WRITE(*,*) 'NORM OF WEIGHTS =',CNORM
      IF (CNORM.LT.0.99D0) THEN
*     .99 is an arbitrary treshold
         WRITE(*,*) ' THE CONFIGURATION WEIGHTS HAVE NOT '
         WRITE(*,*) ' BEEN READ CORRECTLY.'
         STOP
      ENDIF

*     Read the radial wavefunctions and sort
*     them according to the order in IAJCMP.
      CALL READWFN 

*     For debuging the wave function read in ...
      WRITE(*,*) ""
      WRITE(*,*) "ATOM ", ATOM, " TERM ", TERM
      IF( iprintall.EQ.1 ) THEN
          WRITE(*,*) ""
          DO 934 K=1,MAXORB
            WRITE(*,*) ' THE ORBITAL ',ELNAME(K),' HAS ',MAX(K),
     :             ' POINTS, AZ=',AZ(K), '.'
934       CONTINUE
      END IF
      
      IF (NORTH.NE.0) THEN
           WRITE(*,'(A)') 
     :     ' THIS VERSION OF THE PROGRAM ASSUMES ORTHOGONAL ORBITALS' 
           STOP
      ENDIF
      WRITE(*,*) ""
      WRITE( *,* )"ALL WAVEFUNCTIONS EXIST."


********************************************************
*     
*     START THE CALCULATION
*    ========================
* The cut off for the weightproduct is set to EPSILO.
* Weightproducts below this value are not taken into 
* account.
*
* Store the diagonal weightproducts for closed orbitals
* in the CFACTORMATRIX.
*
* The contribution is calculated ( weightproduct x recoupling factor)
* note: use double for offdiagonal elements
*
* Add the contribution of the orbitals in the FACTORMATRIX
* note: we only get the non closed orbitals.
*
*
********************************************************
      WRITE( *,* ) ""
      WRITE( *,* ) ""
      WRITE( *,* ) "START OF THE DENSITY CALCULATION"
      WRITE( *,* ) "================================"
      WRITE( *,* ) ""

*   CALCULATE THE SPIN-ANGULAR PART AND WEIGHT PRODUCT FOR
*   EACH CSF PAIR (CSF1,CSF2) AND ITS ORBITALPAIRS (O1,O2)

*     We print out the Csf contributions from all electrons,
*     given as the recoupling factor and the weightfactor:
      IF( iprintall .eq. 1 ) THEN
          WRITE(*,'(A6,1x,A6,1x,A3,1x,A3,1x,A12,1x,A12)') 
     :          "CSF1","CSF2", "O1","O2","SPIN_ANG","WEIGHTS"
      ENDIF
      
      FIRSTPOINT=0

*     Run through the upper triangle interaction:
*     Calculate values for all CSF with CSF equal or higher
      DO 70 JI=1,NCFG
        DO 80 JF=JI,NCFG

*            Check if we have to calculate anything:
*            If this is a diagonal element: ID=0
*            if not, ID=1, than check with selected
*            typ:
             ID=0
             IF( JI.NE.JF) ID=1

*                WE HAVE TO CALCULATE ...

*                Calculate the weightproduct. Check if it is higher than
*                the cut off EPSILO.
                 WTprod=WT(JI)*WT(JF);
                 IF (DABS(WTprod).GE.EPSILO) THEN
                 
*                    Store the diagonal weightproducts for closed orbitals
*                    in the CFACTORMATRIX.
                     if( ID.eq.0 ) THEN
                        DO 58 I = 1,NCLOSD
                            CFACTORMATRIX(I)=
     :                             CFACTORMATRIX(I)+WTprod;
58                      CONTINUE
                     ENDIF
               
*                   Find out if interaction is between the same
*                   shell or two different ones or if they can not
*                   interact.
                    CALL SETUP(JI,JF,LET)

*                   Do we have no interaction ?					
                    IF (LET .eq. 0) GO TO 80

*                   We have interaction, so get the
*                   decoupling factors for each orbital
                    CALL SPIN_ANGULAR_DENSITY()

                    DO 143  N=1,NHY
*                       For all orbitals, check if we have a contribution
*                       above EPSILO
                        IF (DABS(VHY(N)).GT.EPSILO) THEN
                            IF( iprintall .eq. 1 ) THEN
                               WRITE(*,'(I6,1x,I6,1x,A3,1x,A3,1x,
     :                                   F12.9,1x,F12.9)') 
     :                         JI,JF,ELNAME(IRHY(N)),ELNAME(ISHY(N)),
     :                         VHY(N),WTprod
                            ENDIF
*                           The contribution is calculated,
*                           note: use double for offdiagonal elements
                            CONTRI=VHY(N)*WTprod;
                            IF (ID.EQ.1) CONTRI=2.D0*CONTRI

*                           STORE THE CONTRIBUTION OF THE INTERACTING 
*                           ORBITALS IN THE FACTORMATRIX
*                                WRITE(*,*) "factormatrix [",IRHY(N)
*     :                                      ,"][",ISHY(N),"]"
                                FACTORMATRIX(IRHY(N),ISHY(N))=
     :                            FACTORMATRIX(IRHY(N),ISHY(N))+CONTRI
                        ENDIF
143                 CONTINUE
                 END IF
80       CONTINUE
70    CONTINUE


*     Now we take the radial part.
*     In the same run, we can calculate the modified density at the nucleus
      IF( iprintall .eq. 1 ) THEN
          WRITE( *,* ) ""
          WRITE( *,* ) "MODIFIED DENSITY AT THE NUCLEUS :"
      END IF
*     Calculate electron density from the closed shells.
      IF( NCLOSD.GT.0 ) THEN
          REWIND (IREADW)
          IF( iprintall .eq. 1 ) THEN
              WRITE(*,*) "( CLOSED SHELL CONTRIBUTION )"
          END IF
          DO 55 I = 1,NCLOSD
              READ(IREADW) ATOM,TERM,ELC(I),
     :            MR,Z,ETI,EKI,AZC(I),(PC(JJ,I),JJ=1,MR)
              IF( ELC(I)(1:1).NE.' ') THEN
                  NC(I)=ICHAR( ELC(I)(1:1))-ICHAR('1')
                  LC(I)=LVAL(ELC(I)(2:2))
              ELSE
                  NC(I)=ICHAR(ELC(I)(2:2))-ICHAR('1')
                  LC(I)=LVAL(ELC(I)(3:3))
              END IF    
              VAL=2.*(2.*LC(I)+1.) 
              CFACTORMATRIX(I)=CFACTORMATRIX(I)*VAL

*             For the first density point, add the A0 values only
*             for 's' electrons.
              IF( LC(I) .EQ. 0 ) THEN
                  FIRSTPOINT=FIRSTPOINT+CFACTORMATRIX(I)*AZC(I)*AZC(I)
                  IF( iprintall .eq. 1 ) THEN
                      WRITE(*,*) ELC(I), " : ", 
     :                    CFACTORMATRIX(I)*AZC(I)*AZC(I)  
                  END IF
              END IF

*		TO REMEMBER:
*	    In the Atsp2K <wfn> files, the tabulated points are related to R by:
*                   R_n,l (r) = r^(-1/2) P_n,l(r)
*       The known Legendre function is related to this P by:
*                   P[Legendre] = P[<wfn>] r^1/2
*       Thus the density function is
*       rho =   R^2 = P[Legendre]^2/r^2 = P[<wfn>]^2 / r 
*       Integrating this function gives the number of electrons.

              DO 29 JJ=1,MR
                  DENSITYPLOT(JJ)=DENSITYPLOT(JJ)+
     :              CFACTORMATRIX(I)*PC(JJ,I)*PC(JJ,I)/R(JJ)
29            CONTINUE              
              DENSITYN=DENSITYN+ 
     :           CFACTORMATRIX(I)*
     :             QQUADR(I,I,0)
55        CONTINUE
      ENDIF


*     Calculate electron density from the other shells.
      IF( iprintall .eq. 1 ) THEN
          WRITE(*,*),"( NONCLOSED SHELL CONTRIBUTION )"
      END IF

      DO 18 I=1,NWD
          DO 19 K=1,NWD

              IF( FACTORMATRIX( I ,K ).NE.0 ) THEN
*				  For the first density point, add the A0 values only
*			      for 's' electrons.
                  IF( L(I) .EQ. 0 ) THEN
                      IF( L(I) .EQ. L(K) ) THEN
                          FIRSTPOINT=FIRSTPOINT+FACTORMATRIX(I,K)*
     :                       AZ(I)*AZ(K)
                          IF( iprintall .eq. 1 ) THEN
                              WRITE(*,*) ELNAME(I), " - ",ELNAME(K),
     :                            " : ",FACTORMATRIX(I,K)*AZ(I)*AZ(K)
                          END IF
                      END IF
                  END IF
*		TO REMEMBER:
*	    In the Atsp2K <wfn> files, the tabulated points are related to R by:
*                   R_n,l (r) = r^(-1/2) P_n,l(r)
*       The known Legendre function is related to this P by:
*                   P[Legendre] = P[<wfn>] r^1/2
*       Thus the density function is
*       rho =   R^2 = P[Legendre]^2/r^2 = P[<wfn>]^2 / r 
*       Integrating this function gives the number of electrons.
               DO 30 JJ=1,NOD
                   DENSITYPLOT(JJ)=
     :              DENSITYPLOT(JJ)+P(JJ,I)*P(JJ,K)/R(JJ)*
     :                  FACTORMATRIX( I ,K )
30             CONTINUE
*              Calculate the overlap integral and add it to calculate the
*              number of electrons.
               DENSITYN=DENSITYN+
     :              FACTORMATRIX(I,K)*QQUADR(I,K,0)
            END IF



19        CONTINUE
18    CONTINUE

*     Give the density at r=0 (Firstpoint)
      WRITE (*,*) ""
      WRITE (*,*) "MODIFIED ELECTRON DENSITY AT THE NUCLEUS:"
      WRITE( *,'(A5,1x,F25.17)' ) 'O = ',FIRSTPOINT



********************************************************
*
*   WRITE OUT THE ORBITALS
*  ========================
*
* Write orbitals to NAME(7)
*
********************************************************
      IF ( NOORBITALPRINT .eq. 1 ) THEN
          OPEN(UNIT=IORBITALOUT, FILE=NAME(7), STATUS='UNKNOWN')
      END IF
      IF( iprintall.EQ.1 ) THEN
          WRITE (*,*) ""
          WRITE (*,*) "NORM OF THE ORBITALS "
      END IF    
*     The closed shells.
      IF( NCLOSD.GT.0 ) THEN
          REWIND (IREADW)
          DO 59 I = 1,NCLOSD
              READ(IREADW) ATOM,TERM,ELC(I),
     :            MR,Z,ETI,EKI,AZC(I),(PC(JJ,I),JJ=1,MR)
              IF( NOORBITALPRINT .eq. 1 ) THEN
                  WRITE(IORBITALOUT,*) ""
                  WRITE(IORBITALOUT,*) ""
                  WRITE(IORBITALOUT,*) "# Orbital : ",ELC(I)
                  WRITE(IORBITALOUT,*) "# r  R  P=rR P<wfn>"
              ENDIF
              Do 646 JJ=1,MR
                  ORB(JJ)=PC(JJ,I)*PC(JJ,I)*R(JJ)
                  IF( NOORBITALPRINT .eq. 1 ) THEN
                      write(IORBITALOUT,'(F25.17,4x,F25.17,4x,
     :                         F25.17,4x,F25.17 )') 
     :                         R(JJ),PC(JJ,I)/R2(JJ),
     :                         PC(JJ,I)*R2(JJ),PC(JJ,I)
                  ENDIF
646           CONTINUE
              IF( iprintall.EQ.1 ) THEN
                  WRITE( *,'(A4,A8,F25.17)' ) ELC(I), 
     :                 " Norm = ",QQUADR(I,I,0)
              END IF
59        CONTINUE
      ENDIF
*     The open shells.  
      DO 643 KD=1,MAXORB
          IF( NOORBITALPRINT .eq. 1 ) THEN
              WRITE(IORBITALOUT,*) ""
              WRITE(IORBITALOUT,*) ""
              WRITE(IORBITALOUT,*) "# Orbital : ", ELNAME(KD)
              WRITE(IORBITALOUT,*) "# r  R  P=rR  P<wfn>"
           ENDIF
          Do 644 JJ=1,MAX(KD)
*             Multiply P/r with r^2 from the Jacobian
              ORB(JJ)=P(JJ,KD)*P(JJ,KD)*R(JJ)
              IF( NOORBITALPRINT .eq. 1 ) THEN
                  write(IORBITALOUT,'(F25.17,4x,F25.17,4x,
     :                F25.17,4x,F25.17)') 
     :               R(JJ),P(JJ,KD)/R2(JJ),P(JJ,KD)*R2(JJ),P(JJ,KD)
              ENDIF
644       CONTINUE
          IF( iprintall.EQ.1 ) THEN
              WRITE( *,'(A4,A8,F25.17)' ) ELNAME(KD),
     :            " Norm = ",QQUADR(KD,KD,0)
          END IF
643   CONTINUE


********************************************************
*     
*     CREATE THE DENSITY MATRIX
*    ===========================
* As the FACTORMATRIX holds the diagonal elements of the
* density matrix and the sum of the off diagonal elements
* the DENSITYMATRIX is build from the diagonal elements
* and half of the off diagonal elements.
*
********************************************************
      IF( iprintall .eq. 1 ) THEN
          WRITE( *,*) ""
          WRITE( *,*) "DENSITY MATRIX:"
      END IF
*     Holds the density matrix dimension      
      DENSITYDIM=0
*     Run through all possible orbital pairs ...	  
      DO 181 I=1,MAXORB
          DO 191 K=1,MAXORB
             IF( I .EQ. K ) THEN
*                Take the diagonal elements      			 
                 DENSITYDIM=DENSITYDIM+1
                 DENSITYMATRIX(I,K)=FACTORMATRIX( I ,K )
                 IF( iprintall .eq. 1 ) THEN
                     WRITE (*,*) ELNAME(I), " - ",ELNAME(K),
     :                      " : ",DENSITYMATRIX( I ,K )
                 END IF
             ELSE
*                Take the off diagonal elements  
*                If we do not have this pair, use the
*                reverse pair (i.e. 1s-2s = 2s-1s
                 IF( FACTORMATRIX( I ,K ).NE.0 ) THEN
                    DENSITYMATRIX(I,K)=FACTORMATRIX( I ,K )*.5
                 ELSE
                     DENSITYMATRIX(I,K)=FACTORMATRIX( K ,I )*.5
                 END IF

                IF( FACTORMATRIX(I,K) .NE. 0 ) THEN
                     IF( iprintall .eq. 1 ) THEN
                         WRITE (*,*) ELNAME(I), " - ",ELNAME(K),
     :                      " : ",DENSITYMATRIX( I ,K )
                     END IF
                 END IF

             END IF
191        CONTINUE
181    CONTINUE

********************************************************
*     
*     CALCULATE NATURAL ORBITALS
*    ============================
*
********************************************************
      IF( NONATURALORBITAL.EQ.1 ) THEN
*         We calculate the Eigenvector and Eigenvalues of the density matrix
*         LAPACK Routine
*         dsyev( jobz, uplo, n, a, lda, w, work, lwork, info )
**         NAG Routine:
**         f02faf( JOB, UPLO, N, A, LDA, W, WORK, LWORK, IFAIL )
*         JOB  : 'V" to compute Eigenvalues and Eigenvectors
*         UPLO : 'U' for upper triangle
*         N    : order of the matrix
*         A    : Symmetric matrix, will hold the Eigenvectors
*         LDA  : Number of Eigenvectors
*         W    : Eigenvalues
*         WORK : Workspace
*         LWORK: optimal N*64, size of workspace
*         IFAIL: Errorhandling

          EGVC=DENSITYMATRIX

Cmrg          call f02faf('V','L',DENSITYDIM,EGVC,NWD,
Cmrg     :                EGVL,WORK,3*NWD,IFAIL)
          call dsyev('V','L',DENSITYDIM,EGVC,NWD,
     :                EGVL,WORK,3*NWD,IFAIL)

          IF (IFAIL .NE. 0) THEN 
             print*,' IFAIL = ',IFAIL
             WRITE(6,*) ' DIAGONALIZATION FAILS ...'
             STOP 
          END IF
          
*         Print out the new Matrix

          IF( iprintall.EQ.1 ) THEN  
              WRITE(*,*) ""
              WRITE(*,*) "EIGENVALUES AND EIGENVECTORS"
              DO 78 I=1,DENSITYDIM
                 WRITE(*,*) "Eigenvalue ",I,"=",EGVL(I)
                 DO 79 K=1,DENSITYDIM
                     WRITE(*,*) ELNAME(K)," :", EGVC(K,I)
79               CONTINUE
78            CONTINUE
          END IF


*         INITIALISE NATURAL ORBITALS
          DO 66 I=1,DENSITYDIM
            NAZ(I)=D0
            DO 65 K=1,NO
            NOP(K,I)=D0
65          CONTINUE
66        CONTINUE

*         Print the Eigenvalues, highest first
          WRITE (*,*)
          WRITE (*,*) "EIGENVECTOR:"
*         The sum of the eigenvalues ( should be N)
          SUM=0
*         Counts the number of natural orbitals.
          CO=0
*         We run over all 'l' values, starting with s,p, ...          
          DO 380 LVALUE=1,LEN(ORBITALS)
*             The corresponding lowest n value minus 1:
*             (for l=p -> n=1, l=d -> n=2 )
*             This is added for each following orbital.
              NNR=LVALUE-1
*             We run over all Eigenvectors and search for
*             the ones with the right 'l' symmetry
              DO 381 K=1,DENSITYDIM
*                 We start a new orbital, so we have no contribution
*                 yet. Set this flag to zero. Once we have a contribution,
*                 this is set to 1.

*                 We run over all orbitals in the Eigenvector and find the 
*                 dominant contribution
                  DOMINANT=0
                  CONTRIBUTION=D0
                  DO 332 I=1,DENSITYDIM
                     IF(DABS(EGVC(I,DENSITYDIM-K+1 )).GT.CONTRIBUTION )
     :               THEN
                         DOMINANT=I
                         CONTRIBUTION=DABS(EGVC(I,DENSITYDIM-K+1 ))
                     END IF
332               CONTINUE
*                 Check if dominant term has right 'l' value
                  IF( ELNAME(DOMINANT)(3:3).EQ.ORBITALS(LVALUE:LVALUE) )
     :            THEN
*                     We have a good eigenvector. Print it out
*                     and create its natural orbital.
*                     First, find it a proper name:
*                     Get the n value                                 
                      CO=CO+1
                      SUM=SUM+EGVL(DENSITYDIM-K+1)
                      DO
                          NNR=NNR+1
*                         Create the new orbital name:     
                          NOPNAME(CO)(1:1)=' '
                          NOPNAME(CO)(2:2)=
     :                     CHAR( NNR+ICHAR('0'))
                          NOPNAME(CO)(3:3)=
     :                     ORBITALS(LVALUE:LVALUE)
*                         Look if this orbital is a closed orbital:
                          HAVENAME=1
                          DO NCL=1,NCLOSD
                              IF( ELC(NCL).EQ.NOPNAME(CO)) 
     :                         THEN
                                  HAVENAME=0
                                  EXIT
                              END IF    
                          END DO
                          IF( HAVENAME.EQ.1 ) THEN
                              EXIT
                          END IF
                      END DO
*                     Write out first informations about these new natural
*                     orbital:
                      WRITE (*,*)  
                      WRITE( *,'(I3,A13,I3," : ",E25.17)' ) 
     :                          CO," = Eigenvalue" , DENSITYDIM-K+1, 
     :                               EGVL(DENSITYDIM-K+1)
                      WRITE( *,'(A3,A3)' ) NOPNAME(CO),"'="
                      WRITE( *,'(E25.17,4X,A4,A8,E25.17)') 
     :                           EGVC(DOMINANT,DENSITYDIM-K+1),
     :                           ELNAME(DOMINANT),' AZ=', AZ(DOMINANT)
*                     Now we collect all contributions with
*                     the same 'l' symmetry of the eigenvector
                      DO 382 I=1,DENSITYDIM
                       IF(ELNAME(I)(3:3).EQ.ORBITALS(LVALUE:LVALUE))THEN
*                                 Add the contribution of the orbital to the
*                                 natural orbital (Eigenvector*P
                                  NOP(1:NO,CO)=NOP(1:NO,CO)+
     :                               EGVC(I,DENSITYDIM-K+1)*P(1:NO,I)
*                                 Add the first point AZ:     
                                  NAZ(CO)=NAZ(CO)+
     :                              EGVC(I,DENSITYDIM-K+1)*AZ(I)
                                  IF( I.NE.DOMINANT ) THEN
                                    WRITE( *,'(E25.17,4X,A4,A8,E25.17)')
     :                                EGVC(I,DENSITYDIM-K+1),
     :                                ELNAME(I),' AZ=',AZ(I)
                                  END IF 
                       END IF
382                   CONTINUE        
                  END IF
381           CONTINUE
380       CONTINUE

          WRITE (*,*) 
          WRITE (*,*) "SUM OF EIGENVALUES ", sum

* FIND THE LAST NONEZERO NUMBER FOR EACH NOP
          DO 68 I=1,DENSITYDIM
              MM(I)=NO
              DO 67 K=1,NO
                  IF( NOP(NO+1-K,I).NE.0 ) THEN
*                    We want to have a zero as last entry
*                    therefore NO+2-K
                     MM(I)=NO+2-K
                     GOTO 68
                   END IF
67            CONTINUE
68        CONTINUE


********************************************************
*     
*     WRITE THE NATURAL ORBITALS, FORMATTED
*    =======================================
*
********************************************************
      OPEN(UNIT=IWRITEW, FILE=NAME(9), STATUS='UNKNOWN')
*     Start with the closed shells ...
      IF( NCLOSD.GT.0 ) THEN
          REWIND (IREADW)
          DO 88 I = 1,NCLOSD
              READ(IREADW) ATOM,TERM,ELC(I),
     :            MR,Z,ETI,EKI,AZC(I),(PC(JJ,I),JJ=1,MR)
          
              WRITE(IWRITEW,*) ""
              WRITE(IWRITEW,*) ""
              WRITE(IWRITEW,*) "# Natural Orbital : ",ELC(I)
              WRITE(IWRITEW,*) "# r  R  P=rR P<wfn>"
              Do 746 JJ=1,MR
                  WRITE(IWRITEW,'(F25.17,4x,F25.17,4x,
     :                         F25.17,4x,F25.17 )') 
     :                         R(JJ),PC(JJ,I)/R2(JJ),
     :                         PC(JJ,I)*R2(JJ),PC(JJ,I)
746           CONTINUE
88        CONTINUE
      END IF

*     now the nonclosed shells ...
      DO 481 I=1,DENSITYDIM
          WRITE(IWRITEW,*) ""
          WRITE(IWRITEW,*) ""
          WRITE(IWRITEW,*) "# Natural Orbital : ",NOPNAME(I)
          WRITE(IWRITEW,*) "# r  R  P=rR P<wfn>"
          
          DO 491 JJ=1,MM(I)
              WRITE(IWRITEW,'(F25.17,4x,F25.17,4x,
     :                        F25.17,4x,F25.17 )') 
     :                        R(JJ),NOP(JJ,I)/R2(JJ),
     :                        NOP(JJ,I)*R2(JJ),NOP(JJ,I)
491       CONTINUE
481   CONTINUE
      CLOSE(IWRITEW)
 
********************************************************
*     
*     WRITE THE NATURAL ORBITALS UNFORMATTED
*    ========================================
*
********************************************************
      OPEN(UNIT=IWRITENW, FILE=NAME(10), STATUS='REPLACE',
     :         FORM='UNFORMATTED')
*     Start with the closed shells ...
      IF( NCLOSD.GT.0 ) THEN
          REWIND (IREADW)
          DO 81 I = 1,NCLOSD
              READ(IREADW) ATOM,TERM,ELC(I),
     :            MR,Z,ETI,EKI,AZC(I),(PC(JJ,I),JJ=1,MR)
              WRITE (IWRITENW) ATOM,TERM,ELC(I),MR,
     :          Z,ETI,EKI,AZC(I),(PC(JJ,I),JJ=1,MR)
81        CONTINUE
      END IF
*     now the nonclosed shells ...          
      DO 482 I=1,DENSITYDIM
          WRITE (IWRITENW) ATOM,TERM,NOPNAME(I),MM(I),
     :          ZED,D0,D0,NAZ(I),(NOP(K,I),K=1,MM(I))
482       CONTINUE
       END IF

*      
********************************************************
*
*   WRITE OUT THE DENSITY MATRIX
*  ==============================
*
* Write full matrix to NAME(8)
*
********************************************************
      IF ( NOMATRIXPRINT .eq. 1 ) THEN
          
          OPEN(UNIT=IMATRIXOUT, FILE=NAME(8), STATUS='UNKNOWN')
*         HEADLINE
          WRITE(IMATRIXOUT,*) "# Dimension: ",DENSITYDIM
          DO 282 I=1,MAXORB
            WRITE(IMATRIXOUT,'(A25,$ )')
     :                     ELNAME(I)
282       CONTINUE          
          WRITE(IMATRIXOUT,*) 
          DO 281 I=1,MAXORB
              WRITE(IMATRIXOUT,'(A5,$ )')
     :                     ELNAME(I)
              DO 291 K=1,MAXORB
                  WRITE(IMATRIXOUT,'(F25.17,$ )')
     :                     DENSITYMATRIX(I,K)
291           CONTINUE
              WRITE(IMATRIXOUT,*) 
281       CONTINUE
      END IF


********************************************************
*
*   SOME ADDITIONAL OUTPUT
*  ========================
*
********************************************************
*     Integrate the density function
      WRITE( *,* ) ""
      WRITE( *,* ) "INTEGRAL OF THE DENSITY FUNCTION:"
      WRITE( *,'(A5,1x,F25.17)' ) 'N = ',DENSITYN

*     Write out the density function
      WRITE(IOUT,*) '# R   W(R)   D(R)=W(R)*R^2 '
      write(iout,'(F25.17,4x,F25.17,4x,F25.17)') 
     :           0.0,FIRSTPOINT
*     :           0.0,FIRSTPOINT/(4*PI)
      Do 642 JJ=1,NOD
          IF( DENSITYPLOT(JJ)*R(JJ)*R(JJ).GT.EPSILO ) THEN
              write(iout,'(F25.17,4x,F25.17,4x,F25.17)') 
     :           R(JJ),DENSITYPLOT(JJ),DENSITYPLOT(JJ)*R(JJ)*R(JJ)
          END IF
642   CONTINUE

      WRITE( *,* ) ""
      WRITE( *,* ) "DENSITY FUNCTION IS IN FILE ",NAME(4)
      WRITE( *,* ) "END."
      END


**************************************************************
*     
*     GET THE USER INPUT
*    ====================
* The input files are <name>.c, <name>.w and <name>.l
* The output file name is <name>.d for the density and
* <name>.plt for the orbital data and <name>.matrix for
* the density matrix, natural orbitals are in <name>.n.
*
* The program can calculate only diagonal, only off diagonal
* or both contributions.
*
* The program can write out the radial parts of the orbitals
*
* The program can write out the density matrix
*
* The program can write out the natural orbitals
*
* The program can write out all information to the screen
*
********************************************************
      SUBROUTINE USERINPUT( NAME,NOORBITALPRINT,
     :                      NOMATRIXPRINT,NONATURALORBITAL,
     :                      iprintall )

         CHARACTER*24 NAME
         CHARACTER ANS
         INTEGER NOORBITALPRINT
         INTEGER NOMATRIXPRINT
         INTEGER NONATURALORBITAL
         INTEGER iprintall

*        NAME OF THE CALCULATION	  
10       WRITE(0,*) 
     :     'Give <name> of the <name>.c, <name>.l <name>.w files:'
         READ(5,'(A)') NAME
         K=INDEX(NAME,' ')
         IF (K.EQ.1) THEN
            WRITE(0,*) 'Names may not start with a blank'
            GO TO 10
         ENDIF
         WRITE(*,*) "Files: ",NAME(1:K)
         WRITE(*,*) 
            
*        PRINTING THE ORBITALS ...
*        Ask for printing out the orbital files
         WRITE(*,*) 
21       WRITE(0,*) 'PRINT THE ORBITALS  (*/n) '
         READ(5,'(A1)') ANS
         IF (ANS.EQ.'n' .OR. ANS.EQ.'N') THEN
            NOORBITALPRINT=0
            WRITE(*,*) "Do not printout orbitals"
         ELSEIF (ANS.EQ.'*' .OR. ANS.EQ.'') THEN
            NOORBITALPRINT=1
            WRITE(*,*) "Printout orbitals"
         ELSEIF (ANS.EQ.'y' .OR. ANS.EQ.'Y') THEN
            NOORBITALPRINT=1
            WRITE(*,*) "Printout orbitals"
         ELSE
           WRITE(0,*) 'PLEASE CHOOSE ONE OF N,*'
           GO TO 21
         ENDIF

*        PRINTING THE MATRIX FILE ...
*        Ask for printing out the matrix file
         WRITE(*,*) 
22       WRITE(0,*) 'PRINT THE MATRIX (*/n) '
         READ(5,'(A1)') ANS
         IF (ANS.EQ.'n' .OR. ANS.EQ.'N') THEN
            NOMATRIXPRINT=0
            WRITE(*,*) "Do not printout the matrix"
         ELSEIF (ANS.EQ.'*' .OR. ANS.EQ.'') THEN
            NOMATRIXPRINT=1
            WRITE(*,*) "Printout the matrix"
         ELSEIF (ANS.EQ.'y' .OR. ANS.EQ.'Y') THEN
            NOMATRIXPRINT=1
            WRITE(*,*) "Printout the matrix"
         ELSE
           WRITE(0,*) 'PLEASE CHOOSE ONE OF N,*'
           GO TO 22
         ENDIF

*        CALCULATING NATURAL ORBITALS ...
*        Ask for calculation natural orbitals
         WRITE(*,*) 
23       WRITE(0,*) 'CALCULATE NATURAL ORBITALS (*/n) '
         READ(5,'(A1)') ANS
         IF (ANS.EQ.'n' .OR. ANS.EQ.'N') THEN
            NONATURALORBITAL=0
            WRITE(*,*) "Do not calculate natural orbitals"
         ELSEIF (ANS.EQ.'*' .OR. ANS.EQ.'') THEN
            NONATURALORBITAL=1
            WRITE(*,*) "Calculate natural orbitals"
         ELSEIF (ANS.EQ.'y' .OR. ANS.EQ.'Y') THEN
            WRITE(*,*) "Calculate natural orbitals"
            NONATURALORBITAL=1
         ELSE
           WRITE(0,*) 'PLEASE CHOOSE ONE OF N,*'
           GO TO 23
         ENDIF


*        PRINTING ALL INFORMATIONS TO THE SCREEN
*        Ask for printing out all informations
         WRITE(*,*) 
Cab24       WRITE(0,*) 'PRINT ALL DATA (*/y)'
24       WRITE(0,*) 'PRINT ALL DATA (y/*)'
         READ(5,'(A1)') ANS
         IF (ANS.EQ.'n' .OR. ANS.EQ.'N') THEN
            iprintall=0
            WRITE(*,*) "Do not print all informations"
         ELSEIF (ANS.EQ.'*' .OR. ANS.EQ.'') THEN
            iprintall=0
            WRITE(*,*) "Do not print all informations"
         ELSEIF (ANS.EQ.'y' .OR. ANS.EQ.'Y') THEN
            iprintall=1
            WRITE(*,*) "Print all informations"
         ELSE
           WRITE(0,*) 'PLEASE CHOOSE ONE OF Y,*'
           GO TO 24
         ENDIF

        RETURN
       END
**************************************************************





















































































































































































































