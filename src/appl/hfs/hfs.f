*=======================================================================    
*
*      A General Program for Computing Magnetic Dipole and
*      Electric Quadrupole Hyperfine Constants 
*
*      This version of the program assumes orthogonal orbitals
*
*      P. Jonsson and C.G. Wahlstrom, Department of Physics,  
*                                     Lund Institute of Technology
*                                     P.O. Box 118,S-221 00 Lund
*                                     Sweden 
*      C. Froese Fischer,             Department of Computer Science
*                                     Vanderbilt University,
*                                     Nashville, TN 37235
*                                     USA
*      Comput. Phys. Commun. 74 (1993) 399
*
*      Modified for dynamic memory allocation by P. Jonsson, Nov 91
*
*      Modified to calculate electron densities at the nucleus and
*      level field shifts.                       P. Jonsson, Feb 92
*
*      Modified to calculate f-shells by G. Gaigalas  December 1997
*
*      This program should be linked with libang.a libcom.a librad.a
*                       libdudu.a librecls.a libsqlsf1.a libsqlsf2.a
*
*=======================================================================    
*
*     The current limits are:
*     max nwd=94 different orbitals

      PROGRAM HFS 

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      EXTERNAL INITT
      PARAMETER (EPSILO=1.D-9,NWD=94)
      CHARACTER*24 NAME(6)
      CHARACTER ANS*1,IEM(4)*2
      INTEGER SS,SS1,SS2
      REAL TIME(2)

*                        Dimension statements

      dimension orba1(10,20),dipa1(10,20),conta1(10,20),quadb1(10,20)
      dimension vol1(10,20),nbound(2,10),radint1(nwd,nwd)
      dimension azsqr1(nwd,nwd)
Cper change dimensions January 29,2008 cfgcontri(10000,20,4)->cfgcontri(1000000,20,4)
      dimension cfgcontri(1000000,20,4)

      DIMENSION ORBA(20),DIPA(20),CONTA(20),QUADB(20),VOL(20),VOLF1(20)
      DIMENSION ORBF1(20),DIPF1(20),CONTF1(20),QUADF1(20),WTJIJF(20)

*                        Common blocks

      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/INFORM/IREADC,IWRITE,IOUT,IREADJ,IREADW,ISC(6),ISCW
      COMMON/DIAGNL/IJDIAG,JA,JB
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :J1QN2(31,3),IJFUL(16)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      POINTER(QWT,WT(20,1))
      COMMON /NAN/ QWT
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      POINTER (QIORTH,IORTH(1))
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     1 ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     2 QIORTH
      COMMON/OVRINT/IOVEL(2),IOVER(2),IOVEP(2)
      COMMON/FACT/GAM(100)
      COMMON/NTRM/NTERMS

      COMMON/HENSOR/VVSHELL(50),IIHSH,NNOVLP,IIRHO(50),IISIG(50),
     : MMU(50),MMUP(50),NNU(50),NNUP(50)
      COMMON/COEF/ORBF(20),DIPF(20),CONTF(20),QUADF(20),J(20),JP(20),
     : JJMAX,JJMIN,JJ1MAX,JJ1MIN,JJ2MAX,JJ2MIN,LL1,SS1,LL2,SS2,VOLF(20)
      COMMON/HYPER/VHY(20),NHY,IRHY(20),ISHY(20),JRHY(20),JSHY(20),
     1 IMUHY(20),IMUPHY(20),JMUHY(20),JMUPHY(20),
     2 INUHY(20),INUPHY(20),JNUHY(20),JNUPHY(20)

*                       Data statements

      DATA ORBA,DIPA,CONTA,QUADB,VOL/100*0.D0/
      DATA AMHZ,BMHZ,GS/47.70534D0,234.9647D0,2.002319D0/
      DATA IEM/2HL1,2HCS,2HS1,2HC2/
      data orba1,dipa1,conta1,quadb1/800*0.D0/
Cper change dimensions January 29,2008 cfgcontri(10000,20,4)->cfgcontri(1000000,20,4)
      data cfgcontri/80000000*0.D0/

*                       Output formats

610   FORMAT(/20X,' Hyperfine structure calculation '/)
620   FORMAT(/20X,' The configuration set ')
630   FORMAT(/30X,'**',//20X,17HA factors in MHz ,
     : //,2X,8HJ     J'
     : ,8X,54HOrbital        Spin-dipole       Contact         Total )
640   FORMAT(/,A4,2X,A4,4(2X,F15.7))
650   FORMAT(/20X,17HB factors in MHz ,
     : //,2X,8HJ     J'
     : ,8X,10HQuadrupole )
660   FORMAT(/,A4,2X,A4,2X,F15.7)
670   FORMAT(/30X,'**',//20X,28HHyperfine parameters in a.u.,
     : //20X,47Hal             ad             ac             bq ) 
680   FORMAT(/,9X,4(F17.7)/)
*     IBUG3   debug in recoupling package
*     NBUG6   debug in tensor package

CDBG  WRITE(0,*) ' Input IBUG3, NBUG6 (0/1) '
CDBG  READ(5,*) IBUG3,NBUG6

*     The following section concerns input/output and may be 
*     system dependent. Check allowed unit numbers and file name
*     conventions - modify, if necessary

*     IWRITE  debug file
*     IREADC  configuration file
*     IREADJ  j or l file
*     IREADW  wavefunction file
*     IOUT    output file

      IWRITE=6
      IREADC=4
      IREADJ=2
      IREADW=8
      IOUT=3
      ISC(1)=7
      ISCW=0
      INTGR=0
      I=0
      IBUG1=0
      IBUG2=0
      IBUG3=0
      IBUG6=0
      IBUG7=0
      IFULL=0
      NFACTS=32

CSUN  I=IARGC()
10    IF (I.EQ.0) THEN
         WRITE(0,*) 'Name of state'
         READ(5,'(A)') NAME(1)
      ELSE
CSUN        CALL GETARG(1,NAME)
      ENDIF
      K=INDEX(NAME(1),' ')
      IF (K.EQ.1) THEN
         WRITE(0,*) 'Names may not start with a blank'
         GO TO 10
      ELSE
         NAME(1)(K:K+1)='.c'
         NAME(2)=NAME(1)(1:K-1)//'.j'
         NAME(3)=NAME(1)(1:K-1)//'.w'
         NAME(4)=NAME(1)(1:K-1)//'.h'
Cww         NAME(5)=NAME(1)(1:K-1)//'.m'
         NAME(6)=NAME(1)(1:K-1)//'.l'

      ENDIF
      OPEN(UNIT=IREADC, FILE=NAME(1),STATUS='OLD')
      OPEN(UNIT=IOUT, FILE=NAME(4),STATUS='UNKNOWN')
Cww      OPEN(UNIT=ISC(1), FILE=NAME(5),STATUS='UNKNOWN')

*     Write heading

      WRITE(0,610)
      WRITE(IOUT,610)

15    write(0,'(a)') ' Electron density at the nucleus ? (Y/N)'
      READ(5,'(A)') ANS
      IF (ANS.EQ.'N'.OR.ANS.EQ.'n') THEN
         ndens=0
      ELSEIF (ANS.EQ.'Y'.OR.ANS.EQ.'y') THEN
         ndens=1 
      ELSE 
         GOTO 15
      ENDIF

20    WRITE(0,'(A/A/A)') ' Indicate the type of calculation ',
     : ' 0 => diagonal A and B factors only;',
     : ' 1 => diagonal and off-diagonal A and B factors;'
      READ(5,*) ICALC
30    WRITE(0,*)'Input from an MCHF (M) or CI (C) calculation ?'
      READ(5,'(A)') ANS
      IF (ANS.EQ.'M'.OR.ANS.EQ.'m') THEN
         WRITE(IOUT,*) NAME(1)
         IMCHF=1
         NR=1
      ELSEIF (ANS.EQ.'C'.OR.ANS.EQ.'c') THEN
35       WRITE(0,*)'Is the CI calculation J dependant ? (Y/N)'
         READ(5,'(A)') ANS
         IF (ANS.EQ.'N'.OR.ANS.EQ.'n') THEN
            OPEN(UNIT=IREADJ, FILE=NAME(6),STATUS='OLD')
            WRITE(IOUT,*) NAME(6)
            IMCHF=2
         ELSEIF (ANS.EQ.'Y'.OR.ANS.EQ.'y') THEN
            OPEN(UNIT=IREADJ, FILE=NAME(2),STATUS='OLD')
            WRITE(IOUT,*) NAME(2)
            IMCHF=3
         ELSE
            GO TO 35
         ENDIF
         WRITE(0,*) 'Give the index of the dominant cfg in the CI ' 
         WRITE(0,*) 'expansion for which the hfs is to be calculated ?'
         READ(5,*) NR
      ELSE
         GO TO 30 
      ENDIF
      IF (ICALC.EQ.0) IDIAG=1
      IF (ICALC.EQ.1) IDIAG=0
      OPEN(UNIT=IREADW,FILE=NAME(3),STATUS='OLD',FORM='UNFORMATTED')

*     Give data for summation

      nsubsp = 0

50    WRITE(0,*) 'Give 2*I and nuclear dipole and quadrupole moment
     :s (in n.m. and barns) '
      READ(5,*) II,DIPM,Q
      IF (II.EQ.0) GO TO 50 
      IF (DABS(DIPM).LT.EPSILO) THEN
         IF (DABS(Q).LT.EPSILO) GO TO 50 
      ENDIF
      GI=2.D0*DIPM/DFLOAT(II)
      WRITE(IOUT,'(A,F16.11,A)') ' Nuclear dipole moment     ',DIPM,
     :' n.m.'
      WRITE(IOUT,'(A,F16.11,A)') ' Nuclear quadrupole moment ',Q,
     :' barns'
      WRITE(IOUT,'(A,F16.11)') ' Nuclear spin              ',II/2.d0

*     Set factorials

      CALL FACTRL(NFACTS)

*     Start to perform the calculations.
*     Read the configurationlist and print it.
     
      WRITE(0,620)

      CALL CFGN1(NCLOSD)

*     Determine L,S,Jmax,Jmin and the LSJ dependent factors for the main 
*     configuration. Save the LSJ dependent factors for later use.

      CALL LSJ(QNOC,QJ1,NR,LL,SS,JJMAX,JJMIN)
      CALL LSJ(QNOC,QJ1,NR,LL1,SS1,JJ1MAX,JJ1MIN)
      CALL LSJ(QNOC,QJ1,NR,LL2,SS2,JJ2MAX,JJ2MIN)
      CALL LSJFACT(NR,NR,NJQ,IDIAG,ndens)

      DO 60 K=1,NJQ
         ORBF1(K)=ORBF(K)
         DIPF1(K)=DIPF(K)
         CONTF1(K)=CONTF(K)
         VOLF1(K)=VOLF(K)
         QUADF1(K)=QUADF(K)
60    CONTINUE 

*     Read the J dependent weights of the configurations and store them
*     in a vector WT(JI,K1). Read the radial wavefunctions and sort
*     them according to the order in IAJCMP.

      CALL READWT(MAXORB,IMCHF,NCFG,NR,JJMAX,JJMIN)

* --- Check if the read weights are correct

      cnorm = 0.d0
      do 933 k = 1,ncfg
         cnorm = cnorm + wt(1,k)*wt(1,k)
933   continue
      if (cnorm.lt.0.5d0) then
         write(*,*) ' The configuration weights have not '
         write(*,*) ' been read correctly.'
         write(*,*) ' Norm of the weights =',cnorm
         stop
      endif

      CALL READWFN 

*     Calculate the radial one-particle integrals

      call radial1(radint1,maxorb)
      call radial2(azsqr1,maxorb)

*     Calculate and print the contributions to the A and B factors from
*     every pair of configurations.

      INCFG=NCFG
      DO 70 JI=1,NCFG
	 JA=JI
         IF(MOD(JI,10).EQ.0) WRITE(*,*) '   ja =',JI
         IF (IMCHF.LT.3.OR.IDIAG.EQ.1) INCFG=JI
         DO 80 JF=1,INCFG
	    JB=JF
            ID=0
            IF ((IMCHF.LT.3.OR.IDIAG.EQ.1).AND.(JI.NE.JF)) ID=1
            NOVLPS=0
            JMUP=0
            JNUP=0
            JMU=0
            JNU=0
            IF (NORTH.NE.0) THEN
               write(*,'(a)') 
     :        ' This version of the program assumes orthogonal orbitals' 
               STOP
            ENDIF

ctc b The lines below have moved (see l349) - 12.01.2011

            IF (IMCHF.EQ.3) THEN
               CALL LSJ(QNOC,QJ1,JI,LL1,SS1,JJ1MAX,JJ1MIN)
ctc d            write(*,*) JI, LL1, SS1, JJ1MAX, JJ1MIN
               CALL LSJ(QNOC,QJ1,JF,LL2,SS2,JJ2MAX,JJ2MIN)
ctc d            write(*,*) JI, LL2, SS2, JJ2MAX, JJ2MIN
               IF (LL1.EQ.LL.AND.LL2.EQ.LL.AND.SS1.EQ.SS.AND.SS2.EQ.SS) 
     :         THEN
                  DO 90 K=1,NJQ
                     ORBF(K)=ORBF1(K)
                     DIPF(K)=DIPF1(K)
                     CONTF(K)=CONTF1(K)
                     VOLF(K)=VOLF1(K)
                     QUADF(K)=QUADF1(K)
90                CONTINUE
               ELSE
                  CALL LSJFACT(JI,JF,NJQ,IDIAG,ndens)
               ENDIF
            ENDIF
ctc e
Cww Change 10.6 96
            CALL MULTWT(JI,JF,IMCHF,IDIAG,WTJIJF)
            LSKIP = 0 
            DO 62 K=1,NJQ
               IF (DABS(WTJIJF(K)).GT.1.D-14) LSKIP = LSKIP + 1
62          CONTINUE
            IF (LSKIP.EQ.0) GOTO 80
Cww            
            CALL SETUP(JI,JF,LET)
*           .. skip the rest if LET=0; 
*              the CSFs differ by too many shells or electrons
	    IF (LET .eq. 0) GO TO 80

*  Calculate the LSJ dependent factors. If MCHF cfglist, the factors are
*  identical for every cfg. pair and should not be recalculated.
*  If BREIT cfglist there is no need for recalc. when configurationpair
*  has the same LS terms as the main cfg.

ctc   They were here - 12.01.2011

*           Multiply WT(JI,K1) and WT(JF,K2) and save it in a vector
*           WTJIJF(K)

Cww            CALL MULTWT(JI,JF,IMCHF,IDIAG,WTJIJF) 
            
            IF (SS1.NE.SS2) GOTO 101 
C            CALL ORBITAL(JI,JF)
	    CALL NONHIPER(1)
            DO 100 N=1,NHY
               IF (DABS(VHY(N)).GT.EPSILO) THEN
                  DO 110  K=1,NJQ
                     IF (DABS(ORBF(K)).GT.EPSILO) THEN
                        CONTRI=WTJIJF(K)*ORBF(K)*VHY(N)*
     :                  RADINT1(irhy(n),ishy(n))
                        IF (ID.EQ.1) CONTRI=2.D0*CONTRI

                        do 109 mn = 1,nsubsp
                           if ((ji.ge.nbound(1,mn)).and.
     :                     (ji.le.nbound(2,mn))) then
                              if (((jf.ge.nbound(1,mn)).and.
     :                        (jf.le.nbound(2,mn))).or.(jf.eq.1)) then
                                 orba1(mn,k) = orba1(mn,k)+contri
                              endif
                           endif
109                     continue

                        cfgcontri(JI,K,1) = cfgcontri(JI,K,1)+contri
                        if (ji.ne.jf) then
                          cfgcontri(JF,K,1) = cfgcontri(JF,K,1)+contri
                        endif

                        ORBA(K)=ORBA(K)+CONTRI
                     ENDIF
110               CONTINUE
               ENDIF
100         CONTINUE
101         continue
       
C            CALL DIPOLE(JI,JF)
	    CALL NONHIPER(2)
            DO 120 N=1,NHY
               IF (DABS(VHY(N)).GT.EPSILO) THEN
                  DO 130 K=1,NJQ
                     IF (DABS(DIPF(K)).GT.EPSILO) THEN
                        CONTRI=WTJIJF(K)*DIPF(K)*VHY(N)*
     :                  RADINT1(irhy(n),ishy(n))
                        IF (ID.EQ.1) CONTRI=2.D0*CONTRI

                        do 129 mn = 1,nsubsp
                           if ((ji.ge.nbound(1,mn)).and.
     :                     (ji.le.nbound(2,mn))) then
                              if (((jf.ge.nbound(1,mn)).and.
     :                        (jf.le.nbound(2,mn))).or.(jf.eq.1)) then
                                 dipa1(mn,k) = dipa1(mn,k)+contri
                              endif
                           endif
129                     continue

                        cfgcontri(JI,K,2) = cfgcontri(JI,K,2)+contri
                        if (ji.ne.jf) then
                          cfgcontri(JF,K,2) = cfgcontri(JF,K,2)+contri
                        endif

                        DIPA(K)=DIPA(K)+CONTRI
                     ENDIF
130               CONTINUE
               ENDIF
120         CONTINUE
        
            IF (LL1.NE.LL2) GOTO 141
C            CALL CONTACT(JI,JF)
	    CALL NONHIPER(3)
            DO 140  N=1,NHY
               IF (DABS(VHY(N)).GT.EPSILO) THEN
                  DO 150  K=1,NJQ
                     IF (DABS(CONTF(K)).GT.EPSILO) THEN
                        CONTRI=WTJIJF(K)*CONTF(K)*VHY(N)*
     :                  AZSQR1(irhy(n),ishy(n))
                        IF (ID.EQ.1) CONTRI=2.D0*CONTRI

                        do 149 mn = 1,nsubsp
                           if ((ji.ge.nbound(1,mn)).and.
     :                     (ji.le.nbound(2,mn))) then
                              if (((jf.ge.nbound(1,mn)).and.
     :                        (jf.le.nbound(2,mn))).or.(jf.eq.1)) then
                                  conta1(mn,k) = conta1(mn,k)+contri
                              endif
                           endif
149                     continue

                        cfgcontri(JI,K,3) = cfgcontri(JI,K,3)+contri
                        if (ji.ne.jf) then
                          cfgcontri(JF,K,3) = cfgcontri(JF,K,3)+contri
                        endif

                        CONTA(K)=CONTA(K)+CONTRI
                     ENDIF
150               CONTINUE
               ENDIF
140         CONTINUE
141         continue

*           The electron density at the nucleus
*           To get the density we divide by 4*PI 

            IF (LL1.NE.LL2.OR.SS1.NE.SS2) GOTO 142
            if (ndens.eq.0) goto 142
C            CALL VOLUME(JI,JF)
	    CALL NONHIPER(5)
            DO 143  N=1,NHY
               IF (DABS(VHY(N)).GT.EPSILO) THEN
                  DO 144  K=1,NJQ
                     IF (DABS(VOLF(K)).GT.EPSILO) THEN
                        CONT=VOLF(K)*VHY(N)
                        CONTRI=WTJIJF(K)*CONT*AZSQR1(irhy(n),ishy(n))
                        IF (ID.EQ.1) CONTRI=2.D0*CONTRI

                        do 147 mn = 1,nsubsp
                           if ((ji.ge.nbound(1,mn)).and.
     :                     (ji.le.nbound(2,mn))) then
                              if (((jf.ge.nbound(1,mn)).and.
     :                        (jf.le.nbound(2,mn))).or.(jf.eq.1)) then
                                 vol1(mn,k) = vol1(mn,k)+contri
                              endif
                           endif
147                     continue

                        VOL(K)=VOL(K)+CONTRI
                     ENDIF
144               CONTINUE
               ENDIF
143         CONTINUE
142         continue

            IF (SS1.NE.SS2) GOTO 161
C            CALL QDRPOLE(JI,JF)
	    CALL NONHIPER(4)
            DO 160 N=1,NHY
               IF (DABS(VHY(N)).GT.EPSILO) THEN
                  DO 170 K=1,NJQ
                     IF (DABS(QUADF(K)).GT.EPSILO) THEN
                        CONTRI=WTJIJF(K)*QUADF(K)*VHY(N)*
     :                  RADINT1(irhy(n),ishy(n))
                        IF (ID.EQ.1) CONTRI=2.D0*CONTRI

                        do 169 mn = 1,nsubsp
                           if ((ji.ge.nbound(1,mn)).and.
     :                     (ji.le.nbound(2,mn))) then
                              if (((jf.ge.nbound(1,mn)).and.
     :                        (jf.le.nbound(2,mn))).or.(jf.eq.1)) then
                                 quadb1(mn,k) = quadb1(mn,k)+contri
                              endif
                           endif
169                     continue

                        cfgcontri(JI,K,4) = cfgcontri(JI,K,4)+contri
                        if (ji.ne.jf) then
                          cfgcontri(JF,K,4) = cfgcontri(JF,K,4)+contri
                        endif

                        QUADB(K)=QUADB(K)+CONTRI
                     ENDIF
170               CONTINUE
               ENDIF 
160         CONTINUE
161         continue
80       CONTINUE
70    CONTINUE

*     Print the final results.

      do 179 k = 1,njq
         orba(k) = orba(k)*amhz*gi
         dipa(k) = dipa(k)*amhz*gi
         conta(k) = conta(k)*amhz*gi
         quadb(k) = quadb(k)*bmhz*q
         vol(k) = vol(k)*0.25D0*0.318309886d0
         do 182 n = 1,nsubsp
            orba1(n,k) = orba1(n,k)*amhz*gi
            dipa1(n,k) = dipa1(n,k)*amhz*gi
            conta1(n,k) = conta1(n,k)*amhz*gi
            quadb1(n,k) = quadb1(n,k)*bmhz*q
182      continue
179   continue

      if (nsubsp.eq.0) write(iout,630) 
      DO 180 K=1,NJQ
         if (nsubsp.ne.0) write(iout,630) 
         WRITE(IOUT,640) J(K),JP(K),ORBA(K),
     :   DIPA(K),CONTA(K),ORBA(K)+DIPA(K)+CONTA(K)
         if (nsubsp.ne.0) then
         write(iout,*)
         write(iout,*) '                    Contribution from subspaces'
         do 185 n=1,nsubsp
            write(iout,742) nbound(1,n),nbound(2,n),orba1(n,k),
     :      dipa1(n,k),conta1(n,k),orba1(n,k)+dipa1(n,k)+conta1(n,k) 
185      continue
         endif
180   CONTINUE
742   FORMAT(/,I4,2X,I4,4(2X,F13.5))
      
      if (nsubsp.eq.0) write(iout,650) 
      DO 190 K=1,NJQ
         if (nsubsp.ne.0) write(iout,650) 
         WRITE(IOUT,660) J(K),JP(K),QUADB(K)
         if (nsubsp.ne.0) then
         write(iout,*)
         write(iout,*) '                    Contribution from subspaces'
         do 195 n=1,nsubsp
            write(iout,743) nbound(1,n),nbound(2,n),quadb1(n,k)
195      continue
         endif
190   CONTINUE
743   FORMAT(/,I4,2X,I4,1(2X,F13.5))


      IF (ICALC.EQ.1) THEN
         IF (JJMAX.GT.(JJMIN+2)) K=NJQ-2
         IF (JJMAX.EQ.(JJMIN+2)) K=NJQ-1
         IF (JJMAX.EQ.JJMIN) K=NJQ
      ELSE
         K=NJQ
      ENDIF
      BL=DFLOAT(LL)/2.D0
      BS=DFLOAT(SS)/2.D0
      BJ=DFLOAT(JJMAX)/2.D0
      BLJ=(BJ*(BJ+1.D0)+BL*(BL+1.D0)-BS*(BS+1.D0))/2.D0
      BSJ=(BJ*(BJ+1.D0)-BL*(BL+1.D0)+BS*(BS+1.D0))/2.D0
      BSL=(BJ*(BJ+1.D0)-BL*(BL+1.D0)-BS*(BS+1.D0))/2.D0
      IF(ABS(ORBA(K)).LT.EPSILO) THEN
         AL=0.D0
      ELSE
         AL=ORBA(K)*BL*BJ*(BJ+1.D0)/(2.D0*AMHZ*GI*BLJ)
      ENDIF
      IF(ABS(DIPA(K)).LT.EPSILO) THEN
         AD=0.D0
      ELSE
         AD=DIPA(K)*BS*BL*(2.D0*BL-1.D0)*BJ*(BJ+1.D0)/
     :       (AMHZ*GS*GI*(3.D0*BSL*BLJ-BL*(BL+1.D0)*BSJ))
      ENDIF
      IF(ABS(CONTA(K)).LT.EPSILO) THEN
         AC=0.D0
      ELSE
         AC=3.D0*CONTA(K)*BS*BJ*(BJ+1.D0)/(AMHZ*GS*GI*BSJ)
      ENDIF
      IF(ABS(QUADB(K)).LT.EPSILO) THEN
         BQ=0.D0
      ELSE
         BQ=-QUADB(K)*BL*(2.D0*BL-1.D0)*(BJ+1.D0)*(2.D0*BJ+3.D0)/
     :  (BMHZ*Q*(6.D0*BLJ*BLJ-3.D0*BLJ-2.D0*BL*(BL+1.D0)*BJ*(BJ+1.D0)))
      ENDIF
      WRITE(IOUT,670)
      WRITE(IOUT,680) AL,AD,AC,BQ

*     Calculate the electron density from the closed core

      if (ndens.eq.1) then
         call cdens(ireadw,corden)

*     Output the total electron density at the nucleus

         write(iout,*) '               Electron density at the nucleus'
641      FORMAT(/,A4,2X,A4,2X,F15.8)
         DO 200 K=1,NJQ
            if (j(k).eq.jp(k)) then
               WRITE(IOUT,641) J(K),JP(K),VOL(K)+corden
            endif
200      continue
      endif

*     Output cfgcontribution to hfs <name>.m

Cww643   FORMAT(I5,2X,A4,2X,A4,2X,5(F17.10,2X))
Cww      do 220 n = 1,ncfg
Cww         do 230 k = 1,njq
Cww            write(isc(1),643) n,j(k),jp(k),
Cww     :      cfgcontri(n,k,1)*amhz*gi,cfgcontri(n,k,2)*amhz*gi,
Cww     :      cfgcontri(n,k,3)*amhz*gi,cfgcontri(n,k,4)*bmhz*q,
Cww     :      cfgcontri(n,k,1)*amhz*gi+cfgcontri(n,k,2)*amhz*gi+
Cww     :      cfgcontri(n,k,3)*amhz*gi
Cww230      continue
Cww220   continue

999   WRITE(IOUT,'(30X,A)') '**'
Cww 080123
      close (iout)
Cww 080123
      NWF = MAXORB
      if (north .gt. 0) call dalloc(qiorth,nwf*(nwf-1)/2)
      call dalloc(qlist,ncfg)
      call dalloc(qnoc,ncfg)
      call dalloc(qnelcsh,8*ncfg)
      call dalloc(qnocorb,8*ncfg)
      call dalloc(qj1,15*ncfg)
      call dalloc(qnjcomp,nwf)
      call dalloc(qljcomp,nwf)
      call dalloc(qwt,20*ncfg) 
C      call etime_(time)
      write(*,'(//A/A//A,F8.3,A//)') ' END OF CASE','===========',
     :      ' Total CPU time was ', TIME(1)/60, ' minutes'
      STOP
      END
