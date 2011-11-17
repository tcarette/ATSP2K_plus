***********************************************************************
*
*
* --- This SUBROUTINE computes the Breit-Pauli matrix elements
*     for a variety of operators:
*
*     ONE-ELECTRON OPERATOR (L-integrals)
*     ELECTROSTATIC INTERACTION (F, G, R) integrals)
*     SPIN-ORBIT INTERACTION  (Z-integrals)
*     SPIN-OTHER-ORBIT INTERACTION (N, V integrals)
*     SPIN-SPIN INTERACTION (S-integrals)
*
***********************************************************************
*
      Subroutine BREVALA(MSOO)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NTERMD=31,NWD=60)
*
      PARAMETER (LSDIM=30000)
      POINTER (qcn,cn(lsdim)),(qinptr,inptr(lsdim)),
     :        (qnijptr,nijptr(lsdim)),(qjan,jan(lsdim)),
     :        (qjbn,jbn(lsdim)),(qintptr,intptr(0:2*lmax+1,7)),
     :        (qpackn,ipackn(1)),(qlused,lused(1)),(qico,ico(1))
      COMMON /buffer/qcn,qinptr,qpackn,qlused,qintptr,lmax,qnijptr,
     :               qjan,qjbn,qico
      POINTER  (qjptr, jptr(1))
      COMMON /fout/n,ntot,iflag,nih,nij,qjptr
      CHARACTER ANS*2, string*72
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON/DIAGNL/IDIAG,JA,JB
      COMMON/IMAGNT/ IREL,ISTRICT,IELST
CG 
      COMMON /OPERAT/ ICOLOM,ISOTOP,IORBORB
      COMMON /BREIT/ ISPORB,ISOORB,ISPSPN
CG 
      COMMON /INFORM/IREAD,IWRITE,IOUT,ISC(8),ISCW
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      POINTER(QSIGNFA,SIGNFA(1))
      COMMON/PHASES/QSIGNFA,ICSTAS
      POINTER(QACMULT,ACMULT(1))
      COMMON /SPORB/ QACMULT

      POINTER (IQLSP,LSP(1)) 
      DIMENSION index(NTERMD)
      CHARACTER*3 el(NWD)
*
      LOGICAL skip
      EXTERNAL INITT

    2 FORMAT(1H1////29X,21H---------------------/29X,21HTHE CONFIGURATIO
     :N SET/29X,21H---------------------///)
  105 FORMAT (49H ISPORB=0 AND ISOORB=1 CAUSES THE PROGRAM TO FAIL,
     :  34H BECAUSE THE BLUME WATSON FORMULAE,/
     :  50H CANNOT BE USED FOR CLOSED SUBSHELLS.  TO OVERCOME,
     :  34H THIS, THE CODE HAS SET ISPORB = 1//)
   11 FORMAT(////' THE TYPE OF CALCULATION IS DEFINED BY ',
     :  'THE FOLLOWING PARAMETERS - '/
     : 5X,22H BREIT-PAULI OPERATORS,13X,8HIREL   =,I2/
     : 5X,27H PHASE CONVENTION PARAMETER,8X,8HICSTAS =,I2/)
   13 FORMAT(40H RELATIVISTIC OPERATORS TO BE INCLUDED -/5X,13H SPIN-ORB
     :IT (,I1,22H),  SPIN-OTHER-ORBIT (,I1,15H),  SPIN-SPIN (,I1,1H)/)
   24 FORMAT(36H0INITIAL DEBUG: IN 1-ELECTRON PART =,I2,2H ,,5X,
     : 20H IN 2-ELECTRON PART=,I2/16X,23HIN RECOUPLING PACKAGE =,I2/)
   42 FORMAT(//' MATRIX ELEMENTS CONSTRUCTED USING ',
     :       'THE SPHERICAL HARMONIC PHASE CONVENTION OF'/)
   43 FORMAT(/16X,47HCONDON AND SHORTLEY, THEORY OF ATOMIC STRUCTURE/
     :16X,47H-----------------------------------------------///)
   44 FORMAT(25X,42HFANO AND RACAH, IRREDUCIBLE TENSORIAL SETS/25X,42H--
     :----------------------------------------///)
   50 FORMAT(/20X,'======================='/
     :        20X,' B R E I T - P A U L I '/
     :        20X,'======================='/)
   78 FORMAT(19H DEBUG PARAMETERS -/5X,16H NBUG6(TENSOR) =,I2/5X,
     : ' NBUG7(RELATIVISTIC OPERATORS - SO,SOO,SS) =',I2//)
*
*     .. set rel in COMMON to be rl(input argument)
CG 
      ISOTOP=0
      ICOLOM=0
      IORBORB=0
      MYID = 0 
         write(IWRITE,50)
         write(iscw,'(A/A)')
         write(iscw,*)

        call inp_type(irel,iread,iscw,ISPORB,ISOORB,ISPSPN,
     :     IORBORB,ICOLOM,iwrite,IFULL,ICSTAS,skip)
*
*  --:  Determine debug parameters
*
      IBUG1 = 0
      IBUG2 = 0
      IBUG3 = 0
      NBUG6 = 0
      NBUG7 = 0
*
      WRITE(IWRITE,11) IREL,ICSTAS
*
*     SET FACTORIALS AND LOG OF FACTORIALS
*
      IF(IORBORB.EQ.1) MSOO=MSOO+3
   79 CALL FACTRL(32)
*
* --- READ IN THE SET OF CONFIGURATIONS
*
      WRITE(IWRITE,2)

      CALL ACNFIG(NCLOSD)
*
      IDG = 0
      NZERO = NCFG
      NEW = NCFG
      ISTART = 1
        call inp_interact(new,nzero,istrict,istart,ISCW,NCFG)

*
*     ..  determine the number of electrons
      nnel = 0
      do i = 1,nclosd
        nnel = nnel + 2*(2*ljclsd(i)+1)
      end do
      do i = 1,noccsh(1)
        nnel = nnel + nelcsh(i,1)
      end do

      nwf = maxorb+nclosd
*     .. find maximum l-value
      lmax = 0
      do i=1,maxorb
	lmax = max0(lmax,ljcomp(i))
      end do
*     
*     .. generate list of configurations
      call genintbr(nclosd,maxorb,lmax,qljcomp,qintptr,qpackn,qvalue,
     :                 .false. ,skip,.false.)

*
*     .. allocate the arrays of dimension lsdim and the jptr array
      call alloc(qcn,    lsdim, 8)
      call alloc(qinptr, lsdim, 4)
      call alloc(qjan,   lsdim, 4)
      call alloc(qico,   lsdim, 4)
      call alloc(qjptr,  ncfg,  4)

*     .. start generating the matrix lists
*     n   number of coefficients in the current block
*     nij number of coefficients since the beginning
*     nze number of non-zero matrix elements since beginning
      n   =0
      nij = 0
      nze = 0
      mycol = 0;
      nprocs = 1;
      DO jb = 1, ncfg
         if(mod(jb,100).eq.0) write(ISCW,'(A,I5)') '   jb = ',jb
         if(jb == ncfg) write(ISCW,'(A,I5)') '   jb = ',jb
	 CALL SHELLSJB(JB)
         call BreitGG(NEW,NZERO,IFIRST,idg,skip,nze)
         write(11) nih, (jan(i),i=1,nih);
*        pRINT*, nih, (jan(i),i=1,nih);
         write(12) nih, (ico(i),i=1,nih);
*        pRINT*, nih, (ico(i),i=1,nih);
         mycol = mycol + 1
	 jptr(mycol) = nij
      end do

*     .. finish writing the coefficient data, if non empty arrays
        if (n.eq.lsdim) then
           write(50) lsdim,(cn(j),j=1,lsdim),(inptr(j),j=1,lsdim)
*          write(*,*) lsdim,(cn(j),j=1,lsdim),(inptr(j),j=1,lsdim)
           n=0
           cn(n)=0
           inptr(n)=0
         end if
         write(50) n,(cn(j),j=1,n),(inptr(j),j=1,n)
*        write(*,*) n,(cn(j),j=1,n),(inptr(j),j=1,n)
*
*     .. form list of electrons
*     
      do i = 1,nclosd
        write(el(i),'(A3)') iajcld(i)
      end do
      do i = nclosd+1,nwf
        write(el(i),'(A3)') iajcmp(i-nclosd)
      end do
*     print *, 'el(i)', el(1:nwf)
*     .. write out intial data about the problem
      write(iout) MSOO,nclosd, maxorb, lsdim, ncfg
      write(iout) ljclsd(1:nclosd),ljcomp(1:maxorb)
      write(iout) el(1:nwf)
      write(iout) lmax,nnel,skip

*
      NWFD = MAXORB

      if(qiajcmp.ne.0) call dalloc(qiajcmp,nwfd);
      if(qljcomp.ne.0) call dalloc(qljcomp,nwfd);
      if(qnjcomp.ne.0) call dalloc(qnjcomp,nwfd); 
      if(qiajcld.ne.0) call dalloc(qiajcld,nwcd);
      if(qljclsd.ne.0) call dalloc(qljclsd,nwcd);
      if(qsignfa.ne.0) call dalloc(qsignfa,ncfg);
      if(qacmult.ne.0) call dalloc(qacmult,nwfd);
      if(qnelcsh.ne.0) call dalloc(qnelcsh,8*ncfg);
      if(qnocorb.ne.0) call dalloc(qnocorb,8*ncfg);
      print *, 'finished deallocation'
*
*
*     .. copy information needed by CI before deallocating
* 
      if (.not. skip) then
        call alloc(iqlsp,ncfg,4)
        do i = 1,ncfg
          m = noccsh(i)
          ls = j1qnrd(2*m-1,i)/64
          lli = mod(ls,64) -1
          ksi = (ls/64) - 1
          lsp(i) = (64*lli +ksi)*64
         end do
      end if
      print *, 'finished copying information for CI'
      print *, 'skip =',skip
*     if (.not. skip) then
*
*  *****  Determine the list of TERMS
*
        NTERM = 0
        DO 650 I = 1,NCfg
          lt = mod(lsp(i),64)
          ls = lsp(i)/64
          ksi = mod(ls,64)
          lli = ls/64
*	  print *, i,ncfg,lt
          IF (LT .EQ. 0) THEN
            NTERM = NTERM + 1
            LSP(I) = lsp(i)+ 2*nterm
            INDEX(NTERM) = I
*	    print *, i, nterm, index(nterm)
            DO 660 J = I+1,NCfg
              lt = mod(lsp(j),64)
              ls = lsp(j)/64
              ksj = mod(ls,64)
              llj = ls/64
*	      print *, i,j
              IF (lt .EQ. 0 .AND.  LLi.EQ.LLJ .AND. KSI.EQ.KSJ) THEN
                LSP(J) = lsp(j)+2*nterm
              END IF
 660        CONTINUE
          END IF
 650    CONTINUE
      print *, 'Finished determining list of terms'
*      end if
       write(iout) lsp(1:ncfg), jptr(1:ncfg)
       write(iout) nterm, index(1:nterm)
*      print *, lsp(1:ncfg), jptr(1:ncfg)
*      print *, nterm, index(1:nterm)
*     call dalloc(qnoc,ncfg)
*     call dalloc(qj1,15*ncfg)
      close(unit=IREAD,iostat = ios); 
      close(unit = IOUT,iostat = ios); 
      close(unit = 11,iostat = ios); 
      close(unit = 12,iostat = ios);
      close(unit = 50,iostat = ios);

      return
      END
