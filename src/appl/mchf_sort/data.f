*r
*     ------------------------------------------------------------------
*    3-3       D A T A
*     ------------------------------------------------------------------
*
*       Data concerning the number of configurations (NCFG), the number
*   and type of electrons in each  configuration,  as  well  as  data
*   associated with the energy expression are read and stored.
*
*
      SUBROUTINE DATA(IVAR, eigst_weight,cf_tot)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER (NWD=60,NWC=20,NOD=220,NOFFD=800,MTERM=20,MEIG=20)
      DIMENSION IVAR(NWD)
*
      CHARACTER     EL*3,ATOM*6,TERM*6, couple(15)*3
      LOGICAL       STRONG
      CHARACTER*3   EL1,EL2,ELORT(10,2),
     :              ELI(8),ANS*1,STRING*64,BUFFER*8,LIST*180,HEADER*29
      INTEGER       IQ(8),nterm 
      INTEGER       IN,OUT,PRI,OUC,OUF,OUH,OUD
      COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,IUD,OUC,OUF,OUD,OUH,ISCW
      COMMON/LABEL/ EL(NWD),ATOM,TERM
      LOGICAL       FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
      COMMON/TEST/  FAIL,OMIT,EZERO,REL,ALL,TRACE
      COMMON/WAVE/  EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
      COMMON/PARAM/ H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,
     :              NCFG,IB,IC,ID,
     :              D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,
     :              NSCF,NCLOSD,RMASS
      COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
*
      POINTER (QWT,WT(1)),(QWP,WP(1)),  (IQWPTR,W(1))
      COMMON/CFGS/ ETOTAL,QWT,QWP,IQWPTR
*
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :     (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :     (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :      QMETH,QIEPTR
*
      POINTER (pkval,kval(1)),(pvalue,value(1)),
     :        (pih,ih(1)),(pjh,jh(1)),
     :        (pcoeff,coeff(1)),(pnijptr,nijptr(1)),
     :        (pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PVALUE,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
*
      COMMON ZZ(NWD),IND(NWD),IELI(8)
   
      double precision eigst_weight(meig,mterm)
      logical leigen(meig,mterm), lguess
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm),
     :        nze_max(mterm)
      character term_bl(mterm)*3
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess,nze_max
*
      integer 	ice, ipos_B, ipos_E, cf_tot(mterm)
      character     line*72, t(60)*3, str*72, tmp*20, tmpa*5, tmpb*5
      character*5  who, char*1
      eigst_weight = 0.0
      rewind(29) 
      read(29,'(i4,A64)') nclosd, string
      read(29,'(18(1x,A3))') (el(i), i=1,nclosd)
      read(29,'(i4,A64)') nowf, string
      nwf = nowf + nclosd
      IF(NWF .GT. NWD) THEN
         write(ISCW,'(A,I4)') 'NWF too large for this code: MAX = ',NWD
         stop
      END IF

      nwf1 = nclosd + 1
      read(29,'(18(1x,A3))') (el(i), i=nwf1,nwf)
      read(29,'(I3,2I8,3X,A5)') nblock,idim,ncodim,who
      if (who .ne. 'snonh') then
         write(0,*) 'Data not suitable for simultaneous optimization'
         stop
       end if
   
*    read cfg.h 
      do nb = 1,nblock
          READ(29,'(2X,A3,I6,I10,I15,I10)') 
     :      term_bl(nb),ncfg_bl(nb),NZE_bl(nb),cf_tot(nb),nze_max(nb)
            term_bl(nb) = trim(adjustl(term_bl(nb)))
            !print*, '::', term_bl(nb), '::';
      end do
*     set the TERM value
      if (nblock .eq. 1) then
         TERM = term_bl(1)
      else
         TERM = 'AV_E'
      end if

*     determine eigenvalues
      leigen = .false.
      write(iscw,*) 'cfg.inp has configurations for ',nblock , ' terms'   
      write(iscw,*)
      write(iscw,*) 'Enter eigenvalues and weights: ',
     :        'one line per term, eigenvalues with weights'
      write(iscw,*) 'in parenthesis and separated by commas, ',
     :        'default is 1.0'

      do nb = 1,nblock
         write(iscw,*) term_bl(nb)
         read(in,'(A)') str
         nch = 1
         len = len_trim(str)
*>>>>
         do while (nch <= len)
            ipos = index(str(nch:len),',')
            if (ipos .eq. 0) ipos = len+2-nch
            read (str(nch:nch+ipos-2),*) tmp 
            ipos_B = index(tmp,'(')

            if (ipos_B.eq.0) then
               read (str(nch:nch+ipos-2),*) keigv
               eigst_weight(keigv,nb) = 1.0
            else
               ipos_E = index(tmp,')')
               read (str(nch:nch+ipos_B-2),*) keigv
               read (tmp(ipos_B+1:ipos_E-1),*) eigst_weight(keigv,nb)
            end if 

            if (keigv .gt. meig) then
              write(0,*) 'Too high an eigenvalue requested:',
     :                   'Maximum for current dimensions is',meig
              stop
            end if

            leigen(keigv,nb) = .true.
            nch = nch + ipos 
         end do
*>>>>
         nume(nb) = keigv
      end do

    7 FORMAT(A3,F6.0,I3,I3,F3.1)
*  *****  READ 'ATOM' CARD
*
      NORT = 1
      ID = 0
      icc = 1
      ib = 0
*
      sum_eigst_weights = 0.0
      do ix=1,meig
        do iy =1,mterm
          sum_eigst_weights = sum_eigst_weights + eigst_weight(iy,ix)
        end do
      end do
      do iy = 1, mterm
         do ix = 1, meig
            eigst_weight(ix,iy) = eigst_weight(ix,iy)/sum_eigst_weights
         end do
      end do

      call alloc_mem(nwf)

      do i= 1,nwf
        DPM(I) = D10
        IEPTR(I) = 0
        IND(I) = 0
      end do

*
*  *****  DETERMINE NCLOSD SHELLS
*

         SS = D0
         do i=1,nclosD
             VARIED(I) = .TRUE.
             J = 3
             IF (EL(I)(1:1) .NE. ' ') J = 2
             L(I) = LVAL(EL(I)(J:J))
             N(I) = ICHAR(EL(I)(J-1:J-1)) - ICHAR('1') + 1
             IFULL = 2*(2*L(I)+1)
             S(I) = SS + IFULL/2
             SS = SS + IFULL
             METH(I) = 1
             ACC(I) = D0
             IND(I) = 0
             SUM(I) = 4*L(I)+2
             IF (IUF .NE. 0)  IND(I) = -1
          end do

*
*  *****  DETERMINE THE OTHER ELECTRONS AND ORDER
*

      WRITE(ISCW,19) NWF,(EL(J),J=1,NWF)

13     continue 

   19 FORMAT(/' There are ',I3,' orbitals as follows:'/(1X,18(1X,A3)))
21    WRITE(ISCW,'(A,A)') ' Enter orbitals to be varied:',
     :     ' (ALL,NONE,SOME,NIT=,comma delimited list)'
      READ '(A)', LIST
      IF (LIST(1:3) .EQ. 'ALL' .OR. LIST(1:3) .EQ. 'all') THEN
         NIT = NWF
         do 209 i = 1,nit
            ivar(i) = i
209      continue
      ELSE IF (LIST(1:3) .EQ. 'LLA' .OR. LIST(1:3) .EQ. 'lla') THEN
         NIT = NWF
         do 210 i = 1,nclosd
            ivar(i) = i
210      continue
         do 211 i = nclosd+1,nit
            ivar(i) = nit - i + nclosd +1
211      continue
      ELSE IF (LIST(1:4).EQ.'SOME' .OR. LIST(1:4).EQ.'some') THEN
         NIT = NWF - NCLOSD
         do 212 i = 1,nit
            ivar(i) = nclosd + i
212      continue
      ELSE IF (LIST(1:4).EQ.'NONE' .OR. LIST(1:4).EQ.'none') THEN
         NIT = 0
      ELSE IF (INDEX(LIST,'=') .NE. 0) THEN
         J = INDEX(LIST,'=')
         READ(LIST(J+1:),'(I2)') NIT
         IF (NIT .GT. NWF) THEN
            WRITE(ISCW,*)'NIT greater than NWF: replaced by NIT=NIT/10'
            NIT = NIT/10
         END IF
         do 214 i = 1,nit
            ivar(i) = nwf-nit + i
214      continue
      ELSE
         NIT = 0
         J = 1
22       NEXT = INDEX(LIST(J:),',')

*
*        ***  Search for last electron label which need not be followed
*             by a comma
*
         IF (NEXT .EQ. 0 .AND. LIST(J:J+2) .NE. '   ')
     :       NEXT = INDEX(LIST(J+1:),' ') + 1
         IF (NEXT .NE. 0) THEN
            IF (NEXT .EQ. 4) THEN
               EL1 = LIST(J:J+2)
            ELSE IF (NEXT .EQ. 3) THEN
               EL1 = ' '//LIST(J:J+1)
            ELSE
               WRITE(ISCW,*)'Electron labels must be separated
     :         by commas;'
               WRITE(ISCW,*)'each label must contain 2 or 3 characters'
               GO TO 21
            END IF
            CALL EPTR(EL,EL1,ii,*98)
            NIT = NIT+1
            ivar(nit) = ii
            j = j+next
Cmrg        IF (J .LT. 72) GO TO 22
            IF (J .LT. 180) GO TO 22
Cmrg
         END IF
      END IF
*
      IB = NWF - NIT + 1
      do i = nclosd+1,nwf
            S(I) = SS
            METH(I) = 3
            ACC(I) = D0
            IND(I) = 0
            IF (IUF .NE. 0)  IND(I) = -1
            VARIED(I) = .TRUE.
            J = 2
            IF (EL(I)(1:1) .EQ. ' ') J = 3
            L(I) = LVAL(EL(I)(J:J))
            N(I) = ICHAR(EL(I)(J-1:J-1)) - ICHAR('1') + 1
       end do
32     write(iscw,*) 'Enter those that are spectroscopic'
       read '(A)', list
       j = 1
       do while(list(j:j+2) .ne. '   ')
         NEXT = INDEX(LIST(J:),',')
         if (next == 0) next = index(list(j+1:), ' ')+1 
         if( next .eq. 4) then
            el1 = list(J:J+2)
         else if (next .eq. 3) then
            el1 = ' '//list(j:j+1)
         else
            write(iscw,*) 'Electron labels must be separated',
     :                    ' by commas'
            go to 32
         end if
         call eptr(el,el1, ii, *98)
         meth(ii) = 1
         j = j+next
       end do
      WRITE(ISCW,'(/,A)') ' Default electron parameters ? (Y/N) '
      READ '(A)', ANS

      IF ( ANS.NE.'Y' .AND. ANS.NE.'y' .AND. NIT.NE.0) then
	 write(iscw,'(A)')
     :     ' S, IND, METH, ACC for electrons',
     :     'to be varied (free-format)'
	 do ii = 1,nit
	    i = ivar(ii) 
	    write(iscw,*) el(i) 
	    read(in,*) s(i),ind(i),meth(i),acc(i)
         end do
       end if

* .. define all orbitals of the same l to be orthogonal

      DO 38 I = 2,NWF
        DO 39 J = 1,I-1
           IF (L(I) .EQ. L(J)) THEN
               C = 1.D-5
               IF ((I.LE.NCLOSD .AND. J.LE.NCLOSD))  C = 1.D-10
               call EIJSET(I,J,C)
               CALL EIJSET(J,I,C)
           END IF
  39    CONTINUE
  38  CONTINUE

 50   NO = (NOD)
      ND = NO - 2
      REL = .FALSE.
      STRONG = .FALSE.
      IF (NCFG_bl(1) .GT. 1) STRONG = .TRUE.
      WRITE(PRI,62) ATOM,TERM,Z,(EL(I),4*L(I)+2,I=1,NCLOSD)
62    FORMAT(1H1///9X,33HHARTREE-FOCK WAVE FUNCTIONS FOR  ,2A6,4H Z =,
     1   F5.1//14X,'CORE = ',5(A3,'(',I2,')'))
      WRITE(PRI,'(//11X,A,37X,A//)') 'CONFIGURATION','WEIGHT'
      OMIT = .NOT. STRONG
*
      WRITE(PRI,71)
71    FORMAT(//9X,10HINPUT DATA/9X,10H----- ----//13X,13HWAVE FUNCTION,
     1   11H  PROCEDURE/17X,22HNL  SIGMA METH ACC OPT///)
      DO 79 I = 1,NWF
      WRITE(PRI,78) I,EL(I),N(I),L(I),S(I),METH(I),ACC(I),IND(I)
78    FORMAT(I8, 2X,A3,2I3,F7.1,I4,F4.1,I4)
79    CONTINUE

*
*  .. INITIALIZE radial functions
*
      CALL WAVEFN

* ..  initialize eigenvectors, if available
      if (lguess) then
        ievstart = 0
        do nb = 1,nblock
          write(0,*) 'Reading eigenvectors for term: ',term_bl(nb)
          ncfg = ncfg_bl(nb)
          read(60,'(I5)') niv
          niv_bl(nb) = niv
*         .. read the eigenvectors for the block
          do j = 1,min(niv,nume(nb))
            read(60,'(I8,F15.9/(7F11.7))') 
     :               iv,eiv, (eigvec(ievstart+(j-1)*ncfg+i),i=1,ncfg)
          end do
          write(0,'(A,I10,F15.8)') 'Eigenvalue:',iv ,eiv
          ievstart = ievstart +nume(nb)*ncfg
        end do
	rewind(unit=60)
      end if
      
* ..  initialize integral data
      intptr(1:6) = 0
      call spintgrl(cf_tot)
*
      RETURN
98    WRITE(ISCW,*) ' Case must match as well as position of',
     :                  ' imbedded blanks'
      WRITE(ISCW,*) ' For 3rd character of label to be blank',
     :                    ' follow blank with comma'
      GO TO 13
      END
