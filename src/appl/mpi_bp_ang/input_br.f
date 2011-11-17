
!## read the name of the atom 
      subroutine inp_atom(ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW,
     :      REL,MASS,iscw,in)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220, NWD=60)
      LOGICAL REL, SKIP

      CHARACTER*26 STRING,ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW
      CHARACTER QREL, QMASS

1     WRITE(ISCW,'(/A,A)') ' Enter ATOM, relativistic (Y/N)',
     :               ' with mass correction (Y/N)'
      READ(in,'(A)') STRING
      I = INDEX(STRING,',')
      IF (I .eq. 0)     THEN
         WRITE(ISCW,*) ' Separate with commas'
         GO TO 1
      ELSE IF (I .GT. 22) THEN
         WRITE(ISCW,*) ' ATOM name can have at most 22 characters'
         GO TO 1
      END IF

      QREL  = STRING(I+1:I+1)
      QMASS = STRING(I+3:I+3)
      ATOMC = STRING(1:I-1)//'.c'
      ATOML = STRING(1:I-1)//'.l'
      ATOMJ = STRING(1:I-1)//'.j'
      ATOMW = STRING(1:I-1)//'.w'
      ATOMNEW = STRING(1:I-1)//'.new'
      IF (QREL.EQ.'N' .OR. QREL.EQ.'n') REL = .FALSE.
      MASS = 0
      IF (QMASS.EQ.'Y' .OR. QMASS.EQ.'y') THEN
         WRITE(ISCW,'(A)') ' Gradient or Slater form? (G/S):'
         READ(in,'(A1)') QMASS
         MASS = 1
         IF (QMASS.EQ.'S' .OR. QMASS.EQ.'s') MASS = 2
      end if

      return
      end

!## select Gradient or Slater
      subroutine inp_gs(QREL,REL,MASS,QMASS,in,iscw)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220, NWD=60)
      LOGICAL REL
      CHARACTER QREL, QMASS

      IF (QREL.EQ.'N' .OR. QREL.EQ.'n') REL = .FALSE.
      MASS = 0
      IF (QMASS.EQ.'Y' .OR. QMASS.EQ.'y') THEN
         WRITE(ISCW,'(A)') ' Gradient or Slater form? (G/S):'
         READ(in,'(A1)') QMASS
         MASS = 1
         IF (QMASS.EQ.'S' .OR. QMASS.EQ.'s') MASS = 2
      END IF

      return
      end

!## convert y or n to true or false     
      subroutine inp_yn(yn,iscw,in)
      
      logical yn
      character*3  string 
      
      WRITE(ISCW,*) ' Restarting (Y/y) ?'
      READ(in,'(A)') STRING
      IF(STRING(1:1).EQ.'Y'.OR.STRING(1:1).EQ.'y') THEN
          yn = .true.
      end if 

      return 
      end

!## get the type of calculation: 0, 1 or 2 and S-O,SOO,SS,OO
      subroutine inp_type(irel,in,iscw,ISPORB,ISOORB,ISPSPN,
     :     IORBORB,ICOLOM,iwrite,IFULL,ICSTAS,skip)      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      logical skip
      integer iscw, in, irel 
      character*2 ans

  105 FORMAT (49H ISPORB=0 AND ISOORB=1 CAUSES THE PROGRAM TO FAIL,
     :  34H BECAUSE THE BLUME WATSON FORMULAE,/
     :  50H CANNOT BE USED FOR CLOSED SUBSHELLS.  TO OVERCOME,
     :  34H THIS, THE CODE HAS SET ISPORB = 1//)

      WRITE(ISCW,'(A/A/A/A)') ' Indicate the type of calculation ',
     : ' 0 => non-relativistic Hamiltonian only;',
     : ' 1 => one or more relativistic operators only;',
     : ' 2 => non-relativistic operators and selected relativistic:  '
        READ(5,*) IREL

      IF (IREL.NE.1) ICOLOM=1
      skip = .true.
      if (irel .ne. 0) then
         skip = .false.
      end if
      IFULL = 0; ICSTAS = 1; ISPORB = 0; ISOORB = 0; ISPSPN = 0
      IORBORB=0
      IF (IREL .NE. 0) THEN
         ISPORB = 1
         ISOORB = 1
         ISPSPN = 1
         IORBORB = 1
         WRITE(ISCW,'(A)') ' All relativistic operators ? (Y/N) '
         READ(5,'(A2)') ANS
         IF ( ANS .EQ. 'N ' .OR. ANS .EQ. 'n ') THEN
           WRITE(ISCW,'(A)')
     :      'Spin-orbit,Spin-other-orbit,Spin-spin,Orbit-Orbit (0/1)'
             READ(5,*) ISPORB,ISOORB,ISPSPN,IORBORB
         END IF
         IF(ISPORB.EQ.0.AND.ISOORB.NE.0) THEN
           ISPORB = 1
           WRITE (IWRITE,105)
         END IF
      END IF
  
      return
      end


!##  get other: Interactions, nzero, istrict, istart
      subroutine inp_interact(new,nzero,istrict,istart,ISCW,NCFG)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      logical yn
      character*2 ans
   
      WRITE(ISCW,'(A)') ' All Interactions? (Y/N): '
      READ (5,'(A2)') ANS
      IF (ANS .NE. 'Y' .AND. ANS .NE. 'y') THEN
        WRITE(ISCW,'(A,I8,A/A)') ' Of the ',NCFG,
     :         ' configurations, how many at the end are new? ',
     :         ' How many configurations define the zero-order set?'
        READ (5,*) NEW,NZERO
        IF (NEW .EQ. 0) NEW = NCFG
        ISTART = NCFG - NEW + 1
        IF (NZERO .eq. 0) NZERO = NCFG
        ISTRICT = 1
        WRITE(ISCW,*)  ' Restricted Two-body interactions? (Y/N); '
        READ (5,'(A2)') ANS
        IF (ANS .NE. 'Y' .and. ANS .NE. 'y' ) ISTRICT = 0
      END IF

      return 
      end

!## input reqested eigenvalues
      subroutine inp_eig(njv,MAXJ,MINJ,leigen,nume,iscw,in,onlydvd,
     :     termsh,indx,maxnum)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,NWCD=20,NWD=60,NTERMD=31,LINT=500)
      CHARACTER ANS*2, SYMBOL*11
      DATA SYMBOL/'SPDFGHIKLMN'/
      POINTER (IQLSP, LSP(1))

      character str*72
      logical inputOK, leigen(NTERMD,NTERMD),onlydvd
      integer nume(ntermd),njeig(ntermd),termsh(ntermd),indx(ntermd)

      inputOK = .false.
      do while (.not.inputOK)
         WRITE(ISCW,*) ' Enter Maximum and minimum ',
     :                           'values of 2*J'
         READ(in,*) MAXJ,MINJ
         njv = (maxj-minj)/2 + 1
         if (njv .gt. 20) then
           write(iscw, *) ' Current dimensions ',
     :                    'allow upto 20 J-values: Re-enter'
         else if (maxj < minj) then
           write(iscw, *) ' max must be greater than max: Re-enter'
         else
           inputOK = .true.
         end if
      end do
      ij = 1
      maxnum = 1
      nb = 1
      leigen = .false.
      write(iscw,*)
      write(iscw,*) 'Enter eigenvalues: ',
     :      'one line per term, eigenvalues separated by commas'
      do nj = maxj,minj,-2
        write(iscw,'(A,I3)') '2*J =',nj
        read(in,'(A)') str
        nch = 1
        len = len_trim(str)
        do while (nch <= len)
          ipos = index(str(nch:len),',')
          if (ipos .eq. 0) ipos = len+2-nch
          read (str(nch:nch+ipos-2),*) keigv
          if (keigv .gt. ntermd) then
            write(0,*) 'Too high an eigenvalue requested:',
     :                   'Maximum for current dimensions is',ntermd
            stop
          end if
          leigen(keigv,nb) = .true.
          nch = nch + ipos
        end do
        nume(nb) = keigv
        njeig(nb) = nume(nb)
        if (nume(nb).gt. maxnum) maxnum = nume(nb)
        nb = nb + 1
      end do

      return
      end






