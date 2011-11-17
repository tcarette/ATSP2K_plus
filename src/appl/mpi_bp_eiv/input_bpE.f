      subroutine inp_case(QREL,QMASS,ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW,
     :      REL,MASS,iscw,iread)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220, NWD=60)
      LOGICAL REL, SKIP

      CHARACTER*26 STRING,ATOMC,ATOML,ATOMJ,ATOMW,ATOMNEW
      CHARACTER QREL, QMASS

1     WRITE(ISCW,'(/A,A)') ' Enter ATOM, relativistic (Y/N)',
     :               ' with mass correction (Y/N)'
      READ(IREAD,'(A)') STRING
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
         READ(IREAD,'(A1)') QMASS
         MASS = 1
         IF (QMASS.EQ.'S' .OR. QMASS.EQ.'s') MASS = 2
      END IF

      return
      end



      subroutine inp_eig(njv,MAXJ,MINJ,leigen,nume,iscw,iread,onlydvd,
     :     termsh,indx,maxnum)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,NWCD=20,NWD=60,NTERMD=31,LINT=500)
      CHARACTER ANS*2, SYMBOL*11
      DATA SYMBOL/'SPDFGHIKLMN'/
      POINTER (IQLSP, LSP(1))

      character str*72
      logical inputOK, leigen(NTERMD,NTERMD),onlydvd
      integer nume(ntermd),njeig(ntermd),indx(ntermd)

      inputOK = .false.
      do while (.not.inputOK)
         WRITE(ISCW,*) ' Enter Maximum and minimum ',
     :                           'values of 2*J'
         READ(IREAD,*) MAXJ,MINJ
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
        read(iread,'(A)') str
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

