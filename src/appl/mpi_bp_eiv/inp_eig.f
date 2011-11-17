      subroutine inp_eig(njv,MAXJ,MINJ,leigen,nume,iscw,iread,onlydvd,
     :     maxnum)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,NWCD=20,NWD=60,NTERMD=31,LINT=500)
      CHARACTER ANS*2, SYMBOL*11
      DATA SYMBOL/'SPDFGHIKLMN'/
      POINTER (IQLSP, LSP(1))

      character str*72
      logical inputOK, leigen(NTERMD,NTERMD),onlydvd
      integer nume(ntermd),njeig(ntermd)

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

