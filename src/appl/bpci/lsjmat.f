*
*     ------------------------------------------------------------------
*             L S J M A T
*     ------------------------------------------------------------------
*
*     Read the non-fine structure file, and the fine structure files 
*   (if any), form the interaction matrix, find eigenvalues and 
*   corresponding eigenvectors.  This routine does not use the
*   iselec option of Davidson.  If it does, the number is restricted
*   to NTERMD.  The matrix may be either on disk on in memory
*   as defined by the idisk parameter.  The calculations may be
*   an LS (iLS=1) or LSJ (iLS=0) calculation.
* 
      SUBROUTINE LSJMAT(JJ, nume,skip,nclosd,nwf,ncfg,nze,ec,
     :                  onlydvd,jhere,lhere,iounew,leigen)
*	WITH RESTARTING
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NTERMD=31,NWD=60)
      LOGICAL  SKIP
      logical onlydvd,jhere,lhere
      logical leigen(ntermd)
*
      POINTER (qh,h(ncfg,3)),(qjan,jan(ncfg,3))
      COMMON /buffer/qh, qjan, nrow(3), iflag(3)
      POINTER (IQLSP,LSP(1))
      COMMON /NCFGS/ IQLSP,index(ntermd),IJK,flsj(ntermd,ntermd,2),
     :               termsh(ntermd),nterm
*
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk
*
      INTEGER cptr
      POINTER (iqhmx,hmx(1)),(iqjptr,jptr(1)),(iqimx,imx(1)),
     :        (iqeigvec,eigvec(1)),(iqkval,kval(1)),(iqcptr,cptr(1)),
     :        (qtm,tm(1)),(qtp,tp(1)),(qdiag,hii(1)),(qiwork,iwork(1))
      COMMON /MX/iqhmx, iqjptr, iqimx, iqeigvec, iqkval, iqcptr,
     :           qtm,qtp,qdiag,qiwork
*
      CHARACTER config*66,elc(8)*3,couple(15)*3
      INTEGER   q(8), ipos(3),iselec(NTERMD)
      LOGICAL iupper,hiend
      DATA iupper/.false./
      DATA crite,critc,critr,trhold,ortho/
     :       1.d-14,1.d-8,1.d-8,1.d0,1.d0/
      dimension Ssms(nume),g_J(nume),g_JLS(nume)
*     ------------------------------------------------------------------
*    
*     on entry
*     --------
*
*     N   - dimension of the matrix
*     JJ  - 2J-value (not needed in this restricted application)
*     Nume- number of desired eigenvalues/eigenvectors
*     LS  - type of calculation (LS or LSJ -- not needed here)
*     NZE - size of the non-zero array for the interaction matrix
*
*     ------------------------------------------------------------------
* 
*       .. clear the interaction matrix
        call dfill(nze,0.d0,hmx,1)
	print *, 'Entering LSJMAT with 2J =',jj, ' NUME =',nume
*
      If (iLS .eq. 0) then
*
*  *****  CLEAR TABLE OF J-DEPENDENT FACTORS
*
        DO 3 II = 1,NTERM
	  I = INDEX(II)
	  ls = lsp(i)/64
	  ksi = mod(ls,64)
	  lli = ls/64
	  DO 4 JJJ = 1,NTERM
	    J = INDEX(JJJ)
	    ls = lsp(j)/64
	    ksj = mod(ls,64)
	    llj = ls/64
            PHASE =
     :         (-1)**((LLI+KSJ-JJ+LLI+KSI-JJ+LLJ+KSJ-JJ)/2) 
            CALL GRACAH(LLJ,KSJ,LLI,KSI,JJ,2,W1) 
            CALL GRACAH(LLJ,KSJ,LLI,KSI,JJ,4,W2) 
            FLSJ(II,JJJ,1) = PHASE*W1
	    FLSJ(II,JJJ,2) = PHASE*W2
   4      CONTINUE
   3    CONTINUE
*
*       .. determine which configuration states are "IN" (in=1)
*       
        do 10 i = 1,ncfg
	  ls   = lsp(i)/64
	  ksi   = mod(ls,64)
	  lli   = ls/64
	  if (jj .ge. ABS(LLi-KSi) .and. jj .le. lli+ksi) then
            lsp(i) = 2*(lsp(i)/2) + 1
	  else
	    lsp(i) = 2*(lsp(i)/2)
	  end if
   10   continue
*   
*    ***** Form the interaction matrix
*   
        nij = 0
	istart = 0
	shift = 0.d0
        rewind(iouhn)
        rewind(iouhz)
        rewind(iouhs)
	if (idisk .eq. 1) rewind(iouhm)
        do 30 j = 1,ncfg
	  if (idisk .eq. 1) then
	    nij = 0
            call dfill(nze,0.d0,hmx,1)
	  end if
	  if (mod(lsp(j),2) .eq. 0) then
*            .. csf does not contribute
	    nij = nij+1
	    imx(nij) = j
	    hmx(nij) = 1.d+5
	    read(iouhn) jb
	    read(iouhz) jb
	    read(iouhs) jb
 	  else
            ipos(1) = 0
            ipos(2) = 0
            ipos(3) = 0
	    ntermj = mod(lsp(j),64)/2
	    read(iouhn) jb,m1,(h(i,1),i=1,m1),(jan(i,1),i=1,m1)
	    read(iouhz) jb,m2,(h(i,2),i=1,m2),(jan(i,2),i=1,m2)
	    read(iouhs) jb,m3,(h(i,3),i=1,m3),(jan(i,3),i=1,m3)
            call advance(1,m1,ipos,ncfg)
	    call advance(2,m2,ipos,ncfg)
	    call advance(3,m3,ipos,ncfg)
*           .. while there is more data
*              find which file has lowest row
  100       iipos = min(ipos(1),ipos(2),ipos(3))
	    if (iipos.le. ncfg) then
*             .. there is more data
	      nij = nij + 1
	      irow = min(nrow(1),nrow(2),nrow(3))
	      if (nij .gt. nze) then
	        newlen = nze + nze
	        call realloc(iqhmx,nze,newlen,8)
	        call realloc(iqimx,nze,newlen,4)
*               ..clear new memory
                call dfill(newlen-nze,0.d0,hmx(nze+1),1)
	        nze = newlen
              end if
	      if (nrow(1) .eq. irow) then
	        hmx(nij) = hmx(nij) + h(ipos(1),1)
	        imx(nij) = jan(ipos(1),1)
	        call advance(1,m1,ipos,ncfg)
	      end if
	      if (nrow(2) .eq. irow) then
	        ntermi = mod(lsp(irow),64)/2
	        hmx(nij)= hmx(nij)+h(ipos(2),2)*flsj(ntermi,ntermj,1)
	        imx(nij) = jan(ipos(2),2)
	        call advance(2,m2,ipos,ncfg)
	      end if
	      if (nrow(3) .eq. irow) then
	        ntermi = mod(lsp(irow),64)/2
	        hmx(nij)= hmx(nij)+h(ipos(3),3)*flsj(ntermi,ntermj,2)
	        imx(nij) = jan(ipos(3),3)
	        call advance(3,m3,ipos,ncfg)
	      end if
	      go to 100
	    end if
*
*         .. determine term and correct diagonal
*
   	    ls1 = lsp(j)/64
 	    do ii = 1,nterm
 	      jjj = index(ii)
 	      ls2 = lsp(jjj)/64
 	      if (ls1 .eq. ls2) go to 150
 	    end do
150         if (ii .gt. nterm) then
	      write(iscw,*)  'Term not found', mod(ls1,64),ls1/64
	      write(iscw,*)  'Diagonal not corrected'
	    else
*	      write(iscw,*) 'The term for column ',j,'is ',ii,termsh(ii)
	      hmx(istart+1) = hmx(istart+1) - termsh(ii)
	    end if
	  end if
	  if (shift .eq. 0.d0 .and. hmx(istart+1) .lt. 0.d0)
     :      shift = hmx(istart+1)
	  hmx(istart+1) = hmx(istart+1)-shift
	  hii(j) = hmx(istart+1)
	  if (idisk .eq. 1) then
	    write(iouhm) j,nij,(hmx(i),i=1,nij),(imx(i),i=1,nij)
	  else
	    m1 = nij-istart
	    jptr(j) = nij
	    istart = nij
	  end if
 30     continue
      else
*       .. we have an LS case
	if (idisk .eq. 1 ) then
*     
*      ***** Form the diagonal matrix. search for lowest
*     
          rewind(iouhn)
          do 120 j = 1,ncfg
            ipos(1) = 0
     	    read(iouhn) jb,m1,hii(j)
  120     continue
          shift = hii(1)
*     
*        .. scale the diagonal elements for accuracy
          do 130 j = 1,ncfg
     	    hii(j) = hii(j) - shift
  130     continue
	else
*   
*    ***** Form the interaction matrix
*   
          nij = 0
          rewind(iouhn)
          do 230 j = 1,ncfg
	    read(iouhn) jb,m1,(hmx(nij+i),i=1,m1),(imx(nij+i),i=1,m1)
	    nij = nij + m1
	    jptr(j) = nij
  230     continue
*         .. compute diagonals
	  shift   = hmx(1)
	  hmx(1) = 0.0d0 
          hii(1) = 0.d0  
          do 232 i = 2,ncfg
            hii(i) = hmx(jptr(i-1)+1) - shift
	    hmx(jptr(i-1)+1) = hii(i)
  232     continue
	end if
      end if
*
*  ***** COMPUTE THE EIGENVALUES AND EIGENVECTORS
*   
      if (iLS .eq. 1) then
        iouv = ioul
        jhere = lhere
      else
        iouv = iouj
        lhere = jhere
      end if

!      lim = min(20,ncfg)
      lim = min(2*nume+40,ncfg)
      iworksz = (2*ncfg+lim+nume+10)*lim + nume
      iiwsz = 6*lim + nume
      n = ncfg
      print *, 'LSJMAT with idisk=', idisk, 'Nze =', nze
*     if (idisk .eq. 1) then
*	Nze = jptr(n)
*     else 
*	nze = ncfg
*     end if
c      ilow = 1
c      ihigh = nume
      ilow = -1
      ihigh = -1

      ie = 1
      do ic1 = 1,ntermd
         if (leigen(ic1)) then
           iselec(ie) = ic1
           ie = ie + 1
         end if
      end do
      iselec(ie) = -1

      niv = 0
!      maxiter = max(30*nume,100)
      maxiter = max(30*nume,500)
c      mblock = nume
      mblock = ie -1

*
*      print *, 'hii:', (hii(i),i=1,n)
*      print *, 'hmx:', (hmx(i),i=1,nze)
*      print *, 'jptr:', (jptr(i),i=1,n)
*      print *, 'imx:', (imx(i),i=1,nze)

      CALL dvdson(hmx,imx,jptr,iupper,nze,tm,tp,
     :        Ncfg,lim,hii,ilow,ihigh,iselec,niv,mblock,
     :        crite,critc,critr,ortho,maxiter,
     :        eigvec,iworksz,iwork,iiwsz,hiend,nloops,nmv,ierr)

      if (ierr .ne. 0) then
        write(0,*) ' Dvdson returned ierr=',ierr
	write(0,*) ' ...  continuing'
      end if
      ist = n*nume
      write(iscw,*)
      write(iscw,*) 'Summary of Davidson Performance'
      write(iscw,*) '==============================='
      write(iscw,*) 'Number of Iterations: ', nloops
      write(iscw,*) 'Shifted Eigval''s:', (eigvec(ist+i),i=1,nume)
      write(iscw,*) 'Delta Lambda:    ', (eigvec(ist+nume+i),i=1,nume)
      write(iscw,*) 'Residuals:       ', (eigvec(ist+2*nume+i),i=1,nume)
      WRITE (iscw, '(//I6,A)' ) nume,' Eigenvalues found'
      if (onlydvd) iouv=iounew
      WRITE (iouv, '(//A8,I4,2X,A8,I4)' ) '  2*J = ',JJ,'NUMBER =',
     :           mblock
c     :           nume
*
      rewind(iuc)
         call brci_zeeman(jj,ncfg,nume,eigvec,g_J,g_JLS,iuc)
         rewind (iuc)
         call brci_eig(ncfg,jj,nume,leigen,eigvec,Ssms,g_J,g_JLS,
*    :        eigvec(ncfg*nume+1),EC,shift,iouj,iuc)
     :        eigvec(ncfg*nume+1),EC,shift,iouv,iuc)

      print *, 'Finished with Davidson'
!      DO 40 K = 0,nume-1
!*
!*  *****  SEARCH FOR THE LARGEST COMPONENT IN THE EIGENVECTOR FOR
!*         LABELLING PURPOSES
!*
!      if (leigen(k+1)) then
!      VMAX = 0.D0
!      JMAX = 0
!      istrt = n*nume +1
!      EIGVAL = EIGVec(istrt+k) + EC + shift
!      DO 90 J = 1,N
!        ABSEIG = DABS(EIGVEC(J+K*N))
!        IF ( ABSEIG .LT. VMAX ) GO TO 90
!        JMAX = J
!        VMAX = ABSEIG
!90    CONTINUE
!      IF (EIGVEC(JMAX+K*N) .LT. 0.D0 ) THEN
!         DO 91 J = 1,N
!            EIGVEC(J+K*N) = - EIGVEC(J+K*N)
!91       CONTINUE
!      END IF
!*
!*     ... find configuration JMAX in configuration file
!*
!      print *, 'Finding largest component'
!      READ(iuc,'(A)' ) config
!      READ(iuc,'(20(1x,A3))') (elc(j),j=1,max(1,nclosd))
!      if (config(31:66) .ne. '  ') then
!        lwf = nwf - nclosd
!92      READ(iuc,'(A)' ) elc(1)
!        lwf = lwf -20
!        if (lwf .gt. 0) go to 92
!      end if
!*     .. start reading the configurations
!      do 200 i = 1,jmax
!         READ(iuc,'(8(1X,A3,1X,I2,1X))') (ELC(j),Q(j),j=1,8)
!         READ(iuc,'(15(1X,A3))') (COUPLE(J),J=1,15)
!  200 continue
!      NOCC = 0
!    8 IF (ELC(NOCC+1) .NE. '   ') THEN
!         NOCC = NOCC + 1
!         IF (NOCC .LT. 8) GO TO 8
!      END IF
!      call pack(nocc,elc,q,couple,config)
!      len = 64
!98    if (config(len:len) .eq. ' ') then
!         len = len-1
!         if (len .gt. 1) go to 98
!      end if
!
!*      if (leigen(k+1)) then
!      WRITE (iouv,'(/i6,f16.9,2x,A)')  
!     :       JMAX, EIGVAL, CONFIG(1:LEN)
!      WRITE (iouv,'(7F11.8)') (EIGVEC(J+K*N),J=1,N)
!      WRITE(iscw,'(/1X,F19.12,2X,A/(1X,7F11.8))' ) EIGVAL
!*    :      CONFIG(1:LEN),(EIGVEC(J+K*N),J=1,N)
!      print *, 'Finished Eigenvector ',k+1
!      end if
!      rewind(iuc)
!   40 continue
      print *, 'Leaving LSJMAT'
      END 
