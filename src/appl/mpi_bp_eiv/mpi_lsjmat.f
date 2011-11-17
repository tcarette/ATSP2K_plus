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
     :    onlydvd,jhere,lhere,leigen,pflsj,njv,maxj,ATOMNAME)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NOD=220,NTERMD=31,NWD=60)
*     MPI stuff ***********************************************
*
        INCLUDE 'mpif.h'
        parameter (MAXPROC=100)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)

****************************************************************
      LOGICAL  SKIP
      logical onlydvd,jhere,lhere
      logical leigen(ntermd)
*
      POINTER (qh,h(ncfg,3)),(qjan,jan(ncfg,3))
      COMMON /buffer/qh, qjan, nrow(3), iflag(3)
      POINTER (IQLSP,LSP(ncfg))
      COMMON /NCFGS/ IQLSP,index(ntermd),IJK,nterm,termsh(ntermd)

      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,iout,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk
      COMMON /DSSLSJ/ shift
*
      INTEGER cptr
      POINTER (iqhmx,hmx(30)),(iqjptr,jptr(30)),(iqimx,imx(30)),
     :   (iqeigvec,eigvec(1)),(iqkval,kval(1)),(iqcptr,cptr(1)),
     :        (qtm,tm(1)),(qtp,tp(1)),(qdiag,hii(ncfg)),
     :        (qiwork,iwork(1))
      COMMON /MX/iqhmx, iqjptr, iqimx, iqeigvec, iqkval, iqcptr,
     :           qtm,qtp,qdiag,qiwork
*
      CHARACTER config*66,elc(8)*3,couple(15)*3
      INTEGER   q(8), ipos(3),iselec(NTERMD)
      LOGICAL iupper,hiend
      DATA iupper/.false./
      DATA crite,critc,critr,trhold,ortho/
     :       0.d-17,1.d-8,1.d-8,1.d-10,1.d-15/
!     :       0.d-17,1.d-12,1.d-20,1.d0,1.d0/
      integer req_eigen,ind_jj
      dimension Ssms(nume),g_J(nume),g_JLS(nume)
      pointer(pflsj,flsj(nterm,nterm,2,maxj))
      common/hmx_nze/nze_t,nze_c
      double precision flsj;
      logical lhmx_memory,lhmx_disk
      logical lest

*     ------------------------------------------------------------------
*     on entry
*     --------
*     N   - dimension of the matrix
*     JJ  - J-value (not needed in this restricted application)
*     Nume- number of desired eigenvalues/eigenvectors
*     LS  - type ofm  calculation (LS or LSJ -- not needed here)
*     NZE - size of the non-zero array for the interaction matrix
*     ------------------------------------------------------------------

      ind_jj = 0.5*(maxj - jj + 2); 
      call alcmat(skip,ncfg,nume,nze_c,nze_t,lhmx_memory,lhmx_disk,
     :            jj,njv)
      if (lhmx_memory) nze = nze_t
      if (idisk == 1) nze = nze_c

*       .. clear the interaction matrix
        call dfill(nze,0.d0,hmx,1)
        call dfill(ncfg,0.d0,hii,1)
*
      if (iLS .eq. 0) then
         nij = 0
         istart = 0
         shift = 0.d0
         rewind(iouhn)
         rewind(iouhz)
         rewind(iouhs)
         if (idisk .eq. 1) rewind(iouhm)
         mycol = 0

         do j = myid+1,ncfg,nprocs
           call hmx_lsj(ncfg,j,nze,ind_jj,nij,istart,shift,
     :                  mycol,pflsj,njv)
         end do

!***************************END IM***********************************
*       ..gather all diagonals from processors.
	call mpi_allr_dp(hii,ncfg) !gdsummpi(hii,ncfg,tm)
      end if 

      if (iLS.ne.0) then
*       .. we have an LS case
	 if (idisk .eq. 1 ) then
*      ***** Form the diagonal matrix. search for lowest
            rewind(iouhn)
	    mycol = 0 
            do j = myid+1,ncfg,nprocs
	      mycol = mycol + 1
              ipos(1) = 0
     	      read(iouhn) jb,m1,hii(j)
            end do
*       ..gather all diagonals from processors.
	    call mpi_allr_dp(hii,ncfg) !gdsummpi(hii,ncfg,tm)
            shift = hii(1)
*     
*        .. scale the diagonal elements for accuracy
            do j = 1,ncfg
     	      hii(j) = hii(j) - shift
            end do
	 else
*       ... Form the interaction matrix.. compute diagonals
            nij = 0
	    mycol = 0
            rewind(iouhn)
            do j = myid+1,ncfg,nprocs
	      mycol = mycol + 1
	      read(iouhn) jb,m1,(hmx(nij+i),i=1,m1),(imx(nij+i),i=1,m1)
	      hii(j) = hmx(nij+1)
	      nij = nij + m1
	      jptr(mycol) = nij
            end do
*       ..gather all diagonals from processors.
            call mpi_allr_dp(hii,ncfg) !gdsummpi(hii,ncfg,tm)
	    shift = hii(1)
	    do j = 1,ncfg
              hii(j) = hii(j) - shift
            end do
            hmx(1) = hii(myid+1)
            mycol = 1
	    do j = 1+myid+nprocs,ncfg,nprocs
              mycol = mycol + 1
              hmx(jptr(mycol-1)+1) = hmx(jptr(mycol-1)+1)-shift
            end do
	 end if
      end if
  
      nloc = mycol
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
      lim = min(5*nume,ncfg)
      iworksz = (2*ncfg+lim+nume+10)*lim + nume
      iiwsz = 6*lim + nume
!      print*, 'IWORKSZ = ', iworksz;
      n = ncfg
*     if (idisk .eq. 1) then
*	Nze = jptr(n)
*     else 
*	nze = ncfg
*     end if
***      ilow = 1
***      ihigh = nume

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

      maxiter = max(nume*100,200)
*      mblock = nume
      mblock = ie -1
      niv = 0

      if (onlydvd .and. lhere) then
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
         if (lhere) then
*           .. Read in initial guess else(There is no .l file) Start with
*           .. 0 initial guesses.
            if (myid.eq.0) then
              rewind(iouv)
              read(iouv,*)
              read(iouv,'(//8X,I4,2X,8X,I4)') jj1,number
              niv = number
	      do K=0,number-1
	         read(iouv,*)
	         read(iouv,*)
	         read(iouv,'(7F11.8)') (EIGVEC(J+K*N),J=1,N)
	      enddo
*
*	      ..Brief Modified Gram Scmidt. - - - - - - - - - - - - 
*
	      kcur = 1
	      do k=1,number
                jcur = kcur + N
                call dgemv('T',N,number-k,1.d0,eigvec(jcur),N,
     :                  eigvec(kcur),1, 0.d0,tm,1)
                do j=k+1,number
                   call daxpy(N,-tm(j-k),eigvec(kcur),1,eigvec(jcur),1)
                   jcur = jcur + N
                enddo
                kcur = kcur + N
              enddo
*	      ..Normalize		- - - - - - - - - - - - 
	      kcur = 1
              do k=1,number
                tm(k) = 1./dnrm2(N,eigvec(kcur),1)
                call dscal(n,tm(k),eigvec(kcur),1)
                kcur = kcur + N
 	      enddo
*					- - - - - - - - - - - -
	    endif
	    call MPI_BCAST(niv,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(eigvec,number*N,MPI_DOUBLE_PRECISION,0,
     $			MPI_COMM_WORLD,ierr)
         endif
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      endif

!       write(*,'(A)') 'hii:'
!       write(*,'(5F12.5)') (hii(i),i=1,n)
!       write(*,'(A)') 'hmx:'
!       write(*,'(5F12.5)') (hmx(i),i=1,nze)
!       write(*,'(A)') 'jptr:'
!       write(*,'(5I12)') (jptr(i),i=1,n)
!       write(*,'(A)') 'imx:'
!       write(*,'(5I12)') (imx(i),i=1,nze)
!       write(*,'(A)') 'eigvec:'
!       write(*,'(5F12.9)') (eigvec(i),i=1,nume*N)
!       write(IWRITE,*) iupper,nze,Ncfg,lim,ilow,ihigh,niv,mblock,
!     :		crite,critc,critr,ortho,maxiter
	if (myid .eq. 0) write(iscw,*) 'Starting Davidson'

        if (myid == 0) then
            call brci_est(JJ,eigvec,ncfg,nume,iuc,ATOMNAME,lest)
        end if
        call MPI_BCAST(eigvec,nume*(ncfg+1),MPI_DOUBLE_PRECISION,0,
     $                  MPI_COMM_WORLD,ierr)
        call MPI_BCAST(lest,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        if (lest) niv = nume
        maxiter = 500;
	CALL dvdson(hmx,imx,jptr,iupper,nze,tm,tp,
     :        Ncfg,lim,hii,ilow,ihigh,iselec,niv,mblock,
     :        crite,critc,critr,ortho,maxiter,
     :        eigvec,iworksz,iwork,iiwsz,hiend,nloops,nmv,ierr)
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!       write(0,'(A)') 'eigvec:'
!       write(0,'(5F12.9)') (eigvec(i),i=1,nume*N)

      if (myid == 0) then
         if (ierr .ne. 0) then
            write(0,*) ' Dvdson returned ierr=',ierr
         end if

         ist = n*nume
         req_eigen = 0;
         do icount = 1, nume 
            if (leigen(icount)) req_eigen = req_eigen + 1;
         end do

         write(iscw,*)
         write(iscw,*) 'Summary of Davidson Performance'
         write(iscw,*) '==============================='
         write(iscw,'(A,2I10)') 'Number of Iterations/Matvecs: ', 
     :                         nloops,nmv
         write(iscw,*) 'Shifted Eigval''s:'
         write(iscw,*) (eigvec(ist+i),i=1,nume)
         write(iscw,*) 'Delta Lambda:    '
         write(iscw,*) (eigvec(ist+nume+i),i=1,nume) 
         write(iscw,*) 'Residuals:       '
         write(iscw,*)  (eigvec(ist+2*nume+i),i=1,nume)
         WRITE (iouj, '(//A8,I4,2X,A8,I4)' ) '  2*J = ',JJ,'NUMBER =',
     :           req_eigen; 

         Ssms(1:nume) = 0
!         call prprty_brci(jj,nume,leigen,eigvec,Ssms)
         call brci_zeeman(jj,ncfg,nume,eigvec,g_J,g_JLS,iuc)
         rewind (iuc)
         call brci_eig(ncfg,jj,nume,leigen,eigvec,Ssms,g_J,g_JLS,
     :        eigvec(ncfg*nume+1),EC,shift,iouj,iuc)
      end if

      if (myid == 0) write (iscw,*) ':: Finished Block 2J = ',jj
      call dalloc(iqeigvec,iworksz)
      call dalloc(qtm,ncfg)
      call dalloc(qtp,ncfg)
      call dalloc(qdiag,ncfg)
      call dalloc(qiwork,iiwsz)
      call dalloc(iqjptr,ncfg)
      call dalloc(iqhmx,n_hmx)
      call dalloc(iqimx,n_hmx)
      write (iscw,*)

      end subroutine lsjmat 
