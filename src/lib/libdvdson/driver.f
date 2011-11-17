        progRAM DRIVER

        implicit double precision (A-H,O-z)
*        PARAMETER(Nmax=2149,NZERmax=335416) 
        PARAMETER(Nmax=9000,NZERmax=600000) 
	PARAMETER(NUMEmax=1)
	PARAMETER(LIMmax=NUMEmax+18)

        common/MATRIX/ dMy_Matr(NZERmax),Ind_Col(Nmax),
     :		        Ind_Row(NZERmax),iupper
	common/TEMP/tm(Nmax),tp(Nmax)
	dimension DIAG(Nmax)
        dimension WORK(LIMmax*(2*Nmax+LIMmax+(NUMEmax+10))+NUMEmax)
	dimension IWORK(6*LIMmax+NUMEmax)
	dimension idim(10),MINELEM(NUMEmax),ISELEC(LIMmax)
	real*4 tarray(2),t1,etime

	logical iguess,iupper,hiend
	common /hyperstuff/ iam, nodes
        iam = 0
        nodes=1
	iupper=.true.

******** Reading in the upper matrix.
	IN_UNIT=10
        open(IN_UNIT,file='../RA/Hypercube/SUMMER92/MATRICES/mat.748',
     :		status='unknown',form='formatted')

	read(IN_UNIT,*) N,NZER
** READ THE COLUMN INDICES **
        read(IN_UNIT,*) (WORK(i), i=1,N)
        print*,'read col'

*** READ ROW INDICES AND CREATE IND_COL **
        mycol=0
        iup=0
        localup=0
        do 10 icol=1,N
        idown=iup+1
        iup=WORK(icol)
           if (iam.eq.mod(icol-1,nodes) ) then
              mycol=mycol+1
              localdown=localup + 1
              localup =localdown + iup - idown
              read(IN_UNIT,*) (Ind_Row(irow), irow=localdown,localup)
              Ind_Col(mycol)=localup
           else
              read(IN_UNIT,*) (junk, irow=idown,iup)
           endif
 10     continue
        print*,'read row'

*** READ VALUES **
        mycol=0
        iup=0
        localup=0
        do 20 icol=1,N
        idown=iup+1
        iup=WORK(icol)
           if (iam.eq.mod(icol-1,nodes) ) then
              mycol=mycol+1
              localdown=localup + 1
              localup=Ind_Col(mycol)
              read(IN_UNIT,*) (dMy_Matr(i), i=localdown,localup)
           else
              read(IN_UNIT,*) (rubs, irow=idown,iup)
        endif
 20     continue
        print*,'read val'

        do 133 i=1,N
	   if (iupper) then
    	      DIAG(i)=dmy_matr( Ind_Col(i) )
	   else 
	      DIAG(i)=dmy_matr( Ind_Col(i-1)+1 )
	   endif
	   diamax=max(diamax,ABS(DIAG(i)))
 133 	continue

	print*, nzer
	goto 248
	print*, 'Do you want scaling and shifting ? (1/0)'
	read(*,*) iyes
	if (iyes.ne.1) goto 248
		print*,diamax
*	goto 248
*
* SHIFT
*
	dshift=diag(1)
	DIAG(1)=DIAG(1)-dshift
	dmy_matr(1)=dmy_matr(1)-dshift
	do 243 i=2,N
	   DIAG(i)=DIAG(i)-dshift
	   if (iupper) then
	       dmy_matr(Ind_Col(i))=dmy_matr(Ind_Col(i))-dshift
	   else
 	   dmy_matr(Ind_Col(i-1)+1)=dmy_matr(Ind_Col(i-1)+1)-dshift
	   endif
 243	continue

*
* SCALE
*
	diamax=100000000
	call dscal(N,1/diamax,DIAG,1)
	call dscal(NZER,1/diamax,dmy_matr,1)
*	do 246 i=1,NZER
* 246	   dmy_matr(i)=dmy_matr(i)/diamax
*	print*,diamax
 248	NZER=localup
	ND=N

*************************************************************************
*		TEsting 
************************************************************************
         IRWSZ=LIMmax*(2*Nmax+LIMmax+(NUMEmax+10))+NUMEmax
	 IIWSZ=6*LIMmax+NUMEmax

	N=ND
	do 143 LIM=LIMmax,LIMmax,10
	   HIEND=.false.
	   ILOW=1
	   IHIGH=NUMEmax
	   ISELEC(1)=3
	   ISELEC(2)=5
	   ISELEC(3)=-1
           NIV = 0
	   NUME=NUMEmax
	if (N.lt.ND) then
*           CRITC=1D-5
*           CRITR=1D-5
*	   TRHOLD=1D-2
*	   ORTHO=1D+6
*	   MAXITER=1000
	else
	   CRITE=1.0D-16
	   CRITC=1.D-8
	   MBLOCK=1
           CRITR=1.D-8
	   ORTHO=1D-15
	   MAXITER=150
	endif
	if (NIV.NE.0) then 
*	   call dinit(N,-1.D0,tm,1)
*              do 112 i=1,NUME
*imin= the first not gotten elem( NUME<=N )
*                 do 15 j=1,N
*  15                if (tm(j).lt.0) goto 16
*  16             imin=j
*                 do 222 j=2,N
*  222                if ((tm(j).lt.0).and.
*     :                  (DIAG(j).lt.DIAG(imin))) imin=j
* 
*                 MINELEM(i)=imin
*                 tm(imin)=1.D0
*  112          continue
*	NIV=NUMEmax
*	   call dinit(NIV*N,0.d0,WORK,1)
*	   WORK(2)=0.5
*	   work(11)=0.5 
*	   work(12)=0.5 
*	   work(13)=0.5 
*	   work(4*N+5)=1 
*	   work(5*N+6)=1
*	   work(6*N+7)=1
*	   work(7*N+MINELEM(8))=1
*	   work(8*N+MINELEM(9))=1
*	   work(9*N+MINELEM(10))=1

	endif
*	t1=dclock()
	t1=etime(tarray)
	CALL dvDSon(dMy_Matr,Ind_Row,Ind_Col,iupper,nzermax,tm,tp,
     :             N,LIM,DIAG,
     :		   ILOW,IHIGH,ISELEC,NIV,MBLOCK,
     :             CRITE,CRITC,CRITR,ORTHO,MAXITER,
     :             WORK,IRWSZ,IWORK,IIWSZ,HIEND,nloops,nmv,IERR)
*        t1=(dclock()-t1)
	t1=etime(tarray)-t1

	if (iyes.eq.1) then
	do 135 i=1,NUME
 135	   WORK(NUME*N+i)=diamax*WORK(NUME*N+i)+dshift
	endif

        print*,'Time',t1
	   if (HIEND) NUME=N-ILOW+1
*eigenvalues, differences, residuals, eigenvectors.
	        WRITE(6,4000) IERR,NLOOPS,NMV
        WRITE(6,2000) ((WORK(i), i=NUME*N+j,(N+3)*NUME,NUME),j=1,NUME)
        WRITE(6,3000) ((WORK(I), I=J,J+5), J=1,NUME*N,N)
 4000   FORMAT(/'IERR =',I6,'   Matrix accesses=',I4,
     :          '   Matrix-Vector products=',I4)
 2000   FORMAT(//9X,'Eigenvalues',8X,'Eigval Differences',6X,'Residuals',
     :    //(D25.15,2D20.10))
 3000   FORMAT(//' First six componenents of the eigenvectors'//
     :    (6D25.15))
*     :    (6D12.3))

 143 	continue
	close(IN_UNIT)
        stop 
	end
