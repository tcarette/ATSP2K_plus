*------------------------------------------------------------------------
*        G E N I N T B R
*------------------------------------------------------------------------
*
*       Generate the list of integrals that could arise from a set
*       of orbitals. 

        SUBROUTINE genintbr(nclosd,maxorb,lmax,ql,qintptr,qpackn,qvalue,
     :                                               rel,skip,eval)
        IMPLICIT real*8 (a-h,o-z)
        INCLUDE 'mpif.h'
        parameter (MAXPROC=9999)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
!        POINTER (qtmp,tmp(1))
        POINTER (ql,l(1)),(qintptr,intptr(0:2*lmax+1,7)),
     :          (qpackn,ipackn(21)),(qvalue,value(21))
        POINTER(qtmp_value,tmp_value(1)),(qtmp_ipackn,itmp_ipackn(1))
        INTEGER l,intptr,ipackn,itmp_ipackn
        INTEGER maxorb,lmax,n
        LOGICAL omit, ltriang, gen, rel, skip, eval
	CHARACTER*3 el

C >>>>  The hlc routine requires an array of electron labels 
*	for diagonstic purposes.  The A3 format of MCHF is not 
*	compatible with the Belfast Integer (3H   ) format.
*	No problem will be encountered as long as an error is
*	not detected by hlc.
*       .. sweep through to find dimensions
	gen = .false.
*
*       Generate the list of possible integrals
*
*       Make the FK integrals: i2 <= i4
*
        ntmp = nclosd;  
10      n = 0
        do k = 0,2*lmax+1
          do i2 = 1,maxorb
            l2 = l(i2+ntmp)
            do i4 = i2,maxorb
              l4 = l(i4+ntmp)
              if (ltriang(k,l2,l2) .and. ltriang(k,l4,l4)) then
                n = n+1
		if (gen.and. (mod(n,nprocs) == myid) ) then
                  ipackn(n) = (k*64 + i2)*64 + i4
		  j2 = i2+nclosd
		  j4 = i4+nclosd
CGG		  value(n) = fk(j2,j4,k,rel)
		  if( eval) value(n) = rk(j2,j4,j2,j4,k,rel)
*                 if (abs(value(n)).lt.(1.e-39)) stop
!       write (*,'(A5,I5,I15,3I5,f15.12)') 'Fk',n,ipackn(n),i2,
!     :            i4,k, value(n)
		end if
              end if
            end do
          end do
          if (gen) then
	    intptr(k,1) = n
C           write (*,*) k,intptr(k,1)
	  end if
        end do
*
*       Make the GK integrals: i2 < i4
*
        do k = 0,2*lmax+1
          do i2 = 1,maxorb
            l2 = l(i2+ntmp)
            do i4 = i2+1,maxorb
              l4 = l(i4+ntmp)
              if (ltriang(k,l2,l4) ) then
                n = n+1
		if (gen.and. (mod(n,nprocs) == myid)) then
                  ipackn(n) = (k*64 + i2)*64 + i4
		  j2 = i2+nclosd
		  j4 = i4+nclosd
C		  value(n) = gk(j2,j4,k,rel)
		  if (eval) value(n) = rk(j2,j4,j4,j2,k,rel)
*                 if (abs(value(n)).lt.(1.e-29)) stop
!          write (*,'(A5,I5,I15,3I5,f15.12)') 'Gk',n,ipackn(n),i2,
!     :                        i4,k, value(n)
		end if
              end if
            end do
          end do
          if (gen) then
	    intptr(k,2) = n
C           write (*,*) k,intptr(k,2)
	  end if
        end do
*
*       Make the RK integrals: i1<=i2; i1 <= i3; i1<=i4
*
        do k = 0,2*lmax+1
          do i1 = 1,maxorb
            l1 = l(i1+ntmp)
            do i2 = i1,maxorb
              l2 = l(i2+ntmp)
              do i3=i1,maxorb
                l3 = l(i3+ntmp)
                if (ltriang(k,l1,l3)) then
                  do i4 = i1,maxorb
                    l4 = l(i4+ntmp)
*                   .. omit if FK
                    omit = (i1.eq.i4.and.i2.eq.i3).or.
     :                     (i1.eq.i3.and.i2.eq.i4)
                    if ( .not. omit .and. ltriang(k,l2,l4) ) then
                      n = n+1
		      if (gen.and. (mod(n,nprocs) == myid)) then
                        ipackn(n) = (((k*64+i1)*64+i2)*64+i3)*64+i4
			j1 = i1+nclosd
			j2 = i2+nclosd
			j3 = i3+nclosd
			j4 = i4+nclosd
			if (eval) value(n) = rk(j1,j2,j3,j4,k,rel)
!          write (*,'(A5,i5,i15,5i5,f15.12)') 'Rk',n,ipackn(n),i1,i2,
!     :                        i3,i4,k, value(n)
                      end if
                    end if
                  end do
                end if
              end do
            end do
          end do
          if (gen) then
	    intptr(k,3) = n
C           write (*,*) k,intptr(k,3)
	  end if
        end do
*
*       Make the L integrals: i2 < i4
*
        do i2 = 1,maxorb
          l2 = l(i2+ntmp)
          do i4 = i2,maxorb
            l4 = l(i4+ntmp)
            if (l2 .eq. l4) then
              n = n+1
	      if (gen.and. (mod(n,nprocs) == myid)) then
                ipackn(n) = i2*64 + i4
                j1 = i2 + nclosd
                j2 = i4 + nclosd
		if (eval) value(n) = hlc(el,j1,j2,rel)
!          write (*,'(A5,i5,i15,2i8,f15.12)') 'HL',n,ipackn(n),i2,
!     :                        i4, value(n)
	      end if
            end if
          end do
        end do
        if (gen) then
	  do k = 0,2*lmax+1
	    intptr(k,4) = n
	  end do
C         write (*,*) 0,intptr(0,4)
	end if

*       Omit Breit-Pauli integrals if skip=true
	if (skip) go to 20
*
*       Make the Z integrals: i2 < i4
*
        do i2 = 1,maxorb
          l2 = l(i2+ntmp)
          do i4 = i2,maxorb
            l4 = l(i4+ntmp)
            if (l2 .eq. l4 .and. l2 .gt. 0) then
              n = n+1
	      if (gen.and. (mod(n,nprocs) == myid)) then
                ipackn(n) = i2*64 + i4
		j2 = i2+nclosd
		j4 = i4+nclosd
		if (eval) value(n) = zeta(j2,j4)
!        write (*,'(A5,i5,i15,2i8,f15.12)') 'Zeta',n,ipackn(n),i2,
!     :                        i4, value(n)
	      end if
            end if
          end do
        end do
        if (gen) then
	  do k = 0,2*lmax+1
	    intptr(k,5) = n
	  end do
C         write (*,*) 0,intptr(0,5)
	end if
*
*       Make the NK integrals: i1<=i3; i2 <= i4
*       Note that for these integrals, the packing
*       routine packs K+1, where K=-1,0,1,..  So
*       here k = K+1
*
        do k = 0,2*lmax+1
 	  if (k .eq. 0) then
 	    kk = k+1
 	  else
 	    kk = k-1
 	  end if
	  km = k-1
          do i1 = 1,maxorb
            l1 = l(i1+ntmp)
            do i2 = 1,maxorb
              l2 = l(i2+ntmp)
              do i3=i1,maxorb
                l3 = l(i3+ntmp)
                if (ltriang(kk,l1,l3).or. ltriang(kk+2,l1,l3) ) then
*               if ( ltriang(kk+2,l1,l3) ) then
                  do i4 = i2,maxorb
                    l4 = l(i4+ntmp)
                    if (ltriang(kk,l2,l4) .or.ltriang(kk+2,l2,l4)) then
                      n = n+1
		      if (gen.and. (mod(n,nprocs) == myid)) then
                        ipackn(n) = (((k*64+i1)*64+i2)*64+i3)*64+i4
			j1 = i1+nclosd
			j2 = i2+nclosd
			j3 = i3+nclosd
			j4 = i4+nclosd
			if (eval) value(n) = sn(j1,j2,j3,j4,km)
!         write (*,'(A5,i5,i15,5i5,f15.12)') 'Nk',n,ipackn(n),i1,i2,
!     :                        i3,i4,km,value(n)
                      end if
                    end if
                  end do
                end if
              end do
            end do
          end do
          if (gen) then
	    intptr(k,6) = n
C           write (*,*) k,intptr(k,6)
	  end if
        end do
*
*       Make the VK integrals: i2 < i4
*
        do k = 0,2*lmax+1
	  kp = k+1 
          do i1 = 1,maxorb
            l1 = l(i1+ntmp)
            do i2 = 1,maxorb
              l2 = l(i2+ntmp)
              do i3=1,maxorb
                l3 = l(i3+ntmp)
                if (ltriang(kp,l1,l3)) then
                  do i4 = i2,maxorb
                    l4 = l(i4+ntmp)
                    if ( ltriang(kp,l2,l4) ) then
                      n = n+1
		      if (gen.and. (mod(n,nprocs) == myid)) then
                        ipackn(n) = (((k*64+i1)*64+i2)*64+i3)*64+i4
			j1 = i1+nclosd
			j2 = i2+nclosd
			j3 = i3+nclosd
			j4 = i4+nclosd
			if (eval) value(n) = vk(j1,j2,j3,j4,k)
!          write (*,'(A5,i5,i15,5i5,f15.12)') 'Vk',n,ipackn(n),i1,i2,
!     :                        i3,i4,km,value(n)
                      end if
                    end if
                  end do
                end if
              end do
            end do
          end do
          if (gen) then
	    intptr(k,7) = n
C           write (*,*) k,intptr(k,7)
	  end if
        end do
*
*       Note: SK integrals are the same as Nk

20      continue
*       .. allocate memory for integral book keeping
        if (.not. gen) then
          call alloc(qintptr,((2*lmax+2)*7),4)      
          if(myid ==0) write(*,*) 'Allocating ',n,' integrals'
          call alloc(qpackn,n,4)
          ipackn(1:n) = 0
          if (eval) call alloc(qvalue,n,8)
          if (eval) value(1:n) = 0
          gen = .true.
	  go to 10
	end if

       if (eval) then 
          call alloc(qtmp_value,n,8)  
          tmp_value(1:n) = 0.0
         call MPI_ALLREDUCE(value,tmp_value,n,
     :           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
         value(1:n) = tmp_value(1:n)
         if(qtmp_value.ne.0) call dalloc(qtmp_value,n)
       end if

      call alloc(qtmp_ipackn,n,4)
      itmp_ipackn(1:n) = 0 
      call MPI_ALLREDUCE(ipackn,itmp_ipackn,n,
     :          MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror)
      ipackn(1:n) = itmp_ipackn(1:n)
 
      if (qtmp_ipackn.ne.0) call dalloc(qtmp_ipackn,n)

*        Write(99,*) 'ic, (intptr(ll,ic),ll=0,2*lmax+1)'
! 	do ic = 1,7
! 	   write(99,*) ic, (intptr(ll,ic),ll=0,2*lmax+1)
! 	end do
!        write(99,*) 'n',n
!        write(99,*) ipackn(1:n)

*
        end
