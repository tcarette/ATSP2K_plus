*------------------------------------------------------------------------
*        G E N I N T B R
*------------------------------------------------------------------------
*
*       Generate the list of integrals that could arise from a set
*       of orbitals. 

        SUBROUTINE genintbr(nclosd,maxorb,lmax,ql,qintptr,qpackn,qvalue,
     :                                               rel,skip)
        IMPLICIT real*8 (a-h,o-z)

        POINTER (ql,l(1)),(qintptr,intptr(0:2*lmax+1,7)),
     :          (qpackn,ipackn(1)),(qvalue,value(1))
        INTEGER l,intptr,ipackn
        INTEGER maxorb,lmax,n
        LOGICAL omit, ltriang, gen, rel, skip
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
10      n = 0
        do k = 0,2*lmax+1
          do i2 = 1,maxorb
            l2 = l(i2)
            do i4 = i2,maxorb
              l4 = l(i4)
              if (ltriang(k,l2,l2) .and. ltriang(k,l4,l4)) then
                n = n+1
		if (gen) then
                  ipackn(n) = (k*64 + i2)*64 + i4
		  j2 = i2+nclosd
		  j4 = i4+nclosd
CGG		  value(n) = fk(j2,j4,k,rel)
		  value(n) = rk(j2,j4,j2,j4,k,rel)
                  if (abs(value(n)).lt.(1.e-39)) stop
*                 write (*,*) 'Fk',n,ipackn(n),i2,i4,k, value(n)
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
            l2 = l(i2)
            do i4 = i2+1,maxorb
              l4 = l(i4)
              if (ltriang(k,l2,l4) ) then
                n = n+1
		if (gen) then
                  ipackn(n) = (k*64 + i2)*64 + i4
		  j2 = i2+nclosd
		  j4 = i4+nclosd
C		  value(n) = gk(j2,j4,k,rel)
		  value(n) = rk(j2,j4,j4,j2,k,rel)
                  if (abs(value(n)).lt.(1.e-29)) stop
*                 write (*,*) 'Gk',n,ipackn(n),i2,i4,k, value(n)
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
            l1 = l(i1)
            do i2 = i1,maxorb
              l2 = l(i2)
              do i3=i1,maxorb
                l3 = l(i3)
                if (ltriang(k,l1,l3)) then
                  do i4 = i1,maxorb
                    l4 = l(i4)
*                   .. omit if FK
                    omit = (i1.eq.i4.and.i2.eq.i3).or.
     :                     (i1.eq.i3.and.i2.eq.i4)
                    if ( .not. omit .and. ltriang(k,l2,l4) ) then
                      n = n+1
		      if (gen) then
                        ipackn(n) = (((k*64+i1)*64+i2)*64+i3)*64+i4
			j1 = i1+nclosd
			j2 = i2+nclosd
			j3 = i3+nclosd
			j4 = i4+nclosd
			value(n) = rk(j1,j2,j3,j4,k,rel)
*                       write (*,*) 'Rk',n,ipackn(n),i1,i2,i3,i4,k,
*    :                              value(n)
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
          l2 = l(i2)
          do i4 = i2,maxorb
            l4 = l(i4)
            if (l2 .eq. l4) then
              n = n+1
	      if (gen) then
                ipackn(n) = i2*64 + i4
                j1 = i2 + nclosd
                j2 = i4 + nclosd
		value(n) = hlc(el,j1,j2,rel)
*               write (*,*) 'HL',n,ipackn(n),i2,i4, value(n)
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
          l2 = l(i2)
          do i4 = i2,maxorb
            l4 = l(i4)
            if (l2 .eq. l4 .and. l2 .gt. 0) then
              n = n+1
	      if (gen) then
                ipackn(n) = i2*64 + i4
		j2 = i2+nclosd
		j4 = i4+nclosd
		value(n) = zeta(j2,j4)
*               write (*,*) 'Zeta',n,ipackn(n),i2,i4, value(n)
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
            l1 = l(i1)
            do i2 = 1,maxorb
              l2 = l(i2)
              do i3=i1,maxorb
                l3 = l(i3)
                if (ltriang(kk,l1,l3).or. ltriang(kk+2,l1,l3) ) then
*               if ( ltriang(kk+2,l1,l3) ) then
                  do i4 = i2,maxorb
                    l4 = l(i4)
                    if (ltriang(kk,l2,l4) .or.ltriang(kk+2,l2,l4)) then
                      n = n+1
		      if (gen) then
                        ipackn(n) = (((k*64+i1)*64+i2)*64+i3)*64+i4
			j1 = i1+nclosd
			j2 = i2+nclosd
			j3 = i3+nclosd
			j4 = i4+nclosd
			value(n) = sn(j1,j2,j3,j4,km)
*                       write (*,*) 'Nk',n,ipackn(n),i1,i2,i3,i4,km,
*    :                     value(n)
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
            l1 = l(i1)
            do i2 = 1,maxorb
              l2 = l(i2)
              do i3=1,maxorb
                l3 = l(i3)
                if (ltriang(kp,l1,l3)) then
                  do i4 = i2,maxorb
                    l4 = l(i4)
                    if ( ltriang(kp,l2,l4) ) then
                      n = n+1
		      if (gen) then
                        ipackn(n) = (((k*64+i1)*64+i2)*64+i3)*64+i4
			j1 = i1+nclosd
			j2 = i2+nclosd
			j3 = i3+nclosd
			j4 = i4+nclosd
			value(n) = vk(j1,j2,j3,j4,k)
*                       write (*,*) 'Vk',n,ipackn(n),i1,i2,i3,i4,k,
*    :                        value(n)
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
          write(*,*) 'Allocating space for ',n,' integrals'
          call alloc(qpackn,n,4)
          call alloc(qvalue,n,8)
          gen = .true.
	  go to 10
	end if
*	do ic = 1,7
*	   print *, ic, (intptr(ll,ic),ll=0,2*lmax+1)
*	end do
*
        end
