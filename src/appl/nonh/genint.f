*
*------------------------------------------------------------------------
*        G E N I N T
*------------------------------------------------------------------------
*
        SUBROUTINE genint(maxorb,lmax,ql,qintptr,qpackn,qlused,noint)
        POINTER (ql,l(1)),(qintptr,intptr(0:2*lmax,4)),
     :          (qpackn,ipackn(1)),(qlused,lused(1))
        INTEGER l,intptr,ipackn
        INTEGER maxorb,lmax,n
        LOGICAL omit, ltriang, lused
        DIMENSION noint(4)

*       .. sweep through to find dimensions
*
*       Generate the list of possible integrals
*
*       Make the FK integrals: i2 <= i4
*
        n = 0
        do k = 0,2*lmax
          do i2 = 1,maxorb
            l2 = l(i2)
            do i4 = i2,maxorb
              l4 = l(i4)
              if (ltriang(k,l2,l2) .and. ltriang(k,l4,l4)) then
                n = n+1
              end if
            end do
          end do
        end do
        noint(1)=n
*
*       Make the GK integrals: i2 < i4
*
        do k = 0,2*lmax
          do i2 = 1,maxorb
            l2 = l(i2)
            do i4 = i2+1,maxorb
              l4 = l(i4)
              if (ltriang(k,l2,l4) ) then
                n = n+1
              end if
            end do
          end do
        end do
        noint(2)=n
*
*       Make the RK integrals: i1<=i2; i1 <= i3; i1<=i4
*
        do k = 0,2*lmax
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
                    end if
                  end do
                end if
              end do
            end do
          end do
        end do
        noint(3)=n
*
*       Make the  integrals: i2 < i4
*
        do i2 = 1,maxorb
          l2 = l(i2)
          do i4 = i2,maxorb
            l4 = l(i4)
            if (l2 .eq. l4) then
              n = n+1
            end if
          end do
        end do
        noint(4)=n
*

*       .. allocate memory for integral book keeping

        call alloc(qintptr,((2*lmax+1)*4),4)      
        write(*,*) 'Allocating space for ',n,' integrals'
        call alloc(qpackn,n,4)
        call alloc(qlused,n,4)

*       .. now generate the pointer data
*
*       Generate the list of possible integrals in packed form
*       along with pointer values
*
*       Make the FK integrals: i2 <= i4
*
        n = 0
        do k = 0,2*lmax
          do i2 = 1,maxorb
            l2 = l(i2)
            do i4 = i2,maxorb
              l4 = l(i4)
              if (ltriang(k,l2,l2) .and. ltriang(k,l4,l4)) then
                n = n+1
                ipackn(n) = (k*64 + i2)*64 + i4
              end if
            end do
          end do
          intptr(k,1) = n
        end do
*
*       Make the GK integrals: i2 < i4
*
        do k = 0,2*lmax
          do i2 = 1,maxorb
            l2 = l(i2)
            do i4 = i2+1,maxorb
              l4 = l(i4)
              if (ltriang(k,l2,l4) ) then
                n = n+1
                ipackn(n) = (k*64 + i2)*64 + i4
              end if
            end do
          end do
          intptr(k,2) = n
        end do
*
*       Make the RK integrals: i1<=i2; i1 <= i3; i1<=i4
*
        do k = 0,2*lmax
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
                      ipackn(n) = (((k*64+i1)*64+i2)*64+i3)*64+i4
                    end if
                  end do
                end if
              end do
            end do
          end do
          intptr(k,3) = n
        end do
*
*       Make the  integrals: i2 < i4
*
        do i2 = 1,maxorb
          l2 = l(i2)
          do i4 = i2,maxorb
            l4 = l(i4)
            if (l2 .eq. l4) then
              n = n+1
              ipackn(n) = i2*64 + i4
            end if
          end do
        end do
        intptr(0,4) = n
*

        end
