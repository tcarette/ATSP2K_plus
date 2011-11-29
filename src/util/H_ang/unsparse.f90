SUBROUTINE unsparse(print,save)
  IMPLICIT NONE

  logical, intent(in) :: print,save

  character*32, parameter :: fmt='"<",i3,"|H|",i3,"> ="'
  character*32, parameter :: fmt0='"<",i3,"|H|",i3,"> = 0"'
  ! block properties
  integer ncn,ncfg
  ! global iterators
  integer nz,ncoef,nb
  ! in-block iterators
  integer ja,jb, nijcurr,mei
  integer i,n
  ! integral caracteristics
  integer ipack,icase,iptr,k,l1,l2,r1,r2,ni
  !tmp
  integer j
  character*3 str
  character*50 str2

  type(Hang) cur

!  if(.not.save.and..not.print)then
!    write(0,*) "Nothing to do in unparse"
!    return
!  endif

  if(save) allocate(tot_hang(nblock))

  ncoef = 0
  nz=0
  do nb=1,nblock

    ncfg=ncfg_bl(nb)

    if(save) allocate(tot_hang(nb)%h(ncfg,ncfg))

    ncn=ncn_bl(nb)
    ja = 1
    jb = 1
    nijcurr = nz+1

    ! prepare the first matrix element

    mei=0
    ni=ico(nijcurr)
    if(save)then
      cur%ni=ni
      allocate(cur%int(ni))
      allocate(tot_hang(nb)%h(1,1)%int(ni))
    endif
    if(print)then
      write(*,*)
      write(*,'('//fmt//')') 1,1
    endif

    do i = 1, ncn;

      n = ncoef+i

      !determines the next non-zero element
      !and do what is correspondingly necessary

      if (i.gt.ico(nijcurr))then

        nijcurr = nijcurr + 1

        if(save) tot_hang(nb)%h(ja,jb)=cur

        ! determines row and column of the current matrix elements

        if(ja.lt.jan(nijcurr)) then
          do j=ja+1,jan(nijcurr)-1
            if(save)tot_hang(nb)%h(j,jb)%ni=0
            if(print)then
              write(*,*)
              write(*,'('//fmt0//')') jb,j
            endif
          enddo
        else
          jb=jb+1

          do j=ja+1,ncfg
            if(save) tot_hang(nb)%h(j,jb)%ni=0
            if(print)then
              write(*,*)
              write(*,'('//fmt0//')') jb,j
            endif
          enddo
          do j=jb,jan(nijcurr)-1
            if(save) tot_hang(nb)%h(j,jb)%ni=0
            if(print)then
              write(*,*)
              write(*,'('//fmt0//')') jb,j
            endif
          enddo

        endif
        ja= jan(nijcurr)

        ! prepare the matrix element
        mei=0
        ni=ico(nijcurr)-ico(nijcurr-1)
        if(save)then
          deallocate(cur%int)
          cur%ni=ni
          allocate(cur%int(ni))
          allocate(tot_hang(nb)%h(ja,jb)%int(ni))
        endif
        if(print)then
          write(*,*)
          write(*,'('//fmt//')') jb,ja
        endif
      endif

      !determines integral caracteristics

      mei=mei+1
      iptr  = inptr(n)
      ipack = ipackn(iptr)

        !icase

      if(iptr.le.intptr(1))then
        icase=1
      else if(iptr.le.intptr(2))then
        icase=2
      else if(iptr.le.intptr(3))then
        icase=3
      else
        icase=4
      endif

        !involved electron

      if (icase.eq.3) then
        r2 = mod(ipack,64)
        ipack = ipack/64
        r1 = mod(ipack,64)
        ipack = ipack/64
        l2 = mod(ipack,64)
        ipack = ipack/64
        l1 = mod(ipack,64)
        k = ipack/64
      else
        l2 = mod(ipack,64)
        ipack = ipack/64
        l1 = mod(ipack,64)
        k = ipack/64
      end if

      if(save)then
        cur%int(mei)%c  = cn(n)
        cur%int(mei)%inptr = iptr
        cur%int(mei)%icase = icase
        cur%int(mei)%l1 = k
        cur%int(mei)%l1 = l1
        cur%int(mei)%l2 = l2
        cur%int(mei)%r1 = r1
        cur%int(mei)%r2 = r2
      endif
      if(print)then

        str='no'
        if(mod(mei,1)==0.or.mei==ni) str='yes'

        if (icase==3) then
          write(str2,'('//intfmt(3)//')') &
                   k,el(l1),el(l2),el(r1),el(r2)
        else if (icase==4) then
          write(str2,'('//intfmt(4)//')') el(l1),el(l2)
        else
          write(str2,'('//intfmt(icase)//')') k,el(l1),el(l2)
        endif

        write(*,'(" + ",f14.10,"  ")',advance='no') cn(n)
        write(*,'(a20)',advance=str) str2
      endif

    end do

    ncoef = ncoef + ncn
    nz=nz+nze_bl(nblock)
  enddo
END SUBROUTINE unsparse
