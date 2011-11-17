      subroutine hrestor(ncfg)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      POINTER (qjptr,jptr(1))
      COMMON/COLS/qjptr
      POINTER (qhmx,hmx(1)),(qdiag,hii(1))
      common/spd/qhmx,qtm,qtp,qdiag,qiwork
      POINTER (pih,ih(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PVALUE,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
      COMMON/SORT/PHMATT
      POINTER (phmat,hmat(1))
      POINTER (phmatt,hmatt(1))
      POINTER (phmatc,hmatc(ncfg,ncfg))
      logical flagpr

      
      flagpr = .false.
      if (flagpr) then
      print *, ' coucou from beginning of hrestor'
      print *, 'ncfg:', ncfg
      print *, 'jptr:', jptr(ncfg)
      print *, 'hmx:'
      print '(3F20.16)', (hmx(i),i=1,jptr(ncfg))
      print *, 'ih:'
      print '(4I6)', (ih(i),i=1,jptr(ncfg))
      print *, 'jptr:', (jptr(i),i=1,ncfg)
      print *, 'diag :'
      print '(3F20.16)',(hii(i),i=1,ncfg)
      print *, ' coucou from end of hrestor'
      end if

      ncfg2 = ncfg*ncfg
      nbrel = ncfg*(ncfg+1)/2
      call alloc(phmat,nbrel,8)
      call alloc(phmatc,ncfg2,8)
      call alloc(phmatt,nbrel,8)

      k = 0
      do 1 j = 1,ncfg
        do 2 i = j,ncfg
          k = k+1
          if (j .eq. i) hmat(k) = hii(j)
    2   continue
    1 continue


      do 3 j = 1,ncfg
        do 4 i = 1,ncfg
          if (j .eq. i) hmatc(i,j) = hii(j)
    4   continue
    3 continue

      if(flagpr) then
      print *, 'HMAT:'
      print '(3F20.16)',(hmat(i),i=1,nbrel)
      print *, 'diagonale HMATC:'
      print '(3F20.16)',(hmatc(i,i),i=1,ncfg)
      end if

      k = 0
      nzero = 0
      ncfgm1 = ncfg-1

      do 5 j = 0,ncfgm1
        ncfgmj = ncfg - j
        do 6 i = 1,ncfgmj
          ipj = i + j
          ipk = i + k
          inz = ipk - nzero
          if (ih(inz) .eq. ipj) then
             hmatt(k+i) = hmx(inz)
          else
             hmatt(k+i) = 0.0
             nzero = nzero + 1
          end if
    6   continue
        if(flagpr) then
           print *, ' k = ',k
           print '(3F20.16)',(hmatt(i),i=k+1,k+ncfgmj)
        end if
        k = k + ncfgmj
    5 continue


      return
      end
