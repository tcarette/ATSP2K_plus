C23456789012345678901234567890123456789012345678901234567890123456789012

      subroutine recalcen(ncfg,wt,ec,shift)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      POINTER (phmat,hmat(1))
      COMMON/SORT/PHMAT
      double precision wt(*)

      print*,' coucou: entry in recalcen'
      print*,' wt:',(wt(i),i=1,ncfg)
      ncfg2 = ncfg*(ncfg+1)/2
      print*,' hmat:',(hmat(i),i=1,ncfg2)

      en = 0.0
      k = 0

      print*,' ec    = ', ec
      print*,' shift = ', shift 

      do 1 j = 1,ncfg
        ncfgmjp1 = ncfg - j + 1
      print*, 'j = ',j,' ncfgmjp1 = ',ncfgmjp1
        kp1 = k + 1
        en = en + wt(j)**2 * hmat(kp1)
        do 2 i = 2,ncfgmjp1
          kpi = k + i
          ipjm1 = i + j - 1
          en = en + 2. * wt(ipjm1) * wt(j) * hmat(kpi)
    2   continue
        k = k + ncfgmjp1
    1 continue
      en = en + ec + shift

      print *,' Etot Thomas = ', en 

      en = 0
      do 3 j = 1, ncfg
         do 4 i = j,ncfg
         indx = (j-1)*(2*ncfg-j)/2 + i
           if (i .eq. j ) then
             en = en + wt(i)**2 * hmat(indx) 
           else
             en = en + 2. * wt(i) * wt(j) * hmat(indx)
           end if
    4   continue
    3 continue
      en = en + ec + shift
      print *,' Etot Michel = ', en

      end
