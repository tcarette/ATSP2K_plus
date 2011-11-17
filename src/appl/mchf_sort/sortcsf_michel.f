C23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine sortcsf(ncfg,wt,ec,shift)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      POINTER (phmat,hmat(1))
      COMMON/SORT/PHMAT
      double precision wt(*)

      POINTER (piref,iref(1))
      POINTER (pnoref,noref(1))

      CALL alloc(piref,ncfg,8)
      CALL alloc(pnoref,ncfg,8)


      iref(1) = 1
      ncfgm1 = ncfg - 1
      eps1 = hmat(1) + ec + shift
      print *, ' norm = ', wt(1)**2, 'Total energy = ',eps1
      sumeps = wt(1)**2 * eps1

      ncfgref = 0
      kref = 1

initialiser le vecteur iflag(1,ncfg) a 0
iflag = 0 pas encore pris ; =1 inclus dans la reference

      do 1 i=2,ncfg
        iflag(i) = 0
    1 continue

trouver l'indice pour l''interaction maximale

      ncfgm1 = ncfg - 1
      do 2 i=2,ncfgm1
         do 3 j = 1,ncfg
            contri = 2. wt(iref(i)) * wt(j))
            if (iflag(j). ) 
              
            IF (abscurr .gt. epsmax) THEN
              epsmax = abscurr
              eps = epsapl
              jmax = japl
              print *, ' jmax, epsmax = ',jmax,epsmax
            END IF

      indxref(kref) = jmax
      iflag(kref) = 1

   etendre la reference

      ncfgref = ncfgref + 1
      iflag
 
   






      ncfgref = ncfgref + 1

      do 2 j=


      DO 1 nref = 1,ncfg

        DO 2 japl = 2,ncfg
  
        epsmax = 0.0
        jmax = 0
        k = ncfg


        DO 2 japl = 2,ncfg
          epsapl = 0.0
          IF (noref(japl) .eq. 1) THEN
            kp1 = k + 1
            contri = wt(japl)**2 * (hmat(kp1) + ec + shift)
            epsapl = epsapl + contri

            DO 3 i = 1,nref
              kpiref = k + iref(i)
              contri = 2. * wt(iref(i)) * wt(japl) * hmat(kpiref)
              epsapl = epsapl + contri
    3       CONTINUE

            print *, ' nref, japl, epsapl = ',nref,japl,epsapl
            abscurr = DABS(epsapl)


          END IF
          k = k + ncfg - japl + 1
    2   CONTINUE

        IF ( jmax .eq. 0 ) THEN
          print *, ' erreur debug :'
          print *, ' jmax = 0 (IF B tjs false?) at iteration ', nref
          call exit
        ELSE
          noref(jmax) = 0
          nrefp1 = nref + 1
          iref(nrefp1) = jmax
          print *, ' Slcted config. at iteration ', nref, ': ', jmax
          anorm = 0.0
          do 6 kj = 1,nrefp1
             anorm = anorm + wt(iref(kj))**2
    6     continue
          sumeps = sumeps + eps
c         sumepsn = sumeps/anorm
          print *,' sumeps = ',sumeps
          print *,'norm=',anorm,'Total energy=',sumepsn ,'eps=',eps
           
        END IF

c   1 CONTINUE

      return
      END
