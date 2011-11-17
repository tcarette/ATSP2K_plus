C23456789012345678901234567890123456789012345678901234567890123456789012
      subroutine sortcsf(ncfg,wt,ec,shift)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      logical flagpr

      POINTER (phmat,hmat(1))
      COMMON/SORT/PHMAT
      double precision wt(*)

      POINTER (piref,iref(1))
      POINTER (pnoref,noref(1))

      CALL alloc(piref,ncfg,8)
      CALL alloc(pnoref,ncfg,8)

      DO 10 i=2,ncfg
        noref(i) = 1
   10 CONTINUE

C#######################################################################
C
C  first loop DO the iteration to construct the reordered (same) list 
C
C  other loops : japl span the unselected CSF's
C                iref(i) span the selected CSF's
C
C######################################################################

      flagpr = .false.
      iref(1) = 1
      ncfgm1 = ncfg - 1
      eps1 = hmat(1) + ec + shift
      print *, ' norm = ', wt(1)**2, ' Total energy = ',eps1
      sumeps = wt(1)**2 * eps1

      DO 1 nref = 1,ncfgm1
        epsmax = 0.0
        jmax = 0

c       en = wt(1)**2 * eps1
c       print*, ' enter in sortcsf. en = ',en

        DO 2 japl = 2,ncfg
          epsapl = 0.0
          IF (noref(japl) .eq. 1) THEN
            contri = wt(japl)**2  
     :                   * (hmat(indice(japl,japl,ncfg)) + ec + shift)
            epsapl = epsapl + contri
c           en = en + contri

            DO 3 i = 1,nref
              IF (japl .lt. iref(i)) THEN
                contri = 2. * wt(iref(i)) * wt(japl) 
     :                              * hmat(indice(iref(i),japl,ncfg))
                epsapl = epsapl + contri
c               en = en + contri
              ELSE
                contri = 2. * wt(iref(i)) * wt(japl) 
     :                              * hmat(indice(japl,iref(i),ncfg))
                epsapl = epsapl + contri
c               en = en + contri
              END IF
    3       CONTINUE

            if (flagpr) print *, ' nref, japl, epsapl = ',nref,japl,epsapl
            abscurr = DABS(epsapl)

            IF (abscurr .gt. epsmax) THEN
              epsmax = abscurr
              eps = epsapl
              jmax = japl
              if (flagpr) print *, ' jmax, epsmax = ',jmax,epsmax
            END IF

          END IF
    2   CONTINUE

        IF ( jmax .eq. 0 ) THEN
          print *, ' erreur debug :'
          print *, ' jmax = 0 (IF B tjs false?) at iteration ', nref
          call exit
        ELSE
          noref(jmax) = 0
          nrefp1 = nref + 1
          iref(nrefp1) = jmax
          anorm = 0.0
          do 6 k = 1,nrefp1
             anorm = anorm + wt(iref(k))**2
    6     continue
          sumeps = sumeps + eps
          sumepsn = sumeps/anorm
          print *,' jmax = ', jmax,'  eps = ', eps
          print *,' norm = ',anorm,' Total energy = ',sumepsn
        END IF

    1 CONTINUE
 
      print *,' iref = ', (iref(i),i=1,ncfg)

      return
      END


      FUNCTION indice(i,j,ncfg)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

        indice = (j - 1)*(2 * ncfg - j)/2 + i
        RETURN

      END
