      subroutine coeff_disk(ncoef,ico,nz,hmx,
     :           coeff,value,inptr,jptr,jjh,nijcurr,ncodim)

      implicit double precision(a-h,o-z)
      parameter(MEIG=20,MTERM=20)
      INCLUDE 'mpif.h'
      parameter (MAXPROC=100)
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      common /PVM/ istart,ifinish

      dimension ico(*), hmx(*), coeff(*), value(*), inptr(*)
      dimension jptr(*)


      max_col = 1;
      ncoef = 0;
        do
          read(39) num,(coeff(j),j=1,num),
     :                 (inptr(j),j=1,num)
          do ii = 1, num
*       .. test for next non-zero matrix element
            if (ncoef+ii.gt.ico(nijcurr+nz)) nijcurr = nijcurr + 1
            hmx(nijcurr) = hmx(nijcurr) + coeff(ii)*value(inptr(ii))
            if (nijcurr.gt.jptr(jjh)) then
              jjh = jjh + 1;
              max_col = max_col+1
             end if
*           print '(3I5,3F20.15,I4)',
*     :     ii,nijcurr,inptr(ii),coeff(ii),value(inptr(ii)),
*     :     hmx(nijcurr)
          end do
          ncoef = ncoef + num
          if (num .lt. ncodim) exit
        end do
*     .. compute diagonals

        end
      subroutine coeff_memory(n_cf,ncoef,ico,nz,hmx,
     :           coeff,value,inptr,jptr,jjh,nijcurr)

      implicit double precision(a-h,o-z)
      parameter(MEIG=20,MTERM=20)
      INCLUDE 'mpif.h'
      parameter (MAXPROC=100)
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)
      common /PVM/ istart,ifinish

      dimension ico(*), hmx(*), coeff(*), value(*), inptr(*)
      dimension jptr(*)

      max_col = 1; 
      do ii = 1, n_cf;
*       .. test for next non-zero matrix element
          n_count_tmp = ncoef+ii
          if (ii.gt.ico(nijcurr+nz)) nijcurr = nijcurr + 1
          hmx(nijcurr) = hmx(nijcurr) +
     :            coeff(n_count_tmp)*value(inptr(n_count_tmp))
          if (nijcurr.gt.jptr(jjh)) then
            jjh = jjh + 1;
            max_col = max_col+1
          end if
*            if (myid == 0) then
*          print '(3I5,3F20.15,I4)',
*     :     ii,nijcurr,inptr(ncoef+ii),coeff(ncoef+ii),
*     :                value(inptr(ncoef+ii)),
*     :     hmx(nijcurr)
*             end if
      end do

      ncoef = ncoef + n_cf;

      end 
