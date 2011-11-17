      subroutine diag_memory_all(n_cf,ncoef,ico,nz,hmx,
     :           coeff,value,inptr,jptr,jjh,nijcurr,nze)

      implicit double precision(a-h,o-z)
      parameter(MEIG=20,MTERM=20)

      dimension ico(nze), hmx(nze), coeff(*), value(*), inptr(*)
      dimension jptr(*)
      hmx(1:nze) = 0.0

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
*          print '(3I5,3F20.15,I4)',
*     :     ii,nijcurr,inptr(ncoef+ii),coeff(ncoef+ii),
*     :                value(inptr(ncoef+ii)),
*     :     hmx(nijcurr)
        end do

      ncoef = ncoef + n_cf;

      end 

****  
      subroutine diag_disk_clst(ncoef,ico,hmx,coeff,value,
     :           inptr,nijcurr,ncodim,nze,nz)

      implicit double precision(a-h,o-z)
      parameter(MEIG=20,MTERM=20)

      dimension ico(nz), hmx(nz), coeff(*), value(*), inptr(*)

      hmx(1:nze) = 0
      ncoef = 0;
        do
          read(39) num,(coeff(j),j=1,num),
     :                 (inptr(j),j=1,num)
          do ii = 1, num
*       .. test for next non-zero matrix element
            if (ncoef+ii.gt.ico(nijcurr+nz)) nijcurr = nijcurr + 1
            hmx(nijcurr) = hmx(nijcurr) + coeff(ii)*value(inptr(ii))
*           print '(3I5,3F20.15,I4)',
*     :     ii,nijcurr,inptr(ii),coeff(ii),value(inptr(ii)),
*     :     hmx(nijcurr)
          end do
          ncoef = ncoef + num
          if (num .lt. ncodim) exit
        end do
*     .. compute diagonals

        end

*****
      subroutine diag_disk_ico(ncoef,ico,hmx,coeff,value,inptr,
     :           jptr,jjh,nijcurr,ncodim,ncfg,nze,nze_col)

      implicit double precision(a-h,o-z)
      parameter(MEIG=20,MTERM=20)
      logical write_hmx 
      dimension ico(nze_col),hmx(nze), coeff(*), value(*), inptr(*)
      dimension jptr(ncfg);

      ncoef = 0;
      nijcurr = 1;
      nij_count = 1;
      hmx(1:nze) = 0.0
      read(12) nico, (ico(jj),jj=1,nico)
      itmp = jptr(1)
      write_hmx = .false.
      rewind(13)
 
      do
         read(39) num,(coeff(j),j=1,num), (inptr(j),j=1,num)

         do ii = 1, num

           if (ncoef+ii > ico(nij_count)) then
              nijcurr = nijcurr + 1
              nij_count = nij_count + 1
           end if

           hmx(nijcurr) = hmx(nijcurr) + 
     :                        coeff(ii)*value(inptr(ii))

           if (ncoef+ii > ico(nico)) then 
              if (jjh < ncfg) then
                 read(12) nico ,(ico(jj),jj=1,nico)
                 nij_count = 1
              end if
              if (nijcurr-1 == jptr(jjh)) then
                 jjh = jjh + 1
              end if
           end if

*           print '(3I5,3F20.15,I4)',
*     :     ii,nijcurr,inptr(ii),coeff(ii),value(inptr(ii)),
*     :     hmx(nijcurr)
         end do
         ncoef = ncoef + num
         if (num .lt. ncodim) exit
      end do
      End subroutine diag_disk_ico

*************
      subroutine diag_disk_hmx(ncoef,ico,hmx,coeff,value,inptr,
     :           jptr,jjh,nijcurr,ncodim,ncfg,nze_max,nze,hii,
     :           shift_hmx)

      implicit double precision(a-h,o-z)
      parameter(MEIG=20,MTERM=20)
      logical next_nze,read_ico 
      dimension ico(nze_max), hmx(nze_max), coeff(*), value(*), inptr(*)
      dimension jptr(ncfg), hii(ncfg);

      ncoef = 0;
      nijcurr = 1;
      nij_count = 1;
      hmx(1:nze_max) = 0.0
      ico(1:nze_max) = 0
      read(12) nico, (ico(jj),jj=1,nico)
      itmp = jptr(1)
      read_ico = .false.
      next_nze = .false.
      shift_hmx = 0
      rewind(13)
 

      do
         read(39) num,(coeff(j),j=1,num), (inptr(j),j=1,num)

         do ii = 1, num
           hmx(nij_count) = hmx(nij_count) +
     :                        coeff(ii)*value(inptr(ii))
           if (ncoef+ii == ico(nij_count))  next_nze = .true.

*           print '(3I5,3F20.15,I4)',
*     :     ii,nij_count,inptr(ii),coeff(ii),value(inptr(ii)),
*     :     hmx(nij_count)

           if ((next_nze).and.(nij_count == 1)) then
              if (jjh == 1) shift_hmx = hmx(1)
              hii(jjh)= hmx(nij_count) - shift_hmx
              hmx(nij_count) = hii(jjh)
            end if

           if ((next_nze).and.(nijcurr == jptr(jjh))) then
*              print *, nico, (hmx(ihmx),ihmx=1,nico);
              write(13) nico, (hmx(ihmx),ihmx=1,nico);
              hmx(1:nze_max) = 0.0;
              ico(1:nze_max) = 0.0;
              read_ico = .true.
           end if;

           if (next_nze) then
              nijcurr = nijcurr + 1
              nij_count = nij_count + 1
              next_nze = .false.
           end if

           if (read_ico) then
              if (jjh < ncfg) read(12) nico ,(ico(jj),jj=1,nico)
              nij_count = 1;
              jjh = jjh + 1;
              read_ico = .false.
           end if;

         end do
         ncoef = ncoef + num
         if (num .lt. ncodim) exit
      end do
      hii(1) = 0.0;
      end subroutine diag_disk_hmx

