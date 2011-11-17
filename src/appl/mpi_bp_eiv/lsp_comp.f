      subroutine lsp_comp(ncfg,lsp,jj)
      integer ncfg, lsp(ncfg)

*       .. determine which configuration states are "IN" (in=1)
         do i = 1,ncfg
            ls = lsp(i)/64
            ksi = mod(ls,64)
            lli = ls/64
            if (jj .ge. ABS(LLi-KSi) .and. jj .le. lli+ksi) then
               lsp(i) = 2*(lsp(i)/2) + 1
            else
               lsp(i) = 2*(lsp(i)/2)
            end if
         end do

      return
      end
