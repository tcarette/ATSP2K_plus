      subroutine term_shift(iscw,nterm,ncfg,indx,lsp,termsh)

      integer iscw,nterm,jjj,j,ls,ksj,llj,lsp(ncfg),ncfg
      double precision termsh(nterm)
      CHARACTER SYMBOL*11
      DATA SYMBOL/'SPDFGHIKLMN'/


      Write(iscw,'(A)') 'Enter the term energy shifts (in cm-1)'
      DO jjj = 1,nterm
         j = indx(jjj)
         ls = lsp(j)/64
         ksj = mod(ls,64)+1
         llj = ls/64/2 + 1
         write(iscw,'(I3,A)') ksj,SYMBOL(llj:llj)
         Read(5,*) termsh(jjj)
         termsh(jjj) = termsh(jjj)/2/109737.31534
      END DO

      return
      end
