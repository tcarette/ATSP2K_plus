      subroutine term_shift(ncfg)
      PARAMETER (NOD=220, NWD=60, ntermd = 31)

      integer iscw,nterm,jjj,j,ls,ksj,llj,ncfg, ird

      COMMON /NCFGS/ IQLSP,index(ntermd),IJK,nterm,termsh(ntermd)
      POINTER (IQLSP, LSP(1))

      CHARACTER SYMBOL*11
      DATA SYMBOL/'SPDFGHIKLMN'/


      Write(iscw,'(A)') 'Enter the term energy shifts (in cm-1)'
      DO jjj = 1,nterm
         j = index(jjj)
         ls = lsp(j)/64
         ksj = mod(ls,64)+1
         llj = ls/64/2 + 1
         write(0,'(I3,A)') ksj,SYMBOL(llj:llj)
         Read(5,*) ird;
         termsh(jjj) = ird/2/109737.31534
      END DO

      return
      end

