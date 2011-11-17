      subroutine flsj_comp(ncfg,nterm,index,njv,maxj,minj,
     :               pflsj,lsp)
      PARAMETER (NOD=220, NWD=60, ntermd = 31)
      double precision w1,w2,nze_flsj,flsj;
      integer nterm
      logical lnze(nterm,njv)
      pointer(pflsj,nze_flsj(nterm,nterm,2,njv));
      integer ii,i,ls,ksi,li,jjj,ksj,llj,jj,lsp(ncfg),
     :        index(ntermd);
      integer jc,nt,njv,maxj,minj,nze_term,nze_col,tmp
      pointer(pflsj_save,flsj(1))
      pointer(pind,ind(nterm))
      integer ind_jj
  
      call alloc(plnze,nterm*njv,8)
      lnze(1:nterm,1:njv) = .false.
      call alloc(pflsj,nterm*nterm*2*njv,8);
      nze_flsj(1:nterm,1:nterm,1:2,1:njv) = 0.0;

      call alloc(pind,nterm,4)
      ind(1:nterm) = index(1:nterm);
      ksi = 0; lli = 0; ksj = 0; llj = 0;
      ind_jj = 0
      do jj = maxj,minj,-2
         ind_jj = ind_jj + 1
         DO II = 1,NTERM
            call lsp_cv(lsp(index(ii)),lli,ksi);
            DO JJJ = 1,NTERM
               call lsp_cv(lsp(index(jjj)),llj,ksj);
               PHASE = (-1)**((LLI+KSJ-JJ+LLI+KSI-JJ+LLJ+KSJ-JJ)/2) 
               CALL GRACAH(LLJ,KSJ,LLI,KSI,JJ,2,W1) 
               CALL GRACAH(LLJ,KSJ,LLI,KSI,JJ,4,W2);
               nze_FLSJ(II,JJJ,1,ind_jj) = PHASE*W1;
               nze_FLSJ(II,JJJ,2,ind_jj) = PHASE*W2
            end do
         end do
      end do
      return
      end

      subroutine lsp_cv(ti,lli,ksi)
      integer ti,lli,kli
        ls = ti/64;
        ksi = mod(ls,64);
        lli = (ls/64)
      return
      end 


