*
*     ------------------------------------------------------------------
*             A D V A N C E
*     ------------------------------------------------------------------
*
*     Advance the position of array j to the next relevant position
* 
      SUBROUTINE advance(j,len,ipos,ncfg)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NTERMD=31)
*
      DIMENSION ipos(3)
      POINTER (qh,h(ncfg,3)),(qjan,jan(ncfg,3))
      COMMON /buffer/qh, qjan, nrow(3), iflag(3)
      POINTER (IQLSP,LSP(1))
      COMMON /NCFGS/ IQLSP,index(ntermd),IJK,flsj(ntermd,ntermd,2),
     :               termsh(ntermd),nterm
*
*     ------------------------------------------------------------------
*
      i = ipos(j)
   10 i = i+1
      if ( i .gt. len) then
	ipos(j) = ncfg+1
	nrow(j) = ncfg+1
	return
      else
	irow = jan(i,j)
	if (mod(lsp(irow),2) .eq. 0) then
	  go to 10
	else
	  ipos(j) = i
	end if
      end if
      nrow(j) = jan(i,j)
*     print *, ' File, pos, row:',j,i,nrow(j),len
      end
