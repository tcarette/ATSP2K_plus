*-----------------------------------------------------------------------
*     C O U P E X       written by Thomas Caretten, Stockholm, Nov 2011
*-----------------------------------------------------------------------
*
*
*  Read the CSFs in clist (parent list).
*  Couple the extra/excited electron.
*  Write the result to the new list.
*
      subroutine coupex(lmax,term)

      implicit double precision(a-h,o-z)

      CHARACTER*13 orbm,orbc
      PARAMETER (ORBM='spdfghiklmnoq')
      PARAMETER (ORBC='SPDFGHIKLMNOQ')
      common/inform/iread,iwrite,iout,isc0,isc1,isc2,isc3,
     : iall,jsc(3),iscw

      character*3 term
      character*2 lsp,lsf
      character*1 parity,slc,slm

      character*64 line1,line2
      integer Sf,Sp

9     format(a)

      ireadp=iread+1
      ireadr=iread+2
*
* This is totally unoptimized but it is not supposed to be used for large
* lists
*

      lsf=term(1:2)

      Lf=2*lval(lsf(2:2))
      read(lsf(1:1),'(i1)') Sf
      Sf=Sf-1


   50 read(ireadp,9,end=60) line1
      if (line1(1:1) .ne.'*') then
        read(ireadp,9) line2

        if(parity(line1).eq.term(3:3))then
          l=0
        else
          l=1
        endif
        do 
           l_=l+1
           slc=orbc(l_:l_)
           slm=orbm(l_:l_)

           k2 = LEN_TRIM(line2);
  
           lsp= line2(k2-1:k2)
  
           Lp=2*lval(lsp(2:2))
           read(lsp(1:1),'(i1)') Sp
           Sp=Sp-1
  
           tst=TRITST(sp,sf,1)+TRITST(lp,lf,2*l)
           if(tst.lt.0.5d0)then
             k1 = LEN_TRIM(line1);
             write(iread,9,advance='no') line1(1:k1) 
             write(iread,'(a8)') '  ;'//slm//'( 1)'
             write(iread,9,advance='no') line2(1:k1/2) 
             write(iread,'(a4)',advance='no') ' 2'//slc//'1'
             write(iread,9,advance='no') line2(k1/2+1:k2) 
             write(iread,'(a4)') '  '//lsf//' '
           endif

           if(l_.ge.lmax) goto 700
  
           l=l_+1
         enddo
700      continue
         goto 50
      endif
60    continue

      end subroutine
