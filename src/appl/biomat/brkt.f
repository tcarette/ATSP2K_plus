*
*     ------------------------------------------------------------------
*	B R K T    (for "bra-ket")
*     ------------------------------------------------------------------
*
      subroutine brkt
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=128,NOD=220)
      COMMON /FO/IWF(2),NCLOS(2),NOPEN(2),MAXNFO
      COMMON /RAS/ ninac(11,2),nras1(11,10,2),nras2(11,2),
     :             nras3(11,10,2),NGAS1(2),NGAS3(2),
     :             nl(11,2),elras(2*nwd),elrasi(nwd),elrasf(nwd),
     :             itab(nwd,2),lmax(2)
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON /NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7),el(nwd)
      character*3 elras,elrasi,elrasf,el
Cmrg      dimension braket(nwd,nwd)

      print*,'<LHS|RHS>'

      do 1 i = 1,iwf(1)
        do 2 j = 1,iwf(2)
          braket = quadr(i,j+iwf(1),0)
          if (l(i) .eq. l(j+iwf(1))) then
      print*,'<',elrasi(i),'|',elrasf(j),'> = ',braket
          end if
    2   continue
    1 continue
      return
      end
