*
*     ------------------------------------------------------------------
*	B R K T I N (bra-ket within a space)
*     ------------------------------------------------------------------
*
      subroutine brktin(k)
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
Cmrg  dimension braket(nwd,nwd)
      nwf = iwf(k)
      if (k .eq. 1) then

      print*, ' One-electron overlaps WITHIN SPACE 1 '
      print*,'<LHS|LHS>'

      do 1 i = 1,nwf
        do 2 j = 1,nwf
          braket = quadr(i,j,0)
          if (l(i) .eq. l(j)) then
      print*,'<',elrasi(i),'|',elrasi(j),'> = ',braket
          end if
    2   continue
    1 continue

      else 

      print*, ' One-electron overlaps WITHIN SPACE 2 '
      print*,'<RHS|RHS>'

      do 3 i = 1,nwf
        iiwf = i + iwf(1)
        do 4 j = 1,nwf
          jiwf = j + iwf(1)
          braket = quadr(iiwf,jiwf,0)
          if (l(iiwf) .eq. l(jiwf)) then
      print*,'<',elrasf(i),'|',elrasf(j),'> = ',braket
          end if
    4   continue
    3 continue

      end if
      return
      end
