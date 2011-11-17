*
*     ------------------------------------------------------------------
*	R A S I N
*     ------------------------------------------------------------------
*
* --- Sets up the RAS information of shells. 
*     It is basically the original rasin module, but with all the
*     the active shells in ras2, assuming that both CSF expansions
*     satisfy c.u.d. where de-excitation means n->n' with n' .leq. n, 
*     for a given l-value.
*
      subroutine rasin
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=128)
      character*3 elras,elrasi,elrasf
      character bl*1
      CHARACTER*7 LABEL(2)
      CHARACTER*22 SET
      DATA LABEL/'Initial','Final  '/
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2),iut(2)
      COMMON /RAS/ ninac(11,2),nras1(11,10,2),nras2(11,2),
     :             nras3(11,10,2),NGAS1(2),NGAS3(2),
     :             nl(11,2),elras(2*nwd),elrasi(nwd),elrasf(nwd),
     :             itab(nwd,2),lmax(2)
      COMMON /FO/IWF(2),NCLOS(2),NOPEN(2),MAXNFO
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      DATA SET/'spdfghiklmnSPDFGHIKLMN'/
  10  format(' inactive : ',11I5)
  11  format(' RAS1     : ',11I5)
  14  format(' Number of shells per l for ',A7,' state')
  15  format('          : ',11I5)
  16  format(' Labels of shells for ',A7,' state')
      bl=' '
      write(iscw,*) 
      print*,' RAS-type calculation for initial state: '
      print*,' --------------------------------------: '
      print*,'               s    p    d    f    g    h    i    k    l  
     :  m    n '
      write(iwrite,10) (ninac(i,1),i=1,LMAX(1))
        WRITE(6,'(A,11I5)')
     &  ' GAS      : ',(NRAS2(L,1),L=1,LMAX(1))
      print*
      write(iscw,*) 
      print*,' RAS-type calculation for final   state: '
      print*,' --------------------------------------: '
      print*,'               s    p    d    f    g    h    i    k    l  
     :  m    n '
      write(iwrite,10) (ninac(i,2),i=1,LMAX(2))
        WRITE(6,'(A,11I5)')
     &  ' GAS      : ',(NRAS2(L,2),L=1,LMAX(2))
      print*
*
      do 2 is = 1,2
        write(iwrite,14) label(is)
        do 1 i = 1,lmax(is)
          nl(i,is) = ninac(i,is)+nras2(i,is)
    1   continue
      print*,'               s    p    d    f    g    h    i    k    l  
     :  m    n '
        write(iwrite,15) (nl(i,is),i=1,lmax(is))
    2 continue
      print*
*
      io = 0
      do 7 i = 1,2
        write(iwrite,16) label(i)
        do 8 j = 1,lmax(i)
          n = j
* Inactive
          ic = ninac(j,i)
    3     if (ic.ne.0) then
            ic = ic-1
            io = io+1
            elras(io) = bl//char(ichar('0')+n)//set(j:j)
            print*,' inactive: ',n,set(j:j)
            n = n+1
            go to 3
          end if
*. GAS
	  ic = nras2(j,i)
    5     if (ic.ne.0) then
            ic = ic-1
            io = io+1
            elras(io) = bl//char(ichar('0')+n)//set(j:j)
            print*,'   GAS     : ',n,set(j:j)
            n = n+1
            go to 5
          end if
*
    8  continue
    7  continue
       print*,(elras(k),k=1,maxnfo)
       end
