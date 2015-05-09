*-----------------------------------------------------------------------
*               I N T V A L
*-----------------------------------------------------------------------
*     This routine returns the value of an integral. 
*
      Double Precision Function INTVAL(ICASE,K,i1,I2,I3,I4)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      POINTER (qintptr,intptr(1)),(qpackn,ipackn(1)),(qvalue,value(1))
      COMMON /TABLE/qintptr,qpackn,qvalue,lmax
      LOGICAL found
*
*-----------------------------------------------------------------------
*
*       print *, 'Intval: (icase,i1,i2,i3,i4,k)',
*    :          icase,i1,i2,i3,i4,k
        klocal = k
*       .. pack the data
ctc  NWD limited to 63 > 94 (root 64 -> 95)
        if (icase.le.2) then
          int = (klocal*95 + i2)*95 + i4
        else if (icase.eq.4 .or. icase.eq.5) then
          int = i2*95 + i4
        else
CGG          if (icase .eq. 6 .or. icase .eq. 8) klocal = klocal+1
          if (icase .eq. 6 .or. icase .eq. 8 .or. icase .eq. 9) 
     :      klocal = klocal+1
          int = (((klocal*95+I1)*95+I2)*95+I3)*95+I4
        endif
	n = isearch(icase,int,qpackn,qintptr,lmax)
	intval = value(n)
*        write (*,*) k,i1,i2,i3,i4,int,icase
*        write (*,*) value(n)

	end
