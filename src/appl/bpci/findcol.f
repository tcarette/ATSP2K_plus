************************************************************************
*
*
* --- This SUBROUTINE determines the number of columns already computed
*
*     Three files are considered:
*       Hnr -- non-relativistic integrals 
*       HZ  -- contributions from Z, N, and V integrals
*       HS  -- contributions from S integrals.
*
*
************************************************************************
*
      Subroutine FindCol(skip,jbegin,nze)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      POINTER (qh,h(ncfg,3)),(qjan,jan(ncfg,3))
      COMMON /buffer/qh, qjan, nrow(3), iflag(3)
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,
     :             iouhn,iouhz,iouhs,iouhm,iLS,idisk
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      LOGICAL skip
*
      jcolumn = 0
      if (.not. skip) then
        DO 2 JJ=1,NCFG
          JAN(JJ,2)=1
          JAN(JJ,3)=1
  2     CONTINUE
      end if
100   CONTINUE
      DO 1 JJ=1,NCFG
        JAN(JJ,1)=-1
        if (.not. skip) then
          JAN(JJ,2)=-1
          JAN(JJ,3)=-1
        ENDIF
  1   CONTINUE
*       .. read contents of each record
      read (iouhn,iostat=ion) jb,m,(h(i,1),i=1,m),(jan(i,1),i=1,m)
      nrow(1)=m
      IF(M.EQ.0) THEN
        IT1=1
      ELSE
        IT1=JAN(M,1)
      ENDIF
      if (.not. skip) then
        read(iouhz,iostat=ioz) jb,mz,(h(i,2),i=1,mz),(jan(i,2),i=1,mz)
        nrow(2)=mz
        IF(MZ.EQ.0) THEN
          IT2=1
        ELSE
          IT2=JAN(MZ,2)
        ENDIF
        read(iouhs,iostat=ios) jb,ms,(h(i,3),i=1,ms),(jan(i,3),i=1,ms)
        nrow(3)=ms
        IF(MS.EQ.0) THEN
          IT3=1
        ELSE
          IT3=JAN(MS,3)
        ENDIF
      end if
      IF((IABS(M)+IABS(MZ)+IABS(MS)).EQ.0 .OR. 
     :      MIN(IT1,IT2,IT3).LT.0) THEN
*        .. data for column not read correctly
	 backspace(iouhn)
	 if (.not. skip) then
	   backspace(iouhz)
	   backspace(iouhs)
	 end if
      else
*        ... all records advanced correctly
         if (idisk .eq. 0) then
            if (nze+m .gt. big) then
               nze = ncfg
               idisk = 1
            else
               nze = nze + m
            end if
         end if
	 jcolumn = jb
	 go to 100
       end if
       write(iscw,*) jcolumn, ' columns processed completely'
       jbegin = jcolumn+1
      Return
      END
