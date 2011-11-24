*
*------------------------------------------------------------------------
*        S A V E
*------------------------------------------------------------------------
*
      SUBROUTINE SAVELS(ICASE,C,K,I1,I2,I3,I4,JA,JB,IPTR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER (LSDIM=30000)
      POINTER (qcn,cn(lsdim)),(qinptr,inptr(lsdim)),
     :        (qnijptr,nijptr(lsdim)),(qjan,jan(lsdim)),
     :        (qjbn,jbn(lsdim)),(qintptr,intptr(0:2*lmax,4)),
     :        (qpackn,ipackn(1)),(qlused,lused(1)),(qico,ico(1))
      COMMON /buffer/qcn,qinptr,qpackn,qlused,qintptr,lmax,qnijptr,
     :               qjan,qjbn,qico
      POINTER  (qjptr, jptr(1))
      COMMON /fout/n,ntot,idum(6),nrec(8),iflag,lij,nij,qjptr
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(4),IALL,JSC(3),ISCW

      LOGICAL lused
*
      if (n .eq. LSDIM) then
*        .. write data to disk
         new = n
         write(50) new,(cn(j),j=1,new),(inptr(j),j=1,new)
         n = 1
      else
         n = n + 1
      end if

      nrec(icase) = nrec(icase) + 1
      ntot = ntot + 1

      IF (icase .LE. 2 .or. icase .eq. 4 ) THEN
         IF (I2 .GT. I4) THEN
            II2 = I4
            II4 = I2
         ELSE
            II2 = I2
            II4 = I4
         END IF
         IPACK = (K*64 + II2)*64 + II4
         CN(n) = C
         INPTR(n) = isearch(icase,ipack,qpackn,qlused,qintptr,lmax)
         nijptr(n) = nij + 1
      ELSE 
c@
c@     Rk data
c@
         J = 1
         IMIN = I1
         IF (I2 .LT. IMIN) THEN
            IMIN=I2
            J = 2
         END IF
         IF (I3 .LT. IMIN) THEN
            IMIN = I3
            J = 3
         END IF
         IF (I4 .LT. IMIN) THEN
            IMIN = I4
            J = 4
         END IF
         GO TO (10,20,30,40) J
10       II1 = I1
         II2 = I2
         II3 = I3
         II4 = I4
         Go to 50
        
20       II1 = I2
         II2 = I1
         II3 = I4
         II4 = I3
         GO TO 50

30       II1 = I3
         II2 = I4
         II3 = I1
         II4 = I2
         GO TO 50

40       II1 = I4
         II2 = I3
         II3 = I2
         II4 = I1

50       IPACK = (((K*64+II1)*64+II2)*64+II3)*64+II4
         CN(n) = C
         INPTR(n) = isearch(icase,ipack,qpackn,qlused,qintptr,lmax)
         NIJPTR(n) = nij + 1
      END IF

*
*      ... unpack the electron data
*
            iel1=-1
            iel2=-1
            iel3=-1
            iel4=-1
            kv = ipackn(inptr(n))
            if (icase.le.2.or.icase.eq.4) then
                iel2 = mod(kv,64) + nclosd
                kv = kv/64
                iel1 = mod(kv,64) + nclosd
                kv = kv/64
            else
                iel4 = mod(kv,64) + nclosd
                kv = kv/64
                iel3 = mod(kv,64) + nclosd
                kv = kv/64
                iel2 = mod(kv,64) + nclosd
                kv = kv/64
                iel1 = mod(kv,64) + nclosd
                kv = kv/64
            end if
  101       FORMAT(f15.10,2X,"F",I1,"(",I1,",",I1,")")
  102       FORMAT(f15.10,2X,"G",I1,"(",I1,",",I1,")")
  103       FORMAT(f15.10,2X,"R",I1,"(",I1,",",I1,";",I1,",",I1,")")
  104       FORMAT(f15.10,2X,"L","(",I1,",",I1,")???")
            if(icase.eq.1)then
              write(*,101) c,kv, iel1,iel2
            else if(icase.eq.2)then
              write(*,102) c,kv, iel1,iel2
            else if(icase.eq.3)then
              write(*,103) c,kv, iel1,iel2,iel3,iel4
            else
              write(*,104) c, iel1,iel2
            endif


      IFLAG = 1
      END

