*
*------------------------------------------------------------------------
*        S A V E L S 
*------------------------------------------------------------------------
*
      SUBROUTINE SAVELS(ICASE,C,K,I1,I2,I3,I4,JA,JB,IPTR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER (LSDIM=30000)
      POINTER (qcn,cn(lsdim)), (qpackn,ipackn(lsdim)),
     :        (qjan,jan(lsdim)),(qjbn,jbn(lsdim)) 
      COMMON /buffer/qcn, qpackn, qjan, qjbn
      COMMON /fout/lcount,nrec,iflag,lij,nij
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(4),IALL,JSC(3),ISCW
*
*     L data
*
      
      if (lcount .eq. LSDIM) then
*        .. dump data to scratch disk, opening file if necessary
         iou = isc(4)
         write(iou) (cn(j),j=1,lsdim),(ipackn(j),j=1,lsdim),
     :              (jan(j),j=1,lsdim),(jbn(j),j=1,lsdim)
         nrec = nrec+1
         lcount = 0
       end if

      N = lcount+1
         IF (I2 .GT. I4) THEN
            II2 = I2
            II4 = I4
	    jan(n) = ja
	    jbn(n) = jb
         ELSE
            II2 = I4
            II4 = I2
*. do not forget to permute configuration indices!
	    jan(n) = jb
	    jbn(n) = ja
         END IF
         IPACK = (K*64 + II2)*64 + II4
         CN(n) = C
         IPACKN(n) = IPACK
      
      lcount = lcount + 1

      IFLAG = 1
      END
