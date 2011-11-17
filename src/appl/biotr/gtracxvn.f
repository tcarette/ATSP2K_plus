*
*     ------------------------------------------------------------------
*	G T R A C X V N
*     ------------------------------------------------------------------
*
      SUBROUTINE GTRACXVN(I,LINP,IFIRST,NFOUND,BUF,
     &           IBUF,LBUF,JLEI,IORF,NTESTG)
*
* Obtain coupling coefficients
*   <L! E(jl il) !R>, where L and R labels CSF's
*                           l labels angular momentum of subshells
*                           j and i are principal quantum numbers of subshells
*
* This routine is currently 
* interfaced with the list of one-electron integrals
* produced by YNONH as communicated common block /CSF/
*
* The following is the list in the order : l  j i   R L coef
* i and l-value fixed by input : I, LINP
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BUF(LBUF),IBUF(4,LBUF)
      PARAMETER(NWD =128)
*
* Common block containing all coupling coefficients
*
      POINTER (qintgrl1,intgrl1(1)),(qintptr1,intptr1(1)),
     :        (qcnn1,cnn1(1)),(qjann1,jann1(1)),(qjbnn1,jbnn1(1)),
     :        (qintgrl2,intgrl2(1)),(qintptr2,intptr2(1)),
     :        (qcnn2,cnn2(1)),(qjann2,jann2(1)),(qjbnn2,jbnn2(1))
      COMMON /CSF/qintgrl1,nint1,qintptr1,qcnn1,ncoef1,qjann1,qjbnn1,
     :            qintgrl2,nint2,qintptr2,qcnn2,ncoef2,qjann2,qjbnn2
*. From absolute orbital number to relative (n) number 
      common/nlorb/nq1(nwd),lq1(nwd),nq2(nwd),lq2(nwd)
*. Local common block showing next element to be accessed
      COMMON/GTRACLOC/LISTCUR,IELMCUR
*
      NTESTL = 000000
      NTEST = MAX0(NTESTL,NTESTG)
*
* Loop over lists of coefficients, i.e lists with given l j i
*
C?    IF(NTEST.GE.100) THEN
C?       WRITE(6,*) ' GTRACXVN in service '
C?       WRITE(6,*) ' I LINP JLEI IFIRST ', I,LINP,JLEI,IFIRST
C?    END IF
*
      IF(IFIRST.EQ.1) THEN
*. Start of search
         LISTCUR = 1
         IELMCUR = 1
      END IF
*
      IF(IORF.EQ.1) THEN
        NLIST = NINT1
      ELSE
        NLIST = NINT2
      END IF
*
      LCUR = 0
      LISTOFF = LISTCUR
      DO LIST = LISTOFF,NLIST 
*. l i j in LIST
        IF(IORF.EQ.1) THEN
          LV = INTGRL1(LIST)
          LAST = INTPTR1(LIST)
*. Number of coupling coefficients in list
          IF(LIST.EQ.1) THEN    
            NCOUP = INTPTR1(1)
          ELSE
            NCOUP = INTPTR1(LIST)-INTPTR1(LIST-1)
          END IF
        ELSE IF( IORF.EQ.2) THEN
          LV = INTGRL2(LIST)
          LAST = INTPTR2(LIST)
*. Number of coupling coefficients in list
          IF(LIST.EQ.1) THEN    
            NCOUP = INTPTR2(1)
          ELSE
            NCOUP = INTPTR2(LIST)-INTPTR2(LIST-1)
          END IF
        END IF
        IF(NTEST.GE.10000) WRITE(6,*) ' LIST NCOUP LAST ',
     &  LIST,NCOUP,LAST
*. First element in list
        IOFF = LAST - NCOUP + 1
*
        IV = MOD(LV,64)
        LV = LV/64
        JV = MOD(LV,64)
        LV = LV/64
        IF(NTEST.GE.10000) 
     &  WRITE(6,*) ' IV JV LV ', IV,JV,LV
*. IV and JV runs here from L+1. BIOTRN starts with 1 so 
        if (iorf.eq.1) then
          LIV = lq1(IV)
          IV = nq1(IV) - LIV
          JV = nq1(JV) - LIV
        else
          LIV = lq2(IV)
          IV = nq2(IV) - LIV
          JV = nq2(JV) - LIV 
        endif
C        LIV = L(IV)
C        IV = N(IV) - LIV  
C        JV = N(JV) - LIV
        IF(NTEST.GE.10000) 
     &  WRITE(6,*) ' UPdated IV JV  ', IV,JV
*. Is this list of interest
        IF((JLEI.EQ.1.AND.LV.EQ.LINP.AND.IV.EQ.I.AND.JV.LT.I) .OR.
     &     (JLEI.EQ.2.AND.LV.EQ.LINP.AND.IV.EQ.I.AND.JV.EQ.I)) THEN
*. Congratulations, you belong to the choosen
*. First element in list to be read
          IELMF = IOFF-1+IELMCUR
*. Number of elements to be read
          NELM = MIN0(LBUF-LCUR,NCOUP-IELMCUR+1)
*. Transfer elements
          IF(NTEST.GE.10000) WRITE(6,*) ' IELMF NELM ', IELMF,NELM
Cmrg
      if ((lcur+nelm) .gt. lbuf) then
	print*,' Would you please increase LBUF: '
	print*,'           lcur+nelm in gtracxvn = ', lcur+nelm
	stop
      end if
Cmrg
          DO ICOUP = 1, NELM
            IBUF(1,LCUR+ICOUP) = JV
            IBUF(2,LCUR+ICOUP) = IV
            IF(IORF.EQ.1) THEN
              IBUF(3,LCUR+ICOUP) = JANN1(IELMF+ICOUP-1)
              IBUF(4,LCUR+ICOUP) = JBNN1(IELMF+ICOUP-1)
              BUF(LCUR+ICOUP)    = -2.0*CNN1(IELMF+ICOUP-1)
            ELSE IF (IORF.EQ.2) THEN
              IBUF(3,LCUR+ICOUP) = JANN2(IELMF+ICOUP-1)
              IBUF(4,LCUR+ICOUP) = JBNN2(IELMF+ICOUP-1)
              BUF(LCUR+ICOUP)    = -2.0*CNN2(IELMF+ICOUP-1)
            END IF
          END DO
*. Update pinters (pinters or pointers?)
          LCUR = LCUR + NELM
*. First element to be read the next time around
Cmrg      IF(NELM.EQ.NCOUP) THEN
          IF((ielmcur+nelm-1).EQ.NCOUP) THEN
            IELMCUR = 1
            LISTCUR = LIST+1
          ELSE
            IELMCUR =  IELMCUR + NELM
            LISTCUR = LIST
            GOTO 1001
          END IF
        END IF
      END DO
 1001 CONTINUE
      NFOUND = LCUR

*      
       IF( NTEST.GE.1000) THEN
        WRITE(6,*) ' GTRACXVN to your service '
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' I LINP JLEI IFIRST ', I,LINP,JLEI,IFIRST
        WRITE(6,*) 'Number of elements obtained',NFOUND
        WRITE(6,*) ' IBUF and BUF '
        DO 300 IELMNT = 1, NFOUND
          WRITE(6,'(E14.7,3X,4I4)') 
     &    BUF(IELMNT),(IBUF(K,IELMNT),K=1,4)
  300  CONTINUE
       WRITE(6,*)
*
      END IF
      RETURN
      END
