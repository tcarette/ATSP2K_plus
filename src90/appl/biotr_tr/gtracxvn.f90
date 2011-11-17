!
!     ------------------------------------------------------------------
!     G T R A C X V N
!     ------------------------------------------------------------------
!
      SUBROUTINE GTRACXVN(I, LINP, IFIRST, NFOUND, BUF, IBUF, LBUF, &
           JLEI, IORF, NTESTG) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE GTRACLOC_C 
      use csf_C
      use nlorb_C
!
!
! Obtain coupling coefficients
!   <L! E(jl il) !R>, where L and R labels CSF's
!                           l labels angular momentum of subshells
!                           j and i are principal quantum numbers of subshells
!
! This routine is currently
! interfaced with the list of one-electron integrals
! produced by YNONH as communicated common block /CSF/
!
! The following is the list in the order : l  j i   R L coef
! i and l-value fixed by input : I, LINP
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:23:09  11/20/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I 
      INTEGER , INTENT(IN) :: LINP 
      INTEGER , INTENT(IN) :: IFIRST 
      INTEGER , INTENT(OUT) :: NFOUND 
      INTEGER , INTENT(IN) :: LBUF 
      INTEGER , INTENT(IN) :: JLEI 
      INTEGER , INTENT(IN) :: IORF 
      INTEGER , INTENT(IN) :: NTESTG 
      INTEGER , INTENT(INOUT) :: IBUF(4,LBUF) 
      REAL(DOUBLE) , INTENT(INOUT) :: BUF(LBUF) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NTESTL, NTEST, NLIST, LCUR, LISTOFF, &
                 LIST, LV, LAST, NCOUP, IOFF, IV, JV, LIV, IELMF, &
                 NELM, ICOUP, IELMNT, K 
!-----------------------------------------------
!
!      PARAMETER(NWD =128)
!
! Common block containing all coupling coefficients
!
!      POINTER (qintgrl1,intgrl1(1)),(qintptr1,intptr1(1)),
!     :        (qcnn1,cnn1(1)),(qjann1,jann1(1)),(qjbnn1,jbnn1(1)),
!     :        (qintgrl2,intgrl2(1)),(qintptr2,intptr2(1)),
!     :        (qcnn2,cnn2(1)),(qjann2,jann2(1)),(qjbnn2,jbnn2(1))
!      COMMON /CSF/qintgrl1,nint1,qintptr1,qcnn1,ncoef1,qjann1,qjbnn1,
!     :            qintgrl2,nint2,qintptr2,qcnn2,ncoef2,qjann2,qjbnn2
!. From absolute orbital number to relative (n) number
!      common/nlorb/nq1(nwd),lq1(nwd),nq2(nwd),lq2(nwd)
!. Local common block showing next element to be accessed
!
      NTESTL = 0 
      NTEST = MAX0(NTESTL,NTESTG) 
!
! Loop over lists of coefficients, i.e lists with given l j i
!
!?    IF(NTEST.GE.100) THEN
!?       WRITE(6,*) ' GTRACXVN in service '
!?       WRITE(6,*) ' I LINP JLEI IFIRST ', I,LINP,JLEI,IFIRST
!?    END IF
!
      IF (IFIRST == 1) THEN 
!. Start of search
         LISTCUR = 1 
         IELMCUR = 1 
      ENDIF 
!
      IF (IORF == 1) THEN 
         NLIST = NINT1 
      ELSE 
         NLIST = NINT2 
      ENDIF 
!
      LCUR = 0 
      LISTOFF = LISTCUR 
      DO LIST = LISTOFF, NLIST 
!. l i j in LIST
         IF (IORF == 1) THEN 
            LV = INTGRL1(LIST) 
            LAST = INTPTR1(LIST) 
!. Number of coupling coefficients in list
            IF (LIST == 1) THEN 
               NCOUP = INTPTR1(1) 
            ELSE 
               NCOUP = INTPTR1(LIST) - INTPTR1(LIST-1) 
            ENDIF 
         ELSE IF (IORF == 2) THEN 
            LV = INTGRL2(LIST) 
            LAST = INTPTR2(LIST) 
!. Number of coupling coefficients in list
            IF (LIST == 1) THEN 
               NCOUP = INTPTR2(1) 
            ELSE 
               NCOUP = INTPTR2(LIST) - INTPTR2(LIST-1) 
            ENDIF 
         ENDIF 
         IF (NTEST >= 10000) WRITE (6, *) ' LIST NCOUP LAST ', LIST, NCOUP, &
            LAST 
!. First element in list
         IOFF = LAST - NCOUP + 1 
!
         IV = MOD(LV,64) 
         LV = LV/64 
         JV = MOD(LV,64) 
         LV = LV/64 
         IF (NTEST >= 10000) WRITE (6, *) ' IV JV LV ', IV, JV, LV 
!. IV and JV runs here from L+1. BIOTRN starts with 1 so
         IF (IORF == 1) THEN 
            LIV = LQ1(IV) 
            IV = NQ1(IV) - LIV 
            JV = NQ1(JV) - LIV 
         ELSE 
            LIV = LQ2(IV) 
            IV = NQ2(IV) - LIV 
            JV = NQ2(JV) - LIV 
         ENDIF 
!        LIV = L(IV)
!        IV = N(IV) - LIV
!        JV = N(JV) - LIV
         IF (NTEST >= 10000) WRITE (6, *) ' UPdated IV JV  ', IV, JV 
!. Is this list of interest
         IF (.NOT.(JLEI==1 .AND. LV==LINP .AND. IV==I .AND. JV<I .OR. JLEI==2&
             .AND. LV==LINP .AND. IV==I .AND. JV==I)) CYCLE  
!. Congratulations, you belong to the choosen
!. First element in list to be read
         IELMF = IOFF - 1 + IELMCUR 
!. Number of elements to be read
         NELM = MIN0(LBUF - LCUR,NCOUP - IELMCUR + 1) 
!. Transfer elements
         IF (NTEST >= 10000) WRITE (6, *) ' IELMF NELM ', IELMF, NELM 
!mrg
         IF (LCUR + NELM > LBUF) THEN 
            WRITE (6, *) ' Would you please increase LBUF: ' 
            WRITE (6, *) '           lcur+nelm in gtracxvn = ', LCUR + NELM 
            STOP  
         ENDIF 
!mrg
         IF (IORF == 1) THEN 
            IBUF(1,LCUR+1:NELM+LCUR) = JV 
            IBUF(2,LCUR+1:NELM+LCUR) = IV 
            IBUF(3,LCUR+1:NELM+LCUR) = JANN1(IELMF:NELM-1+IELMF) 
            IBUF(4,LCUR+1:NELM+LCUR) = JBNN1(IELMF:NELM-1+IELMF) 
            BUF(LCUR+1:NELM+LCUR) = -2.0*CNN1(IELMF:NELM-1+IELMF) 
         ELSE 
            IF (IORF == 2) THEN 
               IBUF(1,LCUR+1:NELM+LCUR) = JV 
               IBUF(2,LCUR+1:NELM+LCUR) = IV 
               IBUF(3,LCUR+1:NELM+LCUR) = JANN2(IELMF:NELM-1+IELMF) 
               IBUF(4,LCUR+1:NELM+LCUR) = JBNN2(IELMF:NELM-1+IELMF) 
               BUF(LCUR+1:NELM+LCUR) = -2.0*CNN2(IELMF:NELM-1+IELMF) 
            ELSE 
               IBUF(1,LCUR+1:NELM+LCUR) = JV 
               IBUF(2,LCUR+1:NELM+LCUR) = IV 
            ENDIF 
         ENDIF 
!. Update pinters (pinters or pointers?)
         LCUR = LCUR + NELM 
!. First element to be read the next time around
!mrg      IF(NELM.EQ.NCOUP) THEN
         IF (IELMCUR + NELM - 1 == NCOUP) THEN 
            IELMCUR = 1 
            LISTCUR = LIST + 1 
         ELSE 
            IELMCUR = IELMCUR + NELM 
            LISTCUR = LIST 
            EXIT  
         ENDIF 
      END DO 
      NFOUND = LCUR 
 
!
      IF (NTEST >= 1000) THEN 
         WRITE (6, *) ' GTRACXVN to your service ' 
         WRITE (6, *) ' =======================' 
         WRITE (6, *) ' I LINP JLEI IFIRST ', I, LINP, JLEI, IFIRST 
         WRITE (6, *) 'Number of elements obtained', NFOUND 
         WRITE (6, *) ' IBUF and BUF ' 
         DO IELMNT = 1, NFOUND 
            WRITE (6, '(E14.7,3X,4I4)') BUF(IELMNT), (IBUF(K,IELMNT),K=1,4) 
         END DO 
         WRITE (6, *) 
!
      ENDIF 
      RETURN  
      END SUBROUTINE GTRACXVN 
