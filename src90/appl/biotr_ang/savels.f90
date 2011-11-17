!
!------------------------------------------------------------------------
!        S A V E L S
!------------------------------------------------------------------------
!
      SUBROUTINE SAVELS(ICASE, C, K, I1, I2, I3, I4, JA, JB, IPTR) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      use parameters_biotr_C
      use buffer_C
      use fout_C
      use inform_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  00:10:35  11/17/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: ICASE 
      INTEGER , INTENT(IN) :: K 
      INTEGER  :: I1 
      INTEGER , INTENT(IN) :: I2 
      INTEGER  :: I3 
      INTEGER , INTENT(IN) :: I4 
      INTEGER , INTENT(IN) :: JA 
      INTEGER , INTENT(IN) :: JB 
      INTEGER  :: IPTR 
      REAL(DOUBLE) , INTENT(IN) :: C 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IOU, J, N, II2, II4, IPACK 
!-----------------------------------------------
!
!      PARAMETER (LSDIM=30000)
!      POINTER (qcn,cn(lsdim)), (qpackn,ipackn(lsdim)),
!     :        (qjan,jan(lsdim)),(qjbn,jbn(lsdim))
!      COMMON /buffer/qcn, qpackn, qjan, qjbn
!      COMMON /fout/lcount,nrec,iflag,lij,nij
!      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(4),IALL,JSC(3),ISCW
!
!     L data
!
      IF (LCOUNT == LSDIM) THEN 
!        .. dump data to scratch disk, opening file if necessary
         IOU = ISC(4) 
         WRITE (IOU) (CN(J),J=1,LSDIM), (IPACKN(J),J=1,LSDIM), (JAN(J),J=1,&
            LSDIM), (JBN(J),J=1,LSDIM) 
         NREC = NREC + 1 
         LCOUNT = 0
!         WRITE (50, '(A10,I7)') 'cn(j)', LSDIM 
!         WRITE (50, '(6F12.7)') (CN(J),J=1,LSDIM) 
!         WRITE (50, '(A10,I7)') 'ipackn(j)', LSDIM 
!         WRITE (50, '(6I12)') (IPACKN(J),J=1,LSDIM) 
!         WRITE (50, '(A10,I7)') 'jan(j)', LSDIM 
!         WRITE (50, '(6I12)') (JAN(J),J=1,LSDIM) 
!         WRITE (50, '(A10,I7)') 'jbn(j)', LSDIM 
!         WRITE (50, '(6I12)') (JBN(J),J=1,LSDIM) 
      ENDIF 
 
      N = LCOUNT + 1 
      IF (I2 > I4) THEN 
         II2 = I4 
         II4 = I2 
!. do not forget to permute configuration indices!
         JAN(N) = JB 
         JBN(N) = JA 
      ELSE 
         II2 = I2 
         II4 = I4 
         JAN(N) = JA 
         JBN(N) = JB 
      ENDIF 
!       WRITE(66,'(I1,F17.7,4I10)')K,C,II2,II4,JAN(N),JBN(N)
 
      IPACK = (K*64 + II2)*64 + II4 
      CN(N) = C 
      IPACKN(N) = IPACK 
 
      LCOUNT = LCOUNT + 1 
 
      IFLAG = 1 
      RETURN  
      END SUBROUTINE SAVELS 
