!
!     ------------------------------------------------------------------
!     G T R A C A
!     ------------------------------------------------------------------
!
      SUBROUTINE GTRACA(I, L, IFIRST, NFOUND, BUF, IBUF, LBUF, JLEI, LU) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!
! 1s inactive, 2 electrons in (1s 1p), Singlet S
!
! The following is the list in the order : l  j i   r l coef
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  20:39:47  11/18/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I 
      INTEGER , INTENT(IN) :: L 
      INTEGER , INTENT(IN) :: IFIRST 
      INTEGER , INTENT(OUT) :: NFOUND 
      INTEGER , INTENT(IN) :: LBUF 
      INTEGER , INTENT(IN) :: JLEI 
      INTEGER  :: LU 
      INTEGER , INTENT(INOUT) :: IBUF(4,LBUF) 
      REAL(DOUBLE) , INTENT(INOUT) :: BUF(LBUF) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NCOUP = 2 
      REAL(DOUBLE), PARAMETER :: SQRT2 = 1.414D0 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(5,NCOUP) :: IRACAH 
      INTEGER :: ISKIP, IELMNT, JELMNT, NTEST, K 
      REAL(DOUBLE), DIMENSION(NCOUP) :: RACAH 
!-----------------------------------------------
      DATA IRACAH/ 0, 2, 2, 1, 1, 1, 1, 1, 2, 2/  
!234567
!. Output
!
      DATA RACAH/ 2.0D0, 2.0D0/  
!
!. Position so IFIRST -1 elements have been skipped
!
 
      IF (IFIRST /= 1) THEN 
         ISKIP = 0 
         IF (JLEI == 0) THEN 
            DO IELMNT = 1, NCOUP 
               JELMNT = IELMNT 
               IF (IRACAH(3,IELMNT)/=I .OR. IRACAH(1,IELMNT)/=L) CYCLE  
!
               ISKIP = ISKIP + 1 
               IF (ISKIP /= IFIRST - 1) CYCLE  
               EXIT  
!
            END DO 
         ELSE 
            DO IELMNT = 1, NCOUP 
               JELMNT = IELMNT 
               IF (IRACAH(3,IELMNT)/=I .OR. IRACAH(1,IELMNT)/=L) CYCLE  
!
               IF (JLEI==1 .AND. IRACAH(2,IELMNT)<IRACAH(3,IELMNT)) THEN 
                  ISKIP = ISKIP + 1 
               ELSE 
                  IF (JLEI==2 .AND. IRACAH(2,IELMNT)==IRACAH(3,IELMNT)) THEN 
                     ISKIP = ISKIP + 1 
                  ELSE 
                     IF (JLEI==3 .AND. IRACAH(2,IELMNT)<=IRACAH(3,IELMNT)) &
                        ISKIP = ISKIP + 1 
                  ENDIF 
               ENDIF 
               IF (ISKIP /= IFIRST - 1) CYCLE  
               EXIT  
!
            END DO 
         ENDIF 
      ELSE IF (IFIRST == 1) THEN 
         JELMNT = 0 
      ENDIF 
!
! Obtain the next elements, atmost LBUF
!
      NFOUND = 0 
      IF (JLEI == 0) THEN 
         DO IELMNT = JELMNT + 1, NCOUP 
            IF (IRACAH(3,IELMNT)/=I .OR. IRACAH(1,IELMNT)/=L) CYCLE  
            NFOUND = NFOUND + 1 
            IBUF(1,NFOUND) = IRACAH(2,IELMNT) 
            IBUF(2,NFOUND) = IRACAH(3,IELMNT) 
            IBUF(3,NFOUND) = IRACAH(4,IELMNT) 
            IBUF(4,NFOUND) = IRACAH(5,IELMNT) 
            BUF(NFOUND) = RACAH(IELMNT) 
            IF (NFOUND /= LBUF) CYCLE  
            EXIT  
         END DO 
      ELSE 
         DO IELMNT = JELMNT + 1, NCOUP 
            IF (IRACAH(3,IELMNT)/=I .OR. IRACAH(1,IELMNT)/=L) CYCLE  
            IF (JLEI==1 .AND. IRACAH(2,IELMNT)<IRACAH(3,IELMNT)) THEN 
               NFOUND = NFOUND + 1 
               IBUF(1,NFOUND) = IRACAH(2,IELMNT) 
               IBUF(2,NFOUND) = IRACAH(3,IELMNT) 
               IBUF(3,NFOUND) = IRACAH(4,IELMNT) 
               IBUF(4,NFOUND) = IRACAH(5,IELMNT) 
               BUF(NFOUND) = RACAH(IELMNT) 
            ELSE 
               IF (JLEI==2 .AND. IRACAH(2,IELMNT)==IRACAH(3,IELMNT)) THEN 
                  NFOUND = NFOUND + 1 
                  IBUF(1,NFOUND) = IRACAH(2,IELMNT) 
                  IBUF(2,NFOUND) = IRACAH(3,IELMNT) 
                  IBUF(3,NFOUND) = IRACAH(4,IELMNT) 
                  IBUF(4,NFOUND) = IRACAH(5,IELMNT) 
                  BUF(NFOUND) = RACAH(IELMNT) 
               ELSE 
                  IF (JLEI==3 .AND. IRACAH(2,IELMNT)<=IRACAH(3,IELMNT)) THEN 
                     NFOUND = NFOUND + 1 
                     IBUF(1,NFOUND) = IRACAH(2,IELMNT) 
                     IBUF(2,NFOUND) = IRACAH(3,IELMNT) 
                     IBUF(3,NFOUND) = IRACAH(4,IELMNT) 
                     IBUF(4,NFOUND) = IRACAH(5,IELMNT) 
                     BUF(NFOUND) = RACAH(IELMNT) 
                  ENDIF 
               ENDIF 
            ENDIF 
            IF (NFOUND /= LBUF) CYCLE  
            EXIT  
         END DO 
      ENDIF 
      NTEST = 0 
      IF (NTEST /= 0) THEN 
         WRITE (6, *) ' GTRAC1 to your service ' 
         WRITE (6, *) ' =======================' 
         WRITE (6, *) 'Number of elements obtained', NFOUND 
         WRITE (6, *) ' IBUF and BUF ' 
         DO IELMNT = 1, NFOUND 
            WRITE (6, '(E12.7,3X,4I4)') BUF(IELMNT), (IBUF(K,IELMNT),K=1,4) 
         END DO 
!
      ENDIF 
      RETURN  
      END SUBROUTINE GTRACA 
