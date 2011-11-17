!
!     ---------------------------------------------------------------
!     M E S
!     ---------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                          December 1993   *
!
      SUBROUTINE MES(I) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:52:55  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J 
      CHARACTER, DIMENSION(6) :: STRING5*10 
!-----------------------------------------------
      DATA STRING5/ '   I T L S', ' I T L S 2', ' I T L S 3', 'AW P 1 L S', &
         'WA P 1 L S', '       W 1'/  
      IF (I > 50) THEN 
         J = I - 50 
         WRITE (6, '(A)') ' error in func./sub. ' 
         WRITE (6, '(20X,A10)') STRING5(J) 
         WRITE (6, '(A)') ' susimaise f sluoksnio termu kodavimas  ' 
      ELSE 
         WRITE (6, '(A)') ' yra daugiau nei 2 ele. sluoks. f,g,h,i,k,l,m' 
         WRITE (6, '(3X,I5)') I 
         SELECT CASE (I)  
         CASE (1)  
            WRITE (6, '(A)') ' error in Subroutine   W 1 G ' 
         CASE (2)  
            WRITE (6, '(A)') ' error in Subroutine   W 1 ' 
         CASE (11)  
            WRITE (6, '(A)') ' error in Function     I T L S  ' 
         CASE (12)  
            WRITE (6, '(A)') ' error in Function     I T L S 2  ' 
         CASE (13)  
            WRITE (6, '(A)') ' error in Function     I T L S 3  ' 
         CASE (30)  
            WRITE (6, '(A)') ' error in Subroutine   A 1 A 2 A 3 A 4 L S ' 
         CASE (31)  
            WRITE (6, '(A)') ' error in Subroutine   A 1 A 2 L S ' 
         CASE (32)  
            WRITE (6, '(A)') ' error in Subroutine   A 1 A 2 W 3 L S ' 
         CASE (33)  
            WRITE (6, '(A)') ' error in Subroutine   A 1 A W 2 L S ' 
         CASE (34)  
            WRITE (6, '(A)') ' error in Subroutine   W A 1 A 2 L S ' 
         CASE (35)  
            WRITE (6, '(A)') ' error in Subroutine   W 1 W 2 L S ' 
         CASE DEFAULT 
            WRITE (6, '(A)') ' error in unknown Subroutine  ' 
         END SELECT 
      ENDIF 
      WRITE (6, '(A)') ' Contact to   G. Gediminas please ! ! ! ' 
      STOP  
      END SUBROUTINE MES 
