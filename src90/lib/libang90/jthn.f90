!
!     -------------------------------------------------------------
!      J T H N
!     -------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      INTEGER FUNCTION JTHN (K, N, I) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:53:44  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: K 
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: I 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      SELECT CASE (N)  
      CASE (1)  
         JTHN = MOD(K,I) 
      CASE (2)  
         JTHN = MOD(K/I,I) 
      CASE (3)  
         JTHN = MOD(K/(I*I),I) 
      CASE (4)  
         JTHN = MOD(K/(I*I*I),I) 
      CASE (5)  
         JTHN = MOD(K/(I*I*I*I),I) 
      CASE DEFAULT 
         WRITE (6, '(A)') ' ERROR IN JTHN ' 
         STOP  
      END SELECT 
      RETURN  
      END FUNCTION JTHN 
