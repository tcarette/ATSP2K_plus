!
!     -------------------------------------------------------------
!      E I L E
!     -------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE EILE(JA, JB, JC, JAA, JBB, JCC) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:10:27  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: JA 
      INTEGER , INTENT(IN) :: JB 
      INTEGER , INTENT(IN) :: JC 
      INTEGER , INTENT(OUT) :: JAA 
      INTEGER , INTENT(OUT) :: JBB 
      INTEGER , INTENT(OUT) :: JCC 
!-----------------------------------------------
      JAA = JA 
      JCC = JA 
      JAA = MIN0(JB,JAA) 
      JCC = MAX0(JB,JCC) 
      JAA = MIN0(JC,JAA) 
      JCC = MAX0(JC,JCC) 
      IF (JA>JAA .AND. JA<JCC) JBB = JA 
      IF (JB>JAA .AND. JB<JCC) JBB = JB 
      IF (JC>JAA .AND. JC<JCC) JBB = JC 
      RETURN  
      END SUBROUTINE EILE 
