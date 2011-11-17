!
!     -----------------------------------------------------------------
!     N U M T E R
!     -----------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville             October 1996   *
!
      INTEGER FUNCTION NUMTER (I2Q, I2S, I2L, L, NK, ND) 
      use mt67_C
      use mt_C
      use skmt2_C
      use kron_C
      use terms_C
      use mt15_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:48:12  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE jthn_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I2Q 
      INTEGER , INTENT(IN) :: I2S 
      INTEGER , INTENT(IN) :: I2L 
      INTEGER , INTENT(IN) :: L 
      INTEGER , INTENT(IN) :: NK 
      INTEGER , INTENT(IN) :: ND 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(3) :: LP, LG 
      INTEGER , DIMENSION(6) :: LP3, LG3 
      INTEGER :: IN, IL, JP, JG, J, MF 
!-----------------------------------------------
      !COMMON/MT/MT(40)
      !COMMON/SKMT2/MTF(8),MT3(90)
      !COMMON /MT67/ M76(238)
      !EXTERNAL TRMF,TERMLS,TRMF15
      DATA LP/ 1, 3, 9/  
      DATA LG/ 2, 8, 40/  
      DATA LP3/ 1, 11, 23, 37, 53, 71/  
      DATA LG3/ 10, 22, 36, 52, 70, 90/  
      NUMTER = 0 
      IF (L > 9) RETURN  
      IN = (I2Q*100 + I2S)*100 + I2L 
      IF (L < 3) THEN 
         IL = L + 1 
         JP = LP(IL) 
         JG = LG(IL) 
         J = JP 
    2    CONTINUE 
         IF (IN - MT(J) == 0) GO TO 1 
         J = J + 1 
         IF (J > JG) RETURN  
         GO TO 2 
    1    CONTINUE 
         NUMTER = J 
      ELSE IF (L == 3) THEN 
         IF (MAX0(NK,ND) < 3) THEN 
            JP = 1 
            JG = 8 
            J = JP 
    6       CONTINUE 
            IF (IN - MTF(J) == 0) GO TO 5 
            J = J + 1 
            IF (J > JG) RETURN  
            GO TO 6 
    5       CONTINUE 
            NUMTER = J + 300 
         ELSE 
            IF ((I2Q/2)*2 == I2Q) THEN 
               JP = 1 
               JG = 119 
            ELSE 
               JP = 120 
               JG = 238 
            ENDIF 
            J = JP 
   12       CONTINUE 
            MF = JTHN(M76(J),1,1000000) 
            IF (IN - MF == 0) GO TO 11 
            J = J + 1 
            IF (J > JG) RETURN  
            GO TO 12 
   11       CONTINUE 
            NUMTER = J 
         ENDIF 
      ELSE 
         IL = L - 3 
         JP = LP3(IL) 
         JG = LG3(IL) 
         J = JP 
   22    CONTINUE 
         IF (IN - MT3(J) == 0) GO TO 21 
         J = J + 1 
         IF (J > JG) RETURN  
         GO TO 22 
   21    CONTINUE 
         NUMTER = J 
      ENDIF 
      RETURN  
      END FUNCTION NUMTER 
