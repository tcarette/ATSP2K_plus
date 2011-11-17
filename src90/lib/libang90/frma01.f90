!
!     -----------------------------------------------------------------
!      F R M A 0 1
!     -----------------------------------------------------------------
!
!      JT7 = 1 -- 11                                               *
!
      SUBROUTINE FRMA01(J1, J2, ISKA) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:10:49  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: J1 
      INTEGER , INTENT(IN) :: J2 
      INTEGER , INTENT(OUT) :: ISKA 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(9) :: IS2 
      INTEGER , DIMENSION(12) :: IS3 
      INTEGER , DIMENSION(15) :: IS4 
      INTEGER , DIMENSION(16) :: IS5 
      INTEGER , DIMENSION(17) :: IS6, IS7 
      INTEGER , DIMENSION(8) :: IS10A 
      INTEGER , DIMENSION(13) :: IS10B 
      INTEGER , DIMENSION(8) :: IS11A 
      INTEGER , DIMENSION(19) :: IS11B 
      INTEGER :: I2 
!-----------------------------------------------
!
      DATA IS2/ -6174, 7560, 540, 1980, 1764, -4851, -8316, -594, -5265/  
      DATA IS3/ 10290, 0, 4320, -8910, 0, -8085, 0, -8910, -351, 3780, 3080, &
         14014/  
      DATA IS4/ -294, 80, -70, 0, -336, 0, 144, -126, 0, 54, 198, 0, -72, 208, &
         -182/  
      DATA IS5/ 67914, 0, -32670, -15840, 0, -7056, 0, 20250, -61776, -15246, &
         30030, -51744, 2*0, -34398, 70560/  
      DATA IS6/ -22638, 6160, 440, -2430, 6468, 1323, 16380, 1170, -2673, 2*0, &
         9702, 0, -35672, -2548, 13230, -14994/  
      DATA IS7/ 6006, 4*0, 1911, 0, -2106, -1485, 0, -4004, -220, 2*0, -5616, &
         -4590, -10098/  
      DATA IS10A/ 3360, 61440, -1980, -63504, 4851, -3696, -67584, 5265/  
      DATA IS10B/ 58800, 29645, 0, -21600, 17600, -43200, -12100, -43120, -4704&
         , 23760, -46800, 47520, 32175/  
      DATA IS11A/ -406560, 116160, 572220, 0, -121275, -369600, 105600, 637065&
         /  
      DATA IS11B/ -12005, 56448, 104544, 203456, 0, 56644, -43120, 0, 95040, &
         226512, 0, 63063, -314496, -1345344, -168168, -504504, -310464, -38808&
         , -648648/  
!
      SELECT CASE (J1)  
      CASE (1)  
         IF (J2 == 120) ISKA = -112 
      CASE (2)  
         IF (J2 < 120) RETURN  
         IF (J2 > 128) RETURN  
         I2 = J2 - 119 
         ISKA = IS2(I2) 
      CASE (3)  
         IF (J2 < 120) RETURN  
         IF (J2 > 131) RETURN  
         I2 = J2 - 119 
         ISKA = IS3(I2) 
      CASE (4)  
         IF (J2 < 120) RETURN  
         IF (J2 > 134) RETURN  
         I2 = J2 - 119 
         ISKA = IS4(I2) 
      CASE (5)  
         IF (J2 < 120) RETURN  
         IF (J2 > 135) RETURN  
         I2 = J2 - 119 
         ISKA = IS5(I2) 
      CASE (6)  
         IF (J2 < 120) RETURN  
         IF (J2 > 136) RETURN  
         I2 = J2 - 119 
         ISKA = IS6(I2) 
      CASE (7)  
         IF (J2 < 120) RETURN  
         IF (J2 > 136) RETURN  
         I2 = J2 - 119 
         ISKA = IS7(I2) 
      CASE (8)  
         IF (J2 == 124) ISKA = 100 
         IF (J2 == 137) ISKA = -72 
         IF (J2 == 138) ISKA = 108 
      CASE (9)  
         IF (J2 == 125) ISKA = -220 
         IF (J2 == 139) ISKA = 36 
         IF (J2 == 140) ISKA = -360 
      CASE (10)  
         IF (J2 < 121) RETURN  
         IF (J2 < 129) THEN 
            I2 = J2 - 120 
            ISKA = IS10A(I2) 
         ELSE IF (J2 < 151) THEN 
            IF (J2 < 138) RETURN  
            I2 = J2 - 137 
            ISKA = IS10B(I2) 
         ENDIF 
      CASE (11)  
         IF (J2 < 121) RETURN  
         IF (J2 < 129) THEN 
            I2 = J2 - 120 
            ISKA = IS11A(I2) 
         ELSE IF (J2 < 158) THEN 
            IF (J2 < 139) RETURN  
            I2 = J2 - 138 
            ISKA = IS11B(I2) 
         ENDIF 
      END SELECT 
      RETURN  
      END SUBROUTINE FRMA01 
