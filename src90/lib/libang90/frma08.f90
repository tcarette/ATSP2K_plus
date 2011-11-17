!
!     -----------------------------------------------------------------
!      F R M A 0 8
!     -----------------------------------------------------------------
!
!      JT7 = 47 -- 56                                              *
!
      SUBROUTINE FRMA08(J1, J2, ISKA) 
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
      INTEGER , DIMENSION(12) :: IS47, IS50A 
      INTEGER , DIMENSION(7) :: IS50B 
      INTEGER , DIMENSION(13) :: IS51A 
      INTEGER , DIMENSION(12) :: IS51B 
      INTEGER , DIMENSION(19) :: IS52A 
      INTEGER , DIMENSION(15) :: IS52B 
      INTEGER , DIMENSION(19) :: IS53A 
      INTEGER , DIMENSION(17) :: IS53B 
      INTEGER , DIMENSION(19) :: IS54A 
      INTEGER , DIMENSION(15) :: IS54B 
      INTEGER , DIMENSION(28) :: IS55A 
      INTEGER , DIMENSION(7) :: IS55B 
      INTEGER , DIMENSION(24) :: IS56A 
      INTEGER , DIMENSION(17) :: IS56B 
      INTEGER :: I2 
!-----------------------------------------------
!
      DATA IS47/ -21061215, 0, 69441120, 1492260, 0, 46334673, -29504860, 2*0, &
         -155669888, 132952512, -284898240/  
      DATA IS50A/ 10584, 7056, 2*0, -7200, -3300, 2*0, 8085, 0, 7920, 8775/  
      DATA IS50B/ -5040, 5760, -1100, 2695, 5544, -6336, 2925/  
      DATA IS51A/ 47040, 23716, 0, 17280, -220, 19440, -9680, 539, 5880, -19008&
         , 585, -21384, 25740/  
      DATA IS51B/ -10584, 3*0, -9600, -5940, 14553, 0, 10560, 15795, -18480, &
         20328/  
      DATA IS52A/ -49000, 7056, 0, 1246168, 418176, 231200, -264110, 2*0, &
         1387386, 380160, 257400, -39312, 672672, 2543541, -63063, 155232, &
         586971, -81081/  
      DATA IS52B/ 2*-378378, 0, 464640, -686664, 145530, 0, 422400, -764478, &
         29568, 26880, -432432, -87318, -450450, -389844/  
      DATA IS53A/ -32760, 327600, 0, 72072, 0, -26208, 90090, 2*0, -15246, 0, &
         5544, 155952, -66528, 2079, 33957, 288288, -9009, -72171/  
      DATA IS53B/ -4158, -55902, 2*0, 24024, 30030, 2*0, -5082, 2*0, 528, 18018&
         , 174262, -57596, -228096, 67584/  
      DATA IS54A/ -70560, -7056, 2*0, 104544, 332928, 4*0, 95040, 370656, 39312&
         , 0, -189189, 63063, 0, -43659, 81081/  
      DATA IS54B/ 378378, -42042, 7*0, -66528, -60480, -48048, 87318, -50050, &
         -43316/  
      DATA IS55A/ 2593080, -768320, 2*0, 691200, 1113750, 2*0, 1010625, 0, &
         -1425600, 43875, 3*0, 560560, 2*0, -1164240, 2*0, -408240, 1935360, &
         332640, -332640, 1576960, -1751750, 960960/  
      DATA IS55B/ 1088640, 552960, 1202850, 1091475, -2245320, -1140480, 47385&
         /  
      DATA IS56A/ 118800, -45000, 2*0, -235200, 0, 57600, 137280, 3*0, -458640&
         , 2*0, -294000, 2*0, -221760, -9240, 242760, 272160, 11340, -229320, &
         -136500/  
      DATA IS56B/ -332640, 1485, 15000, 78400, -161280, 720, -45760, 6*0, &
         -111475, 356720, 76440, -192080/  
!
      SELECT CASE (J1)  
      CASE (47)  
         IF (J2 < 135) RETURN  
         IF (J2 == 135) THEN 
            ISKA = 128707425 
         ELSE IF (J2 == 136) THEN 
            ISKA = -283156335 
         ELSE IF (J2 < 193) THEN 
            IF (J2 < 181) RETURN  
            I2 = J2 - 180 
            ISKA = IS47(I2) 
         ENDIF 
      CASE (48)  
         IF (J2 == 146) THEN 
            ISKA = -18 
         ELSE IF (J2 == 193) THEN 
            ISKA = -10 
         ENDIF 
      CASE (49)  
         SELECT CASE (J2)  
         CASE (153)  
            ISKA = -9 
         CASE (154)  
            ISKA = -27 
         CASE (194)  
            ISKA = -2 
         CASE (195)  
            ISKA = -18 
         END SELECT 
      CASE (50)  
         IF (J2 < 137) RETURN  
         IF (J2 < 149) THEN 
            I2 = J2 - 136 
            ISKA = IS50A(I2) 
         ELSE IF (J2 < 203) THEN 
            IF (J2 < 196) RETURN  
            I2 = J2 - 195 
            ISKA = IS50B(I2) 
         ENDIF 
      CASE (51)  
         IF (J2 < 138) RETURN  
         IF (J2 < 151) THEN 
            I2 = J2 - 137 
            ISKA = IS51A(I2) 
         ELSE IF (J2 < 205) THEN 
            IF (J2 < 193) RETURN  
            I2 = J2 - 192 
            ISKA = IS51B(I2) 
         ENDIF 
      CASE (52)  
         IF (J2 < 139) RETURN  
         IF (J2 < 158) THEN 
            I2 = J2 - 138 
            ISKA = IS52A(I2) 
         ELSE IF (J2 < 209) THEN 
            IF (J2 < 194) RETURN  
            I2 = J2 - 193 
            ISKA = IS52B(I2) 
         ENDIF 
      CASE (53)  
         IF (J2 < 139) RETURN  
         IF (J2 < 158) THEN 
            I2 = J2 - 138 
            ISKA = IS53A(I2) 
         ELSE IF (J2 < 211) THEN 
            IF (J2 < 194) RETURN  
            I2 = J2 - 193 
            ISKA = IS53B(I2) 
         ENDIF 
      CASE (54)  
         IF (J2 < 139) RETURN  
         IF (J2 < 158) THEN 
            I2 = J2 - 138 
            ISKA = IS54A(I2) 
         ELSE IF (J2 < 209) THEN 
            IF (J2 < 194) RETURN  
            I2 = J2 - 193 
            ISKA = IS54B(I2) 
         ENDIF 
      CASE (55)  
         IF (J2 < 137) RETURN  
         IF (J2 < 165) THEN 
            I2 = J2 - 136 
            ISKA = IS55A(I2) 
         ELSE IF (J2 < 203) THEN 
            IF (J2 < 196) RETURN  
            I2 = J2 - 195 
            ISKA = IS55B(I2) 
         ELSE IF (J2 == 211) THEN 
            ISKA = -1891890 
         ENDIF 
      CASE (56)  
         IF (J2 < 141) RETURN  
         IF (J2 < 165) THEN 
            I2 = J2 - 140 
            ISKA = IS56A(I2) 
         ELSE IF (J2 < 213) THEN 
            IF (J2 < 196) RETURN  
            I2 = J2 - 195 
            ISKA = IS56B(I2) 
         ENDIF 
      END SELECT 
      RETURN  
      END SUBROUTINE FRMA08 
