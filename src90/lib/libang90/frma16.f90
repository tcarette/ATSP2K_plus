!
!     -----------------------------------------------------------------
!      F R M A 1 6
!     -----------------------------------------------------------------
!
!      JT7 = 96 -- 100                                             *
!
      SUBROUTINE FRMA16(J1, J2, ISKA, IVAR) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:10:49  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: J1 
      INTEGER , INTENT(IN) :: J2 
      INTEGER , INTENT(OUT) :: ISKA 
      INTEGER , INTENT(OUT) :: IVAR 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(18) :: IS96A 
      INTEGER , DIMENSION(69) :: IS96B 
      INTEGER , DIMENSION(19) :: IS97A 
      INTEGER , DIMENSION(69) :: IS97B 
      INTEGER , DIMENSION(26) :: IS98A 
      INTEGER , DIMENSION(30) :: IS98B 
      INTEGER , DIMENSION(9) :: IS99A 
      INTEGER , DIMENSION(14) :: IS99B 
      INTEGER , DIMENSION(3) :: IS99C 
      INTEGER , DIMENSION(22) :: IS99D 
      INTEGER , DIMENSION(28) :: IS100A 
      INTEGER , DIMENSION(37) :: IS100B 
      INTEGER :: I2 
!-----------------------------------------------
!
      DATA IS96A/ 42042, 6*0, 86515, 2*0, 127413, -32670, 2*0, -980, 2*0, -7020&
         /  
      DATA IS96B/ -484484, -4840, 2*0, -9240, 4*0, 339768, -8424, 4*0, -100980&
         , -17280, 3*0, -222156, 2*0, -45144, 3*0, -155727, 4410, 9*0, 81081, 0&
         , 31590, 7*0, 41580, 7*0, 216216, 37908, 4*0, 77760, 5*0, 203148/  
      DATA IS97A/ 183456, -12740, 8*0, 247104, 396000, 2*0, -72765, -2695, 0, &
         -6435, -3465/  
      DATA IS97B/ 361405440, 2*0, 144938640, 0, -127244832, 2*0, -176185152, &
         44963919, -51693213, 41599272, 2*0, 581808000, -368934720, 23889840, &
         88172370, 0, -460791936, 75998670, 0, 44258049, 172114635, 3*0, &
         276701040, 30744560, 8*0, -298980864, 0, 24470160, 13176240, 23272320&
         , 5*0, -289877280, 6901840, 6*0, 112117824, -89927838, 0, 118107990, &
         -813280, 0, 737869440, -6352830, 2*0, 12773970, 201355616, -88516098, &
         -180310570/  
      DATA IS98A/ -24227280, 0, 42301600, 3*0, 34512075, 30274049, -772616, 3*0&
         , 314640, 24249280, 2180250, 2*0, 57246310, 0, -512325, -4684295, 3*0&
         , -2622000, 17724720/  
      DATA IS98B/ 551304, 3991064, -28196608, 5*0, -6530832, 12595176, 7*0, &
         9303255, 0, -22297275, -19953648, 0, 84816, -15903, 2*0, -7194825, &
         45239040, -138105, 17115813/  
      DATA IS99A/ -9542610, 44768185, -8690825, 3*0, -298350, -23585562, &
         14496768/  
      DATA IS99B/ -1415942100, 8388975, -2132697600, 2*0, -86295264, 0, &
         210258180, -1221458700, 3*0, 237437487, -354690567/  
      DATA IS99C/ 2722941, 19712231, -550/  
      DATA IS99D/ -18625761, -214639722, 7*0, -582335, 0, -181384203, 38594556&
         , 0, -53012325, 165419100, 2*0, -295689636, 15046955, 7871985, &
         -215588725/  
      DATA IS100A/ 1081080, 45045, -479160, -1155, 8*0, -198900, 414375, 4*0, &
         42075, -1197900, 4*0, -799425, 2*0, 611325/  
      DATA IS100B/ -8386560, 37440, -506880, 7*0, 987840, 2555520, 560560, 4*0&
         , 8910720, 0, -39780, 4369820, 6*0, -224400, 2*0, 4263600, 8796480, 4*&
         0, 8408400/  
!
      SELECT CASE (J1)  
      CASE (96)  
         IF (J2 < 139) RETURN  
         IVAR = 6534 
         IF (J2 < 157) THEN 
            I2 = J2 - 138 
            ISKA = IS96A(I2) 
         ELSE IF (J2 < 235) THEN 
            IF (J2 < 166) RETURN  
            I2 = J2 - 165 
            ISKA = IS96B(I2) 
         ENDIF 
      CASE (97)  
         IF (J2 < 139) RETURN  
         IF (J2 < 158) THEN 
            I2 = J2 - 138 
            ISKA = IS97A(I2) 
            IVAR = 10164 
         ELSE IF (J2 < 236) THEN 
            IF (J2 < 167) RETURN  
            I2 = J2 - 166 
            ISKA = IS97B(I2) 
            IVAR = 19325152 
         ENDIF 
      CASE (98)  
         IF (J2 == 140) THEN 
            ISKA = -12103 
            IVAR = 279 
         ELSE IF (J2 == 153) THEN 
            ISKA = -12825 
            IVAR = 1364 
         ELSE IF (J2 == 154) THEN 
            ISKA = -10051 
            IVAR = 12276 
         ELSE IF (J2 == 156) THEN 
            ISKA = 2223 
            IVAR = 1364 
         ELSE IF (J2 == 157) THEN 
            ISKA = -3325 
            IVAR = 1364 
         ELSE IF (J2 < 196) THEN 
            IF (J2 < 170) RETURN  
            I2 = J2 - 169 
            ISKA = IS98A(I2) 
            IVAR = 1254880 
         ELSE IF (J2 < 236) THEN 
            IF (J2 < 206) RETURN  
            I2 = J2 - 205 
            ISKA = IS98B(I2) 
            IVAR = 1522224 
         ENDIF 
      CASE (99)  
         IF (J2 == 140) THEN 
            ISKA = -4641 
            IVAR = 341 
         ELSE IF (J2 == 153) THEN 
            ISKA = 780759 
            IVAR = 60016 
         ELSE IF (J2 == 154) THEN 
            ISKA = -357 
            IVAR = 496 
         ELSE IF (J2 == 156) THEN 
            ISKA = 8168823 
            IVAR = 420112 
         ELSE IF (J2 == 157) THEN 
            ISKA = 14739 
            IVAR = 496 
         ELSE IF (J2 < 179) THEN 
            IF (J2 < 170) RETURN  
            I2 = J2 - 169 
            ISKA = IS99A(I2) 
            IVAR = 1575420 
         ELSE IF (J2 < 196) THEN 
            IF (J2 < 182) RETURN  
            I2 = J2 - 181 
            ISKA = IS99B(I2) 
            IF (I2 == 2) THEN 
               IVAR = 266662 
            ELSE 
               IVAR = 82131896 
            ENDIF 
         ELSE IF (J2 < 209) THEN 
            IF (J2 < 206) RETURN  
            I2 = J2 - 205 
            ISKA = IS99C(I2) 
            IVAR = 630168 
         ELSE IF (J2 < 236) THEN 
            IF (J2 < 214) RETURN  
            I2 = J2 - 213 
            ISKA = IS99D(I2) 
            IVAR = 13837439 
         ENDIF 
      CASE (100)  
         IF (J2 == 147) THEN 
            ISKA = 3600 
            IVAR = 77 
         ELSE IF (J2 == 148) THEN 
            ISKA = 2160 
            IVAR = 91 
         ELSE IF (J2 == 155) THEN 
            ISKA = 4860 
            IVAR = 77 
         ELSE IF (J2 < 189) THEN 
            IF (J2 < 161) RETURN  
            I2 = J2 - 160 
            ISKA = IS100A(I2) 
            IVAR = 4004 
         ELSE IF (J2 < 237) THEN 
            IF (J2 < 200) RETURN  
            I2 = J2 - 199 
            ISKA = IS100B(I2) 
            IVAR = 64064 
         ENDIF 
      END SELECT 
      RETURN  
      END SUBROUTINE FRMA16 
