!     ----------------------------------------------------------------
!      C O U P L I N G
!     ----------------------------------------------------------------
!
 
      SUBROUTINE COUPLING(JA, JB) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE MEDEFN_A_C 
      USE OCCUPATION_C 
      use ndims_C
!...Translated by Pacific-Sierra Research 77to90  4.3E  12:43:16  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: JA 
      INTEGER , INTENT(IN) :: JB 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IH, JC, I, IC, NC, I2H, I2H1, JD, I2 
!-----------------------------------------------
!
!      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
!     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
!      COMMON/NDIMS/QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      DO IH = 1, IHSH 
         JC = JA 
         IF (IH == 2) THEN 
            DO I = 1, 2 
               IC = ICG(IH,I) 
               NC = NCG(IH,I) 
               I2H = IHSH + IH - 1 
!
! --- FIRST CONSIDER THE L.H.S. (I=1) OF THE MATRIX ELEMENT. NC=1 MEANS
!     UNOCCUPIED, REPRESENTED BY A DUMMY SINGLET S SHELL, AND THE
!    ADDITIONAL SET OF COUPLING QUANTUM NUMBERS WILL BE THE SAME AS THE
!     LAST SET OF COUPLING QUANTUM NUMBERS ALREADY OBTAINED.
!     NC=2 MEANS OCCUPIED.  THEN ALL THE NEW QUANTUM NUMBERS (BOTH FOR
!     THE SHELL AND FOR THE COUPLING OF THIS SHELL TO THE RESULTANT OF
!     THE PREVIOUS ONES) ARE DEFINED IN THE CORRESPONDING J1QNRD ARRAY.
!     NOSH - THE NUMBER OF ELECTRONS IN THIS SHELL, IS DEFINED BY THE
!     APPROPRIATE ENTRY IN NELCSH .  THE R.H.S. IS THEN CONSIDERED
!     SIMILARLY (I=2)
!
               IF (NC == 1) THEN 
                  J1QN(IH,1,I) = 0 
                  J1QN(IH,2,I) = 1 
                  J1QN(IH,3,I) = 1 
                  J1QN(I2H,1,I) = 0 
                  J1QN(I2H,2,I) = J1QN(1,2,I) 
                  J1QN(I2H,3,I) = J1QN(1,3,I) 
               ELSE 
                  JD = J1QNRD(IC,JC) 
                  J1QN(IH,1,I) = MOD(JD,64) 
                  JD = JD/64 
                  J1QN(IH,2,I) = MOD(JD,64) 
                  J1QN(IH,3,I) = JD/64 
!
!    IS THIS THE FIRST OCCUPIED SHELL OF THIS CONFIGURATION, THOUGH NOT
!     THE FIRST OF THE OTHER CONFIGURATION.  IF SO, THE INTERMEDIATE
!     COUPLING FORMED HAS THE SAME  L,S  VALUES AS THIS OCCUPIED SHELL,
!     SINCE WE COUPLE THE SHELL TO A DUMMY SINGLET S.
!
                  IF (IC <= 1) THEN 
                     I2 = 1 
                  ELSE 
                     I2 = NOCCSH(JC) + IC - 1 
                  ENDIF 
                  JD = J1QNRD(I2,JC) 
                  IF (IC <= 1) THEN 
                     J1QN(I2H,1,I) = 0 
                  ELSE 
                     J1QN(I2H,1,I) = MOD(JD,64) 
                  ENDIF 
                  JD = JD/64 
                  J1QN(I2H,2,I) = MOD(JD,64) 
                  J1QN(I2H,3,I) = JD/64 
               ENDIF 
               JC = JB 
            END DO 
         ELSE 
            IF (IH > 2) THEN 
               DO I = 1, 2 
                  IC = ICG(IH,I) 
                  NC = NCG(IH,I) 
                  I2H = IHSH + IH - 1 
!
! --- FIRST CONSIDER THE L.H.S. (I=1) OF THE MATRIX ELEMENT. NC=1 MEANS
!     UNOCCUPIED, REPRESENTED BY A DUMMY SINGLET S SHELL, AND THE
!    ADDITIONAL SET OF COUPLING QUANTUM NUMBERS WILL BE THE SAME AS THE
!     LAST SET OF COUPLING QUANTUM NUMBERS ALREADY OBTAINED.
!     NC=2 MEANS OCCUPIED.  THEN ALL THE NEW QUANTUM NUMBERS (BOTH FOR
!     THE SHELL AND FOR THE COUPLING OF THIS SHELL TO THE RESULTANT OF
!     THE PREVIOUS ONES) ARE DEFINED IN THE CORRESPONDING J1QNRD ARRAY.
!     NOSH - THE NUMBER OF ELECTRONS IN THIS SHELL, IS DEFINED BY THE
!     APPROPRIATE ENTRY IN NELCSH .  THE R.H.S. IS THEN CONSIDERED
!     SIMILARLY (I=2)
!
                  IF (NC == 1) THEN 
                     J1QN(IH,1,I) = 0 
                     J1QN(IH,2,I) = 1 
                     J1QN(IH,3,I) = 1 
                     I2H1 = I2H - 1 
                     J1QN(I2H,1,I) = J1QN(I2H1,1,I) 
                     J1QN(I2H,2,I) = J1QN(I2H1,2,I) 
                     J1QN(I2H,3,I) = J1QN(I2H1,3,I) 
                  ELSE 
                     JD = J1QNRD(IC,JC) 
                     J1QN(IH,1,I) = MOD(JD,64) 
                     JD = JD/64 
                     J1QN(IH,2,I) = MOD(JD,64) 
                     J1QN(IH,3,I) = JD/64 
!
!    IS THIS THE FIRST OCCUPIED SHELL OF THIS CONFIGURATION, THOUGH NOT
!     THE FIRST OF THE OTHER CONFIGURATION.  IF SO, THE INTERMEDIATE
!     COUPLING FORMED HAS THE SAME  L,S  VALUES AS THIS OCCUPIED SHELL,
!     SINCE WE COUPLE THE SHELL TO A DUMMY SINGLET S.
!
                     IF (IC <= 1) THEN 
                        I2 = 1 
                     ELSE 
                        I2 = NOCCSH(JC) + IC - 1 
                     ENDIF 
                     JD = J1QNRD(I2,JC) 
                     IF (IC <= 1) THEN 
                        J1QN(I2H,1,I) = 0 
                     ELSE 
                        J1QN(I2H,1,I) = MOD(JD,64) 
                     ENDIF 
                     JD = JD/64 
                     J1QN(I2H,2,I) = MOD(JD,64) 
                     J1QN(I2H,3,I) = JD/64 
                  ENDIF 
                  JC = JB 
               END DO 
            ELSE 
               DO I = 1, 2 
                  IC = ICG(IH,I) 
                  NC = NCG(IH,I) 
                  I2H = IHSH + IH - 1 
!
! --- FIRST CONSIDER THE L.H.S. (I=1) OF THE MATRIX ELEMENT. NC=1 MEANS
!     UNOCCUPIED, REPRESENTED BY A DUMMY SINGLET S SHELL, AND THE
!    ADDITIONAL SET OF COUPLING QUANTUM NUMBERS WILL BE THE SAME AS THE
!     LAST SET OF COUPLING QUANTUM NUMBERS ALREADY OBTAINED.
!     NC=2 MEANS OCCUPIED.  THEN ALL THE NEW QUANTUM NUMBERS (BOTH FOR
!     THE SHELL AND FOR THE COUPLING OF THIS SHELL TO THE RESULTANT OF
!     THE PREVIOUS ONES) ARE DEFINED IN THE CORRESPONDING J1QNRD ARRAY.
!     NOSH - THE NUMBER OF ELECTRONS IN THIS SHELL, IS DEFINED BY THE
!     APPROPRIATE ENTRY IN NELCSH .  THE R.H.S. IS THEN CONSIDERED
!     SIMILARLY (I=2)
!
                  IF (NC == 1) THEN 
                     J1QN(IH,1,I) = 0 
                     J1QN(IH,2,I) = 1 
                     J1QN(IH,3,I) = 1 
                  ELSE 
                     JD = J1QNRD(IC,JC) 
                     J1QN(IH,1,I) = MOD(JD,64) 
                     JD = JD/64 
                     J1QN(IH,2,I) = MOD(JD,64) 
                     J1QN(IH,3,I) = JD/64 
!
!     IS THIS THE FIRST OCCUPIED SHELL OF EITHER CONFIGURATION. IF SO,
!    THEN THERE ARE NO INTERMEDIATE COUPLINGS TO CONSIDER AT THIS STAGE
!
                     IF (IH > 1) THEN 
!
!    IS THIS THE FIRST OCCUPIED SHELL OF THIS CONFIGURATION, THOUGH NOT
!     THE FIRST OF THE OTHER CONFIGURATION.  IF SO, THE INTERMEDIATE
!     COUPLING FORMED HAS THE SAME  L,S  VALUES AS THIS OCCUPIED SHELL,
!     SINCE WE COUPLE THE SHELL TO A DUMMY SINGLET S.
!
                        IF (IC <= 1) THEN 
                           I2 = 1 
                        ELSE 
                           I2 = NOCCSH(JC) + IC - 1 
                        ENDIF 
                        JD = J1QNRD(I2,JC) 
                        IF (IC <= 1) THEN 
                           J1QN(I2H,1,I) = 0 
                        ELSE 
                           J1QN(I2H,1,I) = MOD(JD,64) 
                        ENDIF 
                        JD = JD/64 
                        J1QN(I2H,2,I) = MOD(JD,64) 
                        J1QN(I2H,3,I) = JD/64 
                     ENDIF 
                  ENDIF 
                  JC = JB 
               END DO 
            ENDIF 
         ENDIF 
      END DO 
      RETURN  
      END SUBROUTINE COUPLING 
!
