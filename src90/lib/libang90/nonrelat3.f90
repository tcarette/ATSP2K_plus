!
!     -------------------------------------------------------------
!      N O N R E L A T 3
!     -------------------------------------------------------------
!                                                                  *
!     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
!     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
!                                              N'2 = N2 + 1        *
!                                                                  *
!     Written by G. Gaigalas,                                      *
!     Vanderbilt University,  Nashville            February 1994   *
!
      SUBROUTINE NONRELAT3(IA, IB, IC, ID, IIRE) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE OPERAT_C 
      USE MEDEFN_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:17:07  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE nonrelat31_I 
      USE ittk_I 
      USE nonrelat32_I 
      USE nonrelat33_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IA 
      INTEGER  :: IB 
      INTEGER  :: IC 
      INTEGER  :: ID 
      INTEGER  :: IIRE 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      IF (IHSH <= 1) RETURN  
      IF (IB == ID) THEN 
         IF (IA==IB .OR. IC==IB) THEN 
            IF (IA == IC) GO TO 10 
            IF (IC == IB) THEN 
               CALL NONRELAT31 (IC, IA, IA, IB, IC, ID, IIRE) 
            ELSE 
               IF (IIRE == 0) RETURN  
               IF (ISOTOP == 1) THEN 
                  IF (ITTK(LJ(IA),LJ(IA),1) == 0) RETURN  
                  IF (ITTK(LJ(IC),LJ(IA),1) == 0) RETURN  
               ENDIF 
               CALL NONRELAT32 (IC, IA, IA, IB, IC, ID) 
            ENDIF 
         ELSE 
            IF (IIRE == 0) RETURN  
            CALL NONRELAT33 (IC, IA, IB, 1, IA, IB, IC, ID) 
         ENDIF 
      ELSE IF (IA == IC) THEN 
         IF (IB==IA .OR. ID==IA) THEN 
            IF (IB == ID) GO TO 10 
            IF (ID == IA) THEN 
               CALL NONRELAT31 (ID, IB, IA, IB, IC, ID, IIRE) 
            ELSE 
               IF (IIRE == 0) RETURN  
               IF (ISOTOP == 1) THEN 
                  IF (ITTK(LJ(IB),LJ(IB),1) == 0) RETURN  
                  IF (ITTK(LJ(ID),LJ(IB),1) == 0) RETURN  
               ENDIF 
               CALL NONRELAT32 (ID, IB, IA, IB, IC, ID) 
            ENDIF 
         ELSE 
            IF (IIRE == 0) RETURN  
            CALL NONRELAT33 (ID, IB, IA, 1, IA, IB, IC, ID) 
         ENDIF 
      ELSE IF (IA == ID) THEN 
         IF (IB==IA .OR. IC==IA) THEN 
            IF (IB == IC) GO TO 10 
            IF (IC == ID) THEN 
               CALL NONRELAT31 (IC, IB, IA, IB, IC, ID, IIRE) 
            ELSE 
               IF (IIRE == 0) RETURN  
               IF (ISOTOP == 1) THEN 
                  IF (ITTK(LJ(IB),LJ(IB),1) == 0) RETURN  
                  IF (ITTK(LJ(IC),LJ(IB),1) == 0) RETURN  
               ENDIF 
               CALL NONRELAT32 (IC, IB, IA, IB, IC, ID) 
            ENDIF 
         ELSE 
            IF (IIRE == 0) RETURN  
            CALL NONRELAT33 (IC, IB, IA, 2, IA, IB, ID, IC) 
         ENDIF 
      ELSE IF (IB == IC) THEN 
         IF (IA==IB .OR. ID==IB) THEN 
            IF (IA == ID) GO TO 10 
            IF (ID == IB) THEN 
               CALL NONRELAT31 (ID, IA, IA, IB, IC, ID, IIRE) 
            ELSE 
               IF (IIRE == 0) RETURN  
               IF (ISOTOP == 1) THEN 
                  IF (ITTK(LJ(IA),LJ(IA),1) == 0) RETURN  
                  IF (ITTK(LJ(ID),LJ(IA),1) == 0) RETURN  
               ENDIF 
               CALL NONRELAT32 (ID, IA, IA, IB, IC, ID) 
            ENDIF 
         ELSE 
            IF (IIRE == 0) RETURN  
            CALL NONRELAT33 (ID, IA, IB, 2, IA, IB, ID, IC) 
         ENDIF 
      ENDIF 
      RETURN  
   10 CONTINUE 
      WRITE (6, '(A)') ' ERRO IN NONRELAT3' 
      STOP  
      END SUBROUTINE NONRELAT3 
