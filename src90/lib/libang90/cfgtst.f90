!
!     -------------------------------------------------------------
!       C F G T S T
!     -------------------------------------------------------------
!
!     Modified by Gediminas Gaigalas,                September 1997
!
!
!      SUBROUTINE CFGTST(NCFG,QLJCOMP,QNOC,QNELCSH,QNOCORB,QJ1)
      SUBROUTINE CFGTST 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE INFORM_C 
      USE TERMS_C 
      USE MT15_C 
      USE MT67_C 
      use non30_C
      use ndims_C
      use kron_C
      
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  09:19:51  11/16/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
       
      USE jthn_I 
      USE ntab1_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!      INTEGER, PARAMETER :: NCWD = 20 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(14) :: IGMAX 
      INTEGER :: IALLOW, I, NELSUM, N, J, NA, NC, LQU, JD, JA, JB, JC, LQUMAX, &
         IMAXGG, JCB, IA, LASTEL, LS, NR, NROW, I1, I2, I3, J2, J1, JE, JG, JF&
         , NELCS 
!-----------------------------------------------
!      INTEGER, PARAMETER :: NCWD = 20
!
!
!     THIS SUBROUTINE CHECKS ALL THE CONFIGURATION SET TO ENSURE THAT
!     IT SATISFIES ALL THE FOLLOWING CONDITIONS:
!        (1)  EACH CONFIGURATION HAS THE SAME NUMBER OF ELECTRONS
!        (2)  NO SUBSHELL HAS TOO MANY (.GT.2*(2*L+1))  ELECTRONS
!        (3)  THE ELECTRONS IN ANY ONE SUBSHELL ARE COUPLED TO FORM AN
!             ALLOWED TRIAD OF QUANTUM NUMBERS
!        (4)  THE TRIADS COUPLE TOGETHER IN AN ALLOWED WAY
!
!     IN THE EVENT OF AN ERROR, THE PROGRAM HALTS AT THE COMPLETION
!     OF THE CHECKING.  ANY NUMBER OF S, P, D  ELECTRONS ARE ALLOWED,
!     (BUT .LE.2*(2*L+1)), BUT ONLY UP TO TWO ELECTRONS, L >=3.
!     WHEN L>4, THE ONLY ALLOWED TERMS ARE THOSE FOR L=4.
!     A FILLED F-SHELL IS ALSO ALLOWED AS WELL AS A SINGLE ELECTRON
!     WITH L.GT.4
!
!
!     IMPLICIT INTEGER (Q)
!      POINTER(QLJCOMP,LJCOMP(1))
!      POINTER(QNOC,NOCCSH(1))
!      POINTER(QNELCSH,NELCSH(8,1))
!      POINTER(QNOCORB,NOCORB(8,1))
!      POINTER(QJ1,J1QNRD(15,1))
!
!GGf
      DATA IGMAX/ 0, 0, 17, 47, 73, 119, 119, 119, 73, 47, 17, 7, 1, 0/  
!GGf
!
    5 FORMAT(/,' THE TRIAD OF QUANTUM NUMBERS OF SHELL',I3,' IN CONFIGURATION',&
         I8,' IS NOT A RECOGNIZED SET') 
    7 FORMAT(/,' THE COUPLING OF SHELL',I3,' IN CONFIGURATION',I3,&
         ' RESULTS IN AN ILLEGAL COUPLING SCHEME') 
   12 FORMAT(/,/,' CONFIGURATION DATA WRONG, PROGRAM HALTED'/,/) 
   15 FORMAT(/,' IN CONFIGURATION',I3,', SHELL',I3,&
         ' CONTAINS TOO MANY ELECTRONS') 
   17 FORMAT(/,' CONFIGURATION',I3,&
      ' INCLUDES A SHELL OF ANGULAR MOMENTUM L.GE.3 WITH TOO MANY ELECTRONS') 
   18 FORMAT(/,' CONFIGURATION',I3,' HAS AN INCORRECT NUMBER OF ','ELECTRONS') 
      IALLOW = 1 
      DO I = 1, NCFG 
         NELSUM = 0 
         N = NOCCSH(I) 
         DO J = 1, N 
            NA = NOCORB(J,I) 
            NC = NELCSH(J,I) 
            LQU = LJCOMP(NA) 
            NELSUM = NELSUM + NC 
            JD = J1QNRD(J,I) 
            JA = MOD(JD,64) 
            JD = JD/64 
            JB = MOD(JD,64) 
            JC = JD/64 
            LQUMAX = 4*LQU + 2 
            IF (NC > LQUMAX) THEN 
               WRITE (IWRITE, 15) I, J 
               IALLOW = 0 
               CYCLE  
!GGf
            ELSE IF (LQU==3 .AND. NC>2 .AND. NC<14) THEN 
               IMAXGG = IGMAX(NC) 
               JCB = (JC - 1)*100 + JB - 1 
               DO IA = 1, IMAXGG 
                  SELECT CASE (NC)  
                  CASE (3)  
                     LASTEL = M3(IA) 
                  CASE (4)  
                     LASTEL = M4(IA) 
                  CASE (5)  
                     LASTEL = M5(IA) 
                  CASE (6)  
                     LASTEL = M6(IA) 
                  CASE (7)  
                     LASTEL = M7(IA) 
                  CASE (8)  
                     LASTEL = M6(IA) 
                  CASE (9)  
                     LASTEL = M5(IA) 
                  CASE (10)  
                     LASTEL = M4(IA) 
                  CASE (11)  
                     LASTEL = M3(IA) 
                  CASE (12)  
                     LASTEL = M2(IA) 
                  CASE (13)  
                     LASTEL = M1(IA) 
                  CASE DEFAULT 
                     STOP  
                  END SELECT 
                  LS = JTHN(LASTEL,1,10000) 
                  IF (JCB /= LS) CYCLE  
                  NR = JTHN(LASTEL,4,100) 
                  IF (JA == NR) GO TO 21 
               END DO 
               GO TO 24 
            ELSE IF (LQU>4 .AND. NC>2) THEN 
               WRITE (IWRITE, 17) I 
               IALLOW = 0 
               CYCLE  
!GGf
            ELSE IF (NC == 1) THEN 
               IF (JA==1 .AND. JB==2*LQU+1 .AND. JC==2) GO TO 21 
            ELSE 
               IF (LQU>4 .AND. NC==2) LQU = 4 
               IF (NC == LQUMAX) THEN 
                  NROW = 2 
!GGf
               ELSE IF (NC == 0) THEN 
                  NROW = 2 
!GGf
               ELSE 
                  NROW = NTAB1(NC + 1,LQU + 1) 
               ENDIF 
               I1 = ITAB(NROW) 
               I2 = JTAB(NROW) 
               DO IA = 1, I1 
                  I3 = I2 + 3*IA - 1 
                  IF (JB /= NTAB(I3)) CYCLE  
                  I3 = I3 + 1 
                  IF (JC /= NTAB(I3)) CYCLE  
                  I3 = I3 - 2 
                  IF (JA == NTAB(I3)) GO TO 21 
               END DO 
            ENDIF 
!GGf
   24       CONTINUE 
            IALLOW = 0 
!GGf
            WRITE (IWRITE, 5) J, I 
            CYCLE  
!
!     CHECK ON THE COUPLING OF THE TRIADS
!
   21       CONTINUE 
            IF (N<=1 .OR. J<=1) CYCLE  
            J2 = N + J - 1 
            J1 = J2 - 1 
            IF (J == 2) J1 = 1 
            JE = J1QNRD(J1,I)/64 
            JD = MOD(JE,64) 
            JE = JE/64 
            JG = J1QNRD(J2,I)/64 
            JF = MOD(JG,64) 
            JG = JG/64 
            IF (.NOT.(JF>=JB + JD .OR. JF<=IABS(JB-JD) .OR. JG>=JC+JE .OR. JG<=&
               IABS(JC-JE) .OR. MOD(JC+JE-JG,2)==0)) CYCLE  
            WRITE (IWRITE, 7) J, I 
            IALLOW = 0 
         END DO 
         IF (I == 1) THEN 
            NELCS = NELSUM 
         ELSE IF (NELSUM /= NELCS) THEN 
            WRITE (IWRITE, 18) I 
            IALLOW = 0 
         ENDIF 
      END DO 
      IF (IALLOW == 0) THEN 
         WRITE (IWRITE, 12) 
         STOP  
      ENDIF 
      RETURN  
      END SUBROUTINE CFGTST 
