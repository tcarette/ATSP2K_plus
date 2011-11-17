!
!     ------------------------------------------------------------
!     I N I T M 2
!     ------------------------------------------------------------
!
!     Determine the Rydberg constant, RMASS constant, and
!     initialize the radial grid
!
      SUBROUTINE INITM2(IPROB) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE INOUT_C 
      USE PARAM_C 
      USE RADIAL_C 
      USE RYDBERG_C 
      USE COR_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:58:56  11/18/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IPROB 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NWD = 128 
      INTEGER, PARAMETER :: NELEMENT = 50 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J 
      REAL , DIMENSION(NELEMENT) :: MN 
      REAL, DIMENSION(NELEMENT - 9) :: FZ 
      REAL(DOUBLE) :: ME, RY 
      CHARACTER :: YES 
!-----------------------------------------------
!
!
      DATA ME/ 548.579903D-6/  
      DATA RY/ 109737.31534/  
      DATA MN/ 1, 4, 7, 9, 11, 12, 14, 16, 19, 20, 23, 24, 27, 28, 31, 32, 35, &
         40, 39, 40, 45, 48, 51, 52, 55, 56, 59, 58, 63, 64, 69, 74, 75, 80, 79&
         , 84, 85, 88, 89, 90, 93, 98, 99, 102, 103, 108, 107, 114, 115, 120/  
!
      DATA FZ/ 163.5, 198.4, 239.6, 283.9, 332.7, 386.1, 444.2, 507.8, 576.4, &
         650.7, 731.2, 817.9, 911.2, 1011.4, 1119.3, 1234.9, 1359.3, 1491.9, &
         1634.4, 1786.6, 1948.7, 2122.7, 2308.2, 2507.4, 2718.5, 2944.3, 3187.5&
         , 3442.2, 3720.3, 4015.4, 4326.6, 4664.8, 5022.4, 5405.9, 5811.9, &
         6245.0, 6706.8, 7201.3, 7728.7, 8296.5, 8889.3/  
!
 
!
! ****   DETERMINE THE RYDBERG CONSTANT AND F(Z) PARAMETER
!
      WRITE (ISCW, '(/1X,A)') ' Default Rydberg constant (y/n)' 
!      READ(iread,'(A1)') YES
      YES = 'y' 
      IF (YES=='Y' .OR. YES=='y' .AND. INT(Z)<=NELEMENT) THEN 
         ZMU = MN(INT(Z)) 
      ELSE 
         WRITE (ISCW, '(1X,A)') 'Enter the mass of the atom' 
!         READ(iread,*) ZMU
         ZMU = 200 
      ENDIF 
      RYDBRG = RY/(1. + ME/ZMU) 
      WRITE (6, *) ' Rydberg constant used (cm-1) = ', RYDBRG 
      WRITE (6, *) 
      IF (IPROB /= 1) THEN 
         IF (INT(Z) < 10) THEN 
            FZE = 1.566432 
         ELSE IF (INT(Z) <= NELEMENT) THEN 
            FZE = FZ(INT(Z-9)) 
         ELSE 
            WRITE (ISCW, '(1X,A)') 'Enter the f(Z) value in MHz' 
            READ (IREAD, *) FZE 
         ENDIF 
      ENDIF 
      RMASS = ME/ZMU 
!
! ****  INITIALIZE THE RADIAL GRID
!
      RHO = -4.0 
      DO J = 1, NO 
         R(J) = DEXP(RHO)/Z 
         RR(J) = R(J)*R(J) 
         R2(J) = DSQRT(R(J)) 
         RHO = RHO + H 
      END DO 
      RETURN  
      END SUBROUTINE INITM2 
