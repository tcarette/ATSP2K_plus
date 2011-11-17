!     ..........................................................   :
!                                                                  :
!          Block                                                   :
!                  Standard Quantities  -  S Q L S F               :
!                         Part One                                 :
!                                                                  :
!                                       Written by  G. Gaigalas,   :
!                Institute of Theoretical Physics and Astronomy    :
!                Vilnius,  Lithuania                               :
!                                                  December 1993   :
!                                                                  :
!                Laboratoire de Chimie Physique Moleculaire        :
!                Universite Libre de Bruxelles                     :
!                                                  December 1995   :
!                                                                  :
!                Vanderbilt University,  Nashville,  U S A         :
!                                                   October 1996   :
!                                                                  :
!     ..........................................................   :
!
!
!     ---------------------------------------------------------------
!     A 1
!     ---------------------------------------------------------------
!
!     Written by G. Gaigalas,                                      *
!     Vilnius,  Lithuania                          December 1993   *
!
      SUBROUTINE A1(IK, BK, ID, BD, QM1, A) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE CONSTS_C 
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:33:24  11/16/01  
!...Switches:                     
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(IN) :: QM1 
      REAL(DOUBLE) , INTENT(OUT) :: A 
      INTEGER , INTENT(IN) :: IK(7) 
      INTEGER , INTENT(IN) :: ID(7) 
      REAL(DOUBLE)  :: BK(3) 
      REAL(DOUBLE)  :: BD(3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ISUMA, IFAZ 
      REAL(DOUBLE) :: AB 
!-----------------------------------------------
      A = ZERO 
      IF (QM1 < EPS) THEN 
         ISUMA = (ID(5)+1)*(ID(6)+1)*ID(4) 
         AB = DBLE(ISUMA) 
         A = DSQRT(AB) 
         IFAZ = ID(5) + ID(6) + ID(3)*2 + 1 - IK(5) - IK(6) + ID(4)*2 
         IF ((IFAZ/4)*4 /= IFAZ) A = -A 
      ELSE 
         ISUMA = (IK(5)+1)*(IK(6)+1)*IK(4) 
         AB = DBLE(ISUMA) 
         A = DSQRT(AB) 
         IF ((IK(4)/2)*2 /= IK(4)) A = -A 
      ENDIF 
      RETURN  
      END SUBROUTINE A1 
