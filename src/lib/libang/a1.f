*     ..........................................................   :
*                                                                  :
*          Block                                                   :
*                  Standard Quantities  -  S Q L S F               :
*                         Part One                                 :
*                                                                  :
*                                       Written by  G. Gaigalas,   :
*                Institute of Theoretical Physics and Astronomy    :
*                Vilnius,  Lithuania                               :
*                                                  December 1993   :
*                                                                  :
*                Laboratoire de Chimie Physique Moleculaire        :
*                Universite Libre de Bruxelles                     :
*                                                  December 1995   :
*                                                                  :
*                Vanderbilt University,  Nashville,  U S A         :
*                                                   October 1996   :
*                                                                  :
*     ..........................................................   :
*
*
*     ---------------------------------------------------------------
*     A 1
*     ---------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                          December 1993   *
*
      SUBROUTINE A1(IK,BK,ID,BD,QM1,A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IK(7),BK(3),ID(7),BD(3)
      A=ZERO
      IF(QM1.LT.EPS) THEN
        ISUMA=(ID(5)+1)*(ID(6)+1)*ID(4)
        AB=DBLE(ISUMA)
        A=DSQRT(AB)
        IFAZ=ID(5)+ID(6)+ID(3)*2+1-IK(5)-IK(6)+ID(4)*2
        IF((IFAZ/4)*4.NE.IFAZ)A=-A
      ELSE
        ISUMA=(IK(5)+1)*(IK(6)+1)*IK(4)
        AB=DBLE(ISUMA)
        A=DSQRT(AB)
        IF((IK(4)/2)*2.NE.IK(4))A=-A
      ENDIF
      RETURN
      END
