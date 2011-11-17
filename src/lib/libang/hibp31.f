*
*     ------------------------------------------------------------------
*     H I B P 3 1
*     ------------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      SUBROUTINE HIBP31(I,BK,IBK,BD,IBD)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH1(16),NOSH2(16),J1QN1(31,3),
     :     J1QN2(31,3),IJFUL(16)
      DIMENSION BK(3),BD(3),IBK(7),IBD(7)
      IBK(2)=NJ(I)
      IBD(2)=IBK(2)
      IBK(3)=LJ(I)
      IBD(3)=IBK(3)
      IBK(4)=NOSH1(I)
      IBD(4)=NOSH2(I)
      IBK(5)=J1QN1(I,2)-1
      IBD(5)=J1QN2(I,2)-1
      IBK(6)=J1QN1(I,3)-1
      IBD(6)=J1QN2(I,3)-1
      IREZ=0
      IF(IBK(3).LT.3) THEN
        IBK(7)=2*LJ(I)+1-J1QN1(I,1)
      ELSEIF(IBK(4).LT.3) THEN
        IBK(7)=2*LJ(I)+1-J1QN1(I,1)
      ELSEIF(IBK(4).LT.14) THEN
        IREZ=1
        I2NK=J1QN1(I,1)
        IBK(1)=NUMTERF(I2NK,IBK(6),IBK(5),IBK(4),IBK(7))
      ELSE
        IBK(7)=2*LJ(I)+1-J1QN1(I,1)
      ENDIF
      BK(1)=HALF*DBLE(IBK(7))
      BK(2)=HALF*DBLE(IBK(6))
      BK(3)=-HALF*DBLE(2*LJ(I)+1-NOSH1(I))
      IF(IREZ.EQ.0) THEN
        IBK(1)=NUMTER(IBK(7),IBK(6),IBK(5),IBK(3),IBD(4),IBK(4))
      ENDIF
      IREZ=0
      IF(IBD(3).LT.3) THEN
        IBD(7)=2*LJ(I)+1-J1QN2(I,1)
      ELSEIF(IBD(4).LT.3) THEN
        IBD(7)=2*LJ(I)+1-J1QN2(I,1)
      ELSEIF(IBD(4).LT.14) THEN
        IREZ=1
        I2ND=J1QN2(I,1)
        IBD(1)=NUMTERF(I2ND,IBD(6),IBD(5),IBD(4),IBD(7))
      ELSE
        IBD(7)=2*LJ(I)+1-J1QN2(I,1)
      ENDIF
      BD(1)=HALF*DBLE(IBD(7))
      BD(2)=HALF*DBLE(IBD(6))
      BD(3)=-HALF*DBLE(2*LJ(I)+1-NOSH2(I))
      IF(IREZ.EQ.0) THEN
        IBD(1)=NUMTER(IBD(7),IBD(6),IBD(5),IBD(3),IBK(4),IBD(4))
      ENDIF
      RETURN
      END
