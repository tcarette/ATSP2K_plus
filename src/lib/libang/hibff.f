*
*     ------------------------------------------------------------------
*     H I B F F
*     ------------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      *
*     Vilnius,  Lithuania                          December 1993   *
*
      SUBROUTINE HIBFF(JA1,JA2,JA3,JA4,NSLUO)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/TRK/BDS1(3),BDS2(3),BKS1(3),BKS2(3),
     *IBDS1(7),IBDS2(7),IBKS1(7),IBKS2(7)
      COMMON/TRK2/BDS3(3),BDS4(3),BKS3(3),BKS4(3),
     *IBDS3(7),IBDS4(7),IBKS3(7),IBKS4(7)
      EXTERNAL GLCONS,RIBLS
      CALL HIBP31(JA1,BKS1,IBKS1,BDS1,IBDS1)
      IF(NSLUO.EQ.1)RETURN
      CALL HIBP31(JA2,BKS2,IBKS2,BDS2,IBDS2)
      IF(NSLUO.EQ.2)RETURN
      CALL HIBP31(JA3,BKS3,IBKS3,BDS3,IBDS3)
      IF(NSLUO.EQ.3)RETURN
      CALL HIBP31(JA4,BKS4,IBKS4,BDS4,IBDS4)
      RETURN
      END
