*
*    ------------------------------------------------------------------
*             I N I T A
*    ------------------------------------------------------------------
*
*      Initializes basic constants of the program including those
*  which define the average energy of a configuration.
*
*
      SUBROUTINE INITA
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /EAV/CCA(10),CCB(35)
*
* *****  AVERAGE INTERACTIONS FOR EQUIVALENT ELECTRONS
*
* *****  P - P
*
      CCA(1) = 2.D0/25.D0
*
* *****  D - D
*
      CCA(2) = 2.D0/63.D0
      CCA(3) = 2.D0/63.D0
*
* *****  F - F
*
      CCA(4) =   4.D0/ 195.D0
      CCA(5) =   2.D0/ 143.D0
      CCA(6) = 100.D0/5577.D0
*
* *****  G - G
*
      CCA(7) =   20.D0/  1309.D0
      CCA(8) =  162.D0/ 17017.D0
      CCA(9) =   20.D0/  2431.D0
      CCA(10) = 4410.D0/371943.D0
*
*
* ***** AVERAGE INTERACTIONS FOR NON-EQUIVALENT ELECTRONS
*
* *****  S - ( S, P, D, F, G )
*
      CCB(1) = 1.D0/ 2.D0
      CCB(2) = 1.D0/ 6.D0
      CCB(3) = 1.D0/10.D0
      CCB(4) = 1.D0/14.D0
      CCB(5) = 1.D0/18.D0
*
* *****  P - ( P, D, F, G )
*
      CCB(6) = 1.D0/  6.D0
      CCB(7) = 1.D0/ 15.D0
      CCB(8) = 1.D0/ 15.D0
      CCB(9) = 3.D0/ 70.D0
      CCB(10) = 3.D0/ 70.D0
      CCB(11) = 2.D0/ 63.D0
      CCB(12) = 2.D0/ 63.D0
      CCB(13) = 5.D0/198.D0
*
* *****  D - ( D, F, G )
*
      CCB(14) =  1.D0/ 10.D0
      CCB(15) =  1.D0/ 35.D0
      CCB(16) =  1.D0/ 35.D0
      CCB(17) =  3.D0/ 70.D0
      CCB(18) =  2.D0/105.D0
      CCB(19) =  5.D0/231.D0
      CCB(20) =  1.D0/ 35.D0
      CCB(21) = 10.D0/693.D0
      CCB(22) =  5.D0/286.D0
*
* *****  F - ( F, G )
*
      CCB(23) =  1.D0/  14.D0
      CCB(24) =  2.D0/ 105.D0
      CCB(25) =  1.D0/  77.D0
      CCB(26) = 50.D0/3003.D0
      CCB(27) =  2.D0/  63.D0
      CCB(28) =  1.D0/  77.D0
      CCB(29) = 10.D0/1001.D0
      CCB(30) = 35.D0/2574.D0
*
* *****  G - ( G )
*
      CCB(31) =   1.D0/   18.D0
      CCB(32) =  10.D0/  693.D0
      CCB(33) =   9.D0/ 1001.D0
      CCB(34) =  10.D0/ 1287.D0
      CCB(35) = 245.D0/21879.D0
*
* *** Initialize /FACT/
*
      CALL FACTRL(32)
      RETURN
      END
