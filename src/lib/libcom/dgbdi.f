      SUBROUTINE DGBDI(ABD,LDA,N,ML,MU,IPVT,DET)                        00000010
      INTEGER LDA,N,ML,MU,IPVT(1)                                       00000020
      DOUBLE PRECISION ABD(LDA,1),DET(2)                                00000030
C                                                                       00000040
C     DGBDI COMPUTES THE DETERMINANT OF A BAND MATRIX                   00000050
C     USING THE FACTORS COMPUTED BY DGBCO OR DGBFA.                     00000060
C     IF THE INVERSE IS NEEDED, USE DGBSL  N  TIMES.                    00000070
C                                                                       00000080
C     ON ENTRY                                                          00000090
C                                                                       00000100
C        ABD     DOUBLE PRECISION(LDA, N)                               00000110
C                THE OUTPUT FROM DGBCO OR DGBFA.                        00000120
C                                                                       00000130
C        LDA     INTEGER                                                00000140
C                THE LEADING DIMENSION OF THE ARRAY  ABD .              00000150
C                                                                       00000160
C        N       INTEGER                                                00000170
C                THE ORDER OF THE ORIGINAL MATRIX.                      00000180
C                                                                       00000190
C        ML      INTEGER                                                00000200
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.           00000210
C                                                                       00000220
C        MU      INTEGER                                                00000230
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.           00000240
C                                                                       00000250
C        IPVT    INTEGER(N)                                             00000260
C                THE PIVOT VECTOR FROM DGBCO OR DGBFA.                  00000270
C                                                                       00000280
C     ON RETURN                                                         00000290
C                                                                       00000300
C        DET     DOUBLE PRECISION(2)                                    00000310
C                DETERMINANT OF ORIGINAL MATRIX.                        00000320
C                DETERMINANT = DET(1) * 10.0**DET(2)                    00000330
C                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0                  00000340
C                OR  DET(1) = 0.0 .                                     00000350
C                                                                       00000360
C     LINPACK. THIS VERSION DATED 08/14/78 .                            00000370
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      00000380
C                                                                       00000390
C     SUBROUTINES AND FUNCTIONS                                         00000400
C                                                                       00000410
C     FORTRAN DABS                                                      00000420
C                                                                       00000430
C     INTERNAL VARIABLES                                                00000440
C                                                                       00000450
      DOUBLE PRECISION TEN                                              00000460
      INTEGER I,M                                                       00000470
C                                                                       00000480
C                                                                       00000490
      M = ML + MU + 1                                                   00000500
      DET(1) = 1.0D0                                                    00000510
      DET(2) = 0.0D0                                                    00000520
      TEN = 10.0D0                                                      00000530
      DO 50 I = 1, N                                                    00000540
         IF (IPVT(I) .NE. I) DET(1) = -DET(1)                           00000550
         DET(1) = ABD(M,I)*DET(1)                                       00000560
C     ...EXIT                                                           00000570
         IF (DET(1) .EQ. 0.0D0) GO TO 60                                00000580
   10    IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20                          00000590
            DET(1) = TEN*DET(1)                                         00000600
            DET(2) = DET(2) - 1.0D0                                     00000610
         GO TO 10                                                       00000620
   20    CONTINUE                                                       00000630
   30    IF (DABS(DET(1)) .LT. TEN) GO TO 40                            00000640
            DET(1) = DET(1)/TEN                                         00000650
            DET(2) = DET(2) + 1.0D0                                     00000660
         GO TO 30                                                       00000670
   40    CONTINUE                                                       00000680
   50 CONTINUE                                                          00000690
   60 CONTINUE                                                          00000700
      RETURN                                                            00000710
      END                                                               00000720
