      SUBROUTINE DPBDI(ABD,LDA,N,M,DET)                                 00000010
      INTEGER LDA,N,M                                                   00000020
      DOUBLE PRECISION ABD(LDA,1)                                       00000030
      DOUBLE PRECISION DET(2)                                           00000040
C                                                                       00000050
C     DPBDI COMPUTES THE DETERMINANT                                    00000060
C     OF A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE BAND MATRIX     00000070
C     USING THE FACTORS COMPUTED BY DPBCO OR DPBFA.                     00000080
C     IF THE INVERSE IS NEEDED, USE DPBSL  N  TIMES.                    00000090
C                                                                       00000100
C     ON ENTRY                                                          00000110
C                                                                       00000120
C        ABD     DOUBLE PRECISION(LDA, N)                               00000130
C                THE OUTPUT FROM DPBCO OR DPBFA.                        00000140
C                                                                       00000150
C        LDA     INTEGER                                                00000160
C                THE LEADING DIMENSION OF THE ARRAY  ABD .              00000170
C                                                                       00000180
C        N       INTEGER                                                00000190
C                THE ORDER OF THE MATRIX  A .                           00000200
C                                                                       00000210
C        M       INTEGER                                                00000220
C                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.       00000230
C                                                                       00000240
C     ON RETURN                                                         00000250
C                                                                       00000260
C        DET     DOUBLE PRECISION(2)                                    00000270
C                DETERMINANT OF ORIGINAL MATRIX IN THE FORM             00000280
C                DETERMINANT = DET(1) * 10.0**DET(2)                    00000290
C                WITH  1.0 .LE. DET(1) .LT. 10.0                        00000300
C                OR  DET(1) .EQ. 0.0 .                                  00000310
C                                                                       00000320
C     LINPACK.  THIS VERSION DATED 08/14/78 .                           00000330
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      00000340
C                                                                       00000350
C     SUBROUTINES AND FUNCTIONS                                         00000360
C                                                                       00000370
C                                                                       00000380
C     INTERNAL VARIABLES                                                00000390
C                                                                       00000400
      DOUBLE PRECISION S                                                00000410
      INTEGER I                                                         00000420
C                                                                       00000430
C     COMPUTE DETERMINANT                                               00000440
C                                                                       00000450
      DET(1) = 1.0D0                                                    00000460
      DET(2) = 0.0D0                                                    00000470
      S = 10.0D0                                                        00000480
      DO 50 I = 1, N                                                    00000490
         DET(1) = ABD(M+1,I)**2*DET(1)                                  00000500
C     ...EXIT                                                           00000510
         IF (DET(1) .EQ. 0.0D0) GO TO 60                                00000520
   10    IF (DET(1) .GE. 1.0D0) GO TO 20                                00000530
            DET(1) = S*DET(1)                                           00000540
            DET(2) = DET(2) - 1.0D0                                     00000550
         GO TO 10                                                       00000560
   20    CONTINUE                                                       00000570
   30    IF (DET(1) .LT. S) GO TO 40                                    00000580
            DET(1) = DET(1)/S                                           00000590
            DET(2) = DET(2) + 1.0D0                                     00000600
         GO TO 30                                                       00000610
   40    CONTINUE                                                       00000620
   50 CONTINUE                                                          00000630
   60 CONTINUE                                                          00000640
      RETURN                                                            00000650
      END                                                               00000660
