      SUBROUTINE DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)                       00000010
      INTEGER LDA,N,ML,MU,IPVT(1),INFO                                  00000020
      DOUBLE PRECISION ABD(LDA,1)                                       00000030
C                                                                       00000040
C     DGBFA FACTORS A DOUBLE PRECISION BAND MATRIX BY ELIMINATION.      00000050
C                                                                       00000060
C     DGBFA IS USUALLY CALLED BY DGBCO, BUT IT CAN BE CALLED            00000070
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          00000080
C                                                                       00000090
C     ON ENTRY                                                          00000100
C                                                                       00000110
C        ABD     DOUBLE PRECISION(LDA, N)                               00000120
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS      00000130
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND   00000140
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS         00000150
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .                       00000160
C                SEE THE COMMENTS BELOW FOR DETAILS.                    00000170
C                                                                       00000180
C        LDA     INTEGER                                                00000190
C                THE LEADING DIMENSION OF THE ARRAY  ABD .              00000200
C                LDA MUST BE .GE. 2*ML + MU + 1 .                       00000210
C                                                                       00000220
C        N       INTEGER                                                00000230
C                THE ORDER OF THE ORIGINAL MATRIX.                      00000240
C                                                                       00000250
C        ML      INTEGER                                                00000260
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.           00000270
C                0 .LE. ML .LT. N .                                     00000280
C                                                                       00000290
C        MU      INTEGER                                                00000300
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.           00000310
C                0 .LE. MU .LT. N .                                     00000320
C                MORE EFFICIENT IF  ML .LE. MU .                        00000330
C     ON RETURN                                                         00000340
C                                                                       00000350
C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND         00000360
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.          00000370
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       00000380
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          00000390
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       00000400
C                                                                       00000410
C        IPVT    INTEGER(N)                                             00000420
C                AN INTEGER VECTOR OF PIVOT INDICES.                    00000430
C                                                                       00000440
C        INFO    INTEGER                                                00000450
C                = 0  NORMAL VALUE.                                     00000460
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR       00000470
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES        00000480
C                     INDICATE THAT DGBSL WILL DIVIDE BY ZERO IF        00000490
C                     CALLED.  USE  RCOND  IN DGBCO FOR A RELIABLE      00000500
C                     INDICATION OF SINGULARITY.                        00000510
C                                                                       00000520
C     BAND STORAGE                                                      00000530
C                                                                       00000540
C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT      00000550
C           WILL SET UP THE INPUT.                                      00000560
C                                                                       00000570
C                   ML = (BAND WIDTH BELOW THE DIAGONAL)                00000580
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)                00000590
C                   M = ML + MU + 1                                     00000600
C                   DO 20 J = 1, N                                      00000610
C                      I1 = MAX0(1, J-MU)                               00000620
C                      I2 = MIN0(N, J+ML)                               00000630
C                      DO 10 I = I1, I2                                 00000640
C                         K = I - J + M                                 00000650
C                         ABD(K,J) = A(I,J)                             00000660
C                10    CONTINUE                                         00000670
C                20 CONTINUE                                            00000680
C                                                                       00000690
C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .         00000700
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR      00000710
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.            00000720
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .    00000730
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE            00000740
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.          00000750
C                                                                       00000760
C     LINPACK. THIS VERSION DATED 08/14/78 .                            00000770
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      00000780
C                                                                       00000790
C     SUBROUTINES AND FUNCTIONS                                         00000800
C                                                                       00000810
C     BLAS DAXPY,DSCAL,IDAMAX                                           00000820
C     FORTRAN MAX0,MIN0                                                 00000830
C                                                                       00000840
C     INTERNAL VARIABLES                                                00000850
C                                                                       00000860
      DOUBLE PRECISION T                                                00000870
      INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1             00000880
C                                                                       00000890
C                                                                       00000900
      M = ML + MU + 1                                                   00000910
      INFO = 0                                                          00000920
C                                                                       00000930
C     ZERO INITIAL FILL-IN COLUMNS                                      00000940
C                                                                       00000950
      J0 = MU + 2                                                       00000960
      J1 = MIN0(N,M) - 1                                                00000970
      IF (J1 .LT. J0) GO TO 30                                          00000980
      DO 20 JZ = J0, J1                                                 00000990
         I0 = M + 1 - JZ                                                00001000
         DO 10 I = I0, ML                                               00001010
            ABD(I,JZ) = 0.0D0                                           00001020
   10    CONTINUE                                                       00001030
   20 CONTINUE                                                          00001040
   30 CONTINUE                                                          00001050
      JZ = J1                                                           00001060
      JU = 0                                                            00001070
C                                                                       00001080
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                        00001090
C                                                                       00001100
      NM1 = N - 1                                                       00001110
      IF (NM1 .LT. 1) GO TO 130                                         00001120
      DO 120 K = 1, NM1                                                 00001130
         KP1 = K + 1                                                    00001140
C                                                                       00001150
C        ZERO NEXT FILL-IN COLUMN                                       00001160
C                                                                       00001170
         JZ = JZ + 1                                                    00001180
         IF (JZ .GT. N) GO TO 50                                        00001190
         IF (ML .LT. 1) GO TO 50                                        00001200
            DO 40 I = 1, ML                                             00001210
               ABD(I,JZ) = 0.0D0                                        00001220
   40       CONTINUE                                                    00001230
   50    CONTINUE                                                       00001240
C                                                                       00001250
C        FIND L = PIVOT INDEX                                           00001260
C                                                                       00001270
         LM = MIN0(ML,N-K)                                              00001280
         L = IDAMAX(LM+1,ABD(M,K),1) + M - 1                            00001290
         IPVT(K) = L + K - M                                            00001300
C                                                                       00001310
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED          00001320
C                                                                       00001330
         IF (ABD(L,K) .EQ. 0.0D0) GO TO 100                             00001340
C                                                                       00001350
C           INTERCHANGE IF NECESSARY                                    00001360
C                                                                       00001370
            IF (L .EQ. M) GO TO 60                                      00001380
               T = ABD(L,K)                                             00001390
               ABD(L,K) = ABD(M,K)                                      00001400
               ABD(M,K) = T                                             00001410
   60       CONTINUE                                                    00001420
C                                                                       00001430
C           COMPUTE MULTIPLIERS                                         00001440
C                                                                       00001450
            T = -1.0D0/ABD(M,K)                                         00001460
            CALL DSCAL(LM,T,ABD(M+1,K),1)                               00001470
C                                                                       00001480
C           ROW ELIMINATION WITH COLUMN INDEXING                        00001490
C                                                                       00001500
            JU = MIN0(MAX0(JU,MU+IPVT(K)),N)                            00001510
            MM = M                                                      00001520
            IF (JU .LT. KP1) GO TO 90                                   00001530
            DO 80 J = KP1, JU                                           00001540
               L = L - 1                                                00001550
               MM = MM - 1                                              00001560
               T = ABD(L,J)                                             00001570
               IF (L .EQ. MM) GO TO 70                                  00001580
                  ABD(L,J) = ABD(MM,J)                                  00001590
                  ABD(MM,J) = T                                         00001600
   70          CONTINUE                                                 00001610
               CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)              00001620
   80       CONTINUE                                                    00001630
   90       CONTINUE                                                    00001640
         GO TO 110                                                      00001650
  100    CONTINUE                                                       00001660
            INFO = K                                                    00001670
  110    CONTINUE                                                       00001680
  120 CONTINUE                                                          00001690
  130 CONTINUE                                                          00001700
      IPVT(N) = N                                                       00001710
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N                                 00001720
      RETURN                                                            00001730
      END                                                               00001740
