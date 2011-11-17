      SUBROUTINE DGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)                      00000010
      INTEGER LDA,N,ML,MU,IPVT(1),JOB                                   00000020
      DOUBLE PRECISION ABD(LDA,1),B(1)                                  00000030
C                                                                       00000040
C     DGBSL SOLVES THE DOUBLE PRECISION BAND SYSTEM                     00000050
C     A * X = B  OR  TRANS(A) * X = B                                   00000060
C     USING THE FACTORS COMPUTED BY DGBCO OR DGBFA.                     00000070
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
C        B       DOUBLE PRECISION(N)                                    00000290
C                THE RIGHT HAND SIDE VECTOR.                            00000300
C                                                                       00000310
C        JOB     INTEGER                                                00000320
C                = 0         TO SOLVE  A*X = B ,                        00000330
C                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE           00000340
C                            TRANS(A)  IS THE TRANSPOSE.                00000350
C                                                                       00000360
C     ON RETURN                                                         00000370
C                                                                       00000380
C        B       THE SOLUTION VECTOR  X .                               00000390
C                                                                       00000400
C     ERROR CONDITION                                                   00000410
C                                                                       00000420
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A   00000430
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY  00000440
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER       00000450
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE     00000460
C        CALLED CORRECTLY AND IF DGBCO HAS SET RCOND .GT. 0.0           00000470
C        OR DGBFA HAS SET INFO .EQ. 0 .                                 00000480
C                                                                       00000490
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 00000500
C     WITH  P  COLUMNS                                                  00000510
C           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)                    00000520
C           IF (RCOND IS TOO SMALL) GO TO ...                           00000530
C           DO 10 J = 1, P                                              00000540
C              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)                00000550
C        10 CONTINUE                                                    00000560
C                                                                       00000570
C     LINPACK. THIS VERSION DATED 08/14/78 .                            00000580
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      00000590
C                                                                       00000600
C     SUBROUTINES AND FUNCTIONS                                         00000610
C                                                                       00000620
C     BLAS DAXPY,DDOT                                                   00000630
C     FORTRAN MIN0                                                      00000640
C                                                                       00000650
C     INTERNAL VARIABLES                                                00000660
C                                                                       00000670
      DOUBLE PRECISION DDOT,T                                           00000680
      INTEGER K,KB,L,LA,LB,LM,M,NM1                                     00000690
C                                                                       00000700
      M = MU + ML + 1                                                   00000710
      NM1 = N - 1                                                       00000720
      IF (JOB .NE. 0) GO TO 50                                          00000730
C                                                                       00000740
C        JOB = 0 , SOLVE  A * X = B                                     00000750
C        FIRST SOLVE L*Y = B                                            00000760
C                                                                       00000770
         IF (ML .EQ. 0) GO TO 30                                        00000780
         IF (NM1 .LT. 1) GO TO 30                                       00000790
            DO 20 K = 1, NM1                                            00000800
               LM = MIN0(ML,N-K)                                        00000810
               L = IPVT(K)                                              00000820
               T = B(L)                                                 00000830
               IF (L .EQ. K) GO TO 10                                   00000840
                  B(L) = B(K)                                           00000850
                  B(K) = T                                              00000860
   10          CONTINUE                                                 00000870
               CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)                   00000880
   20       CONTINUE                                                    00000890
   30    CONTINUE                                                       00000900
C                                                                       00000910
C        NOW SOLVE  U*X = Y                                             00000920
C                                                                       00000930
         DO 40 KB = 1, N                                                00000940
            K = N + 1 - KB                                              00000950
            B(K) = B(K)/ABD(M,K)                                        00000960
            LM = MIN0(K,M) - 1                                          00000970
            LA = M - LM                                                 00000980
            LB = K - LM                                                 00000990
            T = -B(K)                                                   00001000
            CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)                        00001010
   40    CONTINUE                                                       00001020
      GO TO 100                                                         00001030
   50 CONTINUE                                                          00001040
C                                                                       00001050
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B                         00001060
C        FIRST SOLVE  TRANS(U)*Y = B                                    00001070
C                                                                       00001080
         DO 60 K = 1, N                                                 00001090
            LM = MIN0(K,M) - 1                                          00001100
            LA = M - LM                                                 00001110
            LB = K - LM                                                 00001120
            T = DDOT(LM,ABD(LA,K),1,B(LB),1)                            00001130
            B(K) = (B(K) - T)/ABD(M,K)                                  00001140
   60    CONTINUE                                                       00001150
C                                                                       00001160
C        NOW SOLVE TRANS(L)*X = Y                                       00001170
C                                                                       00001180
         IF (ML .EQ. 0) GO TO 90                                        00001190
         IF (NM1 .LT. 1) GO TO 90                                       00001200
            DO 80 KB = 1, NM1                                           00001210
               K = N - KB                                               00001220
               LM = MIN0(ML,N-K)                                        00001230
               B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)             00001240
               L = IPVT(K)                                              00001250
               IF (L .EQ. K) GO TO 70                                   00001260
                  T = B(L)                                              00001270
                  B(L) = B(K)                                           00001280
                  B(K) = T                                              00001290
   70          CONTINUE                                                 00001300
   80       CONTINUE                                                    00001310
   90    CONTINUE                                                       00001320
  100 CONTINUE                                                          00001330
      RETURN                                                            00001340
      END                                                               00001350
