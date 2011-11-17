      SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)                              00000010
      INTEGER LDA,N,IPVT(1),JOB                                         00000020
      DOUBLE PRECISION A(LDA,1),B(1)                                    00000030
C                                                                       00000040
C     DGESL SOLVES THE DOUBLE PRECISION SYSTEM                          00000050
C     A * X = B  OR  TRANS(A) * X = B                                   00000060
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.                     00000070
C                                                                       00000080
C     ON ENTRY                                                          00000090
C                                                                       00000100
C        A       DOUBLE PRECISION(LDA, N)                               00000110
C                THE OUTPUT FROM DGECO OR DGEFA.                        00000120
C                                                                       00000130
C        LDA     INTEGER                                                00000140
C                THE LEADING DIMENSION OF THE ARRAY  A .                00000150
C                                                                       00000160
C        N       INTEGER                                                00000170
C                THE ORDER OF THE MATRIX  A .                           00000180
C                                                                       00000190
C        IPVT    INTEGER(N)                                             00000200
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.                  00000210
C                                                                       00000220
C        B       DOUBLE PRECISION(N)                                    00000230
C                THE RIGHT HAND SIDE VECTOR.                            00000240
C                                                                       00000250
C        JOB     INTEGER                                                00000260
C                = 0         TO SOLVE  A*X = B ,                        00000270
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE            00000280
C                            TRANS(A)  IS THE TRANSPOSE.                00000290
C                                                                       00000300
C     ON RETURN                                                         00000310
C                                                                       00000320
C        B       THE SOLUTION VECTOR  X .                               00000330
C                                                                       00000340
C     ERROR CONDITION                                                   00000350
C                                                                       00000360
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A   00000370
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY  00000380
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER       00000390
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE     00000400
C        CALLED CORRECTLY AND IF DGECO HAS SET RCOND .GT. 0.0           00000410
C        OR DGEFA HAS SET INFO .EQ. 0 .                                 00000420
C                                                                       00000430
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 00000440
C     WITH  P  COLUMNS                                                  00000450
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)                            00000460
C           IF (RCOND IS TOO SMALL) GO TO ...                           00000470
C           DO 10 J = 1, P                                              00000480
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)                        00000490
C        10 CONTINUE                                                    00000500
C                                                                       00000510
C     LINPACK. THIS VERSION DATED 08/14/78 .                            00000520
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      00000530
C                                                                       00000540
C     SUBROUTINES AND FUNCTIONS                                         00000550
C                                                                       00000560
C     BLAS DAXPY,DDOT                                                   00000570
C                                                                       00000580
C     INTERNAL VARIABLES                                                00000590
C                                                                       00000600
      DOUBLE PRECISION DDOT,T                                           00000610
      INTEGER K,KB,L,NM1                                                00000620
C                                                                       00000630
      NM1 = N - 1                                                       00000640
      IF (JOB .NE. 0) GO TO 50                                          00000650
C                                                                       00000660
C        JOB = 0 , SOLVE  A * X = B                                     00000670
C        FIRST SOLVE  L*Y = B                                           00000680
C                                                                       00000690
         IF (NM1 .LT. 1) GO TO 30                                       00000700
         DO 20 K = 1, NM1                                               00000710
            L = IPVT(K)                                                 00000720
            T = B(L)                                                    00000730
            IF (L .EQ. K) GO TO 10                                      00000740
               B(L) = B(K)                                              00000750
               B(K) = T                                                 00000760
   10       CONTINUE                                                    00000770
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)                       00000780
   20    CONTINUE                                                       00000790
   30    CONTINUE                                                       00000800
C                                                                       00000810
C        NOW SOLVE  U*X = Y                                             00000820
C                                                                       00000830
         DO 40 KB = 1, N                                                00000840
            K = N + 1 - KB                                              00000850
            B(K) = B(K)/A(K,K)                                          00000860
            T = -B(K)                                                   00000870
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)                           00000880
   40    CONTINUE                                                       00000890
      GO TO 100                                                         00000900
   50 CONTINUE                                                          00000910
C                                                                       00000920
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B                         00000930
C        FIRST SOLVE  TRANS(U)*Y = B                                    00000940
C                                                                       00000950
         DO 60 K = 1, N                                                 00000960
            T = DDOT(K-1,A(1,K),1,B(1),1)                               00000970
            B(K) = (B(K) - T)/A(K,K)                                    00000980
   60    CONTINUE                                                       00000990
C                                                                       00001000
C        NOW SOLVE TRANS(L)*X = Y                                       00001010
C                                                                       00001020
         IF (NM1 .LT. 1) GO TO 90                                       00001030
         DO 80 KB = 1, NM1                                              00001040
            K = N - KB                                                  00001050
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)                 00001060
            L = IPVT(K)                                                 00001070
            IF (L .EQ. K) GO TO 70                                      00001080
               T = B(L)                                                 00001090
               B(L) = B(K)                                              00001100
               B(K) = T                                                 00001110
   70       CONTINUE                                                    00001120
   80    CONTINUE                                                       00001130
   90    CONTINUE                                                       00001140
  100 CONTINUE                                                          00001150
      RETURN                                                            00001160
      END                                                               00001170
