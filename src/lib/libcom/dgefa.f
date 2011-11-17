      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)                               00000010
      INTEGER LDA,N,IPVT(1),INFO                                        00000020
      DOUBLE PRECISION A(LDA,1)                                         00000030
C                                                                       00000040
C     DGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION.  00000050
C                                                                       00000060
C     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED            00000070
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          00000080
C     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .                   00000090
C                                                                       00000100
C     ON ENTRY                                                          00000110
C                                                                       00000120
C        A       DOUBLE PRECISION(LDA, N)                               00000130
C                THE MATRIX TO BE FACTORED.                             00000140
C                                                                       00000150
C        LDA     INTEGER                                                00000160
C                THE LEADING DIMENSION OF THE ARRAY  A .                00000170
C                                                                       00000180
C        N       INTEGER                                                00000190
C                THE ORDER OF THE MATRIX  A .                           00000200
C                                                                       00000210
C     ON RETURN                                                         00000220
C                                                                       00000230
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS         00000240
C                WHICH WERE USED TO OBTAIN IT.                          00000250
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       00000260
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          00000270
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       00000280
C                                                                       00000290
C        IPVT    INTEGER(N)                                             00000300
C                AN INTEGER VECTOR OF PIVOT INDICES.                    00000310
C                                                                       00000320
C        INFO    INTEGER                                                00000330
C                = 0  NORMAL VALUE.                                     00000340
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR       00000350
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES        00000360
C                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO  00000370
C                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE   00000380
C                     INDICATION OF SINGULARITY.                        00000390
C                                                                       00000400
C     LINPACK. THIS VERSION DATED 08/14/78 .                            00000410
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      00000420
C                                                                       00000430
C     SUBROUTINES AND FUNCTIONS                                         00000440
C                                                                       00000450
C     BLAS DAXPY,DSCAL,IDAMAX                                           00000460
C                                                                       00000470
C     INTERNAL VARIABLES                                                00000480
C                                                                       00000490
      DOUBLE PRECISION T                                                00000500
      INTEGER IDAMAX,J,K,KP1,L,NM1                                      00000510
C                                                                       00000520
C                                                                       00000530
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING                        00000540
C                                                                       00000550
      INFO = 0                                                          00000560
      NM1 = N - 1                                                       00000570
      IF (NM1 .LT. 1) GO TO 70                                          00000580
      DO 60 K = 1, NM1                                                  00000590
         KP1 = K + 1                                                    00000600
C                                                                       00000610
C        FIND L = PIVOT INDEX                                           00000620
C                                                                       00000630
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1                             00000640
         IPVT(K) = L                                                    00000650
C                                                                       00000660
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED          00000670
C                                                                       00000680
         IF (A(L,K) .EQ. 0.0D0) GO TO 40                                00000690
C                                                                       00000700
C           INTERCHANGE IF NECESSARY                                    00000710
C                                                                       00000720
            IF (L .EQ. K) GO TO 10                                      00000730
               T = A(L,K)                                               00000740
               A(L,K) = A(K,K)                                          00000750
               A(K,K) = T                                               00000760
   10       CONTINUE                                                    00000770
C                                                                       00000780
C           COMPUTE MULTIPLIERS                                         00000790
C                                                                       00000800
            T = -1.0D0/A(K,K)                                           00000810
            CALL DSCAL(N-K,T,A(K+1,K),1)                                00000820
C                                                                       00000830
C           ROW ELIMINATION WITH COLUMN INDEXING                        00000840
C                                                                       00000850
            DO 30 J = KP1, N                                            00000860
               T = A(L,J)                                               00000870
               IF (L .EQ. K) GO TO 20                                   00000880
                  A(L,J) = A(K,J)                                       00000890
                  A(K,J) = T                                            00000900
   20          CONTINUE                                                 00000910
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)                  00000920
   30       CONTINUE                                                    00000930
         GO TO 50                                                       00000940
   40    CONTINUE                                                       00000950
            INFO = K                                                    00000960
   50    CONTINUE                                                       00000970
   60 CONTINUE                                                          00000980
   70 CONTINUE                                                          00000990
      IPVT(N) = N                                                       00001000
      IF (A(N,N) .EQ. 0.0D0) INFO = N                                   00001010
      RETURN                                                            00001020
      END                                                               00001030
