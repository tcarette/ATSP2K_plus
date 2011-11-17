      SUBROUTINE DPBSL(ABD,LDA,N,M,B)                                   00000010
      INTEGER LDA,N,M                                                   00000020
      DOUBLE PRECISION ABD(LDA,1),B(1)                                  00000030
C                                                                       00000040
C     DPBSL SOLVES THE DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE     00000050
C     BAND SYSTEM  A*X = B                                              00000060
C     USING THE FACTORS COMPUTED BY DPBCO OR DPBFA.                     00000070
C                                                                       00000080
C     ON ENTRY                                                          00000090
C                                                                       00000100
C        ABD     DOUBLE PRECISION(LDA, N)                               00000110
C                THE OUTPUT FROM DPBCO OR DPBFA.                        00000120
C                                                                       00000130
C        LDA     INTEGER                                                00000140
C                THE LEADING DIMENSION OF THE ARRAY  ABD .              00000150
C                                                                       00000160
C        N       INTEGER                                                00000170
C                THE ORDER OF THE MATRIX  A .                           00000180
C                                                                       00000190
C        M       INTEGER                                                00000200
C                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.       00000210
C                                                                       00000220
C        B       DOUBLE PRECISION(N)                                    00000230
C                THE RIGHT HAND SIDE VECTOR.                            00000240
C                                                                       00000250
C     ON RETURN                                                         00000260
C                                                                       00000270
C        B       THE SOLUTION VECTOR  X .                               00000280
C                                                                       00000290
C     ERROR CONDITION                                                   00000300
C                                                                       00000310
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS     00000320
C        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES            00000330
C        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE    00000340
C        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED    00000350
C        CORRECTLY AND  INFO .EQ. 0 .                                   00000360
C                                                                       00000370
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX                 00000380
C     WITH  P  COLUMNS                                                  00000390
C           CALL DPBCO(ABD,LDA,N,RCOND,Z,INFO)                          00000400
C           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ...          00000410
C           DO 10 J = 1, P                                              00000420
C              CALL DPBSL(ABD,LDA,N,C(1,J))                             00000430
C        10 CONTINUE                                                    00000440
C                                                                       00000450
C     LINPACK.  THIS VERSION DATED 08/14/78 .                           00000460
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      00000470
C                                                                       00000480
C     SUBROUTINES AND FUNCTIONS                                         00000490
C                                                                       00000500
C     BLAS DAXPY,DDOT                                                   00000510
C     FORTRAN MIN0                                                      00000520
C                                                                       00000530
C     INTERNAL VARIABLES                                                00000540
C                                                                       00000550
      DOUBLE PRECISION DDOT,T                                           00000560
      INTEGER K,KB,LA,LB,LM                                             00000570
C                                                                       00000580
C     SOLVE TRANS(R)*Y = B                                              00000590
C                                                                       00000600
      DO 10 K = 1, N                                                    00000610
         LM = MIN0(K-1,M)                                               00000620
         LA = M + 1 - LM                                                00000630
         LB = K - LM                                                    00000640
         T = DDOT(LM,ABD(LA,K),1,B(LB),1)                               00000650
         B(K) = (B(K) - T)/ABD(M+1,K)                                   00000660
   10 CONTINUE                                                          00000670
C                                                                       00000680
C     SOLVE R*X = Y                                                     00000690
C                                                                       00000700
      DO 20 KB = 1, N                                                   00000710
         K = N + 1 - KB                                                 00000720
         LM = MIN0(K-1,M)                                               00000730
         LA = M + 1 - LM                                                00000740
         LB = K - LM                                                    00000750
         B(K) = B(K)/ABD(M+1,K)                                         00000760
         T = -B(K)                                                      00000770
         CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)                           00000780
   20 CONTINUE                                                          00000790
      RETURN                                                            00000800
      END                                                               00000810
