      SUBROUTINE DGECO(A,LDA,N,IPVT,RCOND,Z)                            00000010
      INTEGER LDA,N,IPVT(1)                                             00000020
      DOUBLE PRECISION A(LDA,1),Z(1)                                    00000030
      DOUBLE PRECISION RCOND                                            00000040
C                                                                       00000050
C     DGECO FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION   00000060
C     AND ESTIMATES THE CONDITION OF THE MATRIX.                        00000070
C                                                                       00000080
C     IF  RCOND  IS NOT NEEDED, DGEFA IS SLIGHTLY FASTER.               00000090
C     TO SOLVE  A*X = B , FOLLOW DGECO BY DGESL.                        00000100
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DGECO BY DGESL.                 00000110
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DGECO BY DGEDI.               00000120
C     TO COMPUTE  INVERSE(A) , FOLLOW DGECO BY DGEDI.                   00000130
C                                                                       00000140
C     ON ENTRY                                                          00000150
C                                                                       00000160
C        A       DOUBLE PRECISION(LDA, N)                               00000170
C                THE MATRIX TO BE FACTORED.                             00000180
C                                                                       00000190
C        LDA     INTEGER                                                00000200
C                THE LEADING DIMENSION OF THE ARRAY  A .                00000210
C                                                                       00000220
C        N       INTEGER                                                00000230
C                THE ORDER OF THE MATRIX  A .                           00000240
C                                                                       00000250
C     ON RETURN                                                         00000260
C                                                                       00000270
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS         00000280
C                WHICH WERE USED TO OBTAIN IT.                          00000290
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       00000300
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          00000310
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       00000320
C                                                                       00000330
C        IPVT    INTEGER(N)                                             00000340
C                AN INTEGER VECTOR OF PIVOT INDICES.                    00000350
C                                                                       00000360
C        RCOND   DOUBLE PRECISION                                       00000370
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .        00000380
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS       00000390
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE             00000400
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . 00000410
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION     00000420
C                           1.0 + RCOND .EQ. 1.0                        00000430
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING           00000440
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF         00000450
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE          00000460
C                UNDERFLOWS.                                            00000470
C                                                                       00000480
C        Z       DOUBLE PRECISION(N)                                    00000490
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.  00000500
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS      00000510
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT           00000520
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    00000530
C                                                                       00000540
C     LINPACK. THIS VERSION DATED 08/14/78 .                            00000550
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      00000560
C                                                                       00000570
C     SUBROUTINES AND FUNCTIONS                                         00000580
C                                                                       00000590
C     LINPACK DGEFA                                                     00000600
C     BLAS DAXPY,DDOT,DSCAL,DASUM                                       00000610
C     FORTRAN DABS,DMAX1,DSIGN                                          00000620
C                                                                       00000630
C     INTERNAL VARIABLES                                                00000640
C                                                                       00000650
      DOUBLE PRECISION DDOT,EK,T,WK,WKM                                 00000660
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM                           00000670
      INTEGER INFO,J,K,KB,KP1,L                                         00000680
C                                                                       00000690
C                                                                       00000700
C     COMPUTE 1-NORM OF A                                               00000710
C                                                                       00000720
      ANORM = 0.0D0                                                     00000730
      DO 10 J = 1, N                                                    00000740
         ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))                         00000750
   10 CONTINUE                                                          00000760
C                                                                       00000770
C     FACTOR                                                            00000780
C                                                                       00000790
      CALL DGEFA(A,LDA,N,IPVT,INFO)                                     00000800
C                                                                       00000810
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .              00000820
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .  00000830
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE      00000840
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE  00000850
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID    00000860
C     OVERFLOW.                                                         00000870
C                                                                       00000880
C     SOLVE TRANS(U)*W = E                                              00000890
C                                                                       00000900
      EK = 1.0D0                                                        00000910
      DO 20 J = 1, N                                                    00000920
         Z(J) = 0.0D0                                                   00000930
   20 CONTINUE                                                          00000940
      DO 100 K = 1, N                                                   00000950
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))                      00000960
         IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 30                  00000970
            S = DABS(A(K,K))/DABS(EK-Z(K))                              00000980
            CALL DSCAL(N,S,Z,1)                                         00000990
            EK = S*EK                                                   00001000
   30    CONTINUE                                                       00001010
         WK = EK - Z(K)                                                 00001020
         WKM = -EK - Z(K)                                               00001030
         S = DABS(WK)                                                   00001040
         SM = DABS(WKM)                                                 00001050
         IF (A(K,K) .EQ. 0.0D0) GO TO 40                                00001060
            WK = WK/A(K,K)                                              00001070
            WKM = WKM/A(K,K)                                            00001080
         GO TO 50                                                       00001090
   40    CONTINUE                                                       00001100
            WK = 1.0D0                                                  00001110
            WKM = 1.0D0                                                 00001120
   50    CONTINUE                                                       00001130
         KP1 = K + 1                                                    00001140
         IF (KP1 .GT. N) GO TO 90                                       00001150
            DO 60 J = KP1, N                                            00001160
               SM = SM + DABS(Z(J)+WKM*A(K,J))                          00001170
               Z(J) = Z(J) + WK*A(K,J)                                  00001180
               S = S + DABS(Z(J))                                       00001190
   60       CONTINUE                                                    00001200
            IF (S .GE. SM) GO TO 80                                     00001210
               T = WKM - WK                                             00001220
               WK = WKM                                                 00001230
               DO 70 J = KP1, N                                         00001240
                  Z(J) = Z(J) + T*A(K,J)                                00001250
   70          CONTINUE                                                 00001260
   80       CONTINUE                                                    00001270
   90    CONTINUE                                                       00001280
         Z(K) = WK                                                      00001290
  100 CONTINUE                                                          00001300
      S = 1.0D0/DASUM(N,Z,1)                                            00001310
      CALL DSCAL(N,S,Z,1)                                               00001320
C                                                                       00001330
C     SOLVE TRANS(L)*Y = W                                              00001340
C                                                                       00001350
      DO 120 KB = 1, N                                                  00001360
         K = N + 1 - KB                                                 00001370
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)      00001380
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110                           00001390
            S = 1.0D0/DABS(Z(K))                                        00001400
            CALL DSCAL(N,S,Z,1)                                         00001410
  110    CONTINUE                                                       00001420
         L = IPVT(K)                                                    00001430
         T = Z(L)                                                       00001440
         Z(L) = Z(K)                                                    00001450
         Z(K) = T                                                       00001460
  120 CONTINUE                                                          00001470
      S = 1.0D0/DASUM(N,Z,1)                                            00001480
      CALL DSCAL(N,S,Z,1)                                               00001490
C                                                                       00001500
      YNORM = 1.0D0                                                     00001510
C                                                                       00001520
C     SOLVE L*V = Y                                                     00001530
C                                                                       00001540
      DO 140 K = 1, N                                                   00001550
         L = IPVT(K)                                                    00001560
         T = Z(L)                                                       00001570
         Z(L) = Z(K)                                                    00001580
         Z(K) = T                                                       00001590
         IF (K .LT. N) CALL DAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)            00001600
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130                           00001610
            S = 1.0D0/DABS(Z(K))                                        00001620
            CALL DSCAL(N,S,Z,1)                                         00001630
            YNORM = S*YNORM                                             00001640
  130    CONTINUE                                                       00001650
  140 CONTINUE                                                          00001660
      S = 1.0D0/DASUM(N,Z,1)                                            00001670
      CALL DSCAL(N,S,Z,1)                                               00001680
      YNORM = S*YNORM                                                   00001690
C                                                                       00001700
C     SOLVE  U*Z = V                                                    00001710
C                                                                       00001720
      DO 160 KB = 1, N                                                  00001730
         K = N + 1 - KB                                                 00001740
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 150                    00001750
            S = DABS(A(K,K))/DABS(Z(K))                                 00001760
            CALL DSCAL(N,S,Z,1)                                         00001770
            YNORM = S*YNORM                                             00001780
  150    CONTINUE                                                       00001790
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)                      00001800
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0                            00001810
         T = -Z(K)                                                      00001820
         CALL DAXPY(K-1,T,A(1,K),1,Z(1),1)                              00001830
  160 CONTINUE                                                          00001840
C     MAKE ZNORM = 1.0                                                  00001850
      S = 1.0D0/DASUM(N,Z,1)                                            00001860
      CALL DSCAL(N,S,Z,1)                                               00001870
      YNORM = S*YNORM                                                   00001880
C                                                                       00001890
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM                         00001900
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0                               00001910
      RETURN                                                            00001920
      END                                                               00001930
