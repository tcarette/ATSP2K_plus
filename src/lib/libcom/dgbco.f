      SUBROUTINE DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)                    00000010
      INTEGER LDA,N,ML,MU,IPVT(1)                                       00000020
      DOUBLE PRECISION ABD(LDA,1),Z(1)                                  00000030
      DOUBLE PRECISION RCOND                                            00000040
C                                                                       00000050
C     DGBCO FACTORS A DOUBLE PRECISION BAND MATRIX BY GAUSSIAN          00000060
C     ELIMINATION AND ESTIMATES THE CONDITION OF THE MATRIX.            00000070
C                                                                       00000080
C     IF  RCOND  IS NOT NEEDED, DGBFA IS SLIGHTLY FASTER.               00000090
C     TO SOLVE  A*X = B , FOLLOW DGBCO BY DGBSL.                        00000100
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DGBCO BY DGBSL.                 00000110
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DGBCO BY DGBDI.               00000120
C                                                                       00000130
C     ON ENTRY                                                          00000140
C                                                                       00000150
C        ABD     DOUBLE PRECISION(LDA, N)                               00000160
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS      00000170
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND   00000180
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS         00000190
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .                       00000200
C                SEE THE COMMENTS BELOW FOR DETAILS.                    00000210
C                                                                       00000220
C        LDA     INTEGER                                                00000230
C                THE LEADING DIMENSION OF THE ARRAY  ABD .              00000240
C                LDA MUST BE .GE. 2*ML + MU + 1 .                       00000250
C                                                                       00000260
C        N       INTEGER                                                00000270
C                THE ORDER OF THE ORIGINAL MATRIX.                      00000280
C                                                                       00000290
C        ML      INTEGER                                                00000300
C                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.           00000310
C                0 .LE. ML .LT. N .                                     00000320
C                                                                       00000330
C        MU      INTEGER                                                00000340
C                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.           00000350
C                0 .LE. MU .LT. N .                                     00000360
C                MORE EFFICIENT IF  ML .LE. MU .                        00000370
C                                                                       00000380
C     ON RETURN                                                         00000390
C                                                                       00000400
C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND         00000410
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.          00000420
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE       00000430
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER          00000440
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.       00000450
C                                                                       00000460
C        IPVT    INTEGER(N)                                             00000470
C                AN INTEGER VECTOR OF PIVOT INDICES.                    00000480
C                                                                       00000490
C        RCOND   DOUBLE PRECISION                                       00000500
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .        00000510
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS       00000520
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE             00000530
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . 00000540
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION     00000550
C                           1.0 + RCOND .EQ. 1.0                        00000560
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING           00000570
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF         00000580
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE          00000590
C                UNDERFLOWS.                                            00000600
C                                                                       00000610
C        Z       DOUBLE PRECISION(N)                                    00000620
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.  00000630
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS      00000640
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT           00000650
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    00000660
C                                                                       00000670
C     BAND STORAGE                                                      00000680
C                                                                       00000690
C           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT      00000700
C           WILL SET UP THE INPUT.                                      00000710
C                                                                       00000720
C                   ML = (BAND WIDTH BELOW THE DIAGONAL)                00000730
C                   MU = (BAND WIDTH ABOVE THE DIAGONAL)                00000740
C                   M = ML + MU + 1                                     00000750
C                   DO 20 J = 1, N                                      00000760
C                      I1 = MAX0(1, J-MU)                               00000770
C                      I2 = MIN0(N, J+ML)                               00000780
C                      DO 10 I = I1, I2                                 00000790
C                         K = I - J + M                                 00000800
C                         ABD(K,J) = A(I,J)                             00000810
C                10    CONTINUE                                         00000820
C                20 CONTINUE                                            00000830
C                                                                       00000840
C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .         00000850
C           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR      00000860
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.            00000870
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .    00000880
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE            00000890
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.          00000900
C                                                                       00000910
C     EXAMPLE..  IF THE ORIGINAL MATRIX IS                              00000920
C                                                                       00000930
C           11 12 13  0  0  0                                           00000940
C           21 22 23 24  0  0                                           00000950
C            0 32 33 34 35  0                                           00000960
C            0  0 43 44 45 46                                           00000970
C            0  0  0 54 55 56                                           00000980
C            0  0  0  0 65 66                                           00000990
C                                                                       00001000
C      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN  00001010
C                                                                       00001020
C            *  *  *  +  +  +  , * = NOT USED                           00001030
C            *  * 13 24 35 46  , + = USED FOR PIVOTING                  00001040
C            * 12 23 34 45 56                                           00001050
C           11 22 33 44 55 66                                           00001060
C           21 32 43 54 65  *                                           00001070
C                                                                       00001080
C     LINPACK. THIS VERSION DATED 08/14/78 .                            00001090
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      00001100
C                                                                       00001110
C     SUBROUTINES AND FUNCTIONS                                         00001120
C                                                                       00001130
C     LINPACK DGBFA                                                     00001140
C     BLAS DAXPY,DDOT,DSCAL,DASUM                                       00001150
C     FORTRAN DABS,DMAX1,MAX0,MIN0,DSIGN                                00001160
C                                                                       00001170
C     INTERNAL VARIABLES                                                00001180
C                                                                       00001190
      DOUBLE PRECISION DDOT,EK,T,WK,WKM                                 00001200
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM                           00001210
      INTEGER IS,INFO,J,JU,K,KB,KP1,L,LA,LM,LZ,M,MM                     00001220
C                                                                       00001230
C                                                                       00001240
C     COMPUTE 1-NORM OF A                                               00001250
C                                                                       00001260
      ANORM = 0.0D0                                                     00001270
      L = ML + 1                                                        00001280
      IS = L + MU                                                       00001290
      DO 10 J = 1, N                                                    00001300
         ANORM = DMAX1(ANORM,DASUM(L,ABD(IS,J),1))                      00001310
         IF (IS .GT. ML + 1) IS = IS - 1                                00001320
         IF (J .LE. MU) L = L + 1                                       00001330
         IF (J .GE. N - ML) L = L - 1                                   00001340
   10 CONTINUE                                                          00001350
C                                                                       00001360
C     FACTOR                                                            00001370
C                                                                       00001380
      CALL DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)                             00001390
C                                                                       00001400
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .              00001410
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .  00001420
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE      00001430
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE  00001440
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID    00001450
C     OVERFLOW.                                                         00001460
C                                                                       00001470
C     SOLVE TRANS(U)*W = E                                              00001480
C                                                                       00001490
      EK = 1.0D0                                                        00001500
      DO 20 J = 1, N                                                    00001510
         Z(J) = 0.0D0                                                   00001520
   20 CONTINUE                                                          00001530
      M = ML + MU + 1                                                   00001540
      JU = 0                                                            00001550
      DO 100 K = 1, N                                                   00001560
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))                      00001570
         IF (DABS(EK-Z(K)) .LE. DABS(ABD(M,K))) GO TO 30                00001580
            S = DABS(ABD(M,K))/DABS(EK-Z(K))                            00001590
            CALL DSCAL(N,S,Z,1)                                         00001600
            EK = S*EK                                                   00001610
   30    CONTINUE                                                       00001620
         WK = EK - Z(K)                                                 00001630
         WKM = -EK - Z(K)                                               00001640
         S = DABS(WK)                                                   00001650
         SM = DABS(WKM)                                                 00001660
         IF (ABD(M,K) .EQ. 0.0D0) GO TO 40                              00001670
            WK = WK/ABD(M,K)                                            00001680
            WKM = WKM/ABD(M,K)                                          00001690
         GO TO 50                                                       00001700
   40    CONTINUE                                                       00001710
            WK = 1.0D0                                                  00001720
            WKM = 1.0D0                                                 00001730
   50    CONTINUE                                                       00001740
         KP1 = K + 1                                                    00001750
         JU = MIN0(MAX0(JU,MU+IPVT(K)),N)                               00001760
         MM = M                                                         00001770
         IF (KP1 .GT. JU) GO TO 90                                      00001780
            DO 60 J = KP1, JU                                           00001790
               MM = MM - 1                                              00001800
               SM = SM + DABS(Z(J)+WKM*ABD(MM,J))                       00001810
               Z(J) = Z(J) + WK*ABD(MM,J)                               00001820
               S = S + DABS(Z(J))                                       00001830
   60       CONTINUE                                                    00001840
            IF (S .GE. SM) GO TO 80                                     00001850
               T = WKM - WK                                             00001860
               WK = WKM                                                 00001870
               MM = M                                                   00001880
               DO 70 J = KP1, JU                                        00001890
                  MM = MM - 1                                           00001900
                  Z(J) = Z(J) + T*ABD(MM,J)                             00001910
   70          CONTINUE                                                 00001920
   80       CONTINUE                                                    00001930
   90    CONTINUE                                                       00001940
         Z(K) = WK                                                      00001950
  100 CONTINUE                                                          00001960
      S = 1.0D0/DASUM(N,Z,1)                                            00001970
      CALL DSCAL(N,S,Z,1)                                               00001980
C                                                                       00001990
C     SOLVE TRANS(L)*Y = W                                              00002000
C                                                                       00002010
      DO 120 KB = 1, N                                                  00002020
         K = N + 1 - KB                                                 00002030
         LM = MIN0(ML,N-K)                                              00002040
         IF (K .LT. N) Z(K) = Z(K) + DDOT(LM,ABD(M+1,K),1,Z(K+1),1)     00002050
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110                           00002060
            S = 1.0D0/DABS(Z(K))                                        00002070
            CALL DSCAL(N,S,Z,1)                                         00002080
  110    CONTINUE                                                       00002090
         L = IPVT(K)                                                    00002100
         T = Z(L)                                                       00002110
         Z(L) = Z(K)                                                    00002120
         Z(K) = T                                                       00002130
  120 CONTINUE                                                          00002140
      S = 1.0D0/DASUM(N,Z,1)                                            00002150
      CALL DSCAL(N,S,Z,1)                                               00002160
C                                                                       00002170
      YNORM = 1.0D0                                                     00002180
C                                                                       00002190
C     SOLVE L*V = Y                                                     00002200
C                                                                       00002210
      DO 140 K = 1, N                                                   00002220
         L = IPVT(K)                                                    00002230
         T = Z(L)                                                       00002240
         Z(L) = Z(K)                                                    00002250
         Z(K) = T                                                       00002260
         LM = MIN0(ML,N-K)                                              00002270
         IF (K .LT. N) CALL DAXPY(LM,T,ABD(M+1,K),1,Z(K+1),1)           00002280
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130                           00002290
            S = 1.0D0/DABS(Z(K))                                        00002300
            CALL DSCAL(N,S,Z,1)                                         00002310
            YNORM = S*YNORM                                             00002320
  130    CONTINUE                                                       00002330
  140 CONTINUE                                                          00002340
      S = 1.0D0/DASUM(N,Z,1)                                            00002350
      CALL DSCAL(N,S,Z,1)                                               00002360
      YNORM = S*YNORM                                                   00002370
C                                                                       00002380
C     SOLVE  U*Z = W                                                    00002390
C                                                                       00002400
      DO 160 KB = 1, N                                                  00002410
         K = N + 1 - KB                                                 00002420
         IF (DABS(Z(K)) .LE. DABS(ABD(M,K))) GO TO 150                  00002430
            S = DABS(ABD(M,K))/DABS(Z(K))                               00002440
            CALL DSCAL(N,S,Z,1)                                         00002450
            YNORM = S*YNORM                                             00002460
  150    CONTINUE                                                       00002470
         IF (ABD(M,K) .NE. 0.0D0) Z(K) = Z(K)/ABD(M,K)                  00002480
         IF (ABD(M,K) .EQ. 0.0D0) Z(K) = 1.0D0                          00002490
         LM = MIN0(K,M) - 1                                             00002500
         LA = M - LM                                                    00002510
         LZ = K - LM                                                    00002520
         T = -Z(K)                                                      00002530
         CALL DAXPY(LM,T,ABD(LA,K),1,Z(LZ),1)                           00002540
  160 CONTINUE                                                          00002550
C     MAKE ZNORM = 1.0                                                  00002560
      S = 1.0D0/DASUM(N,Z,1)                                            00002570
      CALL DSCAL(N,S,Z,1)                                               00002580
      YNORM = S*YNORM                                                   00002590
C                                                                       00002600
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM                         00002610
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0                               00002620
      RETURN                                                            00002630
      END                                                               00002640
