      SUBROUTINE DPBCO(ABD,LDA,N,M,RCOND,Z,INFO)                        00000010
      INTEGER LDA,N,M,INFO                                              00000020
      DOUBLE PRECISION ABD(LDA,1),Z(1)                                  00000030
      DOUBLE PRECISION RCOND                                            00000040
C                                                                       00000050
C     DPBCO FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE      00000060
C     MATRIX STORED IN BAND FORM AND ESTIMATES THE CONDITION OF THE     00000070
C     MATRIX.                                                           00000080
C                                                                       00000090
C     IF  RCOND  IS NOT NEEDED, DPBFA IS SLIGHTLY FASTER.               00000100
C     TO SOLVE  A*X = B , FOLLOW DPBCO BY DPBSL.                        00000110
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DPBCO BY DPBSL.                 00000120
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DPBCO BY DPBDI.               00000130
C                                                                       00000140
C     ON ENTRY                                                          00000150
C                                                                       00000160
C        ABD     DOUBLE PRECISION(LDA, N)                               00000170
C                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER   00000180
C                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE      00000190
C                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE      00000200
C                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS.     00000210
C                                                                       00000220
C        LDA     INTEGER                                                00000230
C                THE LEADING DIMENSION OF THE ARRAY  ABD .              00000240
C                LDA MUST BE .GE. M + 1 .                               00000250
C                                                                       00000260
C        N       INTEGER                                                00000270
C                THE ORDER OF THE MATRIX  A .                           00000280
C                                                                       00000290
C        M       INTEGER                                                00000300
C                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.       00000310
C                0 .LE. M .LT. N .                                      00000320
C                                                                       00000330
C     ON RETURN                                                         00000340
C                                                                       00000350
C        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND         00000360
C                FORM, SO THAT  A = TRANS(R)*R .                        00000370
C                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.   00000380
C                                                                       00000390
C        RCOND   DOUBLE PRECISION                                       00000400
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .        00000410
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS       00000420
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE             00000430
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . 00000440
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION     00000450
C                           1.0 + RCOND .EQ. 1.0                        00000460
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING           00000470
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF         00000480
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE          00000490
C                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED.      00000500
C                                                                       00000510
C        Z       DOUBLE PRECISION(N)                                    00000520
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.  00000530
C                IF  A  IS SINGULAR TO WORKING PRECISION, THEN  Z  IS   00000540
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT           00000550
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .                    00000560
C                IF  INFO .NE. 0 , Z  IS UNCHANGED.                     00000570
C                                                                       00000580
C        INFO    INTEGER                                                00000590
C                = 0  FOR NORMAL RETURN.                                00000600
C                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR    00000610
C                     OF ORDER  K  IS NOT POSITIVE DEFINITE.            00000620
C                                                                       00000630
C     BAND STORAGE                                                      00000640
C                                                                       00000650
C           IF  A  IS A SYMMETRIC POSITIVE DEFINITE BAND MATRIX,        00000660
C           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT.        00000670
C                                                                       00000680
C                   M = (BAND WIDTH ABOVE DIAGONAL)                     00000690
C                   DO 20 J = 1, N                                      00000700
C                      I1 = MAX0(1, J-M)                                00000710
C                      DO 10 I = I1, J                                  00000720
C                         K = I-J+M+1                                   00000730
C                         ABD(K,J) = A(I,J)                             00000740
C                10    CONTINUE                                         00000750
C                20 CONTINUE                                            00000760
C                                                                       00000770
C           THIS USES  M + 1  ROWS OF  A , EXCEPT FOR THE  M BY M       00000780
C           UPPER LEFT TRIANGLE, WHICH IS IGNORED.                      00000790
C                                                                       00000800
C     EXAMPLE..  IF THE ORIGINAL MATRIX IS                              00000810
C                                                                       00000820
C           11 12 13  0  0  0                                           00000830
C           12 22 23 24  0  0                                           00000840
C           13 23 33 34 35  0                                           00000850
C            0 24 34 44 45 46                                           00000860
C            0  0 35 45 55 56                                           00000870
C            0  0  0 46 56 66                                           00000880
C                                                                       00000890
C     THEN  N = 6 , M = 2  AND  ABD  SHOULD CONTAIN                     00000900
C                                                                       00000910
C            *  * 13 24 35 46                                           00000920
C            * 12 23 34 45 56                                           00000930
C           11 22 33 44 55 66                                           00000940
C                                                                       00000950
C     LINPACK.  THIS VERSION DATED 08/14/78 .                           00000960
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      00000970
C                                                                       00000980
C     SUBROUTINES AND FUNCTIONS                                         00000990
C                                                                       00001000
C     LINPACK DPBFA                                                     00001010
C     BLAS DAXPY,DDOT,DSCAL,DASUM                                       00001020
C     FORTRAN DABS,DMAX1,MAX0,MIN0,DREAL,DSIGN                          00001030
C                                                                       00001040
C     INTERNAL VARIABLES                                                00001050
C                                                                       00001060
      DOUBLE PRECISION DDOT,EK,T,WK,WKM                                 00001070
      DOUBLE PRECISION ANORM,S,DASUM,SM,YNORM                           00001080
      INTEGER I,J,J2,K,KB,KP1,L,LA,LB,LM,MU                             00001090
C                                                                       00001100
C                                                                       00001110
C     FIND NORM OF A                                                    00001120
C                                                                       00001130
      DO 30 J = 1, N                                                    00001140
         L = MIN0(J,M+1)                                                00001150
         MU = MAX0(M+2-J,1)                                             00001160
         Z(J) = DASUM(L,ABD(MU,J),1)                                    00001170
         K = J - L                                                      00001180
         IF (M .LT. MU) GO TO 20                                        00001190
         DO 10 I = MU, M                                                00001200
            K = K + 1                                                   00001210
            Z(K) = Z(K) + DABS(ABD(I,J))                                00001220
   10    CONTINUE                                                       00001230
   20    CONTINUE                                                       00001240
   30 CONTINUE                                                          00001250
      ANORM = 0.0D0                                                     00001260
      DO 40 J = 1, N                                                    00001270
         ANORM = DMAX1(ANORM,Z(J))                                      00001280
   40 CONTINUE                                                          00001290
C                                                                       00001300
C     FACTOR                                                            00001310
C                                                                       00001320
      CALL DPBFA(ABD,LDA,N,M,INFO)                                      00001330
      IF (INFO .NE. 0) GO TO 180                                        00001340
C                                                                       00001350
C        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .           00001360
C        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .      00001370
C        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL        00001380
C        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .           00001390
C        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.         00001400
C                                                                       00001410
C        SOLVE TRANS(R)*W = E                                           00001420
C                                                                       00001430
         EK = 1.0D0                                                     00001440
         DO 50 J = 1, N                                                 00001450
            Z(J) = 0.0D0                                                00001460
   50    CONTINUE                                                       00001470
         DO 110 K = 1, N                                                00001480
            IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))                   00001490
            IF (DABS(EK-Z(K)) .LE. ABD(M+1,K)) GO TO 60                 00001500
               S = ABD(M+1,K)/DABS(EK-Z(K))                             00001510
               CALL DSCAL(N,S,Z,1)                                      00001520
               EK = S*EK                                                00001530
   60       CONTINUE                                                    00001540
            WK = EK - Z(K)                                              00001550
            WKM = -EK - Z(K)                                            00001560
            S = DABS(WK)                                                00001570
            SM = DABS(WKM)                                              00001580
            WK = WK/ABD(M+1,K)                                          00001590
            WKM = WKM/ABD(M+1,K)                                        00001600
            KP1 = K + 1                                                 00001610
            J2 = MIN0(K+M,N)                                            00001620
            I = M + 1                                                   00001630
            IF (KP1 .GT. J2) GO TO 100                                  00001640
               DO 70 J = KP1, J2                                        00001650
                  I = I - 1                                             00001660
                  SM = SM + DABS(Z(J)+WKM*ABD(I,J))                     00001670
                  Z(J) = Z(J) + WK*ABD(I,J)                             00001680
                  S = S + DABS(Z(J))                                    00001690
   70          CONTINUE                                                 00001700
               IF (S .GE. SM) GO TO 90                                  00001710
                  T = WKM - WK                                          00001720
                  WK = WKM                                              00001730
                  I = M + 1                                             00001740
                  DO 80 J = KP1, J2                                     00001750
                     I = I - 1                                          00001760
                     Z(J) = Z(J) + T*ABD(I,J)                           00001770
   80             CONTINUE                                              00001780
   90          CONTINUE                                                 00001790
  100       CONTINUE                                                    00001800
            Z(K) = WK                                                   00001810
  110    CONTINUE                                                       00001820
         S = 1.0D0/DASUM(N,Z,1)                                         00001830
         CALL DSCAL(N,S,Z,1)                                            00001840
C                                                                       00001850
C        SOLVE  R*Y = W                                                 00001860
C                                                                       00001870
         DO 130 KB = 1, N                                               00001880
            K = N + 1 - KB                                              00001890
            IF (DABS(Z(K)) .LE. ABD(M+1,K)) GO TO 120                   00001900
               S = ABD(M+1,K)/DABS(Z(K))                                00001910
               CALL DSCAL(N,S,Z,1)                                      00001920
  120       CONTINUE                                                    00001930
            Z(K) = Z(K)/ABD(M+1,K)                                      00001940
            LM = MIN0(K-1,M)                                            00001950
            LA = M + 1 - LM                                             00001960
            LB = K - LM                                                 00001970
            T = -Z(K)                                                   00001980
            CALL DAXPY(LM,T,ABD(LA,K),1,Z(LB),1)                        00001990
  130    CONTINUE                                                       00002000
         S = 1.0D0/DASUM(N,Z,1)                                         00002010
         CALL DSCAL(N,S,Z,1)                                            00002020
C                                                                       00002030
         YNORM = 1.0D0                                                  00002040
C                                                                       00002050
C        SOLVE TRANS(R)*V = Y                                           00002060
C                                                                       00002070
         DO 150 K = 1, N                                                00002080
            LM = MIN0(K-1,M)                                            00002090
            LA = M + 1 - LM                                             00002100
            LB = K - LM                                                 00002110
            Z(K) = Z(K) - DDOT(LM,ABD(LA,K),1,Z(LB),1)                  00002120
            IF (DABS(Z(K)) .LE. ABD(M+1,K)) GO TO 140                   00002130
               S = ABD(M+1,K)/DABS(Z(K))                                00002140
               CALL DSCAL(N,S,Z,1)                                      00002150
               YNORM = S*YNORM                                          00002160
  140       CONTINUE                                                    00002170
            Z(K) = Z(K)/ABD(M+1,K)                                      00002180
  150    CONTINUE                                                       00002190
         S = 1.0D0/DASUM(N,Z,1)                                         00002200
         CALL DSCAL(N,S,Z,1)                                            00002210
         YNORM = S*YNORM                                                00002220
C                                                                       00002230
C        SOLVE  R*Z = W                                                 00002240
C                                                                       00002250
         DO 170 KB = 1, N                                               00002260
            K = N + 1 - KB                                              00002270
            IF (DABS(Z(K)) .LE. ABD(M+1,K)) GO TO 160                   00002280
               S = ABD(M+1,K)/DABS(Z(K))                                00002290
               CALL DSCAL(N,S,Z,1)                                      00002300
               YNORM = S*YNORM                                          00002310
  160       CONTINUE                                                    00002320
            Z(K) = Z(K)/ABD(M+1,K)                                      00002330
            LM = MIN0(K-1,M)                                            00002340
            LA = M + 1 - LM                                             00002350
            LB = K - LM                                                 00002360
            T = -Z(K)                                                   00002370
            CALL DAXPY(LM,T,ABD(LA,K),1,Z(LB),1)                        00002380
  170    CONTINUE                                                       00002390
C        MAKE ZNORM = 1.0                                               00002400
         S = 1.0D0/DASUM(N,Z,1)                                         00002410
         CALL DSCAL(N,S,Z,1)                                            00002420
         YNORM = S*YNORM                                                00002430
C                                                                       00002440
         IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM                      00002450
         IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0                            00002460
  180 CONTINUE                                                          00002470
      RETURN                                                            00002480
      END                                                               00002490
