      SUBROUTINE DGEDI(A,LDA,N,IPVT,DET,WORK,JOB)                       00000010
      INTEGER LDA,N,IPVT(1),JOB                                         00000020
      DOUBLE PRECISION A(LDA,1),DET(2),WORK(1)                          00000030
C                                                                       00000040
C     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX            00000050
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.                     00000060
C                                                                       00000070
C     ON ENTRY                                                          00000080
C                                                                       00000090
C        A       DOUBLE PRECISION(LDA, N)                               00000100
C                THE OUTPUT FROM DGECO OR DGEFA.                        00000110
C                                                                       00000120
C        LDA     INTEGER                                                00000130
C                THE LEADING DIMENSION OF THE ARRAY  A .                00000140
C                                                                       00000150
C        N       INTEGER                                                00000160
C                THE ORDER OF THE MATRIX  A .                           00000170
C                                                                       00000180
C        IPVT    INTEGER(N)                                             00000190
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.                  00000200
C                                                                       00000210
C        WORK    DOUBLE PRECISION(N)                                    00000220
C                WORK VECTOR.  CONTENTS DESTROYED.                      00000230
C                                                                       00000240
C        JOB     INTEGER                                                00000250
C                = 11   BOTH DETERMINANT AND INVERSE.                   00000260
C                = 01   INVERSE ONLY.                                   00000270
C                = 10   DETERMINANT ONLY.                               00000280
C                                                                       00000290
C     ON RETURN                                                         00000300
C                                                                       00000310
C        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.               00000320
C                OTHERWISE UNCHANGED.                                   00000330
C                                                                       00000340
C        DET     DOUBLE PRECISION(2)                                    00000350
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.           00000360
C                OTHERWISE NOT REFERENCED.                              00000370
C                DETERMINANT = DET(1) * 10.0**DET(2)                    00000380
C                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0                  00000390
C                OR  DET(1) .EQ. 0.0 .                                  00000400
C                                                                       00000410
C     ERROR CONDITION                                                   00000420
C                                                                       00000430
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS     00000440
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.           00000450
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY      00000460
C        AND IF DGECO HAS SET RCOND .GT. 0.0 OR DGEFA HAS SET           00000470
C        INFO .EQ. 0 .                                                  00000480
C                                                                       00000490
C     LINPACK. THIS VERSION DATED 08/14/78 .                            00000500
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      00000510
C                                                                       00000520
C     SUBROUTINES AND FUNCTIONS                                         00000530
C                                                                       00000540
C     BLAS DAXPY,DSCAL,DSWAP                                            00000550
C     FORTRAN DABS,MOD                                                  00000560
C                                                                       00000570
C     INTERNAL VARIABLES                                                00000580
C                                                                       00000590
      DOUBLE PRECISION T                                                00000600
      DOUBLE PRECISION TEN                                              00000610
      INTEGER I,J,K,KB,KP1,L,NM1                                        00000620
C                                                                       00000630
C                                                                       00000640
C     COMPUTE DETERMINANT                                               00000650
C                                                                       00000660
      IF (JOB/10 .EQ. 0) GO TO 70                                       00000670
         DET(1) = 1.0D0                                                 00000680
         DET(2) = 0.0D0                                                 00000690
         TEN = 10.0D0                                                   00000700
         DO 50 I = 1, N                                                 00000710
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)                        00000720
            DET(1) = A(I,I)*DET(1)                                      00000730
C        ...EXIT                                                        00000740
            IF (DET(1) .EQ. 0.0D0) GO TO 60                             00000750
   10       IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20                       00000760
               DET(1) = TEN*DET(1)                                      00000770
               DET(2) = DET(2) - 1.0D0                                  00000780
            GO TO 10                                                    00000790
   20       CONTINUE                                                    00000800
   30       IF (DABS(DET(1)) .LT. TEN) GO TO 40                         00000810
               DET(1) = DET(1)/TEN                                      00000820
               DET(2) = DET(2) + 1.0D0                                  00000830
            GO TO 30                                                    00000840
   40       CONTINUE                                                    00000850
   50    CONTINUE                                                       00000860
   60    CONTINUE                                                       00000870
   70 CONTINUE                                                          00000880
C                                                                       00000890
C     COMPUTE INVERSE(U)                                                00000900
C                                                                       00000910
      IF (MOD(JOB,10) .EQ. 0) GO TO 150                                 00000920
         DO 100 K = 1, N                                                00000930
            A(K,K) = 1.0D0/A(K,K)                                       00000940
            T = -A(K,K)                                                 00000950
            CALL DSCAL(K-1,T,A(1,K),1)                                  00000960
            KP1 = K + 1                                                 00000970
            IF (N .LT. KP1) GO TO 90                                    00000980
            DO 80 J = KP1, N                                            00000990
               T = A(K,J)                                               00001000
               A(K,J) = 0.0D0                                           00001010
               CALL DAXPY(K,T,A(1,K),1,A(1,J),1)                        00001020
   80       CONTINUE                                                    00001030
   90       CONTINUE                                                    00001040
  100    CONTINUE                                                       00001050
C                                                                       00001060
C        FORM INVERSE(U)*INVERSE(L)                                     00001070
C                                                                       00001080
         NM1 = N - 1                                                    00001090
         IF (NM1 .LT. 1) GO TO 140                                      00001100
         DO 130 KB = 1, NM1                                             00001110
            K = N - KB                                                  00001120
            KP1 = K + 1                                                 00001130
            DO 110 I = KP1, N                                           00001140
               WORK(I) = A(I,K)                                         00001150
               A(I,K) = 0.0D0                                           00001160
  110       CONTINUE                                                    00001170
            DO 120 J = KP1, N                                           00001180
               T = WORK(J)                                              00001190
               CALL DAXPY(N,T,A(1,J),1,A(1,K),1)                        00001200
  120       CONTINUE                                                    00001210
            L = IPVT(K)                                                 00001220
            IF (L .NE. K) CALL DSWAP(N,A(1,K),1,A(1,L),1)               00001230
  130    CONTINUE                                                       00001240
  140    CONTINUE                                                       00001250
  150 CONTINUE                                                          00001260
      RETURN                                                            00001270
      END                                                               00001280
