      SUBROUTINE DPBFA(ABD,LDA,N,M,INFO)                                00000010
      INTEGER LDA,N,M,INFO                                              00000020
      DOUBLE PRECISION ABD(LDA,1)                                       00000030
C                                                                       00000040
C     DPBFA FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE      00000050
C     MATRIX STORED IN BAND FORM.                                       00000060
C                                                                       00000070
C     DPBFA IS USUALLY CALLED BY DPBCO, BUT IT CAN BE CALLED            00000080
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.          00000090
C                                                                       00000100
C     ON ENTRY                                                          00000110
C                                                                       00000120
C        ABD     DOUBLE PRECISION(LDA, N)                               00000130
C                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER   00000140
C                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE      00000150
C                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE      00000160
C                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS.     00000170
C                                                                       00000180
C        LDA     INTEGER                                                00000190
C                THE LEADING DIMENSION OF THE ARRAY  ABD .              00000200
C                LDA MUST BE .GE. M + 1 .                               00000210
C                                                                       00000220
C        N       INTEGER                                                00000230
C                THE ORDER OF THE MATRIX  A .                           00000240
C                                                                       00000250
C        M       INTEGER                                                00000260
C                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.       00000270
C                0 .LE. M .LT. N .                                      00000280
C                                                                       00000290
C     ON RETURN                                                         00000300
C                                                                       00000310
C        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND         00000320
C                FORM, SO THAT  A = TRANS(R)*R .                        00000330
C                                                                       00000340
C        INFO    INTEGER                                                00000350
C                = 0  FOR NORMAL RETURN.                                00000360
C                = K  IF THE LEADING MINOR OF ORDER  K  IS NOT          00000370
C                     POSITIVE DEFINITE.                                00000380
C                                                                       00000390
C     BAND STORAGE                                                      00000400
C                                                                       00000410
C           IF  A  IS A SYMMETRIC POSITIVE DEFINITE BAND MATRIX,        00000420
C           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT.        00000430
C                                                                       00000440
C                   M = (BAND WIDTH ABOVE DIAGONAL)                     00000450
C                   DO 20 J = 1, N                                      00000460
C                      I1 = MAX0(1, J-M)                                00000470
C                      DO 10 I = I1, J                                  00000480
C                         K = I-J+M+1                                   00000490
C                         ABD(K,J) = A(I,J)                             00000500
C                10    CONTINUE                                         00000510
C                20 CONTINUE                                            00000520
C                                                                       00000530
C     LINPACK.  THIS VERSION DATED 08/14/78 .                           00000540
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.      00000550
C                                                                       00000560
C     SUBROUTINES AND FUNCTIONS                                         00000570
C                                                                       00000580
C     BLAS DDOT                                                         00000590
C     FORTRAN MAX0,DSQRT                                                00000600
C                                                                       00000610
C     INTERNAL VARIABLES                                                00000620
C                                                                       00000630
      DOUBLE PRECISION DDOT,T                                           00000640
      DOUBLE PRECISION S                                                00000650
      INTEGER IK,J,JK,K,MU                                              00000660
C     BEGIN BLOCK WITH ...EXITS TO 40                                   00000670
C                                                                       00000680
C                                                                       00000690
         DO 30 J = 1, N                                                 00000700
            INFO = J                                                    00000710
            S = 0.0D0                                                   00000720
            IK = M + 1                                                  00000730
            JK = MAX0(J-M,1)                                            00000740
            MU = MAX0(M+2-J,1)                                          00000750
            IF (M .LT. MU) GO TO 20                                     00000760
            DO 10 K = MU, M                                             00000770
               T = ABD(K,J) - DDOT(K-MU,ABD(IK,JK),1,ABD(MU,J),1)       00000780
               T = T/ABD(M+1,JK)                                        00000790
               ABD(K,J) = T                                             00000800
               S = S + T*T                                              00000810
               IK = IK - 1                                              00000820
               JK = JK + 1                                              00000830
   10       CONTINUE                                                    00000840
   20       CONTINUE                                                    00000850
            S = ABD(M+1,J) - S                                          00000860
C     ......EXIT                                                        00000870
            IF (S .LE. 0.0D0) GO TO 40                                  00000880
            ABD(M+1,J) = DSQRT(S)                                       00000890
   30    CONTINUE                                                       00000900
         INFO = 0                                                       00000910
   40 CONTINUE                                                          00000920
      RETURN                                                            00000930
      END                                                               00000940
