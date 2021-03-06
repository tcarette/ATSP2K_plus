*
*     ------------------------------------------------------------------
*	C F P D
*     ------------------------------------------------------------------
*
      SUBROUTINE CFPD(N,IVI,LI,ISI,IVJ,LJ,ISJ,COEFP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
*
*     THIS SUBROUTINE EVALUATES THE COEFFICIENTS OF FRACTIONAL PARENTAGE
*     FOR EQUIVALENT D SHELL ELECTRONS FROM TABLES GIVEN IN J.C.SLATER
*     QUANTUM THEORY OF ATOMIC STRUCTURE,VOLUME2,P350(1960)
*     IN THE SUBROUTINE LIST N,THE NO.OF ELECTRONS,V THE SENIORITY QUAN
*     TUM NO.,L THE ANGULAR MOMENTUM QUANTUM NO.,(2S+1) THE SPIN QUANTUM
*     NO. OF BOTH THE STATE IN QUESTION AND ITS PARENT STATE ARE INPUT
*     PARAMETERS  THE RESULT IS OUTPUT AS COEFP
*
*
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(8)
      DIMENSION      IV(5,16),IL(5,16),IS(5,16),
     :          ITAB1(5,1),ITAB2(8,5),ITAB3(16,8),ITAB4(16,16),
     :          NORM1(5),NORM2(8),NORM3(16),NORM4(16)
      DATA IV/1,2,3,4,5,0,2,3,4,5,0,2,3,4,3,0,2,3,2,5,0,0,3,4,3,0,0,1,4,
     :   5,0,0,3,2,3,0,0,3,4,3,0,0,0,4,5,0,0,0,2,3,0,0,0,4,5,0,0,0,4,1,
     :0,0,0,2,3,0,0,0,4,5,0,0,0,0,3,0,0,0,4,5/
      DATA IL/2,3,3,2,0,0,1,1,5,4,0,4,5,4,3,0,2,4,3,2,0,0,3,3,1,0,0,2,2,
     :   6,0,0,2,1,5,0,0,1,1,4,0,0,0,6,4,0,0,0,4,3,0,0,0,4,3,0,0,0,3,2,
     :   0,0,0,2,2,0,0,0,2,2,0,0,0,0,1,0,0,0,0,0/
      DATA IS/2,3,4,5,6,0,3,4,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,0,2,3,
     :   2,0,0,2,3,2,0,0,2,3,2,0,0,0,1,2,0,0,0,1,2,0,0,0,1,2,0,0,0,1,2,
     :   0,0,0,1,2,0,0,0,1,2,0,0,0,1,2,0,0,0,1,2/
      DATA ITAB1/1,1,1,1,1/
      DATA ITAB2/4,-7,-1,21,7,-21,21,-8,-1,-8,0,0,28,-9,-49,7,0,0,1,11,
     :   -25,-9,-25,0,0,0,0,-10,-10,-5,45,15,0,0,0,0,0,16,0,0/
      DATA ITAB3/7,20,-560,224,-112,-21,-56,16,0,0,0,0,0,0,0,0,3,0,0,-56
     :   ,-448,49,-64,-14,0,0,0,0,0,0,0,0,0,26,308,110,220,0,0,0,7,-154,
     :   -28,-132,0,0,0,0,0,-9,297,90,-405,45,0,0,3,66,-507,-3,-60,15,
     :    0,0,0,5,315,-14,-175,-21,-56,-25,0,70,385,-105,28,63,0,0,0,0,
     :   0,315,0,0,135,0,0,189,0,0,105,0,1,0,0,0,200,15,120,60,-35,10,0,
     :   -25,88,200,45,20,0,1,0,0,0,16,-200,-14,-14,25,0,0,0,120,-42,42,
     :    0,0/
      DATA ITAB4/1,-105,-175,-175,-75,12*0,154,-110,0,0,231,286,924,-308
     :   ,220,-396,6*0,-66,-90,180,0,99,-99,891,-5577,-405,-9,0,45,45,0,
     :   0,0,0,224,0,-56,0,-220,1680,0,112,0,-21,21,0,-16,0,0,-70,14,-84
     :   ,56,0,55,945,4235,-175,-315,0,-21,189,-25,0,0,25,-15,-135,35,0,
     :   0,600,968,120,600,0,60,60,10,3,0,0,-56,0,-64,4*0,448,0,-9,-49,
     :   0,14,0,0,0,-16,126,14,4*0,-200,360,0,-14,126,25,0,5*0,-175,182,
     :   -728,-2184,7*0,6*0,220,880,0,-400,0,-9,-25,0,0,0,5*0,-45,-5,845
     :    ,-1215,275,495,0,-11,99,0,0,6*0,33,-7,-2541,105,-525,0,35,35,
     :   -15,0,7*0,-800 ,0,-160,0,-5,45,0,30,0,7*0,-100,1452,180,-100,0,
     :   -10,90,15,-2,11*0,6,16*0,-14,-56,0,0/
      DATA NORM1/1,1,1,1,1/
      DATA NORM2/5,15,2,42,70,60,140,30/
      DATA NORM3/10,60,1680,840,1680,210,360,90,10,504,1008,560,280,140,
     :   1,1/
      DATA NORM4/1,420,700,700,300,550,1100,8400,18480,2800,2800,50,350,
     :   700,150,5/
*
*     READ IN D SHELL PARAMETERS AND TABLES
*     PERIPHERAL 1 IS THE CARD READER
*
*     TEST IF N IS IN THE FIRST HALF OF SHELL
*
99    IF(N-6) 40,103,103
*
*     TEST IF STATE IN QUESTION IS ALLOWED
*     IF IT IS, IDENTIFY THE ROW OF THE TABLE BY J1
*
40    J = 0
101   J = J+1
      IF(J-17) 41,11,11
41    IF(IV(N,J)-IVI) 101,42,101
42    IF(IL(N,J)-LI) 101,43,101
43    IF(IS(N,J)-ISI) 101,44,101
44    J1=J
*
*     TEST IF PARENT STATE IS ALLOWED
*     IF IT IS, IDENTIFY THE COLUMN OF THE TABLE BY J2
*
      IF(N-1) 45,30,45
30    IF(IVJ) 11,31,11
31    IF(LJ) 11,32,11
32    IF(ISJ-1) 11,1,11
45    J = 0
102   J = J+1
      IF(J-17) 46,11,11
46    IF(IV(N-1,J)-IVJ) 102,47,102
47    IF(IL(N-1,J)-LJ)  102,48,102
48    IF(IS(N-1,J)-ISJ) 102,49,102
49    J2=J
      GO TO 100
*
*     SIMILAR SETTING OF J1 AND J2 IF N IS IN SECOND HALF OF SHELL
*
103   M = 10-N
      IF(M) 36,33,36
33    IF(IVI) 11,34,11
34    IF(LI) 11,35,11
35    IF(ISI-1) 11,37,11
36    J = 0
104   J = J+1
      IF(J-17) 50,11,11
50    IF(IV(M,J)-IVI) 104,51,104
51    IF(IL(M,J)-LI) 104,52,104
52    IF(IS(M,J)-ISI) 104,53,104
53    J1=J
37    J = 0
105   J = J+1
      IF(J-17) 54,11,11
54    IF(IV(M+1,J)-IVJ) 105,55,105
55    IF(IL(M+1,J)-LJ)  105,56,105
56    IF(IS(M+1,J)-ISJ) 105,57,105
57    J2=J
*
*     IDENTIFY THE F.P.C AS A UNIQUE ELEMENT OF ITABN(J1,J2)
*
100   GO TO (1,2,3,4,5,12,12,12,12,1),N
1     COEFP = 1.0D0
      GO TO 10
2     COEFP = ITAB1(J1,J2)
      IF(COEFP) 60,10,81
60    COEFP = - DSQRT(-COEFP/NORM1(J1))
      GO TO 10
81    COEFP = DSQRT(COEFP/NORM1(J1))
      GO TO 10
3     COEFP = ITAB2(J1,J2)
      IF(COEFP) 61,10,82
61    COEFP = -DSQRT(-COEFP/NORM2(J1))
      GO TO 10
82    COEFP = DSQRT(COEFP/NORM2(J1))
      GO TO 10
4     COEFP = ITAB3(J1,J2)
      IF(COEFP) 62,10,83
62    COEFP = -DSQRT(-COEFP/NORM3(J1))
      GO TO 10
83    COEFP = DSQRT(COEFP/NORM3(J1))
      GO TO 10
5     COEFP = ITAB4(J1,J2)
      IF(COEFP) 63,10,84
63    COEFP = -DSQRT(-COEFP/NORM4(J1))
      GO TO 10
84    COEFP = DSQRT(COEFP/NORM4(J1))
      GO TO 10
*
*     USE RECURRENCE RELATION EQUATION (19) OF RACAH FOR SECOND HALF OF
*     SHELL
*
12    ISIGN = (-1)**((ISI+ISJ-7)/2 +LI +LJ)
      FACTOR = DSQRT(DFLOAT((11-N)*ISJ*(2*LJ+1))/DFLOAT(N*ISI*(2*LI+1)))
      M1 =N-5
      GO TO(6,7,8,9),M1
6     COEFP = ITAB4(J2,J1)
      IF(COEFP) 64,10,85
64    COEFP = -DSQRT(-COEFP/NORM4(J2))
      GO TO 86
85    COEFP = DSQRT(COEFP/NORM4(J2))
86    COEFP = COEFP*ISIGN*FACTOR
      IF(MOD((IVJ-1)/2,2))  87,10,87
87    COEFP = -COEFP
      GO TO 10
7     COEFP = ITAB3(J2,J1)
      IF(COEFP) 65,10,88
65    COEFP = -DSQRT(-COEFP/NORM3(J2))
      GO TO 89
88    COEFP = DSQRT(COEFP/NORM3(J2))
89    COEFP = COEFP * ISIGN * FACTOR
      GO TO 10
8     COEFP = ITAB2(J2,J1)
      IF(COEFP) 66,10,90
66    COEFP = -DSQRT(-COEFP/NORM2(J2))
      GO TO 91
90    COEFP = DSQRT(COEFP/NORM2(J2))
91    COEFP = COEFP * ISIGN * FACTOR
      GO TO 10
9     COEFP = ITAB1(J2,J1)
      IF(COEFP) 67,10,92
67    COEFP = -DSQRT(-COEFP/NORM1(J2))
      GO TO 93
92    COEFP = DSQRT(COEFP/NORM1(J2))
93    COEFP = COEFP * ISIGN * FACTOR
      GO TO 10
*
*     AN UNALLOWED STATE OR AN UNALLOWED PARENT
*
11    WRITE(IWRITE,1111)
 1111 FORMAT(' ERROR IN SUBROUTINE CFPD - THE STATE OR IT''S PARENT IS N
     :OT ALLOWED')
      CALL EXIT
10    CONTINUE
      RETURN
      END
