 
*     ------------------------------------------------------------------
*       O U T L S J
*     ------------------------------------------------------------------
*
      SUBROUTINE OUTLSJ(NC,C,IPACK,JA,JB,IPTR,IPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON /FOUT/NF,NG,NR,NL,NZ,NN,NV,NS,IFLAG,NIJ,IDIM
      COMMON/INFORM/IREAD,IWRITE,IOUT,ISC(4),JSC(3),IALL,ISCW
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      PARAMETER (LSTACK=128)
      DIMENSION C(NC),IPACK(NC),JA(NC),JB(NC),
     :          IPTR(NC),IPT(NC)
      INTEGER NCOUNT(8),ISTACK(LSTACK),II(4),IEL(4)
      EQUIVALENCE (NCOUNT(1),NF),(II(1),I1),(II(2),I2),(II(3),I3),
     :            (II(4),I4)
      CHARACTER*1 INT(8)
      DATA INT/'F','G','R','L','Z','N','V','S'/
*
*     Formats for Integrals
*
   10 FORMAT(1X,A1,I2,'(',A3,',',A3,')',I6)
   11 FORMAT(1X,A1,I2,'(',2A3,',',2A3,')',I6)
   12 FORMAT(1X,A1,2X,'(',A3,',',A3,')',I6)
*
*     Format for Coefficients
*
   20 FORMAT(F14.8,A1,2I4,I3)
*
*     Format for Coefficient terminator
*
   30 FORMAT(14X,'*',16X,' ')
   31 FORMAT(1X,'*')
*
* --- TERMINATE  INTEGERAL LISTS and Rewind
*
      DO 1 J = 1,8
         ENDFILE(UNIT=ISC(J))
         REWIND(UNIT=ISC(J))
    1 CONTINUE
*
*===  Begin processing data
*
      NINT = 0
      DO 100 ICASE = 1,8
         N = NCOUNT(ICASE)
         IF (ICASE .NE. 3 .AND. ICASE .NE. 4) THEN
            DO 102 J = 1,N
              READ(ISC(ICASE)) C(J),IPACK(J),JA(J),JB(J)
  102       CONTINUE
         ELSE
            DO 104 J = 1,N
               READ(ISC(ICASE)) C(J),IPACK(J),JA(J),JB(J),IPTR(J)
  104       CONTINUE
         END IF
         CALL QSORT(N,IPACK,IPT,ISTACK,LSTACK,IERR)
         IF (IERR .EQ. 1) THEN
            WRITE(ISCW,*) ' Stack dimension not large enough for sort'
            CALL EXIT(1)
         END IF
*
*        Output the list of integrals with pointers to the data
*
         LAST = 0
  110    J = LAST +1
         LAST = J
         IF (J .LE. N) THEN
*
*          Unpack electron data
*
           K = IPACK(IPT(J))
ctc  NWD limited to 63 > 94
ctc           I4 = MOD(K,64)
ctc           K = K/64
           I4 = MOD(K,95)
           K = K/95
ctc
ctc  NWD limited to 63 > 94
           IF (ICASE.LE.2 .OR. ICASE.EQ.4 .OR. ICASE.EQ.5) THEN
ctc             I2 = MOD(K,64)
ctc             K = K/64
             I2 = MOD(K,95)
             K = K/95
           ELSE
ctc              I3 = MOD(K,64)
ctc              K = K/64
ctc              I2 = MOD(K,64)
ctc              K = K/64
ctc              I1 = MOD(K,64)
ctc              K = K/64
              I3 = MOD(K,95)
              K = K/95
              I2 = MOD(K,95)
              K = K/95
              I1 = MOD(K,95)
              K = K/95
              IF (ICASE .GT. 5) K = K - 1
           END IF
ctc
*
*           Find  last item in the list with this integral
*
  120      LAST = LAST + 1
           IF (LAST .LE. N) THEN
             IF (IPACK(IPT(J)) .EQ. IPACK(IPT(LAST))) GO TO 120
           END IF
           LAST = LAST -1
           NINT = NINT + 1
           IF (ICASE .LE. 2) THEN
             WRITE(IOUT,10) INT(ICASE),K,IAJCMP(I2),IAJCMP(I4),LAST
           ELSE IF (ICASE .EQ. 4 .OR. ICASE .EQ. 5) THEN
             WRITE(IOUT,12) INT(ICASE),IAJCMP(I2),IAJCMP(I4),LAST
           ELSE
             DO 140 J = 1,4
               IF (II(J) .GE. 44 .AND. IALL .EQ. 0) THEN
                 IEL(J) = IAJCLD(64-II(J))
               ELSE
                 IEL(J) = IAJCMP(II(J))
               END IF
  140        CONTINUE
             WRITE(IOUT,11) INT(ICASE),K,(IEL(J),J=1,4),LAST
           END IF
           GO TO 110
        END IF
        WRITE(IOUT,31)
*
*             Write out the data for the integrals
*
        DO 150 J = 1,N
          K = IPT(J)
          WRITE(IOUT,20) C(K),INT(ICASE),JA(K),JB(K)
  150   CONTINUE
        WRITE(IOUT,30)
        IF (ICASE .EQ. 2) THEN
          WRITE(IOUT,31)
          WRITE(IOUT,31)
        END IF
  100 CONTINUE
      WRITE(ISCW,*) 'The total number of integrals =',NINT
      IDIM = NINT
      RETURN
      END
