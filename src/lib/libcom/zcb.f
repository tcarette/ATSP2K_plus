C======================================================================
C         Z C B
C======================================================================

      Double precision FUNCTION ZCB(K1,K2,K3)
C
C     CB =  3j(k1,0,k2,0,k3,0)**2
C
      Real *8 A,B
C
      ZCB=0.0
C
      M=K1+K2+K3
      N=M/2
      IF(N+N.NE.M) RETURN
      M=M+1
C
      N1=N-K1
      N2=N-K2
      N3=N-K3
      IF(N1.LT.0.OR.N2.LT.0.OR.N3.LT.0) RETURN
      M1=N1+N1
      M2=N2+N2
      M3=N3+N3
C
      A=1.0
C
      DO 1 I=2,M
      K=-1                    ! the extent of integer I in units of 1/2
      IF(M1.GE.I) K=K+1
      IF(M2.GE.I) K=K+1
      IF(M3.GE.I) K=K+1
      IF(N .GE.I) K=K+2
      IF(N1.GE.I) K=K-2
      IF(N2.GE.I) K=K-2
      IF(N3.GE.I) K=K-2
c
      IK=IABS(K)
      B=DBLE(I)
      IF(K.GT.0) THEN
        DO 3 J=1,IK
         A=A*B
    3   Continue
      END IF
      IF(K.LT.0) THEN
        DO 2 J=1,IK
         A=A/B
    2   Continue
      END IF
c
    1   Continue
c
      ZCB=A
      RETURN
      END
