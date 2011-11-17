*
*     -------------------------------------------------------------
*      S I X J
*     -------------------------------------------------------------
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF 6j COEFFICIENT         *
*                                                                  *
*     | I/2  J/2  K/2 |                                            *
*     | L/2  M/2  N/2 |          (29.1A) [J.B.77]                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
*
      SUBROUTINE SIXJ(I,J,K,L,M,N,ITIK,SI)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      LOGICAL SAVE
      DIMENSION RACA(0:4,0:4,0:4,0:4,0:4,0:4)
      DATA RACA/15625*1.D-20/
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DATA UNDEF/1.D-20/ 
      SI=ZERO
      IF(ITIK.NE.0) THEN
C
C     CHESKED TRIANGULAR CONDITIONS
C
        IF(IXJTIK(I,J,K,L,M,N).EQ.0) RETURN
      ENDIF
      SAVE = .FALSE.
      IF (MAX0(I,J,K,L,M,N).LE.4) THEN
	 SI = RACA(I,J,K,L,M,N)
	 IF (SI .EQ. UNDEF) THEN
	    SAVE = .TRUE.
	 ELSE
	    RETURN
	 END IF
      END IF
C
C     CALCULATED IN CASE WHEN ONE OF PERAMETERS EQUAL 0
C
      IF(I*J*K*L*M*N.EQ.0) THEN
        IF(I.EQ.0) THEN 
          A=DBLE((M+1)*(K+1))
          IFA=L+M+K
        ELSEIF(J.EQ.0) THEN
          A=DBLE((L+1)*(K+1))
          IFA=I+M+N
        ELSEIF(K.EQ.0) THEN
          A=DBLE((I+1)*(L+1))
          IFA=I+M+N
        ELSEIF(L.EQ.0) THEN
          A=DBLE((J+1)*(K+1))
          IFA=I+J+K
        ELSEIF(M.EQ.0) THEN
          A=DBLE((I+1)*(K+1))
          IFA=I+J+K
        ELSE
          A=DBLE((I+1)*(J+1))
          IFA=I+J+K
        ENDIF
        SI=ONE/DSQRT(A)
        IF(MOD(IFA,4).NE.0)SI=-SI
C
C     THE CASE 1/2
C
      ELSEIF(MIN0(I,J,K,L,M,N).EQ.1) THEN
        IF(I.EQ.1) THEN
          CALL SIXJ5(M,K,L,J,N,0,SI)
        ELSEIF(J.EQ.1) THEN
          CALL SIXJ5(I,N,M,L,K,0,SI)
        ELSEIF(K.EQ.1) THEN
          CALL SIXJ5(I,M,N,L,J,0,SI)
        ELSEIF(L.EQ.1) THEN
          CALL SIXJ5(J,K,I,M,N,0,SI)
        ELSEIF(M.EQ.1) THEN
          CALL SIXJ5(I,K,J,L,N,0,SI)
        ELSE
          CALL SIXJ5(I,J,K,L,M,0,SI)
        ENDIF
C
C     THE CASE 1
C
      ELSEIF(MIN0(I,J,K,L,M,N).EQ.2) THEN
        IF(I.EQ.2) THEN
          CALL SIXJ1(M,K,L,J,N,0,SI)
        ELSEIF(J.EQ.2) THEN
          CALL SIXJ1(I,N,M,L,K,0,SI)
        ELSEIF(K.EQ.2) THEN
          CALL SIXJ1(I,M,N,L,J,0,SI)
        ELSEIF(L.EQ.2) THEN
          CALL SIXJ1(J,K,I,M,N,0,SI)
        ELSEIF(M.EQ.2) THEN
          CALL SIXJ1(I,K,J,L,N,0,SI)
        ELSE
          CALL SIXJ1(I,J,K,L,M,0,SI)
        ENDIF
C
C     THE CASE 3/2 
C
      ELSEIF(MIN0(I,J,K,L,M,N).EQ.3) THEN
        IF(I.EQ.3) THEN
          CALL SIXJ35(M,K,L,J,N,0,SI)
        ELSEIF(J.EQ.3) THEN
          CALL SIXJ35(I,N,M,L,K,0,SI)
        ELSEIF(K.EQ.3) THEN
          CALL SIXJ35(I,M,N,L,J,0,SI)
        ELSEIF(L.EQ.3) THEN
          CALL SIXJ35(J,K,I,M,N,0,SI)
        ELSEIF(M.EQ.3) THEN
          CALL SIXJ35(I,K,J,L,N,0,SI)
        ELSE
          CALL SIXJ35(I,J,K,L,M,0,SI)
        ENDIF
C
C     THE CASE 2
C
      ELSEIF(MIN0(I,J,K,L,M,N).EQ.4) THEN
        IF(I.EQ.4) THEN
          CALL SIXJ2(M,K,L,J,N,0,SI)
        ELSEIF(J.EQ.4) THEN
          CALL SIXJ2(I,N,M,L,K,0,SI)
        ELSEIF(K.EQ.4) THEN
          CALL SIXJ2(I,M,N,L,J,0,SI)
        ELSEIF(L.EQ.4) THEN
          CALL SIXJ2(J,K,I,M,N,0,SI)
        ELSEIF(M.EQ.4) THEN
          CALL SIXJ2(I,K,J,L,N,0,SI)
        ELSE
          CALL SIXJ2(I,J,K,L,M,0,SI)
        ENDIF
C
C     THE CASE 5/2
C
      ELSEIF(MIN0(I,J,K,L,M,N).EQ.5) THEN
        CALL GRACAH1(I,J,M,L,K,N,SI)
        IF(MOD(I+J+M+L,4).NE.0)SI=-SI
C
C     CASES 3
C
      ELSEIF(MIN0(I,J,K,L,M,N).EQ.6) THEN
        IF(I.EQ.6) THEN
          CALL SIXJ3(M,K,L,J,N,0,SI)
        ELSEIF(J.EQ.6) THEN
          CALL SIXJ3(I,N,M,L,K,0,SI)
        ELSEIF(K.EQ.6) THEN
          CALL SIXJ3(I,M,N,L,J,0,SI)
        ELSEIF(L.EQ.6) THEN
          CALL SIXJ3(J,K,I,M,N,0,SI)
        ELSEIF(M.EQ.6) THEN
          CALL SIXJ3(I,K,J,L,N,0,SI)
        ELSE
          CALL SIXJ3(I,J,K,L,M,0,SI)
        ENDIF
C
C     THE CASE 7/2
C
      ELSEIF(MIN0(I,J,K,L,M,N).EQ.7) THEN
        CALL GRACAH1(I,J,M,L,K,N,SI)
        IF(MOD(I+J+M+L,4).NE.0)SI=-SI
C
C     CASES 4
C
      ELSEIF(MIN0(I,J,K,L,M,N).EQ.8) THEN
        IF(I.EQ.8) THEN
          CALL SIXJ4(M,K,L,J,N,0,SI)
        ELSEIF(J.EQ.8) THEN
          CALL SIXJ4(I,N,M,L,K,0,SI)
        ELSEIF(K.EQ.8) THEN
          CALL SIXJ4(I,M,N,L,J,0,SI)
        ELSEIF(L.EQ.8) THEN
          CALL SIXJ4(J,K,I,M,N,0,SI)
        ELSEIF(M.EQ.8) THEN
          CALL SIXJ4(I,K,J,L,N,0,SI)
        ELSE
          CALL SIXJ4(I,J,K,L,M,0,SI)
        ENDIF
C
C     THE CASE 9/2
C
      ELSEIF(MIN0(I,J,K,L,M,N).EQ.9) THEN
        CALL GRACAH1(I,J,M,L,K,N,SI)
        IF(MOD(I+J+M+L,4).NE.0)SI=-SI
C
C     CALCULATED OTHER CASES
C
      ELSE
C       CALL DRACAH(I,J,M,L,K,N,SI)
        CALL GRACAH1(I,J,M,L,K,N,SI)
C        CALL GRACAH(I,J,M,L,K,N,SI)
        IF(MOD(I+J+M+L,4).NE.0)SI=-SI
      ENDIF
      IF (SAVE) RACA(I,J,K,L,M,N) =SI
      RETURN
      END
