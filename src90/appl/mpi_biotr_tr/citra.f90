!
      SUBROUTINE CITRA(CIIN, NCSF, NCIV, L, NSHL, T, LBUF, LU, NIN, &
           CIOUT, SCR, BUF, IBUF, NTESTG) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  19:47:46  11/18/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE scalve_I 
      USE copvec_I 
      USE tiini_I 
      USE ti1tv_I 
      USE vecsum_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NCSF 
      INTEGER  :: NCIV 
      INTEGER  :: L 
      INTEGER  :: NSHL 
      INTEGER  :: LBUF 
      INTEGER  :: LU 
      INTEGER , INTENT(IN) :: NIN 
      INTEGER  :: NTESTG 
      INTEGER  :: IBUF(4,LBUF) 
      REAL(DOUBLE)  :: CIIN(NCSF,NCIV) 
      REAL(DOUBLE)  :: T(NSHL,NSHL) 
      REAL(DOUBLE)  :: CIOUT(NCSF,NCIV) 
      REAL(DOUBLE)  :: SCR(NCSF,NCIV) 
      REAL(DOUBLE)  :: BUF(LBUF) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IIN, IPOT, I, N, JLEI 
      REAL(DOUBLE) :: FACTOR, TII, XNFACI 
!-----------------------------------------------
!. Input
!. Output
!. Scratch
!
      IF (NIN /= 0) THEN 
         FACTOR = 1.0D0 
         DO IIN = 1, NIN 
            FACTOR = FACTOR*T(IIN,IIN) 
         END DO 
         IPOT = 2*(2*L + 1) 
         FACTOR = FACTOR**IPOT 
         CALL SCALVE (CIIN, FACTOR, NCIV*NCSF) 
      ENDIF 
      IF (NIN == NSHL) CALL COPVEC (CIIN, CIOUT, NCIV*NCSF) 
!
      DO I = NIN + 1, NSHL 
!
         TII = T(I,I) 
         CALL TIINI (CIIN, NCSF, NCIV, I, L, TII, LBUF, LU, CIOUT, BUF, IBUF, &
            NTESTG) 
         CALL COPVEC (CIOUT, SCR, NCIV*NCSF) 
!
!. Off diagonal contributions
!
         XNFACI = 1.0D0 
         DO N = 1, 4*L + 2 
! T ** (N-1) is supposed to be in SCR, copy to CIIN
! and apply S
            CALL COPVEC (SCR, CIIN, NCIV*NCSF) 
            JLEI = 1 
            CALL TI1TV (CIIN, NCSF, NCIV, I, L, T(1,I), NSHL, JLEI, LBUF, LU, &
               SCR, BUF, IBUF, NTESTG) 
            XNFACI = XNFACI/FLOAT(N) 
            CALL VECSUM (CIOUT, CIOUT, SCR, 1.0D0, XNFACI, NCIV*NCSF) 
         END DO 
         CALL COPVEC (CIOUT, CIIN, NCIV*NCSF) 
!
      END DO 
!
      RETURN  
      END SUBROUTINE CITRA 
