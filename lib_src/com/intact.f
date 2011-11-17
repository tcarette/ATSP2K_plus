*
*     ------------------------------------------------------------------
*	I N T A C T
*     ------------------------------------------------------------------
*
      SUBROUTINE INTACT(L,LP,IEQUIV)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/ENAV/NINTS,KVALUE(15),COEFCT(15)
      COMMON /EAV/CCA(10),CCB(35)
      INTEGER ICBPTR(0:4)
      DATA ICBPTR/1,6,14,23,31/
*
*     THIS SUBROUTINE GIVES THE INTERACTION ENERGY BETWEEN TWO SHELLS,
*     ONE WITH ORBITAL ANGULAR MOMENTUM  L , THE OTHER WITH ORBITAL
*     ANGULAR MOMENTUM  LP .   NOTICE THAT THE FIRST TERM OF THIS
*     INTERACTION ENERGY IS ALWAYS   F0(L,LP)   AND THIS IS NOT GIVEN
*     IN THIS SUBROUTINE.   THUS ONLY THE EXTRA TERMS ARE HERE PRODUCED.
*     FOR EQUIVALENT ELECTRONS (IEQUIV = 1) ,  THERE WILL BE  FK
*     INTEGRALS ONLY.   FOR NON-EQUIVALENT ELECTRONS (IEQUIV = 2) ,
*     THERE WILL BE  GK  INTEGRALS ONLY.
*
*     THE EXPRESSIONS FOR THE INTERACTION ENERGIES ARE GIVEN BY
*     R. D. COWAN, THE THEORY OF ATOMIC SPECTRA, EQUATIONS (6.38)
*      AND (6.39).
*
      I = 0
      IF (IEQUIV .EQ. 1) THEN
          DO 1 K = 2,2*L,2
             I = I+1
             KVALUE(I) = K
             IF (L .LE. 4) THEN
                COEFCT(I) = -CCA((L*(L-1) + K)/2)
             ELSE
               COEFCT(I) = -RME(L,L,K)**2/((2*L+1)*(4*L+1))
             END IF
    1     CONTINUE
      ELSE
          DO 2 K = IABS(L-LP),L+LP,2
             I = I+1
             KVALUE(I) = K
             IF (L .LE. LP) THEN
                L1 = L
                 L2 = LP
              ELSE
                 L1 = LP
                 L2 = L
             END IF
             IF ( L2 .LE. 4) THEN
                COEFCT(I) = -CCB(ICBPTR(L1)+(K+L1-L2)/2+(L1+1)*(L2-L1))
               ELSE
                COEFCT(I) = -RME(L,LP,K)**2/(2*(2*L+1)*(2*LP+1))
             END IF
    2     CONTINUE
      END IF
      NINTS = I
      RETURN
      END
