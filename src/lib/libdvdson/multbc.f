*=======================================================================
        SUBROUTINE MULTBC(N,K,M,C,TEMP,B)
*=======================================================================
*       called by: DVDRVR
*
*       Multiplies B(N,K)*C(K,M) and stores it in B(N,M)
*       Used for collapsing the expanding basis to current estimates,
*       when basis becomes too large, or for returning the results back

*       Subroutines called
*       DINIT, DGEMV, DCOPY
*-----------------------------------------------------------------------
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        DIMENSION B(N*K),C(K*M),TEMP(M)
*-----------------------------------------------------------------------
        DO 10 IROW=1,N
*              CALL DINIT(M,0.d0,TEMP,1)
           CALL DGEMV('Transp',K,M, 1.D0, C,K,B(IROW),N, 0.D0 ,TEMP,1)
           CALL DCOPY(M,TEMP,1,B(IROW),N)
  10    CONTINUE

        RETURN
        END
