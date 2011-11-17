*
*     ---------------------------------------------------------------
*      I T R E X G 2
*     ---------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      FUNCTION ITREXG2(I1,I2,I3,I4,K)
      J=MIN0(IABS(I1-I2),IABS(I3-I4))
      K=MIN0(IABS(I1+I2),IABS(I3+I4))-J+1
      ITREXG2=J
      RETURN
      END
