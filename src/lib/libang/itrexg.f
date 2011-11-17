*
*     ---------------------------------------------------------------
*      I T R E X G 
*     ---------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      FUNCTION ITREXG(I1,I2,I3,I4,K)
      J=MAX0(IABS(I1-I2),IABS(I3-I4))
      K=MIN0(IABS(I1+I2),IABS(I3+I4))-J+1
      ITREXG=J
      RETURN
      END
