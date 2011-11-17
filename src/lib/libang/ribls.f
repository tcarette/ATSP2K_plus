*
*     -------------------------------------------------------------
*      B L O C K   D A T A   R I B L S
*     -------------------------------------------------------------
*
*     Written by G. Gaigalas,                                      * 
*     Vanderbilt University,  Nashville             October 1996   * 
*
      BLOCK DATA RIBLS
      COMMON/RIBOLS/IMPTLS(40),IMGTLS(40),IMPNLS(40),IMGNLS(40)
      COMMON/RIBOLSF/IMPTLSF(8),IMGTLSF(8),IMPNLSF(8),IMGNLSF(8)
      COMMON/RIBOF/IMPTF(238),IMGTF(238),IMPNF(238),IMGNF(238)
      COMMON/RIBOLS3/IMPTLS3(90),IMGTLS3(90),IMPNLS3(90),IMGNLS3(90)
      DATA IMPTLS/1,2,3*3,3*6,16*9,16*25/
      DATA IMGTLS/1,2,3*5,3*8,16*24,16*40/
      DATA IMPNLS/2,1,3*6,3*3,16*25,16*9/
      DATA IMGNLS/2,1,3*8,3*5,16*40,16*24/
      DATA IMPTLSF/301,7*302/
      DATA IMGTLSF/301,7*308/
      DATA IMPNLSF/302,7*301/
      DATA IMGNLSF/308,7*301/
      DATA IMPTF/119*1,119*120/
      DATA IMGTF/119*119,119*238/
      DATA IMPNF/119*120,119*1/
      DATA IMGNF/119*238,119*119/
      DATA IMPTLS3/1,9*2,11,11*12,23,13*24,37,15*38,53,17*54,71,
     *19*72/
      DATA IMGTLS3/1,9*10,11,11*22,23,13*36,37,15*52,53,17*70,71,
     *19*90/
      DATA IMPNLS3/2,9*1,12,11*11,24,13*23,38,15*37,54,17*53,72,
     *19*71/
      DATA IMGNLS3/10,9*1,22,11*11,36,13*23,52,15*37,70,17*53,90,
     *19*71/
      END
