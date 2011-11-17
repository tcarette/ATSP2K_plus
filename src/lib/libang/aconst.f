*
*     -------------------------------------------------------------
*      A C O N S T
*     -------------------------------------------------------------
*                                                                  *
*                                                                  *
*     Written by G. Gaigalas,                                      * 
*     Vilnius, LITHUANIA                              January 1997 *
*
      SUBROUTINE ACONST(IQ,IP,IR,IS,IT,A)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      AI=DBLE(IR+IQ-IP)*DBLE(IP-IR+IQ)*DBLE(IP+IR-IQ+2)*
     :DBLE(IP+IR+IQ+2)*DBLE(IT+IQ-IS)*DBLE(IS-IT+IQ)*
     :DBLE(IS+IT-IQ+2)*DBLE(IS+IT+IQ+2)
      A=DSQRT(AI/DBLE(256))
      RETURN
      END
