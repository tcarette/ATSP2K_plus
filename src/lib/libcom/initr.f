*
*     ---------------------------------------------------------------
*        I N I T R
*     ---------------------------------------------------------------
*
*
      SUBROUTINE INITR
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220)
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
* *****  SET THE COMMONLY USED DOUBLE PRECISION CONSTANTS
*
      D0  =  0.D0
      D1  =  1.D0
      D2  =  2.D0
      D3  =  3.D0
      D4  =  4.D0
      D5  =  1.D0/2.D0
      D6  =  6.D0
      D8  =  8.D0
      D10 = 10.D0
      D12 = 12.D0
      D16 = 16.D0
      D30 = 30.D0
*
* *****  SET THE STARTING POINT, STEP SIZE, AND RELATED PARAMETERS
*
      RHO = -4.D0
      H   = 1./16.D0
      H1 = H/1.5
      H3 = H/3.
      CH = H*H/12.
      EH = DEXP(-H)
      NO=NOD
      ND = NO - 2
*
* *****  SET THE FINE-STRUCTURE CONSTANT
*
      FINE = 0.25D0/(137.036D0)**2
      RETURN
      END
