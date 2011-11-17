* 
*     ------------------------------------------------------------------ 
*             E P T R
*     ------------------------------------------------------------------ 
* 
*       Determines the position of the electron in the electron list 
* 
      SUBROUTINE EPTR(EL,ELSYMB, IEL, *) 
      IMPLICIT DOUBLE PRECISION(A-H,O-Z) 
      PARAMETER(IWRITE=6)
      CHARACTER CONFIG*66, COUPLE*3, EL(*)*3, ELSYMB*3, BL*3 
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      DATA BL/'   '/ 
* 
* ***** SEARCH ELECTRON LIST FOR LSYMB 
* 
      IF ( ELSYMB .EQ. BL ) THEN 
         IEL = 0 
         RETURN 
      ENDIF 
      DO 10 I=1,NWF 
         IF (EL(I) .EQ. ELSYMB ) THEN 
            IEL = I 
            RETURN 
         ENDIF 
10    CONTINUE 
      IEL = -1 
      WRITE (IWRITE,20) ELSYMB 
20    FORMAT(/10X,A3,' NOT FOUND IN ELECTRON LIST') 
      RETURN 1 
      END 
