*     ------------------------------------------------------------
*     I N I T M
*     ------------------------------------------------------------
*
*     Determine the Rydberg constant, RMASS constant, and
*     initialize the radial grid
*
      SUBROUTINE INITM(IPROB)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=60,NOD=220)
      PARAMETER (NELEMENT=50)


****************************************************************
*
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc,iuw,ioul,iouj,
     ;             iouhn,iouhz,iouhs,iouhm,iLS,idisk
*
      COMMON /PARAM/ H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,
     :ID,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
*
      COMMON /RYDBERG/ RYDBRG,ZMU,ENERGY
*
      LOGICAL CORE,YM,YF
      COMMON /COR/COREF,COREG,CORER,CORES,ZM1,ZM2,FZE,FZA,DR2,CORE,YM,YF
*
      CHARACTER YES
      DOUBLE PRECISION ME
      DOUBLE PRECISION MN(nelement), FZ(NELEMENT-9)
*
      DATA ME/548.579903D-6/,RY/109737.31534/ MN/1,4,7,9,11,12,14,16,
     : 19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,58,
     : 63,64,69,74,75,80,79,84,85,88,89,90,93,98,99,102,103,108,
     : 107,114,115,120/
*
      DATA FZ/163.5,198.4,239.6,283.9,332.7,386.1,444.2,507.8,
     : 576.4,650.7,731.2,817.9,911.2,1011.4,1119.3,1234.9,1359.3,
     : 1491.9,1634.4,1786.6,1948.7,2122.7,2308.2,2507.4,2718.5,
     : 2944.3,3187.5,3442.2,3720.3,4015.4,4326.6,4664.8,5022.4,
     : 5405.9,5811.9,6245.0,6706.8,7201.3,7728.7,8296.5,8889.3/
*

*
* ****   DETERMINE THE RYDBERG CONSTANT AND F(Z) PARAMETER 
*
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      myid = 0; nprocs = 1;
	if (myid.eq.0) then

      WRITE(iscw,'(/1X,A)') 'Default Rydberg constant (y/n)'
      READ(iread,'(A1)') YES
      IF (YES.EQ.'Y'.OR.YES.EQ.'y'.AND.(INT(Z).LE.NELEMENT)) THEN
         ZMU=MN(INT(Z))
      ELSE
         WRITE(iscw,'(1X,A)') 'Enter the mass of the atom'
         READ(iread,*) ZMU
      END IF
      RYDBRG = RY/(1.+ME/ZMU)
      IF (IPROB.NE.1) THEN
         IF (INT(Z).LT.10) THEN
            FZE=1.566432
         ELSEIF (INT(Z).LE.NELEMENT) THEN
            FZE=FZ(INT(Z-9))
         ELSE
            WRITE(iscw,'(1X,A)') 'Enter the f(Z) value in MHz'
            READ(iread,*) FZE
         ENDIF
       ENDIF
      ENDIF
	  
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

       RMASS = ME/ZMU
*
* ****  INITIALIZE THE RADIAL GRID
*
      RHO=-4.0
      DO 10 J=1,NO
         R(J)=DEXP(RHO)/Z
         RR(J)=R(J)*R(J)
         R2(J)=DSQRT(R(J))
         RHO=RHO+H
10    CONTINUE
      RETURN
      END
