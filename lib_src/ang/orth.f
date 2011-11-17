*
*     ------------------------------------------------------------------
*	O R T H
*
*     ------------------------------------------------------------------
*
      SUBROUTINE ORTH
*
*	Determine the orthogonality between initial and final state
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=30,NCD=2000,NCD4=4*NCD)
      COMMON/STATES/NCFG,MAXORB,IAJCMP(2*NWD),LJCOMP(2*NWD),
     :NJCOMP(2*NWD),NOCCSH(NCD4),NELCSH(5,NCD4),NOCORB(5,NCD4),
     :J1QNRD(9,NCD4)
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     1 ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     2 IORTH(NWD*NWD)

      COMMON/NOR/NCOM,NORBI,NORBF,IWAR
*
*   SET UP OF IORTH VECTOR
*   THE COMMON SET (NCOM) IS ASSUMED TO BE ORTHOGONAL TO BOTH
*   NORBI AND NORBF SETS
*
      IF (NORBI .EQ. 0) RETURN
      M1 = NORBI*NORBF
      DO 70 I = 1,M1
      IORTH(I) = 0
   70 CONTINUE
      NORTH = 0
*
*  | 1....NCOM | NCOM1.....NOR1 | NOR11.....NOR2|
*  |   NCOM    |     NORBI      |     NORBF     |
*  |   <= NWD  |     <= NWD     |     <=NWD     |
*
*  THIS LIMITATION IS LINKED TO THE DIMENSION OF BUFFER(NWD) IN
*  ANALYSE SUBROUTINE, where <= stands for LESS THAN or EQUAL.
*
      NCOM1 = NCOM+1
      NOR1  = NCOM + NORBI
      NOR11 = NOR1 + 1
      NOR2  = NOR1 + NORBF
      DO 78 J = NCOM1,NOR1
      DO 79 I = NOR11,NOR2
         IJ = NORBF*(J-NCOM1) + I - NOR1
         IF (LJCOMP(I) .EQ. LJCOMP(J)) THEN
            NORTH = NORTH + 1
            IORTH(IJ) = 1
         ENDIF
   79 CONTINUE
   78 CONTINUE
      RETURN
      END
