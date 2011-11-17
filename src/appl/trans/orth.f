*
*     ------------------------------------------------------------------
*	O R T H
*
*     ------------------------------------------------------------------
*
      SUBROUTINE orth
*
*	Determine the orthogonality between initial and final state
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NWD=80)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,MAXORB
      POINTER (QIORTH,IORTH(1))
      COMMON/OVRLAP/MU,NU,MUP,NUP,NONORT,NOVLPS,IROWMU,IROWNU,ICOLMU,
     1 ICOLNU,NORTH,IORDER,NCALLS,LMU,LNU,LMUP,LNUP,JMU,JNU,JMUP,JNUP,
     2 QIORTH

      COMMON /NOR/NCOM,NCLOSI,NCLOSF,NORBI,NORBF,IWAR
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
      IF (NORTH .NE. 0) THEN
        WRITE(ISCW,'(A)') ' this prog. com. with orthogonal orbitals'
        STOP
      ENDIF
      RETURN
      END
