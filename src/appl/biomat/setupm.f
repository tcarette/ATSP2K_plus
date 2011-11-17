*
*     ------------------------------------------------------------------
*       S E T U P m
*     ------------------------------------------------------------------
*
      SUBROUTINE SETUPm(ish,j1,j2,JA,JB,na,nb)
*
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH(16,2),J1QN(31,3,2),IJFUL(16)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      COMMON/RED/INDL(20),INDR(20),IORST(20),NCI1,KPL,KPR,NC,LIS(16),
     : JIST,JFST,NCONTR
*
*     Inserts the current subshell i1 into the left and right
*     coupling tree at position ish. 
*     j1, j2 : location in the configurations of the shell
*     na, nb : occupied or not: =1 not occupied, =2 occupied
*
      if (na .eq. 2) then
	I1 = nocorb(j1,ja)
      else
	I1 = nocorb(j2,jb)
      end if
*
*     This code allows a maximum of 8 open shells but at
*     most two shells may differ in the present case.  IHSH is an
*     upper bound to the number of shells.
*
      I2HSH = ihsh + ISH - 1
*
* --- FIRST CONSIDER THE L.H.S. (I=1) OF THE MATRIX ELEMENT.  NCC=1 
*     MEANS UNOCCUPIED, REPRESENTED BY A DUMMY SINGLET S SHELL, AND THE
*     ADDITIONAL SET OF COUPLING QUANTUM NUMBERS WILL BE THE SAME AS THE
*     LAST SET OF COUPLING QUANTUM NUMBERS ALREADY OBTAINED.
*     NCC=2 MEANS OCCUPIED.  THEN ALL THE NEW QUANTUM NUMBERS (BOTH FOR
*     THE SHELL AND FOR THE COUPLING OF THIS SHELL TO THE RESULTANT OF
*     THE PREVIOUS ONES) ARE DEFINED IN THE CORRESPONDING J1QNRD ARRAY.
*     NOSH - THE NUMBER OF ELECTRONS IN THIS SHELL, IS DEFINED BY THE
*     APPROPRIATE ENTRY IN NELCSH .  THE R.H.S. IS THEN CONSIDERED
*     SIMILARLY (I=2)
*
      I=1
25    if (i .eq. 1) then
	ic = j1
	jc = ja
	NCC = na
      else
	ic = j2
	jc = jb
	NCC = nb
      end if
      if (NCC .eq. 1 ) then
*       .. shell is not occupied
        NOSH(ISH,I)=0
        J1QN(ISH,1,I)=0
        J1QN(ISH,2,I)=1
        J1QN(ISH,3,I)=1
        IF(ISH .eq. 2) then
          J1QN(I2HSH,1,I)=0
          J1QN(I2HSH,2,I)=J1QN(1,2,I)
          J1QN(I2HSH,3,I)=J1QN(1,3,I)
	else if (ish .gt. 2) then
          DO 27 K=1,3
            J1QN(I2HSH,K,I)=J1QN(I2HSH-1,K,I)
   27     CONTINUE
	end if
      else
*       .. shell is occupied
        NOSH(ISH,I)=NELCSH(IC,JC)
        JD = J1QNRD(IC,JC)
        J1QN(ISH,1,I)=MOD(JD,64)
        JD = JD/64
        J1QN(ISH,2,I) = MOD(JD,64)
        J1QN(ISH,3,I) = JD/64
*
        IF (ISH .gt. 1) THEN
*
*         .. a resultant coupling is present
*
*     IS THIS THE FIRST OCCUPIED SHELL OF THIS CONFIGURATION, THOUGH NOT
*     THE FIRST OF THE OTHER CONFIGURATION.  IF SO, THE INTERMEDIATE
*     COUPLING FORMED HAS THE SAME  L,S  VALUES AS THIS OCCUPIED SHELL,
*     SINCE WE COUPLE THE SHELL TO A DUMMY SINGLET S.
*
          IF(IC .eq.1) THEN
            I2 = 1
          ELSE
            I2 = NOCCSH(JC)+IC-1
          END IF
          JD = J1QNRD(I2,JC)
          IF (IC .eq. 1) THEN
             J1QN(I2HSH,1,I) = 0
          ELSE
             J1QN(I2HSH,1,I) = MOD(JD,64)
          END IF
          JD = JD/64
          J1QN(I2HSH,2,I) = MOD(JD, 64)
          J1QN(I2HSH,3,I) = JD/64
        END IF
      end if
*
      I=I+1
      If (i .eq. 2) GO TO 25
*
* --- SET THE NJ AND LJ VALUES OF THE OCCUPIED SHELLS
*
   24 NJ(ISH)=NJCOMP(I1)
      IJFUL(ISH)=I1
      LJ(ISH)=LJCOMP(I1)
      IF (NC .NE. 0) then
        IF (INDL(NCI1).EQ.I1) KPL=ISH
	IF (INDR(NCI1).EQ.I1) KPR=ISH
      ENDIF
*
      END
