*
*     ------------------------------------------------------------------
*	S E T U P  
*     ------------------------------------------------------------------
*
*	In order to keep NC in COMMON /RED/ as in other
*   routines, the role of NCC and NC has been reversed from
*   the earlier version of this code.
*
      SUBROUTINE setup(JA,JB,LET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH(16,2),J1QN(31,3,2),IJFUL(16)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      LOGICAL finisha, finishb
*
*     NOTICE THE DIFFERENT NAMES IN THE COMMON BLOCK MEDEFN  -  WE
*      STORE NOSH1(I=1,10) AS NOSH((I=1,10),1) AND NOSH2(I=1,10) AS
*     NOSH((I=1,10),2)   AND USE THE FACT THAT NOSH1 AND NOSH2 WILL THEN
*     BE EQUIVALENT TO THE SINGLE 2-DIMENSIONAL ARRAY NOSH.  SIMILARLY
*     FOR J1QN
*
* === GENERATES THE ARRAYS  NJ,LJ - DEFINING THE QUANTUM NUMBERS OF THE
*     SHELLS,   NOSH - DEFINING THE OCCUPATION OF THE SHELLS,  J1QN -
*     DEFINING THE COUPLING OF THE SHELLS,   FOR EACH OF THE TWO
*     CONFIGURATIONS CONSIDERED.    ONLY THOSE SHELLS OCCURRING IN AT
*     LEAST ONE CONFIGURATION ARE INCLUDED.
*                   AT LEAST TWO SHELLS MUST BE CONSIDERED OCCUPIED.
*     THUS (1S)**2    HELIUM  MUST BE TREATED AS ,E.G., (1S)**2(2S)**0
*     THE SIZE OF THE ARRAYS HERE CALCULATED IS ARRANGED TO BE NO
*     GREATER THAN IS NECESSARY TO INCLUDE ALL ORBITALS WHICH ARE
*     DEEMED TO BE OCCUPIED IN EITHER OR BOTH OF THE CONFIGURATIONS
*     JA,JB
*
* --- INITIALIZE BASIC QUANTITIES - (I1+1) RUNS OVER 1,MAXORB,  IHSH IS
*     THE CURRENT VALUE OF THE HIGHEST OCCUPIED SHELL YET CONSIDERED,
*     WHILE I2HSH=2*IHSH-1
*
      Let = 1
      Diff = 0
      IA=NOCCSH(JA)
      IB=NOCCSH(JB)
      if (abs(ia-ib) .gt. 1) then
        Let =0
        Return
      end if
      ipa = 1
      ipb = 1
      j1 = nocorb(ipa,ja)
      j2 = nocorb(ipb,jb)
      finisha = .false.
      finishb = .false.

      if (JA .ne. JB) then
*       ..determine which shells to advance
*      icase=1 left: =2 both; =3 right
*
10      if (finisha) then 
	  icase = 3
        else if (finishb) then
	  icase = 1
        else
	  if (j1 .lt. j2) then
	    icase = 1
	  else if (j1 .eq. j2) then
	    icase = 2
	  else
	    icase = 3
	  end if
        end if
        if (icase .eq. 1 ) then
*         .. shell occupied only in left
	  diff = diff+nelcsh(ipa,ja)
	  if (ipa .lt. ia) then
            ipa = ipa +1
	    j1 = nocorb(ipa,ja)
	  else 
	    finisha = .true.
	  end if
        else if (icase .eq. 2) then
*         .. shell occupied in both
	  diff = diff+iabs(nelcsh(ipa,ja)-nelcsh(ipb,jb))
	  if (ipa .lt. ia) then
            ipa = ipa +1
	    j1 = nocorb(ipa,ja)
	  else 
	    finisha = .true.
	  end if
	  if (ipb .lt. ib) then
            ipb = ipb +1
	    j2 = nocorb(ipb,jb)
	  else 
	    finishb = .true.
	  end if
        else if (icase .eq. 3 ) then
*         .. shells occupied only in right
	  diff = diff+nelcsh(ipb,jb)
	  if (ipb .lt. ib) then
            ipb = ipb +1
	    j2 = nocorb(ipb,jb)
	  else 
	    finishb = .true.
	  end if
        end if
        if (diff .gt. 2 ) then
 	  Let = 0
	  return
        end if
        if (.not. (finisha .and. finishb)) go to 10
      end if
*
*     We have two configurations differing by one electron
*     same angular configuration.
*
      ipa = 1
      ipb = 1
      j1 = nocorb(ipa,ja)
      j2 = nocorb(ipb,jb)
      finisha = .false.
      finishb = .false.
      ish = 0
100   if (finisha) then 
	  icase = 3
       else if (finishb) then
	  icase = 1
       else
	  if (j1 .lt. j2) then
	    icase = 1
	  else if (j1 .eq. j2) then
	    icase = 2
	  else
	    icase = 3
	  end if
       end if
       ish = ish + 1
       if (icase .eq. 1) then
*       .. insert shell on left; dummy for right
*       .. in the call to setupm, 1=.not.occupied; 2=occupied
	call setupm(ish,ipa,ipb,ja,jb,2,1)
	if (ipa .lt. ia) then
	  ipa = ipa+1
	  j1 = nocorb(ipa,ja)
	else
	  finisha = .true.
	end if
      else if (icase .eq. 2) then
*       .. insert same shell for left and right
	call setupm(ish,ipa,ipb,ja,jb,2,2)
	if (ipa .lt. ia) then
	  ipa = ipa+1
	  j1 = nocorb(ipa,ja)
	else
	  finisha = .true.
	end if
	if (ipb .lt. ib) then
	  ipb = ipb+1
	  j2 = nocorb(ipb,jb)
	else
	  finishb = .true.
	end if
      else
*       .. insert dummy shell on left, shell on right
	call setupm(ish,ipa,ipb,ja,jb,1,2)
	if (ipb .lt. ib) then
	  ipb = ipb+1
	  j2 = nocorb(ipb,jb)
	else
	  finishb = .true.
	end if
      end if
      if (.not. (finisha .and. finishb)) go to 100
      ihsh = ish
*
*     Shift the resultant coupling to be contiguous
*     to the shell coupling.  Note that for one-electron operators
*     the two configurations may differ by no more than two shells
*     so the maximum number of shells in the initial and final state
*     is IHSH=10.
*
      if (ihsh .gt. 1 .and. ihsh .lt. 10) then
	do 200 i =1,2
	  do 210 k = 1,3
	    do 220 ish = 1,ihsh-1
	      J1QN(ihsh+ish,k,i) = J1QN(10+ish,k,i)
220         continue
210       continue
200     continue
      end if
      end
