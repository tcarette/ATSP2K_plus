*     ------------------------------------------------------------------

*       PROGRAM  LEVELS -- PRINT THE LEVELS IN A NAME.J FILE
*
*                By C. Froese Fischer
*                   Bin Liu
*                   Vanderbilt University
*                   Nashville, TN 37235 USA
*
*                May, 1983
*
*     ------------------------------------------------------------------
*
*       Output Energy Level Table from NAME.J directly
*
 
****************
*               Declarations
*
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	PARAMETER (NCD2=4000)
        CHARACTER    CH1*1,NAME5*8, NAME*24
        CHARACTER*40 LAB(NCD2),LAB1,LAB2,COF1,COF2
        CHARACTER*2  CH2,TERM(NCD2)
        CHARACTER*1 YES
        CHARACTER    TERM1*2,TERM2*2
        DOUBLE PRECISION  Z5,COE(7)
        DOUBLE PRECISION  EN(NCD2),LOW,EN1,T
        INTEGER      JJ(NCD2),PS(NCD2),icount
	PARAMETER (NELEMENT=50)
        character*80 rd_line
*
      DOUBLE PRECISION ME,  RY, RZ, ZMU
      INTEGER MN(NELEMENT)
      DATA ME/548.579903D-6/,RY/109737.31534/ MN/1,4,7,9,11,12,14,16,
     : 19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,58,
     : 63,64,69,74,75,80,79,84,85,88,89,90,93,98,99,102,103,108,
     : 107,114,115,120/
*

 
*
*      LAB  =  Array for labels with packed form
*      TERM  =  Array for final Term
*       JJ  =  Array for the Values of J*2
*       EN  =  Array for energy
*       PS  =  Pointer to the location of STATE
*
 
***************
*               Read the header of NAME.J
*
*
        iarg = iargc()
        if (iarg .eq. 0) then
	   WRITE(0,*) 'Enter name and type (.l or .j) of file'
99         read(5,'(A)') NAME
	   j = index(NAME, '.')
	   if (j .eq. 0) then
	      print *, 'Enter name along with type'
	      go to 99
	   end if
*        else
*           call getarg(1,NAME)
        end if
        OPEN(UNIT=15,FILE=NAME,STATUS='OLD')

        NS = 0
100     icount = 1;
        write(*,*)
        read(15,'(A)',END=200) rd_line        
        if(rd_line(11:11) == "Z") then
           READ (rd_line,50,END=200) NAME5,Z5,N5,NCFG 
        end if
50      FORMAT (2X,A6,6X,F5.1,8X,I3,10X,I6)
        NZ = INT(Z5)
        NEL = N5
 
*
*       Decide the number of records for the coefficients
*
        LINE = INT(NCFG/7)
        IF (LINE*7 .LT. NCFG) LINE = LINE+1
 
 
 
****************
*               Read NAME.J
*
            READ(15,'(1X,A1)',END=200) CH1
*            IF (CH1 .EQ. '*') GOTO 100
*            if (icount-1== LINE*number) go to 100
102         if (rd_line(11:11) == 'Z') read (15,51,end=200)
            READ (15,51,END=200) J1,NUMBER
51          FORMAT (9X,I3,11X,I3)
            DO 104 I=1,NUMBER
                READ (15,52,END=200) MAX,EN1,LAB1
                write(6,'(A30,F16.9)') LAB1,EN1
52              FORMAT (/I6,F16.9,2X,A30)
                DO 106 J=1,LINE
                    icount = icount + 1
                    READ (15,53,END=200) (COE(K), K=1,7)
53                  FORMAT (7F11.8)
106             CONTINUE
 
*
*       Delete the extra characters
*
                N = INDEX(LAB1,' ')-1
                IF (INDEX(LAB1,'(0)') .GT. 0) THEN
                    J = INDEX(LAB1,'.')
                    LAB2 = LAB1(J+1:)
                    LAB1 = LAB2
                ENDIF
 
*
*       1. Separate the final term
                IF (LAB1(N:N) .GE. '0' .AND.
     :              LAB1(N:N) .LE. '9'      ) THEN
                    N = N-4
		    K = 4
                ELSE
                    N = N-3
		    K = 3
                ENDIF
                TERM1 = LAB1(N+2:N+3)
                LAB1(N+1:N+K) = '   '
 
*
*       2. Delete set subscript
*
                CH1 = LAB1(N:N)
                IF (CH1.GE.'0' .AND. CH1.LE.'9') LAB1(N:N) = ' '
 
*
*       3. If after removing the final term, there are no other
*   intermediate couplings prefaced by '_' and the last coupling
*   is the same as the final term, then the coupling for the TERMal
*   term is omitted .
*
                J = INDEX(LAB1,'_')
                IF (J .EQ. 0) THEN
                    CH2 = LAB1(N-2:N-1)
                    IF (CH2 .EQ. TERM1) LAB1(N-2:N-1) = '  '
                ENDIF
 
****************
*       Assign the values into arrays
*
                NS = NS+1
                LAB(NS) = LAB1
                TERM(NS) = TERM1
                JJ(NS) = J1
                EN(NS) = EN1
                PS(NS) = NS
104         CONTINUE
            if (icount -1 == line*number) go to 100
        GO TO 102
 
 
***************
*       Sort the energy(a.u.) in increasing order
*
200     DO 202 I=1,NS-1
            MIN = I
            T = EN(I)
            DO 204 J=I+1,NS
                IF (EN(J) .LT. T) THEN
                    MIN = J
                    T = EN(J)
                ENDIF
204         CONTINUE
            EN(MIN) = EN(I)
            EN(I) = T
            K = PS(I)
            PS(I) = PS(MIN)
            PS(MIN) = K
202     CONTINUE
        LOW = EN(1)
 
***************
*       Compute   Rz = 109737.31534/(1.+548.579903D-6/ZMU)
*
      WRITE(iscw,'(/1X,A)') ' Default Rydberg constant (y/n)'
      READ(5,'(A1)') YES
      IF (YES.EQ.'Y'.OR.YES.EQ.'y'.AND.(INT(Z5).LE.NELEMENT)) THEN
         ZMU=MN(INT(Z5))
      ELSE
         WRITE(iscw,'(1X,A)') 'Enter the mass of the atom'
         READ(5,*) ZMU
      END IF

        RZ = RY/(1.+ME/ZMU)
*	print *, ry, me, zmu
 
 
***************
*       Output header of the table from Line Printer
*
 
        WRITE (6,212)
212     FORMAT (//' ENERGY LEVELS')
        WRITE (6,214) NZ,NEL
214     FORMAT (//' Z = ',I3,5X,I3,' electrons'///)
        WRITE(6,'(A, F15.5)') 'Rydberg Constant Used = ', RZ
        WRITE (6,216)
216     FORMAT (64('-'))
        WRITE (6,*) 'Configuration         Term  J     Total ',
     :'Energy  Energy Level '
        WRITE (6,*) '                                     (a.u.)',
     :'        (cm-1)'
        WRITE (6,216)
 
 
****************
*       Compute Energy(cm)i by (Ei(a,u)-Eground(a,u))*2*Rz, and output
*   Energy Level Table  from line printer
*
        COF1 = ' '
        TERM1 = ' '
        DO 222 I=1,NS
            EN1 = EN(I)
            ECM = (EN1-LOW)*2*RZ
 
*
*       Assign values to configuration, final term, and J
*
            J = PS(I)
            COF2 = LAB(J)
            TERM2 = TERM(J)
            PJ = JJ(J)/2.
 
*
*       Omit printing of repeated configurations
*
            IF (COF1 .EQ. COF2) THEN
                COF2 = ' '
                IF (TERM1 .EQ. TERM2) THEN
                        TERM2 = ' '
                ELSE
                        TERM1 = TERM2
                        WRITE (6,*)
                ENDIF
            ELSE
                COF1 = COF2
                TERM1 = TERM2
                WRITE (6,*)
                WRITE (6,*)
            ENDIF
 
*
*       Omit printing of repeated final terms
*
 
            WRITE (6,224) COF2,TERM2,PJ,EN1,ECM
224         FORMAT (1X,A20,2X,A2,2X,F4.1,2X,F14.9,2X,F11.4,5X,1PE9.2)
222     CONTINUE
 
        WRITE (6,216)
        write (*,*) 
        END
