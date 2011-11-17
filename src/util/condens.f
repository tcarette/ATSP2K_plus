*     ------------------------------------------------------------------
*
*       A GENERAL PROGRAM TO CONDENSE THE LIST OF CONFIGURATIONS
*
*       by C. Froese Fischer
*          Vanderbilt University
*          Nashville, TN 37235 USA
*
*       March, 1984
*       Modified August, 1989 for Dynamic Memory allocation
*     ------------------------------------------------------------------
*
*
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
Ctc	PARAMETER (NWD = 60,NCD=250000)      !!!!! Change 6/07/2009 for large calculations
	PARAMETER (NWD = 60,NCD=4000000)
        CHARACTER S1*1,S4*4,S7*7,S19*19,S21*21,S23*23,S24*24
	CHARACTER*3 EL(NWD)
        CHARACTER*24 INFILE,OUTFILE,NAME
        CHARACTER HEADER*30,CLOSD*72,COUPLE*64
        DIMENSION INDX(8)
	POINTER (QCFG,CFG(8,1)),(QU,U(1)),(QV,V(1)),(QIP,IP(1)),
     :          (QCOUP,COUP(8,1)),(QKEEP,KEEP(1))
	LOGICAL KEEP
*
* ...   Machine dependent unit numbers for Standard input/output/error
*       ISCW = 5
        ISCW = 0
	IREAD = 5
	IWRITE = 6
 
        WRITE(ISCW,'(//A/A/A/A/A)')' Selection Process:','   1  name.c',
     :        '   2  name.l ','   3  name.j','   4  User Defined'
        WRITE(ISCW,'(//A)') 
     :  ' Enter selection along with ISORT (0=NO/1=YES)'
        READ(IREAD, *) ICASE,ISORT
        iarg = iargc()
        if ( iarg .eq. 0) then
           WRITE(ISCW,*) 'Enter name of atom'
           READ(IREAD,'(A)') NAME
        else
           call getarg(1,name)
        end if
	j = index(name,' ')
	infile = name(1:j-1)//'.c'
 
        OPEN(UNIT=1,FILE=infile,STATUS='OLD')
 
        OUTFILE = 'cfg.out'
        OPEN(UNIT=3,FILE=OUTFILE,STATUS='UNKNOWN')
 
        READ(1,1)  HEADER,NCLOSD,NWF,NCFG,IDIM,NCDIM,NIJ
  1     FORMAT(A30,I3,I4,I7,3I8)
	print *, 'NCLOSD,NWF,NCFG',nclosd,nwf,ncfg
100	READ(1,'(A)') CLOSD
	if (closd(1:30) .eq. header) go to 100
	IF (NWF .GT. NCLOSD) 
     :      READ (1,'(20(1x,A3))') (EL(k),k=nclosd+1,nwf)
*
* ...   allocate memory
*
	
        ncfgd= max(ncfg,NCD)
	Print *, 'Allocating memory for ',ncfgd, 'CFS'
	call alloc(qcfg,8*ncfgd,8)
	call alloc(qcoup,8*ncfgd,8)
        call alloc(qu,ncfgd,8)
	call alloc(qv,ncfgd,8)
	call alloc(qip,ncfgd,4)
	call alloc(qkeep,ncfgd,4)
*
* ...   Initialize arrays
*
	do 2 i = 1,ncfgd
	   v(i) = 0.d0
	   keep(i) = .true.
  2     continue
 
	DO 10 I = 1,NCFGd
           READ(1,'(8A8,F10.8)',END=20) (CFG(k,I),k=1,8),U(I)
	   READ(1,'(A60)',END=20) COUPLE
           READ(COUPLE,'(8A8)') (COUP(k,I),k=1,8)
           IF (ICASE .EQ. 1) V(I) = ABS(U(I))
 10     CONTINUE
 
 20     N = i-1
	print *, 'ncfg=',n
        IF (ICASE .EQ. 2 .OR. ICASE .EQ. 3) THEN
              IF (ICASE .EQ. 2) THEN
              INFILE =NAME(1:j-1)//'.l'
           ELSE
              INFILE =NAME(1:j-1)//'.j'
           END IF
           OPEN(UNIT=2,FILE=INFILE,STATUS='OLD')
	   READ(2,'(A)') HEADER
25         READ(2,'(//22X,I4)',END=45) NUMBER
           DO 30 J = 1,NUMBER
              READ(2,'(//(7F11.8))') (U(I),I=1,N)
              DO 40 I = 1,N
                 V(I) = MAX(V(I),ABS(U(I)))
 40           CONTINUE
 30        CONTINUE
           GO TO 25
 45        CONTINUE
        ELSE IF (ICASE .EQ. 4) THEN
           DO 32 I = 1,N
              V(I) = 1.D0
 32        CONTINUE
           WRITE(ISCW,'(A/A)')
     :          ' Enter the index of the configurations to',
     :          ' be deleted  (five at a time, terminate with a zero)'
           M = 1
 34        WRITE(ISCW,'(I3,A2)') M,': '
           READ(IREAD,*) INDX
           DO 35 K = 1,8
              I = INDX(K)
              IF (I .NE. 0) THEN
                 V(I) = 0.
              ELSE
                 GO TO 36
              END IF
              GO TO 34
 35        CONTINUE
        END IF
 
 36     TOL = 1.D-8
        IF (.NOT.(ICASE .EQ. 4)) THEN
           WRITE(ISCW,'(//A)')' Tolerance for acceptance -- FORMAT(F): '
           READ(IREAD,*) TOL
        END IF
	CALL SORT(N,QV,QIP,ISORT)
	DO 60 II = 1,N
           IF (V(II) .Lt. TOL) THEN
	      KEEP(II) = .FALSE.
	      ncfgold=ncfg
	      ncfg = ncfg-1
	   END IF
 60     CONTINUE

*       WRITE(3,1) HEADER,NCLOSD,NWF,NCFG,IDIM,NCDIM,NIJ
        WRITE(3,1) HEADER
	WRITE(3,'(A)') CLOSD
*	IF (NWF .GT. NCLOSD) 
*     :      WRITE(3,'(20(1x,A3))') (EL(k),k=nclosd+1,nwf)
        DO 70 II = 1,N
	   I = IP(II)
	   IF (KEEP(I)) THEN
		write(couple,'(8A8)') (coup(k,i),k=1,8)
                KK = 64
72              IF (COUPLE(KK:KK) .EQ. ' ') THEN
                   KK = KK-1
                   GO TO 72
                END IF
		IF (ICASE .EQ. 1 .OR. ICASE .EQ. 4) THEN
                  WRITE(3,71) (CFG(k,I),k=1,8),U(I),COUPLE(1:KK)
		ELSE
                  WRITE(3,71) (CFG(k,I),k=1,8),V(I),COUPLE(1:KK)
		END IF
           END IF
 71        FORMAT(8A8,F11.8/A)
 70     CONTINUE
*       Terminate with an asterisk in column 1
*       write(3,'(A1)') '*'
        CLOSE(UNIT=2)
        CLOSE(UNIT=3)
*
* ...   deallocate memory
*
	call dalloc(qcfg,8*ncfgold)
	call dalloc(qcoup,8*ncfgold)
        call dalloc(qu,ncfgold)
	call dalloc(qv,ncfgold)
	call dalloc(qip,ncfgold)
	call dalloc(qkeep,ncfgold)
*
        END
*--------------------------------------------------------------
*       S O R T
*--------------------------------------------------------------
*
	SUBROUTINE SORT(N,QV,QIP,ISORT)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	POINTER (QV,V(1)), (QIP,IP(1))

	DO 1 I = 1,N
	   IP(I) = I
  1     CONTINUE

	IF (ISORT .EQ. 0) RETURN

	DO 10 I = 1,N-1
          JP = I
          DO 12 J = I+1,N
             IF (V(IP(J)) .GT. V(IP(JP))) JP = J
 12       CONTINUE
	  ITEMP = IP(I)
	  IP(I) = IP(JP)
	  IP(JP) = ITEMP
 10     CONTINUE
	END
