      PROGRAM DRIVER 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:22  11/16/01  
!...Switches:                     
!        PARAMETER(Nmax=2149,NZERmax=335416)
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dscal_I 
      USE dvdson_I 
      IMPLICIT NONE
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NMAX = 9000 
      INTEGER, PARAMETER :: NZERMAX = 600000 
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
!...  /MATRIX/ 
      COMMON /MATRIX/ DMY_MATR(NZERMAX), IND_COL(NMAX), IND_ROW(NZERMAX), &
         IUPPER 
      INTEGER   IND_COL, IND_ROW 
      REAL(DOUBLE) :: DMY_MATR 
      LOGICAL   IUPPER 
!...  /TEMP/ 
      COMMON /TEMP/ TM(NMAX), TP(NMAX) 
      REAL(DOUBLE) :: TM, TP 
!...  /HYPERSTUFF/ 
      COMMON /HYPERSTUFF/ IAM, NODES 
      INTEGER   IAM, NODES 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NUMEMAX = 1 
      INTEGER, PARAMETER :: LIMMAX = NUMEMAX + 18 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(6*LIMMAX + NUMEMAX) :: IWORK 
      INTEGER , DIMENSION(10) :: IDIM 
      INTEGER , DIMENSION(NUMEMAX) :: MINELEM 
      INTEGER , DIMENSION(LIMMAX) :: ISELEC 
      INTEGER :: IN_UNIT, N, NZER, I, MYCOL, IUP, LOCALUP, ICOL, IDOWN, &
         LOCALDOWN, IROW, JUNK, IYES, ND, IRWSZ, IIWSZ, LIM, ILOW, IHIGH, NIV, &
         NUME, MBLOCK, MAXITER, NLOOPS, NMV, IERR, J 
      REAL , DIMENSION(2) :: TARRAY 
      REAL :: T1 
      REAL(DOUBLE), DIMENSION(NMAX) :: DIAG 
      REAL(DOUBLE), DIMENSION(LIMMAX*(2*NMAX + LIMMAX + (NUMEMAX + 10)) + &
         NUMEMAX) :: WORK 
      REAL(DOUBLE) :: RUBS, DIAMAX, DSHIFT, CRITE, CRITC, CRITR, ORTHO 
      LOGICAL :: IGUESS, HIEND 
!-----------------------------------------------
!
!
      IAM = 0 
      NODES = 1 
      IUPPER = .TRUE. 
 
!******* Reading in the upper matrix.
      IN_UNIT = 10 
      OPEN(IN_UNIT, FILE='../RA/Hypercube/SUMMER92/MATRICES/mat.748', STATUS=&
         'unknown', FORM='formatted', POSITION='asis') 
 
      READ (IN_UNIT, *) N, NZER 
!* READ THE COLUMN INDICES **
      READ (IN_UNIT, *) (WORK(I),I=1,N) 
      WRITE (6, *) 'read col' 
 
!** READ ROW INDICES AND CREATE IND_COL **
      MYCOL = 0 
      IUP = 0 
      LOCALUP = 0 
      DO ICOL = 1, N 
         IDOWN = IUP + 1 
         IUP = WORK(ICOL) 
         IF (IAM == MOD(ICOL - 1,NODES)) THEN 
            MYCOL = MYCOL + 1 
            LOCALDOWN = LOCALUP + 1 
            LOCALUP = LOCALDOWN + IUP - IDOWN 
            READ (IN_UNIT, *) (IND_ROW(IROW),IROW=LOCALDOWN,LOCALUP) 
            IND_COL(MYCOL) = LOCALUP 
         ELSE 
            READ (IN_UNIT, *) (JUNK,IROW=IDOWN,IUP) 
         ENDIF 
      END DO 
      WRITE (6, *) 'read row' 
 
!** READ VALUES **
      MYCOL = 0 
      IUP = 0 
      LOCALUP = 0 
      DO ICOL = 1, N 
         IDOWN = IUP + 1 
         IUP = WORK(ICOL) 
         IF (IAM == MOD(ICOL - 1,NODES)) THEN 
            MYCOL = MYCOL + 1 
            LOCALDOWN = LOCALUP + 1 
            LOCALUP = IND_COL(MYCOL) 
            READ (IN_UNIT, *) (DMY_MATR(I),I=LOCALDOWN,LOCALUP) 
         ELSE 
            READ (IN_UNIT, *) (RUBS,IROW=IDOWN,IUP) 
         ENDIF 
      END DO 
      WRITE (6, *) 'read val' 
 
      IF (IUPPER) THEN 
         DO I = 1, N 
            DIAG(I) = DMY_MATR(IND_COL(I)) 
            DIAMAX = MAX(DIAMAX,ABS(DIAG(I))) 
         END DO 
      ELSE 
         DO I = 1, N 
            DIAG(I) = DMY_MATR(IND_COL(I-1)+1) 
            DIAMAX = MAX(DIAMAX,ABS(DIAG(I))) 
         END DO 
      ENDIF 
 
      WRITE (6, *) NZER 
!     do 246 i=1,NZER
! 246    dmy_matr(i)=dmy_matr(i)/diamax
!     print*,diamax
      NZER = LOCALUP 
      ND = N 
 
!************************************************************************
!      TEsting
!***********************************************************************
      IRWSZ = LIMMAX*(2*NMAX + LIMMAX + (NUMEMAX + 10)) + NUMEMAX 
      IIWSZ = 6*LIMMAX + NUMEMAX 
 
      N = ND 
      DO LIM = LIMMAX, LIMMAX, 10 
         HIEND = .FALSE. 
         ILOW = 1 
         IHIGH = NUMEMAX 
         ISELEC(1) = 3 
         ISELEC(2) = 5 
         ISELEC(3) = -1 
         NIV = 0 
         NUME = NUMEMAX 
!           CRITC=1D-5
!           CRITR=1D-5
!        TRHOLD=1D-2
!        ORTHO=1D+6
!        MAXITER=1000
         IF (.NOT.N<ND) THEN 
            CRITE = 1.0D-16 
            CRITC = 1.D-8 
            MBLOCK = 1 
            CRITR = 1.D-8 
            ORTHO = 1D-15 
            MAXITER = 150 
         ENDIF 
!     t1=dclock()
         T1 = ETIME(TARRAY) 
         CALL DVDSON (DMY_MATR, IND_ROW, IND_COL, IUPPER, NZERMAX, TM, TP, N, &
            LIM, DIAG, ILOW, IHIGH, ISELEC, NIV, MBLOCK, CRITE, CRITC, CRITR, &
            ORTHO, MAXITER, WORK, IRWSZ, IWORK, IIWSZ, HIEND, NLOOPS, NMV, IERR&
            ) 
!        t1=(dclock()-t1)
         T1 = ETIME(TARRAY) - T1 
 
         IF (IYES == 1) THEN 
            WORK(NUME*N+1:NUME*(N+1)) = DIAMAX*WORK(NUME*N+1:NUME*(N+1)) + &
               DSHIFT 
         ENDIF 
 
         WRITE (6, *) 'Time', T1 
         IF (HIEND) NUME = N - ILOW + 1 
!eigenvalues, differences, residuals, eigenvectors.
         WRITE (6, 4000) IERR, NLOOPS, NMV 
         WRITE (6, 2000) ((WORK(I),I=NUME*N + J,(N + 3)*NUME,NUME),J=1,NUME) 
         WRITE (6, 3000) ((WORK(I),I=J,J + 5),J=1,NUME*N,N) 
 4000    FORMAT(/,'IERR =',I6,'   Matrix accesses=',I4,&
            '   Matrix-Vector products=',I4) 
 2000    FORMAT(/,/,9X,'Eigenvalues',8X,'Eigval Differences',6X,'Residuals'/,/(&
            D25.15,2D20.10)) 
 3000    FORMAT(/,/,' First six componenents of the eigenvectors'/,/(6D25.15)) 
!     :    (6D12.3))
 
      END DO 
      CLOSE(IN_UNIT) 
      STOP  
      END PROGRAM DRIVER 
