      subroutine   droti   ( nz, x, indx, y, c, s )
C
C     ==================================================================
C     ==================================================================
C     ====  DROTI  --  APPLY INDEXED DOUBLE PRECISION GIVENS        ====
C     ====             ROTATION                                     ====
C     ==================================================================
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
C     PURPOSE
C     -------
C
C         DROTI APPLIES A GIVENS ROTATION TO 
C             A SPARSE VECTOR   X  STORED IN COMPRESSED FORM  (X,INDX) 
C         AND 
C             ANOTHER VECTOR  Y  IN FULL STORAGE FORM.
C
C         DROTI DOES NOT HANDLE FILL-IN IN  X  AND THEREFORE, IT IS
C         ASSUMED THAT ALL NONZERO COMPONENTS OF  Y  ARE LISTED IN 
C         INDX.  ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN
C         INDX  ARE REFERENCED OR MODIFIED.  THE VALUES IN  INDX  MUST
C         BE DISTINCT TO ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
C
C         ALTHOUGH DISTINCT INDICES WILL ALLOW VECTOR OR PARALLEL
C         EXECUTION, MOST COMPILERS FOR HIGH-PERFORMANCE MACHINES WILL
C         BE UNABLE TO GENERATE BEST POSSIBLE CODE WITHOUT SOME 
C         MODIFICATION, SUCH AS COMPILER DIRECTIVES, TO THIS CODE.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICIES OF THE
C                             COMPRESSED FORM.  IT IS ASSUMED THAT
C                             THE ELEMENTS IN  INDX  ARE DISTINCT.
C         C,S     DOUBLE      THE TWO SCALARS DEFINING THE GIVENS
C                             ROTATION.
C
C     UPDATED ...
C
C         X       DOUBLE      ARRAY CONTAINING THE VALUES OF THE 
C                             SPARSE VECTOR IN COMPRESSED FORM.
C         Y       DOUBLE      ARRAY WHICH CONTAINS THE VECTOR  Y
C                             IN FULL STORAGE FORM.  ONLY THE
C                             ELEMENTS WHOSE INDICIES ARE LISTED IN
C                             INDX  HAVE BEEN REFERENCED OR MODIFIED.
C
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C
      DOUBLE PRECISION    X (*), Y (*), C, S
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
      DOUBLE PRECISION    TEMP
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      IF  ( ( C .EQ. 1.0D0 ) .AND. ( S .EQ. 0.0D0 ) )  RETURN
C
      DO 10 I = 1, NZ
          TEMP        = - S * X (I)  +  C * Y (INDX(I))
          X (I)       =   C * X (I)  +  S * Y (INDX(I))
          Y (INDX(I)) = TEMP
   10 CONTINUE
C
      RETURN
      END
