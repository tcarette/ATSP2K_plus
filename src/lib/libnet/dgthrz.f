      subroutine dgthrz ( nz, y, x, indx )
C
C     ==================================================================
C     ==================================================================
C     ====  DGTHRZ -- DOUBLE PRECISION GATHER AND ZERO              ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         DGTHRZ GATHERS THE SPECIFIED ELEMENTS FROM 
C             A DOUBLE PRECISION VECTOR  Y  IN FULL STORAGE FORM 
C         INTO 
C             A DOUBLE PRECISION VECTOR  X  IN COMPRESSED FORM  (X,INDX).  
C         FURTHERMORE THE GATHERED ELEMENTS OF  Y  ARE SET TO ZERO.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE REFERENCED OR MODIFIED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS TO BE GATHERED INTO 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             VALUES TO BE GATHERED INTO COMPRESSED FORM.  
C
C     UPDATED ...
C
C         Y       DOUBLE      ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  THE 
C                             GATHERED COMPONENTS IN  Y  ARE SET TO ZERO.  
C                             ONLY THE ELEMENTS CORRESPONDING TO THE
C                             INDICES IN  INDX  HAVE BEEN ACCESSED.
C
C     OUTPUT ...
C
C         X       DOUBLE      ARRAY CONTAINING THE VALUES GATHERED INTO
C                             THE COMPRESSED FORM.
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
      DOUBLE PRECISION    Y (*), X (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          X(I)       = Y(INDX(I))
          Y(INDX(I)) = 0.0D0
   10 CONTINUE
C
      RETURN
      END
