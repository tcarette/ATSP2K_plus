*
*------------------------------------------------------------------------
*          I N V E R T
*------------------------------------------------------------------------
*
      SUBROUTINE INVERT(NC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION C(NC),ID(NC,5)
      POINTER (QC,C),(QID,ID)
      call alloc(qc,nc,8)
      call alloc(qid,5*nc,4)
      CALL OUTLSJ(NC,C,ID(1,1),ID(1,2),ID(1,3),ID(1,4),ID(1,5))
      RETURN
      END
