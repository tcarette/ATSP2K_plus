*
*     ------------------------------------------------------------------
*      cfg_inf
*     ------------------------------------------------------------------
*
*        This routine writes the cfg.inp file with data needed
*

      SUBROUTINE cfg_inf(i,ncoff,itotal,NCLOSD,NCDIM,IDIM,NIJ,who,
     :                   term)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
     :       (QLJCLSD,LJCLSD(1))
      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
      character term*3, who*5
      integer i 
*
  1   FORMAT(A12,i10,I8,i3,I4,I6,I8,I8,I8,2x,a5)
      NWF = MAXORB + NCLOSD
      WRITE(25,1) term,ncoff,itotal,NCLOSD,NWF,NCFG,IDIM,NCDIM,NIJ,who
      RETURN
      END

