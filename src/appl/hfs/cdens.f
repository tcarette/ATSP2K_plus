*-----------------------------------------------------------------------
*     C D E N S
*-----------------------------------------------------------------------
*
*     This subroutine calculates the electron density
*     contribution from the core.

      SUBROUTINE cdens(ireadw,corden)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,NWD=70)
      CHARACTER EL(NWD)*3,ATOM*6,TERM*6
      dimension az(nwd),p(nod,nwd),n(nwd),l(nwd)

      common /closed/b1elc(4),nclosd,ibk 
      
      rewind (ireadw)
*
*  *****  READ THE RADIAL FUNCTIONS FOR THE CORE
*
      do 5 I = 1,nclosd
12       READ(ireadw)ATOM,TERM,EL(I),MR,Z,ETI,EKI,AZ(I),
     :       (P(J,I),J=1,MR)
         IF (EL(I)(1:1).NE.' ') THEN
            N(I)=ICHAR(EL(I)(1:1))-ICHAR('1')
            L(I)=LVAL(EL(I)(2:2))
         ELSE
            N(I)=ICHAR(EL(I)(2:2))-ICHAR('1')
            L(I)=LVAL(EL(I)(3:3))
         ENDIF
5     continue

*  Calculate electron density from the core

      corden = 0.d0
      do 10 k = 1,nclosd
         if (l(k).eq.0) then
            write(*,*) 'l(K)',l(k)
            contri = 2.d0*az(k)*az(k)
            corden = contri + corden
         endif
10    continue
      corden = corden*0.25d0*0.318309886d0
      return
      end
