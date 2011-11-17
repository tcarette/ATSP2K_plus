*
*     ------------------------------------------------------------------
*	S L O P E  
*     ------------------------------------------------------------------
*
      subroutine slope
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=128,NOD=220)
      character*3 el
      COMMON /RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON /NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7),el(nwd)
      COMMON /FO/IWF(2),NCLOS(2),NOPEN(2),MAXNFO
      dimension xa(4),ya(4)
C     Determine the starting parameter az =  p(r)/r^(l+1)|r=0
      do 10 j = 1,maxnfo
        do 20 i = 4,1,-1
          xa(i) = r(i)
          ya(i) = p(i,j)*dsqrt(r(i))/(dexp((l(j)+1.d0)*dlog(r(i))))
20      continue
        call polint(xa,ya,4,0.d0,az(j),daz)
10    continue
      return
      end
