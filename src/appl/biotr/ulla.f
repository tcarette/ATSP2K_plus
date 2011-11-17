*
*     ------------------------------------------------------------------
*	U L L A    
*     ------------------------------------------------------------------
*
      Subroutine Ulla(A,U,L,Ndim,scr)
*
* Obtain U L decomposition of matrix A
*   A = U L
*
* Note that L and U are returned in full matrix form
*
* Quick and dirty routine, Jeppe Olsen, November 1991
*
      Implicit double precision  (a-h,o-z)
*
*. Input
      Dimension A(Ndim,Ndim)
*. Output
      Real * 8  U(Ndim,Ndim),L(Ndim,Ndim)
*. Scratch
      Dimension Scr(*)
*
* In order to change into LU form introduce the orthogonal matrix
*    p(i,j) = delta(i,ndim-i)
* and rewrite
* A = P P A P P  = P L U P = PLP PUP
* where LU is an LU decomposition of PAP, since PLP is upper
* tringular AND PUP is lower traingular we have obtained the goal
*
* 1 : PAP in scr(klPAP)
*
      KLFREE = 1
      KLPAP = KLFREE
      KLFREE = KLFREE + NDIM ** 2
*
      Do 100 I = 1, Ndim
        Do 90 J = 1, Ndim
         SCR(KlPAP-1+(J-1)*NDIM+I) = A(Ndim+1-i,Ndim+1-j)
   90   Continue
 100  Continue
* 2 : Lu decompose PAP
      klL = Klfree
      KLfree = Klfree + ndim*(Ndim+1)/2
*
      klU = Klfree
      KLfree = Klfree + ndim*(Ndim+1)/2
      Call Lulu(SCR(KLPAP),Scr(klL),Scr(klU),Ndim)
*          LULU(A,L,U,NDIM)
* Storage modes
*   L(I,J) = L(I*(I-1)/2 + J ) ( I .GE. J )
*   U(I,J) = U(J*(J-1)/2 + I ) ( J .GE. I )
*
*. 3 : Obtain U as PLP and L as PUP
      Call Setvec(U,0.0D0,Ndim ** 2 )
      Call Setvec(L,0.0D0,Ndim ** 2 )
*
      Do 200 Imax = 1, Ndim
        Do 150 Imin = 1, Imax
          U(Imin,Imax) =
     &    Scr(klL-1+(Ndim+1-Imin)*(Ndim+1-Imin-1)/2+(Ndim+1-Imax) )
          L(Imax,Imin) =
     &    Scr(klU-1+(Ndim+1-Imin)*(Ndim+1-Imin-1)/2+(Ndim+1-Imax) )
  150  Continue
  200 Continue
*
      Ntest = 0
      If(Ntest.Ne.0) then
        Write(6,*) ' U and L from Ulla '
        Call Wrtmat(U,Ndim,Ndim,Ndim,Ndim)
        Call Wrtmat(L,Ndim,Ndim,Ndim,Ndim)
      End if
*

      Return
      End
