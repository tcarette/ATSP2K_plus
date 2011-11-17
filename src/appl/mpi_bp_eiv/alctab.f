*****:******************************************************************
*
*
* --- This SUBROUTINE allocates memory for the tables of integrals
*     for the radial arrays, and for the local arrays for the three
*     components of the Breit-Pauli interaction matrix.
*
*****:******************************************************************
*
      Subroutine ALCTAB(ncfg,nwf,skip)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NOD=220,LINT=500)
      LOGICAL skip
*
      POINTER (qh,h(ncfg,3)),(qjan,jan(ncfg,3))
      COMMON /buffer/qh, qjan, nrow(3), iflag(3)
*
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :       (IQMAX,MAX(1))
      COMMON /NEL/IQP,IQN,IQL,IQAZ,IQMAX,IQ(7)
*
      LOGICAL rel
*
      call alloc(qh,3*ncfg,8)
      call alloc(qjan,3*ncfg,4)
*
      call alloc(iqp,nod*nwf,8)
      call alloc(iqn,nwf,4)
      call alloc(iql,nwf,4)
      call alloc(iqaz,nwf,8)
      call alloc(iqmax,nwf,4)
*
      end
