        SUBROUTINE GDSUMMPI(x,n,y)
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)

        INCLUDE 'mpif.h'
        parameter (MAXPROC=100)
        parameter (ITUNE=4096)
        common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)

        dimension x(*),y(*)

*        POINTER(P,tmp1(1))
        INTEGER itmp,I,ierror
*       In order to make life easier and also have a tuning possibility
*       let's split x a little bit
*       should not create as much swapping as before ?
*
        DO I=1,N,ITUNE
*            p=loc(x(I))
            itmp = min(ITUNE,N-I+1)
            call MPI_ALLREDUCE(x(i),y,itmp,MPI_DOUBLE_PRECISION,MPI_SUM,
     $                          MPI_COMM_WORLD,ierror)
            call dcopy(itmp,y,1,x(i),1)
        ENDDO

        return
        end

