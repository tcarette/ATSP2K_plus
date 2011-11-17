*----------------------------------------------------------------------
*     A L C S T S
*----------------------------------------------------------------------
*    This routine allocates arrays associated with states.  For the 
*    yint.lst arrays, memory needs to be allcoated for all blocks,
*    but c.lst arrays are read in groups of size lsdim.

      SUBROUTINE ALCSTS(cf_tot);
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'mpif.h'
      parameter (MAXPROC=100)
      common /MPI/ myid, nprocs, ierr, istat(MPI_STATUS_SIZE)

      PARAMETER(MEIG=20,MTERM=20)
      PARAMETER(MAXCLST=250000000);!keep clst on disk, memory swapping, 2GB
      POINTER(pkval,kval(1)),(pvalue,value(1)),
     :       (pih,ih(1)),(pjh,jh(1)),
     :       (pcoeff,coeff(1)),(pnijptr,nijptr(1)),
     :       (pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PVALUE,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR

*     .. Logical array passed from nonh to give the used integrals
      POINTER (plused,lused(1))
      COMMON/USED/plused
      LOGICAL lused
      POINTER (pico,ico(1))
      COMMON/ICOFF/pico
      logical ::clst_disk, clst_memory,clst_memory_t
      common/use_disk/clst_disk, clst_memory
      integer ierr_mem_cf,ierr_mem_inp
      INTEGER       IN,OUT,PRI,IUC,IUF,IUD,OUC,OUF,OUD,OUH,ISCW
      COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,IUD,OUC,OUF,OUD,OUH,ISCW

      logical leigen(meig,mterm), lguess
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm)
      character term_bl(mterm)*3
      integer 	cf_tot(nblock)
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess

*

      clst_memory_t = .true.
        !IDIM = IDIM + 1
        nze_tot = sum(nze_bl(1:nblock))
        ncfg_tot = sum(ncfg_bl(1:nblock))
        nnn = sum(cf_tot(1:nblock));
*
        call alloc(pkval,idim,4)
        call alloc(pvalue,idim,8)
        call alloc(plused,idim,4)
        call alloc(pih,nze_tot,4)
        call alloc(pico,nze_tot,4)

        if (clst_memory) then
            pcoeff = malloc(%val(nnn*8));
            pinptr = malloc(%val(nnn*4));
!           call hpalloc(pcoeff,nnn*8,ierr_mem_cf,0);
!           call hpalloc(pinptr,nnn*4,ierr_mem_inp,0);
        end if

        if ((pcoeff==0).or.(pinptr==0).or.(nnn > MAXCLST)) then
          clst_memory_t = .false.
        end if

        call MPI_ALLREDUCE(clst_memory_t,clst_memory,1,MPI_LOGICAL,
     :                     MPI_LAND,MPI_COMM_WORLD,ierr) 
*          call MPI_BCAST(clst_memory,1,MPI_LOGICAL,0,
*     :                  MPI_COMM_WORLD,ierr)

        if ((.not.clst_memory).and.(clst_disk)) then
          if (myid == 0) then
        write(0,'(A)') 'Insufficient memory for loading c.lst: ',
     :                 'On each iteration read c.lst from disk'
          end if
           call free(%val(pcoeff));
           call free(%val(pinptr));
           call disclaim(%val(pcoeff));
           call disclaim(%val(pinptr));
!          call hpdeallc(pcoeff,ierr_mem_cf,0)
!          call hpdeallc(pinptr,ierr_mem_inp,0)
          call alloc(pcoeff,ncodim,8)
          call alloc(pinptr,ncodim,4)
        else if (.not.clst_disk) then
          write(0,'(A)') 'Insufficient memory for loading c.lst: ', 
     :                   'Stop in spalcsts'
          Stop;
        end if

*
        RETURN
        END
