*----------------------------------------------------------------------
*     A L C S T S
*----------------------------------------------------------------------
*    This routine allocates arrays associated with states.  For the 
*    yint.lst arrays, memory needs to be allcoated for all blocks,
*    but c.lst arrays are read in groups of size lsdim.

      SUBROUTINE ALCSTS(cf_tot);
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

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
      integer ierr_mem_ih,ierr_mem_ico,ierr_mem_cf,ierr_mem_inp
      logical   :: hmx_memory, ico_memory, ih_memory, clst_memory
      common/memory_use/hmx_memory, ico_memory, ih_memory, clst_memory


      INTEGER       IN,OUT,PRI,IUC,IUF,IUD,OUC,OUF,OUD,OUH,ISCW
      COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,IUD,OUC,OUF,OUD,OUH,ISCW
*
      logical leigen(meig,mterm), lguess
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm),
     :        nze_max(mterm)
      character term_bl(mterm)*3
      integer 	cf_tot(nblock)
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess,nze_max
      POINTER (qhmx,hmx(1)),(qtm,tm(1)),(qtp,tp(1))
      common/spd/qhmx,qtm,qtp,qdiag,qiwork 
cgd   integer*8 new_size
*

!      IDIM = IDIM + 1
      nze_tot = sum(nze_bl(1:nblock))
      ncfg_tot = sum(ncfg_bl(1:nblock))
      nnn = sum(cf_tot(1:nblock));
*     nze_max_bl is maximum nze of a block 
      nze_max_bl = maxval(nze_bl(1:nblock))
*     nze_max_col is maximum nze of a cloumn of the matrix 
      nze_max_col = maxval(nze_max(1:nblock))

      call alloc(pkval,idim,4)
      call alloc(pvalue,idim,8)
      call alloc(plused,idim,4)
      
      hmx_memory  = .true.
      ico_memory  = .true. 
      ih_memory   = .true.
      clst_memory = .true.

cgd   qhmx = malloc(nze_max_bl*8)
      call alloc(qhmx,nze_max_bl,8)
cgd   pih = malloc(nze_tot*4)
      call alloc(pih,nze_tot,4)
!      qhmx = 0  ! DBG code
      if((qhmx == 0).or.(pih == 0).or.(nze_tot>MAXCLST)) then
        write(0,*) ' Interaction Matrix on Disk '
        call free(qhmx)
        qhmx = 0
        call free(pih)
        pih = 0
        hmx_memory = .false.;
        ico_memory = .false.
        ih_memory = .false.;
        clst_memory = .false.;
cgd     pico = malloc(nze_max_col*4)
        call alloc(pico,nze_max_col,4)
        call alloc(pcoeff,ncodim,8)
        call alloc(pinptr,ncodim,4)
      end if

      if((ih_memory).and.(hmx_memory)) then
cgd     pico = malloc(nze_tot*4)
        call alloc(pico,nze_tot,4)
!      pico = 0 ! DBG code
        if(pico == 0) then
          write(0,*) ' ico.lst on disk '
          call free(pico)
          pico = 0
cgd       pico = malloc(nze_max_col*4)
          call alloc(pico,nze_max_col,4)
          clst_memory = .false.;
          ico_memory = .false.;
          call alloc(pcoeff,ncodim,8)
          call alloc(pinptr,ncodim,4)
        end if
      end if

      if(ico_memory) then
        new_size = nnn*1.1
cgd     pcoeff = malloc(new_size*8)
        call alloc(pcoeff,new_size,8)
cgd     pinptr = malloc(new_size*4)
        call alloc(pinptr,new_size,4)
        pcoeff = 0;  ! code to make all c.lst on disk
        if ((pcoeff == 0).or.(pinptr == 0)) then
          write(0,*) ' Coefficients on Disk ';
          call free(pcoeff)
          call free(pinptr)
          pcoeff = 0
          pinptr = 0
          call alloc(pcoeff,ncodim,8)
          call alloc(pinptr,ncodim,4)
          clst_memory = .false.;
        else 
          call free(pcoeff)
          call free(pinptr)
          pcoeff = 0
          pinptr = 0
cgd     pcoeff = malloc(nnn*8)
        call alloc(pcoeff,nnn,8)
cgd     pinptr = malloc(nnn*4)
        call alloc(pinptr,nnn,4)
        end if
      end if

        RETURN
        END

