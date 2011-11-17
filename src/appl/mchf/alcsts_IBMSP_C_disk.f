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
      integer MB


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
      integer new_size
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


      MB = 1024*1024

      write(0,*)
      write(0,'(A40,A10,A10)') 'Requesting memory for the matrix',
     : 'Elements', '       MB'  
      write(0,'(A40,i10,i10)') 'Largest Matrix block',
     :   nze_max_bl, (nze_max_bl*8)/MB
      qhmx = malloc(%val(nze_max_bl*8))
      write(0,'(A40,i10,i10)') 'Total ih pointers',
     :   nze_tot, (nze_tot*4)/MB
      pih = malloc(%val(nze_tot*4))
!      qhmx = 0  ! DBG code
      if((qhmx == 0).or.(pih == 0).or.(nze_tot>MAXCLST)) then
        write(0,'(A40)') 'Insufficient memory, I.Matrix on Disk '
        call free(%val(qhmx))
        call disclaim(%val(qhmx))
        qhmx = 0
        call free(%val(pih))
        call disclaim(%val(pih))
        pih = 0
        hmx_memory = .false.;
        ico_memory = .false.
        ih_memory = .false.;
        clst_memory = .false.;
        pico = malloc(%val(nze_max_col*4))
        write(0,'(A40,i10,i10)') 'Alocating memory for pcoeff',
     :                            ncodim
        write(0,'(A40,i10,i10)') 'Requesting memory for inptr',
     :                            ncodim 
        call alloc(pcoeff,ncodim,8)
        call alloc(pinptr,ncodim,4)
      end if

      if((ih_memory).and.(hmx_memory)) then
        pico = malloc(%val(nze_tot*4))
      write(0,'(A40,A10,A10)') 'Requesting memory for ico indices ',
     : 'Elements', '   MB   '
      write(0,'(A40,i10,i10)') 'Total pointers to nonzero elements',
     :   nze_tot, (nze_tot*4)/MB

!      pico = 0 ! DBG code
        if(pico == 0) then
          write(0,'(A40)') 'Insufficient memory, ico.lst on disk '
          call free(%val(pico))
          call disclaim(%val(pico))
          pico = 0
          pico = malloc(%val(nze_max_col*4))
          clst_memory = .false.;
          ico_memory = .false.;
          call alloc(pcoeff,ncodim,8)
          call alloc(pinptr,ncodim,4)
        end if
      end if

      if(ico_memory) then
        new_size = nnn*1.05
        write(0,'(A40,A10,A10)') 'Requesting memory for coef ',
     : 'Elements', '   MB   '
      write(0,'(A40,i10,i10)') 'Total coef from c.lst',
     :   new_size, (new_size*8)/MB 
        pcoeff = malloc(%val(new_size*8))
      write(0,'(A40,i10,i10)') 'Total inptr to coef from c.lst',
     :   new_size, (new_size*4)/MB
        pinptr = malloc(%val(new_size*4))
        pcoeff = 0;  ! code to make all c.lst on disk
        if ((pcoeff == 0).or.(pinptr == 0)) then
          write(0,'(A40)') ' Coefficients on Disk ';
          call free(%val(pcoeff))
          call disclaim(%val(pcoeff))
          call free(%val(pinptr))
          call disclaim(%val(pinptr))
          pcoeff = 0
          pinptr = 0
          call alloc(pcoeff,ncodim,8)
          call alloc(pinptr,ncodim,4)
          clst_memory = .false.;
        else 
          call free(%val(pcoeff))
          call disclaim(%val(pcoeff))
          call free(%val(pinptr))
          call disclaim(%val(pinptr))
          pcoeff = 0
          pinptr = 0
        pcoeff = malloc(%val(nnn*8))
        pinptr = malloc(%val(nnn*4))
        end if
      end if

        write(0,*)
        RETURN
        END

