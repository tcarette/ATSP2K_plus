*-----------------------------------------------------------------------
*      E I G _ O U T
*-----------------------------------------------------------------------

*     This routine writes first the eigenvectors, and then
*   the eigenvalues to Term.l 

      subroutine eig_out(iblock,wt,e,z,isom_shift)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
        PARAMETER (NWD=70,NOD=220,NOFFD=800)
        parameter (MTERM=20, MEIG=20)
*
      character string*64
        CHARACTER EL*3,ATOM*6,TERMs(20)*6
        COMMON/LABEL/ EL(NWD),ATOM,TERMs
      logical leigen(meig,mterm), lguess
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm)
      character term_bl(mterm)*3
      integer q(8)
      CHARACTER elc(8)*3, couple(15)*3
      character (len = 64)	:: config_tmp
      character (len = 66), allocatable, dimension(:) :: config0
      character (len = 66)      :: config(MEIG)
      integer, dimension(2)	:: max_comp(MEIG)
      integer, dimension(1)	:: mx
      integer		 	:: max_read, lstring, unit_n
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess
     
      double precision wt(*),e(*)
      double precision isom_shift(MEIG,nblock)

      unit_n = 40 + iblock
      rewind(unit_n)

      ncfg = ncfg_bl(iblock) 
      neig = 0
      max_comp = 0
      com_order = 0

      do i = 1,nume(iblock)
         if (leigen(i,iblock)) neig = neig+1
      end do

      write(60,'(I5,3X,A3)') neig, term_bl(iblock)

      do i = 1,nume(iblock)
         if (leigen(i,iblock)) write(60,'(I8,F15.9/(7F11.7))') 
     :      i, e(i), (wt((i-1)*ncfg+j), j=1,ncfg)
      end do

*     ... read 2 blank lines from file 'cfg.inp'
      Read (30,'(A64)') string
      Read (30,'(A64)') string

*<<< find the configurations with max comp of the eigen vector
      ics = 1
      ice = ncfg 
*     ... find the locations of extrema for all eigenvalues
      do ic1 = 1, nume(iblock)
	mx = MaxLoc(abs(wt(ics:ice)))
        max_comp(ic1) = mx(1)
	ics = ics + ncfg
        ice = ice + ncfg 
      end do

*     ... find the maximum of max_comp to read max_read lines  
      max_read = MaxVal(max_comp)
*
*     allocate enough space to read the cfg
      allocate(config0(max_read+1));

       !print*, max_read , '::::::::'

*     ... read all lines.le.max_read of file 'cfg.inp'
      do i = 1, max_read 
        READ(30,'(8(1X,A3,1X,I2,1X))') (ELC(j),Q(j),j=1,8)
        READ(30,'(15(1X,A3))') (COUPLE(J),J=1,15)
        
*     ... find NOCC
        nocc = 0
        do while(ELC(nocc+1).ne. '   ')
           NOCC = NOCC + 1
        end do

*     ... pack config 
        call pack(nocc,elc,q,couple,config_tmp)
*        lstring = Len_Trim(config_tmp)
        config0(i) = trim(config_tmp)
      end do
*>>>
*     ... order each configurations to correspond to max comp of the eig vector
      do i = 1, nume(iblock)
        config(i) = config0(max_comp(i))
      end do
      
      deallocate(config0);
      jj = 0
      nnel = 0
      write(unit_n,'(2X,A6,A,F5.1,A,I3,A,I6)' ) ATOM,'  Z = ',Z,
     :        '  NEL = ', nnel, '   NCFG = ',NCFG
      WRITE (unit_n, '(//A8,I4,2X,A8,I4)' ) '  2*J = ',JJ,
     :     'NUMBER =', neig;
*     ... write $TERM.l file for each block
      do i = 1, nume(iblock);
         if (leigen(i,iblock)) then;
           config_tmp = config(i) 
           lstring = Len_Trim(config_tmp)
           write (unit_n,'(4x,A8,f16.9)') 'Ssms = ',
     :           isom_shift(i,iblock);
           WRITE (unit_n,'(i6,f16.9,2x,A)') 
     :                 max_comp(i), e(i), config_tmp(1:lstring) 
           WRITE (unit_n,'(7F11.8)') (wt((i-1)*ncfg+j), j=1,ncfg)
*          WRITE(iscw,'(/1X,F19.12,2X,A/(1X,7F11.8))' ) EIGVAL
         end if
      end do
*#####
      do while (string(1:1).ne.'*')
         Read (30,'(A64)') string
      end do

      return
      end
