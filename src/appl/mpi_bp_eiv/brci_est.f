*-----------------------------------------------------------------------
*     brci_est
*-----------------------------------------------------------------------
*     This routine reads .j and .c vectors from lower n and generates
*     an estimate. The program needs <name>.j.est and <name>.c.est to be 
*     present in the working directory 
*     THis routin has been tested with the parallel code only
*     and it was foiund that for very large cases (100K), it may result 
*     in 30% better performance 
*
      subroutine brci_est(JJ,A,ncfg,nume,iuc,file_base,lest)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      parameter (MTERM=31, MEIG=31)
      character*64 string_c, string_e, s_tmp, string
      CHARACTER EL*3,ATOM*6,TERMs(MTERM)*6
!      COMMON/LABEL/ EL(NWD),ATOM,TERMs
      logical leigen(MEIG)
      character term_bl(mterm)*3
      integer q(8),iouj,iuc
      character (len = 64)      :: config_tmp
      character (len = 66)      :: config(MEIG)
      integer, dimension(2)     :: max_comp(MEIG)
      integer, dimension(1)     :: mx
      integer                   :: MAX_read, lstring
      integer                   :: ui_c_est, ui_j_est
      double precision A(nume*ncfg+nume)
      integer JJ, ncfg, nume
      character*128 file_c,file_c_est,file_j_est,file_base
      logical  l_match, lest 
      integer len_atom, len_str

      integer Z_NUMBER, NCFG_EST
      integer J_EST, NUME_EST 
      double precision, allocatable, dimension(:) :: VECT_EST
      character (len = 66), allocatable, dimension(:) :: config_est 
      character (len = 66), allocatable, dimension(:) :: coupl_est
      character*64 string_co, string_cu
     
      l_match = .false.
      ui_c_est = 44
      ui_j_est = 45

      len_atom = LEN_TRIM(file_base)
      
      file_c_est = file_base(1:len_atom)//'.c.est'
      file_j_est = file_base(1:len_atom)//'.j.est'

      inquire(FILE=file_c_est, exist = lest)
      if (.not.lest) then
         print*, ' No initial estimates: ', file_c_est
         return
      end if

      inquire(FILE=file_j_est, exist = lest)
      if (.not.lest) then
         print*, ' No initial estimates: ', file_j_est
         return
      endif

      OPEN(UNIT=ui_c_est,FILE=file_c_est,STATUS='OLD')
      OPEN(UNIT=ui_j_est,FILE=file_j_est,STATUS='OLD')

      !write 
!     ... read estimate: *.j.est
      rewind (ui_j_est)
      read(ui_j_est,'(A4,i2,A34,i7)') s_tmp, Z_NUMBER, s_tmp, NCFG_EST
      read(ui_j_est,'(A64)') s_tmp
      read(ui_j_est,'(A64)') s_tmp 

*     ... find the corresponding J 
      J_EST = -1;
      do while (J_EST.NE.JJ) 
         read(ui_j_est,'(A64)') string
         len_str = len_trim(string)
         if (len_str == 26) then  
            read(string,'(A8,i4,A10,i4)') s_tmp,J_EST,s_tmp,NUME_EST 
         end if
      end do

      if(J_EST==JJ) then
          A(1:nume*ncfg) = 0;
      else 
          print*, 'No 2*J = ', JJ, ' found'
          lest = .false.
          return
      end if
      
*     allocate memory for VECT_EST = NUME_EST*NCFG_EST
      allocate(VECT_EST(NUME_EST*NCFG_EST+NUME_EST));
     
*     ... allocate memory for conf+coupl from the estimate
      allocate(config_est(ncfg));
      allocate(coupl_est(ncfg));
        
*     ... read all eigenvectors (NUME_EST) for the required J
      iin = min(NUME_EST,nume)
      do i = 1,iin
         read(ui_j_est,'(A64)') string
         read(ui_j_est,'(A6,F16.9,A42)') 
     :              string,VECT_EST(NUME_EST*NCFG_EST+i),string
         read(ui_j_est,'(7F11.8)') 
     :                (VECT_EST((i-1)*NCFG_EST+j), j=1,NCFG_EST)    
      end do
      
*     ... read *.c.est (configurations of the estimates)
      read(ui_c_est,'(A64)') string 
      read(ui_c_est,'(A64)') string 
      do i = 1, NCFG_EST
         read(ui_c_est,'(A64)') config_est(i)
         read(ui_c_est,'(A64)') coupl_est(i)
      end do

*     ... read 2 blank lines from file 'clist'
      rewind (iuc);
      Read (iuc,'(A64)') string_c
      Read (iuc,'(A64)') string_c
*
      ir = 0;
*     ... map all estimates 
      do i = 1,NCFG_EST 
        do while (.not.l_match) 
           READ(iuc,'(A64)',end = 1) string_co
           READ(iuc,'(A64)',end = 1) string_cu
           !print*, string_co, 'string_co'
           !print*, string_cu, 'string_cu'
           !print*,config_est(i), 'config_est(i)'
           !print*, coupl_est(i), 'coupl_est(i)'
           ir = ir + 1;
           if(lle(string_co,config_est(i))
     :           .and.lge(string_co,config_est(i))) then
              if(lle(string_cu,coupl_est(i)) 
     :           .and.lge(string_cu,coupl_est(i))) then
                 !print*, 'match'
                 l_match = .true.
                 do in = 1,NUME_EST
                    A(ncfg*(in-1)+ir) = 
     :                   VECT_EST(NCFG_EST*(in-1)+i)
                    !print*, VECT_EST(NCFG_EST*(in-1)+i)
                 end do 
              end if
           end if
        end do
        l_match = .false.
1       continue 
      end do

      do i = 1, NCFG_EST
         A(NCFG*nume+i) = VECT_EST(NCFG_EST*NUME_EST+i)
      end do 

      close(ui_c_est);
      close(ui_j_est);

      deallocate(config_est);
      deallocate(coupl_est);
      deallocate(VECT_EST);

      return
      end

