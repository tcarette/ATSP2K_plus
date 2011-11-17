	Program  Select
	CHARACTER*74  LIST
	CHARACTER*72 COUPL
	CHARACTER*2 el(30)
	CHARACTER*24 name
        
	nwf = 1
	write(0,*) 'Enter electrons, one per line: blank terminates'
 1	read(5,'(A2)') el(nwf)
	if (el(nwf) .ne. ' ') then
	  nwf = nwf+1
	  go to 1
	end if
	nwf = nwf -1
	write(0,*) 'Searching for electrons' ,(el(i),i=1,nwf)
	write(0,*) 'Enter name of file to be searched'
	read(5,'(A)') name
	open(unit=20,file=name, status='old')
10      read(20,'(A74)',END=99)  list
	read(20,'(A72)',END=99) coupl
	do 20 i = 1,nwf
	   if ( index(list,el(i)) .ne. 0) then
	      write(6,'(A74)') list
	      k = 72
30            if (coupl(k:k) .eq. ' ') then
		 k = k-1
	         go to 30
	      end if
	      write(6,'(A)') coupl(1:k)
	      go to 10
	   end if
20      continue
	go to 10
99      stop
	end
