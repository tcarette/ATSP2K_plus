!Program:       FLS_trends
!Author:        Andrei Irimia
!Affiliation:   Department of Electrical Engineering and Computer Science
!               Vanderbilt University
!               Nashville, Tennessee, USA
!------------------------------PROGRAM FLS_TRENDS-------------------------------
program FLS_trends

implicit none
!Variable declarations
logical      		 :: isthere, lenOne = .false.
character (len = 12) :: filename,   filename12
character (len = 11) :: filename11, filenames
character (len = 10) :: filename10
character (len =  9) :: filename9
integer 		 :: ELcntr, Zcntr, cntr = 0, k = 0
double precision, dimension (109) :: masses
double precision :: mass
character (len =  1), dimension (10) :: digits
character (len =  1), dimension (14) :: elements1
character (len =  2), dimension (95) :: elements2

digits    = (/"0","1","2","3","4","5","6","7","8","9"/)
elements1 = (/"H","B","C","N","O","F","P","S","K","V","Y","I","W","U"/)
elements2 = (/"He","Li", "Be", "Ne", "Na", "Mg", "Al", "Si", "Cl", "Ar", &
	     "Ca", "Sc", "Ti", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", &
	     "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Zr", "Nb", &
	     "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", &
	     "Te", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", &
	     "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", &
	     "Ta", "Rh", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", &
	     "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "Np", "Pu", &
	     "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", &
	     "Db", "Sg", "Bh", "Hs", "Mt"/)

masses = (/1.00794,4.00260,6.941,9.01218,10.811,12.0107,14.00674,15.9994,& 
	18.99840,20.1797,22.9897,24.3050,26.98154,28.0855,30.97,32.066,  &
	35.4527,39.0983,39.948,40.078,44.95591,47.867,50.9415,51.9961,   &
	54.93805,55.845,58.93320,58.6934,63.546,65.39,69.723,72.61,      &
	74.92160,78.96,79.904,83.80,85.4678,87.62,88.90585,91.224,       &
	92.90638,95.94,98,101.07,102.90550,106.42,107.8682,112.411,      &
        114.818,118.710,121.760,126.90447,127.60,131.29,132.90545,       &
	137.327,138.9055,140.116,140.90765,144.24,145,150.36,151.964,    &
	157.25,158.92534,162.50,164.93032,167.26,168.93421,173.04,       &
	174.967,178.49,180.9479,183.84,186.207,190.23,192.217,195.078,   &
	196.96655,200.59,204.3833,207.2,208.98038,209,210,222,223,226,   &
	227,232.0381,231.03588,237,238.0289,244,243,247,247,251,252,257, &
        258,259,262,261,262,263,264,265,268/)

lenOne = .true.

do
	ELcntr = ELcntr + 1
	if (lenOne) then
		if (ELcntr > 14) then
			lenOne = .false.	
			ELcntr = 1
		end if
	else 
		if (ELcntr > 95) exit
	end if
	Zcntr = 0
	do
		Zcntr = Zcntr + 1
		if (Zcntr > 112) exit
		filename (1:3) = "LS_"
		mass = masses (Zcntr)
		if (lenOne) then
			filename (4:5) = elements1 (ELcntr)
			filename (5:7) = "_Z_"
			if (Zcntr < 10) then
				filename ( 8: 9) = digits (Zcntr + 1)
			else if (Zcntr > 9 .and. Zcntr < 100) then
				filename ( 8: 9) = digits (Zcntr /  10 + 1)
				filename ( 9:10) = digits (mod (Zcntr, 10) + 1)
			else if (Zcntr > 99) then
				filename ( 8: 9) = digits (Zcntr / 100  + 1)
				filename ( 9:10) = digits (Zcntr /  10  + 1)
				filename (10:11) = digits (mod (Zcntr, 10) + 1)
			end if
		else
			filename (4:6) = elements2 (ELcntr)
			filename (6:9) = "_Z_"
			if (Zcntr < 10) then
				filename ( 9:10) = digits (Zcntr + 1)
			else if (Zcntr > 9 .and. Zcntr < 100) then
				filename ( 9:10) = digits (Zcntr /  10 + 1)
				filename (10:11) = digits (mod (Zcntr, 10) + 1)
			else if (Zcntr > 99) then
				filename ( 9:10) = digits (Zcntr / 100  + 1)
				filename (10:11) = digits (Zcntr /  10  + 1)
				filename (11:12) = digits (mod (Zcntr, 100) + 1)
			end if			
		end if

		cntr = 13
		do
			cntr = cntr - 1
			if (cntr < 1) exit
			if (filename (cntr:cntr) == " " .and. &
			    filename (cntr - 1:cntr - 1) /= " ") then
				if (cntr ==  9) filename9  (1: 9) = &
					filename (1: 9)
				if (cntr == 10) filename10 (1:10) = &
					filename (1:10)
				if (cntr == 11) filename11 (1:11) = &
					filename (1:11)
				if (cntr == 12) filename12 (1:12) = &
					filename (1:12)

				if (lenOne) then
					do k = 1, 11
						filenames (k:k+1) = filename(k:k+1)
					end do
				end if

				!see if the file is there
				if (lenOne) then
					inquire (file = filenames, exist = isthere)
				else
					inquire (file = filename,  exist = isthere)
				end if
				if (isthere) then
					if (lenOne) then
						call process (filenames, mass)
					else
						call process (filename, mass)
					end if
				end if
				exit
			end if
		end do
	end do
end do
end program FLS_trends 

!------------------------------SUBROUTINE PROCESS------------------------------
subroutine process (filename, mass)
implicit none

character (len = *), intent (inout) :: filename
character (len = len(filename) + 1) :: anewfile

logical			  		 :: outisthere
integer               :: in = 3, out = 6, ndx1 = 0, ndx2 = 0
integer,dimension(2,4):: limits
double precision, intent (in) :: mass
double precision :: line, lower, upper, lowerT, upperT,lowest1,lowest2, &
				 		  tr_en = 0
integer :: ios, n, z, nT, zT, zmax, curn, cntr, ELcntr, Zcntr, lencntr, & 			
			  cntrunit, i, format, ind, length
character (len = 80)  :: input, inputT, cfg_lbl1, cfg_lbl2
character (len = 12)  :: EL, EU, ELm, EUm
character (len = 10)  :: SL, SV, GFL, GFV, SLm, SVm, GFLm, GFVm
character (len =  6)  :: ERR, ERRm
character (len =  2)  :: trm1 = "  ", trm2 = "  "

!assign name of new file
anewfile = "F" // filename

!Ensure that configuration labels are blank
cfg_lbl1 = ' '
cfg_lbl2 = ' '

!assign positive values to the energy levels
lowest1 = 1
lowest2 = 1
tr_en   = 0

open (unit = in, pad = "YES", file = filename, status = "OLD", iostat = ios)

!see if the file is there
inquire (file = anewfile, exist = outisthere)

if (outisthere) then
	open (unit = out, file = anewfile, status = "OLD", iostat = ios)
else
	open (unit = out, file = anewfile, status = "NEW", iostat = ios)
end if

101 format (i3,i2,f14.7,f13.7,a10,a10,a10,a10,a6)

201 format (i3,i2,f12.7,f11.7,a10,a10,a10,a10,a6)
202 format (i3,i2,f13.7,f12.7,a10,a10,a10,a10,a6)
203 format (i3,i2,f14.7,f13.7,a10,a10,a10,a10,a6)
204 format (i3,i2,f15.7,f14.7,a10,a10,a10,a10,a6)

301 format (1xi2,1xi1,f13.7,1xf12.7,1xa9,1xa9,1xa9,1xa9,1xa5)
302 format (t1a79)
303 format (a80)
304 format (t1a33)

do
	read (unit = in, iostat = ios, fmt = 303) input
	ind = index (input, "transitions")

	if (input (2:2) /= "-" .and. input (2:2) /= "Z" .and. &
		 input (3:3) /= " " .and. input (4:4) == " " .and. &
		 ind == 0)  then
		backspace 3
		read (in, *) z, n, lower, upper, SL, SV, GFL, GFV, ERR
		backspace 3

		!see how long the first number is
		if 	  (lower <     0 .and. lower >     -10) then
			length = 1
			read (in, *) z, n, lower, upper, SL, SV, GFL, GFV, ERR
		else if (lower <   -10 .and. lower >    -100) then
			length = 2
			read (in, *) z, n, lower, upper, SL, SV, GFL, GFV, ERR
		else if (lower <  -100 .and. lower >   -1000) then
			length = 3
			read (in, *) z, n, lower, upper, SL, SV, GFL, GFV, ERR
		else if (lower < -1000 .and. lower > -100000) then
			length = 4
			read (in, *) z, n, lower, upper, SL, SV, GFL, GFV, ERR
		end if
	else
		n = 0
	end if

	if (ind /= 0) then
		write (unit = out, fmt = 302) input
		read (unit = in, iostat = ios, fmt = '(a80)') inputT
		write (unit = out, advance = "NO", fmt = 302) inputT
		write (unit = out, fmt = '(a)') "--------------"
		exit
	end if

	if (n == 0) then
		if (curn > 0) then
			!print current values
			call extract_energy (tr_en, trm1, trm2, zmax, curn, cfg_lbl1, &
					cfg_lbl2, lowest1, lowest2, mass)

			call println (zmax, curn, lowest1, lowest2, tr_en, SLm, & 
				SVm, GFLm, GFVm, ERRm, length)
			tr_en   = 0
			zmax	  = 0
			curn    = 0
			lowest1 = 0
			lowest2 = 0
			SL		  = ''
			SV		  = ''
			GFL	  = ''
			GFV	  = ''
			ERR	  = ''
		end if

		if (input (2:2) == "-") then
			write (unit = out, advance = "NO", fmt = 302) input
			write (unit = out, fmt = '(a)') "--------------"
		else if (input (2:2) == "Z") then
			write (unit = out, advance = "NO", fmt = 304) input (1:33)
			write (unit = out, advance = "NO", &
				fmt = '(a15)') "       TE      "
			write (unit = out, fmt = '(a)') input (34:80)
		else
			write (unit = out, fmt = 303) input

			!find what the terms for the current configuration are
			if (input (4:4) /= " ") then
				!find the first term
				ndx1 = index (input (4:80), " ")
				trm1 = input ((ndx1 + 1):(ndx1 + 2))

				if (input (1:1) == " ")	then
					cfg_lbl1 = input (3:(ndx1 + 2))
				else 
					cfg_lbl1 = input (1:(ndx1 + 2))
				end if

				!find the second term
				call swap (input)
				ndx2 = index (input, "_")
				call swap (input)
				trm2 = input (80 - ndx2 + 1:80 - ndx2 + 2)
				cfg_lbl2 = input (ndx1 + 4:80 - ndx2 + 2)
			end if
		end if

		curn = 0
	else
		if (n < curn .and. curn /= 0) then
			call extract_energy (tr_en, trm1, trm2, zmax, curn, &
				cfg_lbl1, cfg_lbl2, lowest1, lowest2, mass)
			call println (zmax, curn, lowest1, lowest2, tr_en, SLm, &
				SVm, GFLm, GFVm, ERRm, length)
			tr_en = 0
			!set max values to current values
			curn    = n
			zmax	= z
			lowest1 = lower
			lowest2 = upper
			SLm 	= SL
			SVm 	= SV
			GFLm	= GFL
			GFVm	= GFV
			ERRm	= ERR
		end if
		if (curn == 0) then
			!set max values to current values
			curn    = n
			zmax	= z
			lowest1 = lower
			lowest2 = upper
			SLm 	= SL
			SVm 	= SV
			GFLm	= GFL
			GFVm	= GFV
			ERRm	= ERR
		end if
		if (n > curn .and. curn /= 0) then
			!print current values
			call extract_energy (tr_en, trm1, trm2, zmax, curn, &
				cfg_lbl1, cfg_lbl2, lowest1, lowest2, mass)
			call println (zmax, curn, lowest1, lowest2, tr_en, SLm, &
				SVm, GFLm, GFVm, ERRm, length)
			tr_en = 0
			!set curn to n
			curn    = n
			zmax	= z
			lowest1 = lower
			lowest2 = upper
			SLm 	= SL
			SVm 	= SV
			GFLm	= GFL
			GFVm	= GFV
			ERRm	= ERR
		end if
		if (n == curn .and. n /= 0) then
			if (lower <= lowest1 .and. upper <= lowest2) then
				lowest1 = lower
				lowest2 = upper
				SLm 	= SL
				SVm 	= SV
				GFLm	= GFL
				GFVm	= GFV
				ERRm	= ERR
			end if
		end if
	end if

	! Negative ios means end of file
	if (ios  <  0) exit
end do
endfile out
end subroutine process
!-------------------------SUBROUTINE EXTRACT_ENERGY-----------------------------
subroutine extract_energy (tr_en, trm1, trm2, zmax, curn, cfg_lbl1, cfg_lbl2, &
	E1, E2, mass)
implicit none

double precision, intent (inout) :: tr_en
double precision, intent (in)    :: E1, E2, mass
double precision 						:: Ry, Me
character (len = 2), intent (in) :: trm1, trm2
character (len = *), intent (in) :: cfg_lbl1, cfg_lbl2
character (len = 2)  :: zstr = "  "
integer, intent (in) :: zmax, curn
integer :: i = 0, j = 0, k = 0, finished = 0, numcfg = 0
logical :: doesexist = .false., smcfg = .false., fromf = .false.
character (len = 19) :: name
character (len = 1), dimension (9) :: digits

Ry       = 109737.31534
Me       = 548.579903E-6
tr_en    = 0
finished = 0
i 	 = 0
j        = 0
k        = 0
fromf    = .true.

Ry = Ry / (1 + Me / mass)
tr_en = abs (E1 - E2) * 2 * Ry

end subroutine extract_energy

!---------------------------SUBROUTINE PRINTLN----------------------------------
subroutine println (zmax, curn, lowest1, lowest2, tr_en, SLm, SVm, GFLm, &
	GFVm, ERRm, length)
implicit none
double precision, intent (in)	   	:: lowest1,lowest2, tr_en
integer, intent (in)		  		:: zmax, curn, length
character (len = 10), intent (in)  	:: SLm, SVm, GFLm, GFVm
character (len =  6), intent (in)  	:: ERRm

!see how long the first number is
if 	(length == 1) then
	write (6, 201) zmax, curn, lowest1, lowest2, tr_en, SLm, &
		SVm, GFLm, GFVm, ERRm 
else if (length == 2) then
	write (6, 202) zmax, curn, lowest1, lowest2, tr_en, SLm, &
		SVm, GFLm, GFVm, ERRm
else if (length == 3) then
	write (6, 203) zmax, curn, lowest1, lowest2, tr_en, SLm, &
		SVm, GFLm, GFVm, ERRm
else if (length == 4) then
	write (6, 204) zmax, curn, lowest1, lowest2, tr_en, SLm, &
		SVm, GFLm, GFVm, ERRm
end if

201 format (i3,i3,f11.7,f16.7,f14.2,1xa10,a10,a10,a10,a6)
202 format (i3,i3,f12.7,f15.7,f14.2,1xa10,a10,a10,a10,a6)
203 format (i3,i3,f13.7,f14.7,f14.2,1xa10,a10,a10,a10,a6)
204 format (i3,i3,f14.7,f14.7,f14.2,1xa10,a10,a10,a10,a6)

end subroutine println
!-------------------------------SUBROUTINE SWAP---------------------------------
subroutine swap (input)
implicit none
integer :: i = 0
character (len = 80), intent (inout) :: input
character (len = 80) line2
line2 = input
do i = 0, 79
	input (i:i) = line2 (80 - i:80 - i)
end do
end subroutine swap
!-------------------------------------------------------------------------------
