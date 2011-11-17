!Program name:		T_dependence
!File name:			T_dependence.f90
!Author:				Andrei Irimia
!Affiliation:		Department of Electrical Engineering and Computer Science
!						Vanderbilt University
!						Nashville, Tennessee, 37235
!Language:			ForTran 90
!Purpose:			This program determines the term dependence between a number
!						of LS terms for a given set of configurations. A dynamic
!						array of linked lists is created for the leading percentages
!						associated with each term.
!-------------------------------------------------------------------------------
!--------------------------------MODULE TERM------------------------------------
module term
	implicit none
	!type definitions
	type LS_term
		double precision	:: energy, coef
		character (len = 2) :: term_label
		type (LS_term), pointer :: next
	end type LS_term
contains
	!------------------SUBROUTINE MAKE FOR TYPE LS_TERM--------------------
	!this subroutine is for making a term for the linked list of dependants
	subroutine make (energy, coef, term_label, term_obj)
		implicit none
		!Arguments
		integer :: error
		character (len = 2), intent (in) :: term_label
		double precision, intent (in) 	:: energy, coef
		type (LS_term), pointer 			:: term_obj

		!allocate memory for this new term object
		allocate (term_obj, stat = error)
		if (error /= 0) then
			print *, "Machine out of memory - exiting."
			stop
		end if

		!Note: the coefficient value must be squared before it is 
		!added, see below
		term_obj % coef 	  	 = coef
		term_obj % energy     = energy
		term_obj % term_label = term_label

		!no successor to current term
		nullify (term_obj % next)
	end subroutine make
end module term
!---------------------------MODULE DP_LINKED_LIST-------------------------------
module DP_linked_list
	use term
	implicit none
contains
	!-------------------------SUBROUTINE INITIALIZE-------------------------
	subroutine initialize (head, tail)
		implicit none
		!initialize the empty list
		type (LS_term), pointer :: head, tail
		!there exists no successor
		nullify (head, tail)
	end subroutine initialize
	!---------------------------SUBROUTINE APPEND---------------------------
	subroutine append (head, tail, term_obj)
		implicit none
		!Add a new item to the end of the list
		type (LS_term), pointer :: head, tail, term_obj, iterator

		if (associated (head)) then
			!first, see if this term has already been added
			iterator => head
			do
				!if we reached the end of the segment, the term
				!is not in the segment so we need to add it
				!Note: it is not necessary to add the coefficient
				!because this is the first value
				if (.not. associated (iterator)) then
					!add a term, but first allocate the iterator
					tail % next => term_obj
					nullify (term_obj % next)
					tail => term_obj 
					exit
				end if
				!if there is a value that has been already added to 
				!the desired term

				if (iterator % term_label == term_obj % term_label) then
					!add the coefficient value squared to the value
					!already there
					iterator % coef = iterator % coef + term_obj % coef
					!!!!!!!!!!!!!!!
					if (iterator % term_label == "3P") print *, iterator % coef
					!get rid of the memory since we don't need this 
					!anymore
					nullify (term_obj)
					exit
				end if
				!continue with the next item
				iterator => iterator % next
			end do
		else
			!Add a new item to the end of the list, i.e., add an LS_term
			head => term_obj
			nullify (term_obj % next)
			tail => term_obj
		end if
	end subroutine append
	!--------------------------SUBROUTINE DISPLAY---------------------------
	subroutine display (head)
		implicit none
		!List the contents of the list
		type (LS_term), pointer :: head, ptr
		integer 				:: count = 0, fact = 0
		double precision	:: energy, swapcoef
		character (len = 2):: label = "  "

		if (associated (head)) then
			!Set current pointer to head of list
			ptr => head
			do
				!print details of current request item
				if (.not. associated (ptr % next)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					write (*, advance = "YES", '(f16.2,a1,a2)'), &
						100 * ptr % coef, "-", ptr % term_label
				else
					write (*, advance = "NO", &
						'(f16.2,a1,a2,x1a1)'), 100 * ptr % coef, &
						"-", ptr % term_label, " "
				end if
				ptr => ptr % next
				if (.not. associated (ptr)) exit
			end do
		end if
	end subroutine display
end module DP_linked_list
!--------------------------------MODULE LS_SEGMENT------------------------------
module LS_segment
	use DP_linked_list
	implicit none
	type segment
		character (len = 2)	:: LS_label
		integer					:: twoJ
		character (len = 80)	:: cnf
		!A new linked list for the terms on which this term depends
		type (LS_term), pointer	:: head, tail
		type (segment), pointer	:: next
	end type segment
end module LS_segment
!------------------------------MODULE LS_LINKED_LIST----------------------------
!A linked list of LS terms and their mixing coefficients
module LS_linked_list
	use LS_segment
	use DP_linked_list
	implicit none
contains
	!---------------------SUBROUTINE INITIALIZE-----------------------------
	subroutine initialize_segment (head, tail)
		implicit none
		!initialize the empty list
		type (segment), pointer :: head, tail, next
		!there exists no successor
		nullify (head, tail, next)
	end subroutine initialize_segment
	!-----------------------SUBROUTINE MAKE_SEGMENT-------------------------
	!this subroutine is for making a term for the linked list of 
	!actual terms
	subroutine make_segment (LS_label, twoJ, term_obj, wrapper, cnf)
		implicit none

		!Variable declarations
		character (len = 2), intent (in):: LS_label
		integer, intent (in)			:: twoJ
		character (len = 80),intent (in):: cnf
		type (LS_term), pointer			:: term_obj
		type (segment), pointer			:: wrapper
		integer							:: error

		allocate (wrapper, stat = error)

		if (error /= 0) then
			print *, "Machine out of memory - exiting."
			stop
		end if

		!Copy the data from the argument to the data structure
		!A new linked list for the terms on which this term depends
		wrapper % head 		=> term_obj
		wrapper % tail 		=> term_obj
		wrapper % LS_label 	= LS_label
		wrapper % twoJ 		= twoJ
		wrapper % cnf		= cnf
	
		!no next element
		nullify (wrapper % next)
		nullify (wrapper % tail % next)
	end subroutine make_segment
	!-----------------------SUBROUTINE APPEND_SEGMENT-----------------------
	subroutine append_segment (head, tail, LS_label, twoJ, wrapper)
		implicit none
		type (segment), pointer :: head, tail, wrapper, iterator
		character (len = *), intent (in) :: LS_label
		integer, intent (in) :: twoJ
		integer :: error
		type (LS_term), pointer :: newhead, newtail, newobj

		if (associated (head)) then
			!first, seek an LS_segment with the same LS_label
			iterator => head

			do
				!if we reached the end of the list by now, the label
				!is not in the list so we need to add it
				if (.not. associated (iterator)) then
					iterator => wrapper
					exit
				else
					!if there is a list already with the same label
					if (iterator % LS_label == LS_label .and. &
						iterator % twoJ     == twoJ) then
						!tell the segment to add term_obj to its tail, 
						!i.e. we call the append function of the small 
						!LS_terms the append function receives as 
						!arguments only pointers to objects of type 
						!LS_term
						call make (wrapper % head % energy, 	& 
								   wrapper % head % coef,		&
								   wrapper % head % term_label,	&
								   newobj)
						call append (iterator % head, iterator % tail, &
							newobj)
						exit
					end if					
				end if

				!continue with the next item
				iterator => iterator % next
			end do
		else
			!Add a new item to the end of the list, i.e., add a segment
			head => wrapper
			nullify (wrapper % next)
			head => wrapper
		end if
	end subroutine append_segment
	!----------------------SUBROUTINE DISPLAY_SEGMENTS----------------------
	subroutine display_list (head, cfgl, maxcfgl)
		implicit none
		!tell each linked list in the set of segments to display itself
		!List the contents of the list
		type (segment), pointer :: head, ptr
		integer, intent (in)	:: cfgl, maxcfgl
		
		if (associated (head)) then
			!Set current pointer to head of list
			ptr => head
			do
				if (cfgl) then
					write (*, advance="NO", fmt = '(a5,i4,a2)'), &
						ptr % LS_label, ptr % twoJ, "  "
					write (*, advance="NO", '(a)') ptr % cnf &
						(1:maxcfgl -2)
					write (*, advance = "NO", '(f16.8,a3)') &
						ptr % head % energy, '   '
				else
					write (*, advance="NO", fmt = '(a5,i4)'), &
						ptr % LS_label, ptr % twoJ
					write (*, advance = "NO", '(f16.8,a3)') &
						ptr % head % energy, '   '
				end if
				call display (ptr % head)
				ptr => ptr % next
				if (.not. associated (ptr)) exit
			end do
		end if
	end subroutine display_list
end module LS_linked_list
!---------------------------PROGRAM T_DEPENDENCE-------------------------------
program T_dependence
use LS_linked_list
implicit none

!Variable declarations

integer :: ios1, ios2, Cfile = 1, Jfile = 2, ncfg, twoJ = 0, ind, status, &
	i = 1, frmt, Jnum, j = 0, k = 0, limit = 0, cfgread = 0, end = 0, &
	Jcntr = 0, eof = 0, cnf = 0, cfgl = 0, maxcfgl = 0, begin = 0, m = 0
character (len = 30) :: Cfilename, Jfilename, filename
character (len = 80) :: line1 = ' ', line2 = ' ', line3 = ' ', line4 = ' ', &
	line5 = ' ', message = ' ', config = ' '
character (len = 2)  :: LSterm, LS_label = '  '
character :: cfgans = ' '
double precision	 :: energy, coef, all = 0
type (segment), pointer :: head, tail, wrapper
type (LS_term), pointer :: term_obj

write (0, '(a21)') "Enter name of .c file"
read  *, filename
write (0, '(a33)') "Print configuration labels? (y/n)"
read  *, cfgans

cfgl = 0
if (cfgans == 'y' .or. cfgans == 'Y') cfgl = 1

ind	  = index (filename, " ")
Cfilename = filename (1:ind - 1) // ".c"
Jfilename = filename (1:ind - 1) // ".j"

open (unit = Cfile, file = Cfilename, status = "OLD", iostat = ios1)
open (unit = Jfile, file = Jfilename, status = "OLD", iostat = ios2)

if (ios1 /=0) then
	write (0, '(a39)') "Input file with .c extension not found."
	stop
end if

if (ios2 /=0) then
	write (0, '(a39)') "Input file with .j extension not found."
	stop
end if

!read the number of configurations
read (Jfile,101) line1, ncfg

!Initialize the linked list of LS terms
call initialize_segment (head, tail)
nullify (wrapper)

!while the end of file has not been reached
!do loop for each set of values with a certain 2*J
do
	if (eof == 1) exit
	!read current value of 2 * J and number of terms for 2 * J
	!Ensure that eof has not been reached
	read (2,100) line1

	if (line1 (2:2) == "*" .or. line1 (2:3) == "EN") then
		eof = 1
		exit
	end if

	read (2,100) line1
	read (2,102) line1, twoJ, line2, Jnum

	do Jcntr = 0, Jnum - 1
		if (eof == 1) exit
		
		!nullify objects for current iteration
		nullify (head, tail, term_obj, wrapper)
		call initialize_segment (head, tail)

		rewind Cfile
		!Skip first 2 lines of configuration file
		read (Cfile,100) line1
		read (Cfile,100) line1

		!read header lines
		read (2,100) line1
		read (2,103) line2, energy, line1

		end 	 = 0
		call swap (line1)
		ind 	 = index (line1 (:80), "_")
		cnf 	 = index (line1 (ind:80), " ")
		call swap (line1)
		LS_label = line1 (80 - ind + 1:80 - ind + 2)
		config   = line1 (3:cnf + 1)
		
		!find maximum length of configuration string
		if (cnf + 1 > maxcfgl) maxcfgl = cnf + 1

		do j = 0, ((ncfg + 7) / 7) !was + 1
			if (end == 1 .or. eof == 1) exit

			if (j < (ncfg  + 7) / 7) then !was + 1 with < 
				if (end == 1 .or. eof == 1) exit

				!test next line
				line4 = ''
				read (2, 100) line4
				backspace 2

				do k = 0, 6
					if (end == 1) exit

					!read data from configuration file
					read (1, 100) line1
					read (1, 100) line1

					ind = index (line1 (:80), "    ")
					LSterm = line1 (ind - 2:ind - 1)

					!read coefficients from .j file
					!Some machines may not display 0.00, but  .00, in which case the
					!test below fails unless three " "'s are there
					if (index (line4, "   ") > 67) then !was >=
						if (mod (k, 6) == 0 .and. k /= 0) then
							read (2, 104, advance = "YES") coef
							
							!in some cases, the eigenvector ends at the
							!end of the line
							read (2, 100) line5
							backspace 2
							if (index (line5, "   ") < 11) then
								end = 1
								exit
							end if
						else
							read (2, 104, advance =  "NO") coef
						end if
					else
						if (mod (k, 6) <= (index (line4, "  ") + 10) &
										/ 11) then 
							read (2, 104, advance = "YES") coef
						else
							read (2, 104, advance =  "NO") coef
						end if
						end = 1
						exit
					end if

					if (coef /= 0) then
						!square the coefficient

						coef = coef * coef
						!append the term created to the end of the list
						call make(energy, coef, LSterm, term_obj)
						call make_segment(LS_label, twoJ, term_obj, &
							wrapper, config)
						call append_segment (head, tail, LS_label, &
							twoJ, wrapper)
						nullify (wrapper, term_obj)
					end if
				end do
			end if
		end do
		if (Jcntr == 0 .and. begin == 0) then
			!print headings
			begin = 1
			if (cfgl == 0) then
				print *, "Term  2*J      Energy      ", &
					"Leading percentages"
			else
				write (*, advance = "NO", fmt = '(a25)') &
					"Term  2*J   Configuration"
				do m = 0, maxcfgl - 11
					write (*, advance = "NO", '(a1)') " "
				end do
				print *, "Energy      Leading percentages"
			end if			
		end if
		call display_list (head, cfgl, maxcfgl)
	end do
end do

100 format (a80)
101 format (a39,i7)
102 format (a7,i5,a10,i4)
103 format (a6,f15.9,a50)
104 format (1xf10.8)

end program T_dependence
!-----------------------------SUBROUTINE ALPHABETIZE---------------------------
subroutine alphabetize (dep, LSterm, CLSterm, TwoJ, E)
implicit none

integer, intent (in) :: TwoJ
integer :: i = 0, added = 0
character (len = 2), intent (in) :: LSterm, CLSterm
real, intent (in) 	 :: E
real, dimension (100, 50, 2), intent (inout) :: dep

character (len = 1), dimension (10) :: terms = (/"S", "P", "D", "F", "G", &
	"H", "I", "K", "L", "M"/)

!do i = 1, 100
	!if we hit a blank space, add current and quit
	!if (sum (i, 1) == "  ") then
		!add the current
		!sum (i, 1) = LSterm
		!sum (i, 2) = E
		!dep (i, 1) = CLSterm !need to have a class
		!exit
	!end if

	!search for the wanted term
	!if (sum (i, 1) == LSterm) then
		!we have found the term TO which we add
		!added = 1
	!end if
!end do
end subroutine alphabetize
!-------------------------------SUBROUTINE SWAP---------------------------------
subroutine swap (line1)
implicit none
integer :: i = 0
character (len = 80), intent (inout) :: line1
character (len = 80) line2
line2 = line1
do i = 0, 79
	line1 (i:i) = line2 (80 - i:80 - i)
end do
end subroutine swap
!-------------------------------------------------------------------------------
