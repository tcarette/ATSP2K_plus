*
* Read the lists produced by nonh and show the angularly integrated 
* Hamiltonian to a human
*
* Written by Thomas Carette
*                   November, 2011
*
*

      Program show_H_ang
        IMPLICIT NONE


        OPEN(UNIT=8, FILE='yint.lst',STATUS='UNKNOWN',
     :     FORM='unformatted')
        OPEN(UNIT=11, FILE='ih.'//ih_file//'.lst',STATUS='UNKNOWN',
     :     FORM='unformatted')
        OPEN(UNIT=12, FILE='ico.lst',STATUS='UNKNOWN',
     :     FORM='unformatted')
        OPEN (UNIT=50,FILE='c.lst',STATUS='UNKNOWN',
     :     FORM='UNFORMATTED')




      end program show_H_and
