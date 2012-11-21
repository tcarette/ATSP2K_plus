!From Luca's hmk_standalone

module SafeIO

  !Parameters
  integer, parameter, private :: MAX_UNIT = 65535 

  integer, parameter :: Open_Error = 1 
  integer, parameter :: Read_Error = 2
  integer, parameter :: Read_End_Of_File = 3
  integer, parameter :: Read_End_Of_Record = 4
  integer, parameter :: Write_Error = 5
  integer, parameter :: Close_Error = 6
  integer, parameter :: Data_Consistency_Error = 7

  integer, parameter :: STDERR = 0
  integer, parameter :: STDOUT = 6
  integer, parameter :: STDINP = 5

contains 

  integer function Free_Unit()
    implicit none
    logical :: Unit_is_open
    Free_Unit=6
    Unit_is_open=.TRUE.
    do while(Unit_is_Open.and.Free_Unit<MAX_UNIT)
       Free_Unit = Free_Unit+1
       inquire(UNIT=Free_Unit,OPENED=Unit_is_Open)
    enddo
    if(Free_Unit==MAX_UNIT.and.Unit_is_Open)then
       Free_Unit=-1
       return
    endif
    return
  end function Free_Unit


  integer function Safe_Open(FI,AC,FO,ST,PO)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Description:
    ! Safe_Open returns the first available
    ! Logical Unit to be associated with file
    ! FI with the required Options
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !#########################################
    ! Debug:
    ! 11/6/2002
    ! Sembrerebbe funzionare ma non ho fatto
    ! test estesi
    !#########################################


    implicit none

    !Dummy Variables ======================
    ! -*- FI -*-
    ! File Name
    character(len=*), intent(in) :: FI
    ! -*- AC -*-
    ! ACTION
    ! "R" => ACTION = "READ"
    ! "W" => ACTION = "WRITE"
    character, intent(in), optional :: AC
    ! -*- FO -*-
    ! "F" => FORM = "FORMATTED"
    ! "U" => FORM = "UNFORMATTED"
    character, intent(in), optional :: FO
    ! -*- ST -*-
    ! "O" => STATUS = "OLD"
    ! "N" => STATUS = "NEW"
    ! "U" => STATUS = "UNKNOWN"
    ! "R" => STATUS = "REPLACE"
    character, intent(in), optional :: ST
    ! -*- PO -*-
    ! "A" => APPEND
    ! "I" => ASIS
    ! "R" => REWIND
    character, intent(in), optional :: PO 
    !======================================

    !Local Copies of optional arguments ===
    character :: LAC
    character :: LFO
    character :: LST
    character :: LPO
    !======================================

    !External Functions ===================
    integer :: ACCESS
    !======================================

    !Local Variables ======================
    logical :: File_Exists
    logical :: Unit_is_open
    logical :: File_is_already_open
    logical :: File_cannot_be_read
    logical :: File_cannot_be_writ
    integer :: IOS
    character(len=16):: ACSTRN
    character(len=16):: FOSTRN
    character(len=16):: STSTRN
    character(len=16):: POSTRN
    character(len=4) :: nstrn
    character(len=1024) :: dir
    integer :: i
    !======================================

    !Assegnamento delle copie locali
    !delle variabili opzionali
    LAC="B";if(present(AC))LAC=AC
    LFO="F";if(present(FO))LFO=FO
    LST="U";if(present(ST))LST=ST
    LPO="I";if(present(PO))LPO=PO

    !Estrae la directory dal nome del file e verifica che esista. 
    !dir=trim(adjustl(FI))
    !i=index(dir,"/",back=.TRUE.)
    !if(i>0)then
    !   dir=dir(1:i-1)
    !   if(.not.chkdir(dir))then
    !      Safe_Open=-1
    !      write(STDERR,"(a)") "(EEE) - Directory '"//trim(dir)//"' doesn't exist."
    !      return
    !   endif
    !endif

    !Verifica che il file esista, se richiesto da ST,
    !e che non sia gia' aperto 
    inquire(FILE   = trim(FI)            ,&
         EXIST  = File_Exists         ,&
         OPENED = File_is_already_open,&
         ACTION = ACSTRN              ,&
         NUMBER = Safe_Open            )
    write(nstrn,"(i4)") Safe_Open
    nstrn = adjustl(nstrn)
    if(.not.File_Exists.and.LST=="O")then
       write(STDERR,"(a)") "(EEE) - File '"//trim(FI)//"' doesn't exist."
       Safe_Open=-1
       return 
    endif
    if(File_Exists.and.LST=="N")then
       write(STDERR,"(a)") "(EEE) - File '"//trim(FI)//"' already exists" 
       Safe_Open=-1
       return 
    endif
    if(File_is_already_open)then
       write(STDERR,"(a)") "(WWW) - File '"//trim(FI)//"' is already open in "//trim(ACSTRN)//" mode."
       write(STDERR,"(a)") "        on unit "//trim(nstrn)
       write(STDERR,"(a)") "        Current logical unit will be returned"
       return
    endif

    !Verifica l'accessibilita` in lettura, in scrittura, o in entrambe
    !le modalita`, a seconda dell'azione AC specificata
    File_cannot_be_read = ACCESS(trim(FI),'r').ne.0
    File_cannot_be_writ = ACCESS(trim(FI),'w').ne.0

    select case(LAC)
    case("W")
       if(File_Exists.and.File_cannot_be_writ)then
          write(STDERR,"(a)") "(EEE) - File '"//trim(FI)//"' cannot be accessed in write mode."
          Safe_Open=-1
          return
       endif
       ACSTRN="WRITE"
    case("R")
       if(File_Exists.and.File_cannot_be_read)then
          write(STDERR,"(a)") "(EEE) - File '"//trim(FI)//"' cannot be accessed in read mode."
          Safe_Open=-1
          return
       endif
       ACSTRN="READ"
    case("B")
       if(File_Exists.and.(File_cannot_be_read.or.File_cannot_be_writ))then
          write(STDERR,"(a)") "(EEE) - File '"//trim(FI)//"' cannot be accessed in read/write mode."
          Safe_Open=-1
          return
       endif
       ACSTRN="READWRITE"
    case default
       write(STDERR,"(a)") "(EEE) - In subroutine Safe_Open :"
       write(STDERR,"(a)") "        Invalid Option '"//LAC//"' for AC"
       write(STDERR,"(a)") "        Available options are 'R', 'W', or 'B'"
       Safe_Open=-1
       return
    end select

    select case(LFO)
    case("F")
       FOSTRN="FORMATTED"
    case("U")
       FOSTRN="UNFORMATTED"
    case default
       write(STDERR,"(a)") "(EEE) - In subroutine Safe_Open :"
       write(STDERR,"(a)") "        Invalid Option '"//LFO//"' for FO"
       write(STDERR,"(a)") "        Available options are 'F' or 'U'"
       Safe_Open=-1
       return
    end select

    select case(LST)
    case("O")
       STSTRN="OLD"
    case("N")
       STSTRN="NEW"
    case("U")
       STSTRN="UNKNOWN"
    case("R")
       STSTRN="REPLACE"
    case default
       write(STDERR,"(a)") "(EEE) - In subroutine Safe_Open :"
       write(STDERR,"(a)") "        Invalid Option '"//LST//"' for ST"
       write(STDERR,"(a)") "        Available options are 'O', 'N', or 'U'"
       Safe_Open=-1
       return
    end select

    select case(LPO)
    case("A")
       POSTRN="APPEND"
    case("I")
       POSTRN="ASIS"
    case("R")
       POSTRN="REWIND"
    case default
       write(STDERR,"(a)") "(EEE) - In subroutine Safe_Open :"
       write(STDERR,"(a)") "        Invalid Option '"//LPO//"' for PO"
       write(STDERR,"(a)") "        Available options are 'A', 'I', or 'R'"
       Safe_Open=-1
       return
    end select

    !Seleziona un'unita` fra quelle disponibili tra 7 e 255
    Safe_Open   = 6
    Unit_is_open = .TRUE.
    do while( Unit_is_Open .and. Safe_Open < MAX_UNIT )
       Safe_Open = Safe_Open+1
       inquire(UNIT   = Safe_Open   ,&
            OPENED = Unit_is_Open  )
    enddo

    !Se non ha trovato unita` libere produce un messaggio d'errore
    !e stoppa tutto
    if(Safe_Open==MAX_UNIT.and.Unit_is_Open)then
       write(STDERR,"(a)") "(EEE) - Cannot find more free units"
       Safe_Open=-1
       return
    endif

    !Finalmente Apre il File
    open(UNIT   = Safe_Open   ,&
         FILE   = trim(FI)    ,&
         ACTION = trim(ACSTRN),&
         FORM   = trim(FOSTRN),&
         STATUS = trim(STSTRN),&
         POSITION=trim(POSTRN),&
         IOSTAT = IOS          )

    !Se ci sono errori di apertura stoppa tutto
    if(IOS.ne.0)then
       if(IOS>0)then
          write(STDERR,"(a)")   "(EEE) - An Error Occurred While Opening File '"//trim(FI)//"'"
          write(STDERR,"(a,i3)") "        IOSTAT = ",IOS
          !integer, allocatable :: vettore(:)
          !vettore=0
       endif
       if(IOS<0)then
          write(STDERR,"(a)")   "(EEE) - EOF or EOR Opening File '"//trim(FI)//"'"
          write(STDERR,"(a,i3)") "        IOSTAT = ",IOS
       endif
       Safe_Open=-1
       return
    endif


    return
  end function Safe_Open

  integer function Safe_Close(fp)
    implicit none
    integer, intent(in) :: fp
    character (len=4) :: fpstrn
    character (len=4) :: iostrn
    Close(fp,IOSTAT=Safe_Close)
    if(Safe_Close.ne.0)then
       write(fpstrn,"(i4)") fp
       write(iostrn,"(i4)") Safe_Close 
       write(STDERR,"(a)") "(EEE) - Close Error on unit "//fpstrn
       write(STDERR,"(a)") "        IOSTAT = "//iostrn
    endif
    return
  end function Safe_Close

  logical function chkdir(dir)
    character(len=*),intent(in) :: dir
    integer, parameter :: MAXLEN = 512
    character(len=MAXLEN):: strn
    integer :: n
    n=len_trim(adjustl(dir))
    if(n>MAXLEN)then
       write(*,*) "Directory name '"//trim(adjustl(dir))//"' is ",n," characters long,"
       write(*,*) "hence exceding maximal lenght, equal to ",MAXLEN
       return
    endif
    strn=adjustl(dir)
    if(strn(n:n)=="/")strn(n:n)=" "
    if(len_trim(strn)==0)return
    inquire(FILE=trim(strn),EXIST=chkdir)
    return
  end function chkdir

  subroutine mkdir(dir)
    character(len=*),intent(in) :: dir
    integer, parameter :: MAXLEN = 512
    character(len=MAXLEN):: strn
    integer :: n
    integer :: system 
    logical :: bool
    n=len_trim(adjustl(dir))
    if(n>MAXLEN)then
       write(*,*) "Directory name '"//trim(adjustl(dir))//"' is ",n," characters long,"
       write(*,*) "hence exceding maximal lenght, equal to ",MAXLEN
       return
    endif
    strn=adjustl(dir)
    if(strn(n:n)=="/")strn(n:n)=" "
    !se la directory e` vuota ossia e` quella corrente, non prova nemmeno a 
    !generarla.
    if(len_trim(strn)==0)return
    inquire(FILE=trim(strn),EXIST=bool)
    if(.not.bool)then
       n=system("mkdir -p '"//trim(adjustl(strn))//"'")
       if(n/=0)then
          write(*,*) "(EEE) - Cannot create dir '"//trim(strn)//"'"
          Stop
       endif
       write(0,*) "Directory "//trim(strn)//" has been created"
    endif
    return
  end subroutine mkdir

  integer function file_len(strn)
    implicit none
    character(len=*), intent(in) :: strn
    integer :: fp,stat
    fp=Safe_Open(strn,AC="R",FO="F",ST="O")
    file_len=-1;if(fp<0)return
    file_len=0
    do
       read(fp,*,iostat=stat)
       if(stat/=0)exit
       file_len=file_len+1
    enddo
    if(Safe_Close(fp)/=0)call crash
    return
  end function file_len

end module SafeIO
