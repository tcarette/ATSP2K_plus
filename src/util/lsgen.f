*------------------------------------------------------------------
*    L S G EN
*    By  L. Sturreson and C. Froese Fischre
*    
*    Computer Physics Communication 74 (1993) 432
*------------------------------------------------------------------
*     last edited April 26, 1993
      subroutine Adder(closed,slut,resS,resL,anel,par,dyn,expand)
      integer pop(1:15,0:10),resS,resL,skal,anel,par,i,j,stopp
      logical closed(1:15,0:10),slut,finns,dyn,expand
      character rad1*64,rad2*72,L2(0:20)
      data (L2(i),i=0,20)/'S','P','D','F','G','H','I','K','L','M','N'
     :                         ,'O','Q','R','T','U','V','W','X','Y','Z'/
      if (dyn) then
         inquire(file='clist.inp',exist=finns)
      else
         inquire(file='cfg.old',exist=finns)
      endif
      if (finns) then
         slut = .FALSE.
         if (expand) then
            if (dyn) then
               open(unit=8,file='clist.inp',status='old')
            else
               open(unit=8,file='cfg.old',status='old')
            endif
         else
            if (dyn) then
               open(unit=7,file='clist.inp',status='old')
            else
               open(unit=7,file='cfg.old',status='old')
            endif
         endif
         call Lockad(closed,slut,expand)
         if (.NOT.slut) then
            if (expand) then
               call Lasa1(8,rad1,pop,skal,slut,dyn) 
            else
               call Lasa1(7,rad1,pop,skal,slut,dyn)
            endif
         endif
         if (.NOT.slut) then
            anel = 0
            par  = 0
            do 10 i=1,15
               do 10 j=0,min(10,i-1)
                  if (closed(i,j)) then
                     anel = anel + 2 + 4*j
                  else
                     anel = anel + pop(i,j)
                     par  = mod(par+j*pop(i,j),2)
                  endif
   10       continue
            if (dyn) then
               stopp =  8*skal - 4
            else
               stopp = 16*skal - 8
            endif
            if (expand) then
               read(8,100,end=99) rad2(1:stopp)
            else
               read(7,100,end=99) rad2(1:stopp)
            endif
            resS = ichar(rad2(stopp-2:stopp-2)) - ichar('0')
            do 20 j=0,20
               if (rad2(stopp-1:stopp-1).EQ.L2(j)) then
                  resL = j
                  if (expand) then
                     rewind(8)
                  else
                     rewind(7)
                  endif
                  return
               endif
   20       continue
         endif
   99    slut = .TRUE.
         if (expand) then
            close(8)
         else
            close(7)
         endif
      else
         slut = .TRUE.
      endif
      return
  100 format(A)
      end
*     last edited April 26, 1993
      subroutine Blanda(org,varmax,virmax,lock,resS,resL,skal,nmax,dyn,
     :                  first,low,posn,posl,breit,J2min,J2max,minS,minL,
     :                  virtu,vir,lim,dubbel)
      integer org(1:15,0:10),antel(1:15,0:10),start(1:15,0:10),skal,cf
      integer ansats(1:15,0:10),varupp(1:15,0:10),varned(1:15,0:10)
      integer an10,an20,an21,an30,an31,an32,an40,an41,an42,an43,minS
      integer an50,an51,an52,an53,an54,an60,an61,an62,an63,an64,an65
      integer an70,an71,an72,an73,an74,an75,an76,stopp(1:15,0:10),minL
      integer an80,an81,an82,an83,an84,an85,an86,an87,low(1:15,0:10)
      integer an90,an91,an92,an93,an94,an95,an96,an97,an98,J2min,J2max
      integer anA0,anA1,anA2,anA3,anA4,anA5,anA6,anA7,anA8,anA9,nmax
      integer anB0,anB1,anB2,anB3,anB4,anB5,anB6,anB7,anB8,anB9,anBA
      integer anC0,anC1,anC2,anC3,anC4,anC5,anC6,anC7,anC8,anC9,anCA
      integer anD0,anD1,anD2,anD3,anD4,anD5,anD6,anD7,anD8,anD9,anDA
      integer anE0,anE1,anE2,anE3,anE4,anE5,anE6,anE7,anE8,anE9,anEA
      integer anF0,anF1,anF2,anF3,anF4,anF5,anF6,anF7,anF8,anF9,anFA
      integer varmax,par0,par,resS,resL,i,j,antal,posn(110),posl(110)
      integer virmax,lim(15),steg(1:15,0:10),dum 
      logical lock(1:15,0:10),dyn,first,breit,virtu(1:15,0:10),vir,finns
      logical dubbel(1:15,0:10),napp
      cf    = 0
      antal = 0
      par0  = 0
      finns = .FALSE.
      do 1 i=1,nmax
         do 1 j=0,min(10,i-1)
            if (dubbel(i,j)) then
               steg(i,j) = -2
            else
               steg(i,j) = -1
            endif
            antal = antal + org(i,j)
   1        par0  = mod(par0+j*org(i,j),2)
      if (nmax.LT.15) then
         do 2 i=nmax+1,15
            do 2 j=0,min(10,i-1)
   2           steg(i,j) = -1
      endif
*     1s
      call Slug(1,0,varmax,varupp,varned,ansats,org,lock(1,0),       
     :                      dubbel,low,start(1,0),stopp(1,0))
      do 10 an10 = start(1,0),stopp(1,0),steg(1,0)
         antel(1,0) = an10
         if (antel(1,0).GT.antal .OR. antel(1,0).LT.lim(1)) goto 10
         ansats(1,0) = an10
*     2s
      call Slug(2,0,varmax,varupp,varned,ansats,org,lock(2,0),       
     :                      dubbel,low,start(2,0),stopp(2,0))
      do 20 an20 = start(2,0),stopp(2,0),steg(2,0)
         antel(2,0) = an20 + antel(1,0)
         if (antel(2,0).GT.antal) goto 20
         ansats(2,0) = an20
*     2p
      call Slug(2,1,varmax,varupp,varned,ansats,org,lock(2,1),       
     :                      dubbel,low,start(2,1),stopp(2,1))
      do 21 an21 = start(2,1),stopp(2,1),steg(2,1)
         antel(2,1) = an21 + antel(2,0)
         if (antel(2,1).GT.antal .OR. antel(2,1).LT.lim(2)) goto 21
         ansats(2,1) = an21
*     3s
      call Slug(3,0,varmax,varupp,varned,ansats,org,lock(3,0),       
     :                      dubbel,low,start(3,0),stopp(3,0))
      do 30 an30 = start(3,0),stopp(3,0),steg(3,0)
         antel(3,0) = an30 + antel(2,1)
         if (antel(3,0).GT.antal) goto 30
         ansats(3,0) = an30
*     3p
      call Slug(3,1,varmax,varupp,varned,ansats,org,lock(3,1),       
     :                      dubbel,low,start(3,1),stopp(3,1))
      do 31 an31 = start(3,1),stopp(3,1),steg(3,1)
         antel(3,1) = an31 + antel(3,0)
         if (antel(3,1).GT.antal) goto 31
         ansats(3,1) = an31
*     3d
      call Slug(3,2,varmax,varupp,varned,ansats,org,lock(3,2),       
     :                      dubbel,low,start(3,2),stopp(3,2))
      do 32 an32 = start(3,2),stopp(3,2),steg(3,2)
         antel(3,2) = an32 + antel(3,1)
         if (antel(3,2).GT.antal .OR. antel(3,2).LT.lim(3)) goto 32
         ansats(3,2) = an32
*     4s
      call Slug(4,0,varmax,varupp,varned,ansats,org,lock(4,0),       
     :                      dubbel,low,start(4,0),stopp(4,0))
      do 40 an40 = start(4,0),stopp(4,0),steg(4,0)
         antel(4,0) = an40 + antel(3,2)
         if (antel(4,0).GT.antal) goto 40
         ansats(4,0) = an40
*     4p
      call Slug(4,1,varmax,varupp,varned,ansats,org,lock(4,1),       
     :                      dubbel,low,start(4,1),stopp(4,1))
      do 41 an41 = start(4,1),stopp(4,1),steg(4,1)
         antel(4,1) = an41 + antel(4,0)
         if (antel(4,1).GT.antal) goto 41
         ansats(4,1) = an41
*     4d
      call Slug(4,2,varmax,varupp,varned,ansats,org,lock(4,2),       
     :                      dubbel,low,start(4,2),stopp(4,2))
      do 42 an42 = start(4,2),stopp(4,2),steg(4,2)
         antel(4,2) = an42 + antel(4,1)
         if (antel(4,2).GT.antal) goto 42
         ansats(4,2) = an42
*     4f
      call Slug(4,3,varmax,varupp,varned,ansats,org,lock(4,3),       
     :                      dubbel,low,start(4,3),stopp(4,3))
      do 43 an43 = start(4,3),stopp(4,3),steg(4,3)
         antel(4,3) = an43 + antel(4,2)
         if (antel(4,3).GT.antal .OR. antel(4,3).LT.lim(4)) goto 43
         ansats(4,3) = an43
*     5s
      call Slug(5,0,varmax,varupp,varned,ansats,org,lock(5,0),       
     :                      dubbel,low,start(5,0),stopp(5,0))
      do 50 an50 = start(5,0),stopp(5,0),steg(5,0)
         antel(5,0) = an50 + antel(4,3)
         if (antel(5,0).GT.antal) goto 50
         ansats(5,0) = an50
*     5p
      call Slug(5,1,varmax,varupp,varned,ansats,org,lock(5,1),       
     :                      dubbel,low,start(5,1),stopp(5,1))
      do 51 an51 = start(5,1),stopp(5,1),steg(5,1)
         antel(5,1) = an51 + antel(5,0)
         if (antel(5,1).GT.antal) goto 51
         ansats(5,1) = an51
*     5d
      call Slug(5,2,varmax,varupp,varned,ansats,org,lock(5,2),       
     :                      dubbel,low,start(5,2),stopp(5,2))
      do 52 an52 = start(5,2),stopp(5,2),steg(5,2)
         antel(5,2) = an52 + antel(5,1)
         if (antel(5,2).GT.antal) goto 52
         ansats(5,2) = an52
*     5f
      call Slug(5,3,varmax,varupp,varned,ansats,org,lock(5,3),       
     :                      dubbel,low,start(5,3),stopp(5,3))
      do 53 an53 = start(5,3),stopp(5,3),steg(5,3)
         antel(5,3) = an53 + antel(5,2)
         if (antel(5,3).GT.antal) goto 53
         ansats(5,3) = an53
*     5g
      call Slug(5,4,varmax,varupp,varned,ansats,org,lock(5,4),       
     :                      dubbel,low,start(5,4),stopp(5,4))
      do 54 an54 = start(5,4),stopp(5,4),steg(5,4)
         antel(5,4) = an54 + antel(5,3)
         if (antel(5,4).GT.antal .OR. antel(5,4).LT.lim(5)) goto 54
         ansats(5,4) = an54
*     6s
      call Slug(6,0,varmax,varupp,varned,ansats,org,lock(6,0),       
     :                      dubbel,low,start(6,0),stopp(6,0))
      do 60 an60 = start(6,0),stopp(6,0),steg(6,0)
         antel(6,0) = an60 + antel(5,4)
         if (antel(6,0).GT.antal) goto 60
         ansats(6,0) = an60
*     6p
      call Slug(6,1,varmax,varupp,varned,ansats,org,lock(6,1),       
     :                      dubbel,low,start(6,1),stopp(6,1))
      do 61 an61 = start(6,1),stopp(6,1),steg(6,1)
         antel(6,1) = an61 + antel(6,0)
         if (antel(6,1).GT.antal) goto 61
         ansats(6,1) = an61
*     6d
      call Slug(6,2,varmax,varupp,varned,ansats,org,lock(6,2),       
     :                      dubbel,low,start(6,2),stopp(6,2))
      do 62 an62 = start(6,2),stopp(6,2),steg(6,2)
         antel(6,2) = an62 + antel(6,1)
         if (antel(6,2).GT.antal) goto 62
         ansats(6,2) = an62
*     6f
      call Slug(6,3,varmax,varupp,varned,ansats,org,lock(6,3),       
     :                      dubbel,low,start(6,3),stopp(6,3))
      do 63 an63 = start(6,3),stopp(6,3),steg(6,3)
         antel(6,3) = an63 + antel(6,2)
         if (antel(6,3).GT.antal) goto 63
         ansats(6,3) = an63
*     6g
      call Slug(6,4,varmax,varupp,varned,ansats,org,lock(6,4),       
     :                      dubbel,low,start(6,4),stopp(6,4))
      do 64 an64 = start(6,4),stopp(6,4),steg(6,4)
         antel(6,4) = an64 + antel(6,3)
         if (antel(6,4).GT.antal) goto 64
         ansats(6,4) = an64
*     6h
      call Slug(6,5,varmax,varupp,varned,ansats,org,lock(6,5),       
     :                      dubbel,low,start(6,5),stopp(6,5))
      do 65 an65 = start(6,5),stopp(6,5),steg(6,5)
         antel(6,5) = an65 + antel(6,4)
         if (antel(6,5).GT.antal .OR. antel(6,5).LT.lim(6)) goto 65
         ansats(6,5) = an65
*     7s
      call Slug(7,0,varmax,varupp,varned,ansats,org,lock(7,0),       
     :                      dubbel,low,start(7,0),stopp(7,0))
      do 70 an70 = start(7,0),stopp(7,0),steg(7,0)
         antel(7,0) = an70 + antel(6,5)
         if (antel(7,0).GT.antal) goto 70
         ansats(7,0) = an70
*     7p
      call Slug(7,1,varmax,varupp,varned,ansats,org,lock(7,1),       
     :                      dubbel,low,start(7,1),stopp(7,1))
      do 71 an71 = start(7,1),stopp(7,1),steg(7,1)
         antel(7,1) = an71 + antel(7,0)
         if (antel(7,1).GT.antal) goto 71
         ansats(7,1) = an71
*     7d
      call Slug(7,2,varmax,varupp,varned,ansats,org,lock(7,2),       
     :                      dubbel,low,start(7,2),stopp(7,2))
      do 72 an72 = start(7,2),stopp(7,2),steg(7,2)
         antel(7,2) = an72 + antel(7,1)
         if (antel(7,2).GT.antal) goto 72
         ansats(7,2) = an72
*     7f
      call Slug(7,3,varmax,varupp,varned,ansats,org,lock(7,3),       
     :                      dubbel,low,start(7,3),stopp(7,3))
      do 73 an73 = start(7,3),stopp(7,3),steg(7,3)
         antel(7,3) = an73 + antel(7,2)
         if (antel(7,3).GT.antal) goto 73
         ansats(7,3) = an73
*     7g
      call Slug(7,4,varmax,varupp,varned,ansats,org,lock(7,4),       
     :                      dubbel,low,start(7,4),stopp(7,4))
      do 74 an74 = start(7,4),stopp(7,4),steg(7,4)
         antel(7,4) = an74 + antel(7,3)
         if (antel(7,4).GT.antal) goto 74
         ansats(7,4) = an74
*     7h
      call Slug(7,5,varmax,varupp,varned,ansats,org,lock(7,5),       
     :                      dubbel,low,start(7,5),stopp(7,5))
      do 75 an75 = start(7,5),stopp(7,5),steg(7,5)
         antel(7,5) = an75 + antel(7,4)
         if (antel(7,5).GT.antal) goto 75
         ansats(7,5) = an75
*     7i
      call Slug(7,6,varmax,varupp,varned,ansats,org,lock(7,6),       
     :                      dubbel,low,start(7,6),stopp(7,6))
      do 76 an76 = start(7,6),stopp(7,6),steg(7,6)
         antel(7,6) = an76 + antel(7,5)
         if (antel(7,6).GT.antal .OR. antel(7,6).LT.lim(7)) goto 76
         ansats(7,6) = an76
*     8s
      call Slug(8,0,varmax,varupp,varned,ansats,org,lock(8,0),       
     :                      dubbel,low,start(8,0),stopp(8,0))
      do 80 an80 = start(8,0),stopp(8,0),steg(8,0)
         antel(8,0) = an80 + antel(7,6)
         if (antel(8,0).GT.antal) goto 80
         ansats(8,0) = an80
*     8p
      call Slug(8,1,varmax,varupp,varned,ansats,org,lock(8,1),       
     :                      dubbel,low,start(8,1),stopp(8,1))
      do 81 an81 = start(8,1),stopp(8,1),steg(8,1)
         antel(8,1) = an81 + antel(8,0)
         if (antel(8,1).GT.antal) goto 81
         ansats(8,1) = an81
*     8d
      call Slug(8,2,varmax,varupp,varned,ansats,org,lock(8,2),       
     :                      dubbel,low,start(8,2),stopp(8,2))
      do 82 an82 = start(8,2),stopp(8,2),steg(8,2)
         antel(8,2) = an82 + antel(8,1)
         if (antel(8,2).GT.antal) goto 82
         ansats(8,2) = an82
*     8f
      call Slug(8,3,varmax,varupp,varned,ansats,org,lock(8,3),       
     :                      dubbel,low,start(8,3),stopp(8,3))
      do 83 an83 = start(8,3),stopp(8,3),steg(8,3)
         antel(8,3) = an83 + antel(8,2)
         if (antel(8,3).GT.antal) goto 83
         ansats(8,3) = an83
*     8g
      call Slug(8,4,varmax,varupp,varned,ansats,org,lock(8,4),       
     :                      dubbel,low,start(8,4),stopp(8,4))
 
      do 84 an84 = start(8,4),stopp(8,4),steg(8,4)
         antel(8,4) = an84 + antel(8,3)
         if (antel(8,4).GT.antal) goto 84
         ansats(8,4) = an84
*     8h
      call Slug(8,5,varmax,varupp,varned,ansats,org,lock(8,5),       
     :                      dubbel,low,start(8,5),stopp(8,5))
      do 85 an85 = start(8,5),stopp(8,5),steg(8,5)
         antel(8,5) = an85 + antel(8,4)
         if (antel(8,5).GT.antal) goto 85
         ansats(8,5) = an85
*     8i
      call Slug(8,6,varmax,varupp,varned,ansats,org,lock(8,6),       
     :                      dubbel,low,start(8,6),stopp(8,6))
      do 86 an86 = start(8,6),stopp(8,6),steg(8,6)
         antel(8,6) = an86 + antel(8,5)
         if (antel(8,6).GT.antal) goto 86
         ansats(8,6) = an86
*     8k
      call Slug(8,7,varmax,varupp,varned,ansats,org,lock(8,7),       
     :                      dubbel,low,start(8,7),stopp(8,7))
      do 87 an87 = start(8,7),stopp(8,7),steg(8,7)
         antel(8,7) = an87 + antel(8,6)
         if (antel(8,7).GT.antal .OR. antel(8,7).LT.lim(8)) goto 87
         ansats(8,7) = an87
*     9s
      call Slug(9,0,varmax,varupp,varned,ansats,org,lock(9,0),       
     :                      dubbel,low,start(9,0),stopp(9,0))
      do 90 an90 = start(9,0),stopp(9,0),steg(9,0)
         antel(9,0) = an90 + antel(8,7)
         if (antel(9,0).GT.antal) goto 90
         ansats(9,0) = an90
*     9p
      call Slug(9,1,varmax,varupp,varned,ansats,org,lock(9,1),       
     :                      dubbel,low,start(9,1),stopp(9,1))
      do 91 an91 = start(9,1),stopp(9,1),steg(9,1)
         antel(9,1) = an91 + antel(9,0)
         if (antel(9,1).GT.antal) goto 91
         ansats(9,1) = an91
*     9d
      call Slug(9,2,varmax,varupp,varned,ansats,org,lock(9,2),       
     :                      dubbel,low,start(9,2),stopp(9,2))
      do 92 an92 = start(9,2),stopp(9,2),steg(9,2)
         antel(9,2) = an92 + antel(9,1)
         if (antel(9,2).GT.antal) goto 92
         ansats(9,2) = an92
*     9f
      call Slug(9,3,varmax,varupp,varned,ansats,org,lock(9,3),       
     :                      dubbel,low,start(9,3),stopp(9,3))
      do 93 an93 = start(9,3),stopp(9,3),steg(9,3)
         antel(9,3) = an93 + antel(9,2)
         if (antel(9,3).GT.antal) goto 93
         ansats(9,3) = an93
*     9g
      call Slug(9,4,varmax,varupp,varned,ansats,org,lock(9,4),       
     :                      dubbel,low,start(9,4),stopp(9,4))
      do 94 an94 = start(9,4),stopp(9,4),steg(9,4)
         antel(9,4) = an94 + antel(9,3)
         if (antel(9,4).GT.antal) goto 94
         ansats(9,4) = an94
*     9h
      call Slug(9,5,varmax,varupp,varned,ansats,org,lock(9,5),       
     :                      dubbel,low,start(9,5),stopp(9,5))
      do 95 an95 = start(9,5),stopp(9,5),steg(9,5)
         antel(9,5) = an95 + antel(9,4)
         if (antel(9,5).GT.antal) goto 95
         ansats(9,5) = an95
*     9i
      call Slug(9,6,varmax,varupp,varned,ansats,org,lock(9,6),       
     :                      dubbel,low,start(9,6),stopp(9,6))
      do 96 an96 = start(9,6),stopp(9,6),steg(9,6)
         antel(9,6) = an96 + antel(9,5)
         if (antel(9,6).GT.antal) goto 96
         ansats(9,6) = an96
*     9k
      call Slug(9,7,varmax,varupp,varned,ansats,org,lock(9,7),       
     :                      dubbel,low,start(9,7),stopp(9,7))
      do 97 an97 = start(9,7),stopp(9,7),steg(9,7)
         antel(9,7) = an97 + antel(9,6)
         if (antel(9,7).GT.antal) goto 97
         ansats(9,7) = an97
*     9l
      call Slug(9,8,varmax,varupp,varned,ansats,org,lock(9,8),       
     :                      dubbel,low,start(9,8),stopp(9,8))
      do 98 an98 = start(9,8),stopp(9,8),steg(9,8)
         antel(9,8) = an98 + antel(9,7)
         if (antel(9,8).GT.antal .OR. antel(9,8).LT.lim(9)) goto 98
         ansats(9,8) = an98
*     10s
      call Slug(10,0,varmax,varupp,varned,ansats,org,lock(10,0),     
     :                      dubbel,low,start(10,0),stopp(10,0))
      do 100 anA0 = start(10,0),stopp(10,0),steg(10,0)
         antel(10,0) = anA0 + antel(9,8)
         if (antel(10,0).GT.antal) goto 100
         ansats(10,0) = anA0
*     10p
      call Slug(10,1,varmax,varupp,varned,ansats,org,lock(10,1),     
     :                      dubbel,low,start(10,1),stopp(10,1))
      do 101 anA1 = start(10,1),stopp(10,1),steg(10,1)
         antel(10,1) = anA1 + antel(10,0)
         if (antel(10,1).GT.antal) goto 101
         ansats(10,1) = anA1
*     10d
      call Slug(10,2,varmax,varupp,varned,ansats,org,lock(10,2),     
     :                      dubbel,low,start(10,2),stopp(10,2))
      do 102 anA2 = start(10,2),stopp(10,2),steg(10,2)
         antel(10,2) = anA2 + antel(10,1)
         if (antel(10,2).GT.antal) goto 102
         ansats(10,2) = anA2
*     10f
      call Slug(10,3,varmax,varupp,varned,ansats,org,lock(10,3),     
     :                      dubbel,low,start(10,3),stopp(10,3))
      do 103 anA3 = start(10,3),stopp(10,3),steg(10,3)
         antel(10,3) = anA3 + antel(10,2)
         if (antel(10,3).GT.antal) goto 103
         ansats(10,3) = anA3
*     10g
      call Slug(10,4,varmax,varupp,varned,ansats,org,lock(10,4),     
     :                      dubbel,low,start(10,4),stopp(10,4))
      do 104 anA4 = start(10,4),stopp(10,4),steg(10,4)
         antel(10,4) = anA4 + antel(10,3)
         if (antel(10,4).GT.antal) goto 104
         ansats(10,4) = anA4
*     10h
      call Slug(10,5,varmax,varupp,varned,ansats,org,lock(10,5),     
     :                      dubbel,low,start(10,5),stopp(10,5))
      do 105 anA5 = start(10,5),stopp(10,5),steg(10,5)
         antel(10,5) = anA5 + antel(10,4)
         if (antel(10,5).GT.antal) goto 105
         ansats(10,5) = anA5
*     10i
      call Slug(10,6,varmax,varupp,varned,ansats,org,lock(10,6),     
     :                      dubbel,low,start(10,6),stopp(10,6))
      do 106 anA6 = start(10,6),stopp(10,6),steg(10,6)
         antel(10,6) = anA6 + antel(10,5)
         if (antel(10,6).GT.antal) goto 106
         ansats(10,6) = anA6
*     10k
      call Slug(10,7,varmax,varupp,varned,ansats,org,lock(10,7),     
     :                      dubbel,low,start(10,7),stopp(10,7))
      do 107 anA7 = start(10,7),stopp(10,7),steg(10,7)
         antel(10,7) = anA7 + antel(10,6)
         if (antel(10,7).GT.antal) goto 107
         ansats(10,7) = anA7
*     10l
      call Slug(10,8,varmax,varupp,varned,ansats,org,lock(10,8),     
     :                      dubbel,low,start(10,8),stopp(10,8))
      do 108 anA8 = start(10,8),stopp(10,8),steg(10,8)
         antel(10,8) = anA8 + antel(10,7)
         if (antel(10,8).GT.antal) goto 108
         ansats(10,8) = anA8
*     10m
      call Slug(10,9,varmax,varupp,varned,ansats,org,lock(10,9),     
     :                      dubbel,low,start(10,9),stopp(10,9))
      do 109 anA9 = start(10,9),stopp(10,9),steg(10,9)
         antel(10,9) = anA9 + antel(10,8)
         if (antel(10,9).GT.antal .OR. antel(10,9).LT.lim(10)) goto 109
         ansats(10,9) = anA9
*     11s
      call Slug(11,0,varmax,varupp,varned,ansats,org,lock(11,0),     
     :                      dubbel,low,start(11,0),stopp(11,0))
      do 110 anB0 = start(11,0),stopp(11,0),steg(11,0)
         antel(11,0) = anB0 + antel(10,9)
         if (antel(11,0).GT.antal) goto 110
         ansats(11,0) = anB0
*     11p
      call Slug(11,1,varmax,varupp,varned,ansats,org,lock(11,1),     
     :                      dubbel,low,start(11,1),stopp(11,1))
      do 111 anB1 = start(11,1),stopp(11,1),steg(11,1)
         antel(11,1) = anB1 + antel(11,0)
         if (antel(11,1).GT.antal) goto 111
         ansats(11,1) = anB1
*     11d
      call Slug(11,2,varmax,varupp,varned,ansats,org,lock(11,2),     
     :                      dubbel,low,start(11,2),stopp(11,2))
      do 112 anB2 = start(11,2),stopp(11,2),steg(11,2)
         antel(11,2) = anB2 + antel(11,1)
         if (antel(11,2).GT.antal) goto 112
         ansats(11,2) = anB2
*     11f
      call Slug(11,3,varmax,varupp,varned,ansats,org,lock(11,3),     
     :                      dubbel,low,start(11,3),stopp(11,3))
      do 113 anB3 = start(11,3),stopp(11,3),steg(11,3)
         antel(11,3) = anB3 + antel(11,2)
         if (antel(11,3).GT.antal) goto 113
         ansats(11,3) = anB3
*     11g
      call Slug(11,4,varmax,varupp,varned,ansats,org,lock(11,4),     
     :                      dubbel,low,start(11,4),stopp(11,4))
      do 114 anB4 = start(11,4),stopp(11,4),steg(11,4)
         antel(11,4) = anB4 + antel(11,3)
         if (antel(11,4).GT.antal) goto 114
         ansats(11,4) = anB4
*     11h
      call Slug(11,5,varmax,varupp,varned,ansats,org,lock(11,5),     
     :                      dubbel,low,start(11,5),stopp(11,5))
      do 115 anB5 = start(11,5),stopp(11,5),steg(11,5)
         antel(11,5) = anB5 + antel(11,4)
         if (antel(11,5).GT.antal) goto 115
         ansats(11,5) = anB5
*     11i
      call Slug(11,6,varmax,varupp,varned,ansats,org,lock(11,6),     
     :                      dubbel,low,start(11,6),stopp(11,6))
      do 116 anB6 = start(11,6),stopp(11,6),steg(11,6)
         antel(11,6) = anB6 + antel(11,5)
         if (antel(11,6).GT.antal) goto 116
         ansats(11,6) = anB6
*     11k
      call Slug(11,7,varmax,varupp,varned,ansats,org,lock(11,7),     
     :                      dubbel,low,start(11,7),stopp(11,7))
      do 117 anB7 = start(11,7),stopp(11,7),steg(11,7)
         antel(11,7) = anB7 + antel(11,6)
         if (antel(11,7).GT.antal) goto 117
         ansats(11,7) = anB7
*     11l
      call Slug(11,8,varmax,varupp,varned,ansats,org,lock(11,8),     
     :                      dubbel,low,start(11,8),stopp(11,8))
      do 118 anB8 = start(11,8),stopp(11,8),steg(11,8)
         antel(11,8) = anB8 + antel(11,7)
         if (antel(11,8).GT.antal) goto 118
         ansats(11,8) = anB8
*     11m
      call Slug(11,9,varmax,varupp,varned,ansats,org,lock(11,9),     
     :                      dubbel,low,start(11,9),stopp(11,9))
      do 119 anB9 = start(11,9),stopp(11,9),steg(11,9)
         antel(11,9) = anB9 + antel(11,8)
         if (antel(11,9).GT.antal) goto 119
         ansats(11,9) = anB9
*     11n
      call Slug(11,10,varmax,varupp,varned,ansats,org,lock(11,10),   
     :                      dubbel,low,start(11,10),stopp(11,10))
      do 1110 anBA = start(11,10),stopp(11,10),steg(11,10)
         antel(11,10) = anBA + antel(11,9)
         if (antel(11,10).GT.antal .OR. antel(11,10).LT.lim(11))
     :                                                         goto 1110
         ansats(11,10) = anBA
*     12s
      call Slug(12,0,varmax,varupp,varned,ansats,org,lock(12,0),     
     :                      dubbel,low,start(12,0),stopp(12,0))
      do 120 anC0 = start(12,0),stopp(12,0),steg(12,0)
         antel(12,0) = anC0 + antel(11,10)
         if (antel(12,0).GT.antal) goto 120
         ansats(12,0) = anC0
*     12p
      call Slug(12,1,varmax,varupp,varned,ansats,org,lock(12,1),     
     :                      dubbel,low,start(12,1),stopp(12,1))
      do 121 anC1 = start(12,1),stopp(12,1),steg(12,1)
         antel(12,1) = anC1 + antel(12,0)
         if (antel(12,1).GT.antal) goto 121
         ansats(12,1) = anC1
*     12d
      call Slug(12,2,varmax,varupp,varned,ansats,org,lock(12,2),     
     :                      dubbel,low,start(12,2),stopp(12,2))
      do 122 anC2 = start(12,2),stopp(12,2),steg(12,2)
         antel(12,2) = anC2 + antel(12,1)
         if (antel(12,2).GT.antal) goto 122
         ansats(12,2) = anC2
*     12f
      call Slug(12,3,varmax,varupp,varned,ansats,org,lock(12,3),     
     :                      dubbel,low,start(12,3),stopp(12,3))
      do 123 anC3 = start(12,3),stopp(12,3),steg(12,3)
         antel(12,3) = anC3 + antel(12,2)
         if (antel(12,3).GT.antal) goto 123
         ansats(12,3) = anC3
*     12g
      call Slug(12,4,varmax,varupp,varned,ansats,org,lock(12,4),     
     :                      dubbel,low,start(12,4),stopp(12,4))
      do 124 anC4 = start(12,4),stopp(12,4),steg(12,4)
         antel(12,4) = anC4 + antel(12,3)
         if (antel(12,4).GT.antal) goto 124
         ansats(12,4) = anC4
*     12h
      call Slug(12,5,varmax,varupp,varned,ansats,org,lock(12,5),     
     :                      dubbel,low,start(12,5),stopp(12,5))
      do 125 anC5 = start(12,5),stopp(12,5),steg(12,5)
         antel(12,5) = anC5 + antel(12,4)
         if (antel(12,5).GT.antal) goto 125
         ansats(12,5) = anC5
*     12i
      call Slug(12,6,varmax,varupp,varned,ansats,org,lock(12,6),     
     :                      dubbel,low,start(12,6),stopp(12,6))
      do 126 anC6 = start(12,6),stopp(12,6),steg(12,6)
         antel(12,6) = anC6 + antel(12,5)
         if (antel(12,6).GT.antal) goto 126
         ansats(12,6) = anC6
*     12k
      call Slug(12,7,varmax,varupp,varned,ansats,org,lock(12,7),     
     :                      dubbel,low,start(12,7),stopp(12,7))
      do 127 anC7 = start(12,7),stopp(12,7),steg(12,7)
         antel(12,7) = anC7 + antel(12,6)
         if (antel(12,7).GT.antal) goto 127
         ansats(12,7) = anC7
*     12l
      call Slug(12,8,varmax,varupp,varned,ansats,org,lock(12,8),     
     :                      dubbel,low,start(12,8),stopp(12,8))
      do 128 anC8 = start(12,8),stopp(12,8),steg(12,8)
         antel(12,8) = anC8 + antel(12,7)
         if (antel(12,8).GT.antal) goto 128
         ansats(12,8) = anC8
*     12m
      call Slug(12,9,varmax,varupp,varned,ansats,org,lock(12,9),     
     :                      dubbel,low,start(12,9),stopp(12,9))
      do 129 anC9 = start(12,9),stopp(12,9),steg(12,9)
         antel(12,9) = anC9 + antel(12,8)
         if (antel(12,9).GT.antal) goto 129
         ansats(12,9) = anC9
*     12n
      call Slug(12,10,varmax,varupp,varned,ansats,org,lock(12,10),   
     :                      dubbel,low,start(12,10),stopp(12,10))
      do 1210 anCA = start(12,10),stopp(12,10),steg(12,10)
         antel(12,10) = anCA + antel(12,9)
         if (antel(12,10).GT.antal .OR. antel(12,10).LT.lim(12))
     :                                                         goto 1210
         ansats(12,10) = anCA
*     13s
      call Slug(13,0,varmax,varupp,varned,ansats,org,lock(13,0),     
     :                      dubbel,low,start(13,0),stopp(13,0))
      do 130 anD0 = start(13,0),stopp(13,0),steg(13,0)
         antel(13,0) = anD0 + antel(12,10)
         if (antel(13,0).GT.antal) goto 130
         ansats(13,0) = anD0
*     13p
      call Slug(13,1,varmax,varupp,varned,ansats,org,lock(13,1),     
     :                      dubbel,low,start(13,1),stopp(13,1))
      do 131 anD1 = start(13,1),stopp(13,1),steg(13,1)
         antel(13,1) = anD1 + antel(13,0)
         if (antel(13,1).GT.antal) goto 131
         ansats(13,1) = anD1
*     13d
      call Slug(13,2,varmax,varupp,varned,ansats,org,lock(13,2),     
     :                      dubbel,low,start(13,2),stopp(13,2))
      do 132 anD2 = start(13,2),stopp(13,2),steg(13,2)
         antel(13,2) = anD2 + antel(13,1)
         if (antel(13,2).GT.antal) goto 132
         ansats(13,2) = anD2
*     13f
      call Slug(13,3,varmax,varupp,varned,ansats,org,lock(13,3),     
     :                      dubbel,low,start(13,3),stopp(13,3))
      do 133 anD3 = start(13,3),stopp(13,3),steg(13,3)
         antel(13,3) = anD3 + antel(13,2)
         if (antel(13,3).GT.antal) goto 133
         ansats(13,3) = anD3
*     13g
      call Slug(13,4,varmax,varupp,varned,ansats,org,lock(13,4),     
     :                      dubbel,low,start(13,4),stopp(13,4))
      do 134 anD4 = start(13,4),stopp(13,4),steg(13,4)
         antel(13,4) = anD4 + antel(13,3)
         if (antel(13,4).GT.antal) goto 134
         ansats(13,4) = anD4
*     13h
      call Slug(13,5,varmax,varupp,varned,ansats,org,lock(13,5),     
     :                      dubbel,low,start(13,5),stopp(13,5))
      do 135 anD5 = start(13,5),stopp(13,5),steg(13,5)
         antel(13,5) = anD5 + antel(13,4)
         if (antel(13,5).GT.antal) goto 135
         ansats(13,5) = anD5
*     13i
      call Slug(13,6,varmax,varupp,varned,ansats,org,lock(13,6),     
     :                      dubbel,low,start(13,6),stopp(13,6))
      do 136 anD6 = start(13,6),stopp(13,6),steg(13,6)
         antel(13,6) = anD6 + antel(13,5)
         if (antel(13,6).GT.antal) goto 136
         ansats(13,6) = anD6
*     13k
      call Slug(13,7,varmax,varupp,varned,ansats,org,lock(13,7),     
     :                      dubbel,low,start(13,7),stopp(13,7))
      do 137 anD7 = start(13,7),stopp(13,7),steg(13,7)
         antel(13,7) = anD7 + antel(13,6)
         if (antel(13,7).GT.antal) goto 137
         ansats(13,7) = anD7
*     13l
      call Slug(13,8,varmax,varupp,varned,ansats,org,lock(13,8),     
     :                      dubbel,low,start(13,8),stopp(13,8))
      do 138 anD8 = start(13,8),stopp(13,8),steg(13,8)
         antel(13,8) = anD8 + antel(13,7)
         if (antel(13,8).GT.antal) goto 138
         ansats(13,8) = anD8
*     13m
      call Slug(13,9,varmax,varupp,varned,ansats,org,lock(13,9),     
     :                      dubbel,low,start(13,9),stopp(13,9))
      do 139 anD9 = start(13,9),stopp(13,9),steg(13,9)
         antel(13,9) = anD9 + antel(13,8)
         if (antel(13,9).GT.antal) goto 139
         ansats(13,9) = anD9
*     13n
      call Slug(13,10,varmax,varupp,varned,ansats,org,lock(13,10),   
     :                      dubbel,low,start(13,10),stopp(13,10))
      do 1310 anDA = start(13,10),stopp(13,10),steg(13,10)
         antel(13,10) = anDA + antel(13,9)
         if (antel(13,10).GT.antal .OR. antel(13,10).LT.lim(13))
     :                                                         goto 1310
         ansats(13,10) = anDA
*     14s
      call Slug(14,0,varmax,varupp,varned,ansats,org,lock(14,0),     
     :                      dubbel,low,start(14,0),stopp(14,0))
      do 140 anE0 = start(14,0),stopp(14,0),steg(14,0)
         antel(14,0) = anE0 + antel(13,10)
         if (antel(14,0).GT.antal) goto 140
         ansats(14,0) = anE0
*     14p
      call Slug(14,1,varmax,varupp,varned,ansats,org,lock(14,1),     
     :                      dubbel,low,start(14,1),stopp(14,1))
      do 141 anE1 = start(14,1),stopp(14,1),steg(14,1)
         antel(14,1) = anE1 + antel(14,0)
         if (antel(14,1).GT.antal) goto 141
         ansats(14,1) = anE1
*     14d
      call Slug(14,2,varmax,varupp,varned,ansats,org,lock(14,2),     
     :                      dubbel,low,start(14,2),stopp(14,2))
      do 142 anE2 = start(14,2),stopp(14,2),steg(14,2)
         antel(14,2) = anE2 + antel(14,1)
         if (antel(14,2).GT.antal) goto 142
         ansats(14,2) = anE2
*     14f
      call Slug(14,3,varmax,varupp,varned,ansats,org,lock(14,3),     
     :                      dubbel,low,start(14,3),stopp(14,3))
      do 143 anE3 = start(14,3),stopp(14,3),steg(14,3)
         antel(14,3) = anE3 + antel(14,2)
         if (antel(14,3).GT.antal) goto 143
         ansats(14,3) = anE3
*     14g
      call Slug(14,4,varmax,varupp,varned,ansats,org,lock(14,4),     
     :                      dubbel,low,start(14,4),stopp(14,4))
      do 144 anE4 = start(14,4),stopp(14,4),steg(14,4)
         antel(14,4) = anE4 + antel(14,3)
         if (antel(14,4).GT.antal) goto 144
         ansats(14,4) = anE4
*     14h
      call Slug(14,5,varmax,varupp,varned,ansats,org,lock(14,5),     
     :                      dubbel,low,start(14,5),stopp(14,5))
      do 145 anE5 = start(14,5),stopp(14,5),steg(14,5)
         antel(14,5) = anE5 + antel(14,4)
         if (antel(14,5).GT.antal) goto 145
         ansats(14,5) = anE5
*     14i
      call Slug(14,6,varmax,varupp,varned,ansats,org,lock(14,6),     
     :                      dubbel,low,start(14,6),stopp(14,6))
      do 146 anE6 = start(14,6),stopp(14,6),steg(14,6)
         antel(14,6) = anE6 + antel(14,5)
         if (antel(14,6).GT.antal) goto 146
         ansats(14,6) = anE6
*     14k
      call Slug(14,7,varmax,varupp,varned,ansats,org,lock(14,7),     
     :                      dubbel,low,start(14,7),stopp(14,7))
      do 147 anE7 = start(14,7),stopp(14,7),steg(14,7)
         antel(14,7) = anE7 + antel(14,6)
         if (antel(14,7).GT.antal) goto 147
         ansats(14,7) = anE7
*     14l
      call Slug(14,8,varmax,varupp,varned,ansats,org,lock(14,8),     
     :                      dubbel,low,start(14,8),stopp(14,8))
      do 148 anE8 = start(14,8),stopp(14,8),steg(14,8)
         antel(14,8) = anE8 + antel(14,7)
         if (antel(14,8).GT.antal) goto 148
         ansats(14,8) = anE8
*     14m
      call Slug(14,9,varmax,varupp,varned,ansats,org,lock(14,9),     
     :                      dubbel,low,start(14,9),stopp(14,9))
      do 149 anE9 = start(14,9),stopp(14,9),steg(14,9)
         antel(14,9) = anE9 + antel(14,8)
         if (antel(14,9).GT.antal) goto 149
         ansats(14,9) = anE9
*     14n
      call Slug(14,10,varmax,varupp,varned,ansats,org,lock(14,10),   
     :                      dubbel,low,start(14,10),stopp(14,10))
      do 1410 anEA = start(14,10),stopp(14,10),steg(14,10)
         antel(14,10) = anEA + antel(14,9)
         if (antel(14,10).GT.antal .OR. antel(14,10).LT.lim(14))
     :                                                         goto 1410
         ansats(14,10) = anEA
*     15s
      call Slug(15,0,varmax,varupp,varned,ansats,org,lock(15,0),     
     :                      dubbel,low,start(15,0),stopp(15,0))
      do 150 anF0 = start(15,0),stopp(15,0),steg(15,0)
         antel(15,0) = anF0 + antel(14,10)
         if (antel(15,0).GT.antal) goto 150
         ansats(15,0) = anF0
*     15p
      call Slug(15,1,varmax,varupp,varned,ansats,org,lock(15,1),     
     :                      dubbel,low,start(15,1),stopp(15,1))
      do 151 anF1 = start(15,1),stopp(15,1),steg(15,1)
         antel(15,1) = anF1 + antel(15,0)
         if (antel(15,1).GT.antal) goto 151
         ansats(15,1) = anF1
*     15d
      call Slug(15,2,varmax,varupp,varned,ansats,org,lock(15,2),     
     :                      dubbel,low,start(15,2),stopp(15,2))
      do 152 anF2 = start(15,2),stopp(15,2),steg(15,2)
         antel(15,2) = anF2 + antel(15,1)
         if (antel(15,2).GT.antal) goto 152
         ansats(15,2) = anF2
*     15f
      call Slug(15,3,varmax,varupp,varned,ansats,org,lock(15,3),     
     :                      dubbel,low,start(15,3),stopp(15,3))
      do 153 anF3 = start(15,3),stopp(15,3),steg(15,3)
         antel(15,3) = anF3 + antel(15,2)
         if (antel(15,3).GT.antal) goto 153
         ansats(15,3) = anF3
*     15g
      call Slug(15,4,varmax,varupp,varned,ansats,org,lock(15,4),     
     :                      dubbel,low,start(15,4),stopp(15,4))
      do 154 anF4 = start(15,4),stopp(15,4),steg(15,4)
         antel(15,4) = anF4 + antel(15,3)
         if (antel(15,4).GT.antal) goto 154
         ansats(15,4) = anF4
*     15h
      call Slug(15,5,varmax,varupp,varned,ansats,org,lock(15,5),     
     :                      dubbel,low,start(15,5),stopp(15,5))
      do 155 anF5 = start(15,5),stopp(15,5),steg(15,5)
         antel(15,5) = anF5 + antel(15,4)
         if (antel(15,5).GT.antal) goto 155
         ansats(15,5) = anF5
*     15i
      call Slug(15,6,varmax,varupp,varned,ansats,org,lock(15,6),     
     :                      dubbel,low,start(15,6),stopp(15,6))
      do 156 anF6 = start(15,6),stopp(15,6),steg(15,6)
         antel(15,6) = anF6 + antel(15,5)
         if (antel(15,6).GT.antal) goto 156
         ansats(15,6) = anF6
*     15k
      call Slug(15,7,varmax,varupp,varned,ansats,org,lock(15,7),     
     :                      dubbel,low,start(15,7),stopp(15,7))
      do 157 anF7 = start(15,7),stopp(15,7),steg(15,7)
         antel(15,7) = anF7 + antel(15,6)
         if (antel(15,7).GT.antal) goto 157
         ansats(15,7) = anF7
*     15l
      call Slug(15,8,varmax,varupp,varned,ansats,org,lock(15,8),     
     :                      dubbel,low,start(15,8),stopp(15,8))
      do 158 anF8 = start(15,8),stopp(15,8),steg(15,8)
         antel(15,8) = anF8 + antel(15,7)
         if (antel(15,8).GT.antal) goto 158
         ansats(15,8) = anF8
*     15m
      call Slug(15,9,varmax,varupp,varned,ansats,org,lock(15,9),     
     :                      dubbel,low,start(15,9),stopp(15,9))
      do 159 anF9 = start(15,9),stopp(15,9),steg(15,9)
         antel(15,9) = anF9 + antel(15,8)
         if (antel(15,9).GT.antal) goto 159
         ansats(15,9) = anF9
*     15n
      call Slug(15,10,varmax,varupp,varned,ansats,org,lock(15,10),   
     :                      dubbel,low,start(15,10),stopp(15,10))
      do 1510 anFA = start(15,10),stopp(15,10),steg(15,10)
         antel(15,10) = anFA + antel(15,9)
         if (antel(15,10).NE.antal) goto 1510
         ansats(15,10) = anFA
         par = 0
         do 6 i=1,15
            do 6 j=0,min(10,i-1)
    6          par = mod(par+j*ansats(i,j),2)
         if (par.EQ.par0) then
            if (.NOT. vir) then
               napp = .FALSE.
               call Gen(ansats,posn,posl,resS,resL,skal,
     :                    cf,dyn,first,breit,J2min,J2max,minS,minL,napp)
            else
               napp = .TRUE.
               call Gen(ansats,posn,posl,resS,resL,skal,
     :                    cf,dyn,first,breit,J2min,J2max,minS,minL,napp)
               if (napp) then
                  if (finns) then
                     open(unit=21,status='scratch')
                     call Blandb(ansats,nmax,virmax,lock,21,low,virtu)
                     rewind(21)
                     call Mergeb(nmax,dum)
                  else
                     open(unit=20,status='scratch')
                     call Blandb(ansats,nmax,virmax,lock,20,low,virtu)
                     rewind(20)
                     finns = .TRUE.
                  endif
               endif
            endif
         endif
 1510 continue
  159 continue
  158 continue
  157 continue
  156 continue
  155 continue
  154 continue
  153 continue
  152 continue
  151 continue
  150 continue
 1410 continue
  149 continue
  148 continue
  147 continue
  146 continue
  145 continue
  144 continue
  143 continue
  142 continue
  141 continue
  140 continue
 1310 continue
  139 continue
  138 continue
  137 continue
  136 continue
  135 continue
  134 continue
  133 continue
  132 continue
  131 continue
  130 continue
 1210 continue
  129 continue
  128 continue
  127 continue
  126 continue
  125 continue
  124 continue
  123 continue
  122 continue
  121 continue
  120 continue
 1110 continue
  119 continue
  118 continue
  117 continue
  116 continue
  115 continue
  114 continue
  113 continue
  112 continue
  111 continue
  110 continue
  109 continue
  108 continue
  107 continue
  106 continue
  105 continue
  104 continue
  103 continue
  102 continue
  101 continue
  100 continue
   98 continue
   97 continue
   96 continue
   95 continue
   94 continue
   93 continue
   92 continue
   91 continue
   90 continue
   87 continue
   86 continue
   85 continue
   84 continue
   83 continue
   82 continue
   81 continue
   80 continue
   76 continue
   75 continue
   74 continue
   73 continue
   72 continue
   71 continue
   70 continue
   65 continue
   64 continue
   63 continue
   62 continue
   61 continue
   60 continue
   54 continue
   53 continue
   52 continue
   51 continue
   50 continue
   43 continue
   42 continue
   41 continue
   40 continue
   32 continue
   31 continue
   30 continue
   21 continue
   20 continue
   10 continue
      if (vir) then
         if (nmax.LT.15) then
            do 391 i=nmax+1,15
               do 391 j=0,min(10,i-1)
  391             ansats(i,j) = 0
         endif
         cf   = 0
         napp = .FALSE.
  490    do 491 i=1,nmax
  491       read(20,5000,end=492) (ansats(i,j),j=0,min(10,i-1))
         call Gen(ansats,posn,posl,resS,resL,
     :             skal,cf,dyn,first,breit,J2min,J2max,minS,minL,napp)
         goto 490
      endif
  492 continue
      if (first) then
         write(7,1000) '*'
         rewind(7)
      else
         write(8,1000) '*'
         rewind(8)
      endif
      if (cf.EQ.0) then
         write(*,1005) 'No configuration state has been generated.'
      elseif (cf.EQ.1) then
         write(*,1005) 'One configuration state has been generated.'
      elseif (cf.LT.10) then
         write(*,1001) cf,' configuration states have been generated.'
      elseif (cf.LT.100) then
         write(*,1002) cf,' configuration states have been generated.'
      elseif (cf.LT.1000) then
         write(*,1003) cf,' configuration states have been generated.'
      elseif (cf.LT.10000) then
         write(*,1004) cf,' configuration states have been generated.'
      elseif (cf.LT.100000) then
         write(*,1006) cf,' configuration states have been generated.'
      else
         write(*,*) cf,' configuration states have been generated.'
      endif
 1000 format(A)
 1001 format(' ',I1,A)
 1002 format(' ',I2,A)
 1003 format(' ',I3,A)
 1004 format(' ',I4,A)
 1005 format(' ',A)
 1006 format(' ',I5,A)
 5000 format(11I2)
      return
      end
*     last edited September 26, 1993
      subroutine Blandb(org,nmax,varmax,lock,fil,low,virtu)
      integer org(1:15,0:10),antel(1:15,0:10),start(1:15,0:10)
      integer ansats(1:15,0:10),varupp(1:15,0:10),varned(1:15,0:10)
      integer an10,an20,an21,an30,an31,an32,an40,an41,an42,an43
      integer an50,an51,an52,an53,an54,an60,an61,an62,an63,an64,an65
      integer an70,an71,an72,an73,an74,an75,an76,stopp(1:15,0:10)
      integer an80,an81,an82,an83,an84,an85,an86,an87,low(1:15,0:10)
      integer an90,an91,an92,an93,an94,an95,an96,an97,an98
      integer anA0,anA1,anA2,anA3,anA4,anA5,anA6,anA7,anA8,anA9,fil
      integer anB0,anB1,anB2,anB3,anB4,anB5,anB6,anB7,anB8,anB9,anBA
      integer anC0,anC1,anC2,anC3,anC4,anC5,anC6,anC7,anC8,anC9,anCA
      integer anD0,anD1,anD2,anD3,anD4,anD5,anD6,anD7,anD8,anD9,anDA
      integer anE0,anE1,anE2,anE3,anE4,anE5,anE6,anE7,anE8,anE9,anEA
      integer anF0,anF1,anF2,anF3,anF4,anF5,anF6,anF7,anF8,anF9,anFA
      integer varmax,par0,par,i,j,antal,cf
      logical lock(1:15,0:10),virtu(1:15,0:10)
C     write(*,*) 'nmax,varmax,fil =',nmax,varmax,fil
      cf    = 0
      antal = 0
      par0  = 0
C     do 1 i=1,nmax
      do 1 i=1,15
         do 1 j=0,min(10,i-1)
            antal = antal + org(i,j)
   1        par0  = mod(par0+j*org(i,j),2)
C     write(*,*) 'antal,par0 =',antal,par0
C     write(*,*) 'org'
C     do 2 i=1,15
C  2     write(*,*) (org(i,j),j=0,min(10,i-1))
C     write(*,*) 'lock'
C     do 3 i=1,15
C  3     write(*,*) (lock(i,j),j=0,min(10,i-1))
C     write(*,*) 'low'
C     do 4 i=1,15
C  4     write(*,*) (low(i,j),j=0,min(10,i-1))
C     write(*,*) 'low'
C     do 5 i=1,15
C  5     write(*,*) (virtu(i,j),j=0,min(10,i-1))
*     1s
      call Sluggo(1,0,varmax,varupp,varned,ansats,org,lock(1,0),       
     :                             low,start(1,0),stopp(1,0),virtu)
      do 10 an10 = start(1,0),stopp(1,0),-1
         antel(1,0)  = an10
         ansats(1,0) = an10
*     2s
      call Sluggo(2,0,varmax,varupp,varned,ansats,org,lock(2,0),       
     :                             low,start(2,0),stopp(2,0),virtu)
      do 20 an20 = start(2,0),stopp(2,0),-1
         antel(2,0) = an20 + antel(1,0)
         if (antel(2,0).GT.antal) goto 20
         ansats(2,0) = an20
*     2p
      call Sluggo(2,1,varmax,varupp,varned,ansats,org,lock(2,1),       
     :                             low,start(2,1),stopp(2,1),virtu)
      do 21 an21 = start(2,1),stopp(2,1),-1
         antel(2,1) = an21 + antel(2,0)
         if (antel(2,1).GT.antal) goto 21
         ansats(2,1) = an21
*     3s
      call Sluggo(3,0,varmax,varupp,varned,ansats,org,lock(3,0),       
     :                             low,start(3,0),stopp(3,0),virtu)
      do 30 an30 = start(3,0),stopp(3,0),-1
         antel(3,0) = an30 + antel(2,1)
         if (antel(3,0).GT.antal) goto 30
         ansats(3,0) = an30
*     3p
      call Sluggo(3,1,varmax,varupp,varned,ansats,org,lock(3,1),       
     :                             low,start(3,1),stopp(3,1),virtu)
      do 31 an31 = start(3,1),stopp(3,1),-1
         antel(3,1) = an31 + antel(3,0)
         if (antel(3,1).GT.antal) goto 31
         ansats(3,1) = an31
*     3d
      call Sluggo(3,2,varmax,varupp,varned,ansats,org,lock(3,2),       
     :                             low,start(3,2),stopp(3,2),virtu)
      do 32 an32 = start(3,2),stopp(3,2),-1
         antel(3,2) = an32 + antel(3,1)
         if (antel(3,2).GT.antal) goto 32
         ansats(3,2) = an32
*     4s
      call Sluggo(4,0,varmax,varupp,varned,ansats,org,lock(4,0),       
     :                             low,start(4,0),stopp(4,0),virtu)
      do 40 an40 = start(4,0),stopp(4,0),-1
         antel(4,0) = an40 + antel(3,2)
         if (antel(4,0).GT.antal) goto 40
         ansats(4,0) = an40
*     4p
      call Sluggo(4,1,varmax,varupp,varned,ansats,org,lock(4,1),       
     :                             low,start(4,1),stopp(4,1),virtu)
      do 41 an41 = start(4,1),stopp(4,1),-1
         antel(4,1) = an41 + antel(4,0)
         if (antel(4,1).GT.antal) goto 41
         ansats(4,1) = an41
*     4d
      call Sluggo(4,2,varmax,varupp,varned,ansats,org,lock(4,2),       
     :                             low,start(4,2),stopp(4,2),virtu)
      do 42 an42 = start(4,2),stopp(4,2),-1
         antel(4,2) = an42 + antel(4,1)
         if (antel(4,2).GT.antal) goto 42
         ansats(4,2) = an42
*     4f
      call Sluggo(4,3,varmax,varupp,varned,ansats,org,lock(4,3),       
     :                             low,start(4,3),stopp(4,3),virtu)
      do 43 an43 = start(4,3),stopp(4,3),-1
         antel(4,3) = an43 + antel(4,2)
         if (antel(4,3).GT.antal) goto 43
         ansats(4,3) = an43
*     5s
      call Sluggo(5,0,varmax,varupp,varned,ansats,org,lock(5,0),       
     :                             low,start(5,0),stopp(5,0),virtu)
      do 50 an50 = start(5,0),stopp(5,0),-1
         antel(5,0) = an50 + antel(4,3)
         if (antel(5,0).GT.antal) goto 50
         ansats(5,0) = an50
*     5p
      call Sluggo(5,1,varmax,varupp,varned,ansats,org,lock(5,1),       
     :                             low,start(5,1),stopp(5,1),virtu)
      do 51 an51 = start(5,1),stopp(5,1),-1
         antel(5,1) = an51 + antel(5,0)
         if (antel(5,1).GT.antal) goto 51
         ansats(5,1) = an51
*     5d
      call Sluggo(5,2,varmax,varupp,varned,ansats,org,lock(5,2),       
     :                             low,start(5,2),stopp(5,2),virtu)
      do 52 an52 = start(5,2),stopp(5,2),-1
         antel(5,2) = an52 + antel(5,1)
         if (antel(5,2).GT.antal) goto 52
         ansats(5,2) = an52
*     5f
      call Sluggo(5,3,varmax,varupp,varned,ansats,org,lock(5,3),       
     :                             low,start(5,3),stopp(5,3),virtu)
      do 53 an53 = start(5,3),stopp(5,3),-1
         antel(5,3) = an53 + antel(5,2)
         if (antel(5,3).GT.antal) goto 53
         ansats(5,3) = an53
*     5g
      call Sluggo(5,4,varmax,varupp,varned,ansats,org,lock(5,4),       
     :                             low,start(5,4),stopp(5,4),virtu)
      do 54 an54 = start(5,4),stopp(5,4),-1
         antel(5,4) = an54 + antel(5,3)
         if (antel(5,4).GT.antal) goto 54
         ansats(5,4) = an54
*     6s
      call Sluggo(6,0,varmax,varupp,varned,ansats,org,lock(6,0),       
     :                             low,start(6,0),stopp(6,0),virtu)
      do 60 an60 = start(6,0),stopp(6,0),-1
         antel(6,0) = an60 + antel(5,4)
         if (antel(6,0).GT.antal) goto 60
         ansats(6,0) = an60
*     6p
      call Sluggo(6,1,varmax,varupp,varned,ansats,org,lock(6,1),       
     :                             low,start(6,1),stopp(6,1),virtu)
      do 61 an61 = start(6,1),stopp(6,1),-1
         antel(6,1) = an61 + antel(6,0)
         if (antel(6,1).GT.antal) goto 61
         ansats(6,1) = an61
*     6d
      call Sluggo(6,2,varmax,varupp,varned,ansats,org,lock(6,2),       
     :                             low,start(6,2),stopp(6,2),virtu)
      do 62 an62 = start(6,2),stopp(6,2),-1
         antel(6,2) = an62 + antel(6,1)
         if (antel(6,2).GT.antal) goto 62
         ansats(6,2) = an62
*     6f
      call Sluggo(6,3,varmax,varupp,varned,ansats,org,lock(6,3),       
     :                             low,start(6,3),stopp(6,3),virtu)
      do 63 an63 = start(6,3),stopp(6,3),-1
         antel(6,3) = an63 + antel(6,2)
         if (antel(6,3).GT.antal) goto 63
         ansats(6,3) = an63
*     6g
      call Sluggo(6,4,varmax,varupp,varned,ansats,org,lock(6,4),       
     :                             low,start(6,4),stopp(6,4),virtu)
      do 64 an64 = start(6,4),stopp(6,4),-1
         antel(6,4) = an64 + antel(6,3)
         if (antel(6,4).GT.antal) goto 64
         ansats(6,4) = an64
*     6h
      call Sluggo(6,5,varmax,varupp,varned,ansats,org,lock(6,5),       
     :                             low,start(6,5),stopp(6,5),virtu)
      do 65 an65 = start(6,5),stopp(6,5),-1
         antel(6,5) = an65 + antel(6,4)
         if (antel(6,5).GT.antal) goto 65
         ansats(6,5) = an65
*     7s
      call Sluggo(7,0,varmax,varupp,varned,ansats,org,lock(7,0),       
     :                             low,start(7,0),stopp(7,0),virtu)
      do 70 an70 = start(7,0),stopp(7,0),-1
         antel(7,0) = an70 + antel(6,5)
         if (antel(7,0).GT.antal) goto 70
         ansats(7,0) = an70
*     7p
      call Sluggo(7,1,varmax,varupp,varned,ansats,org,lock(7,1),       
     :                             low,start(7,1),stopp(7,1),virtu)
      do 71 an71 = start(7,1),stopp(7,1),-1
         antel(7,1) = an71 + antel(7,0)
         if (antel(7,1).GT.antal) goto 71
         ansats(7,1) = an71
*     7d
      call Sluggo(7,2,varmax,varupp,varned,ansats,org,lock(7,2),       
     :                             low,start(7,2),stopp(7,2),virtu)
      do 72 an72 = start(7,2),stopp(7,2),-1
         antel(7,2) = an72 + antel(7,1)
         if (antel(7,2).GT.antal) goto 72
         ansats(7,2) = an72
*     7f
      call Sluggo(7,3,varmax,varupp,varned,ansats,org,lock(7,3),       
     :                             low,start(7,3),stopp(7,3),virtu)
      do 73 an73 = start(7,3),stopp(7,3),-1
         antel(7,3) = an73 + antel(7,2)
         if (antel(7,3).GT.antal) goto 73
         ansats(7,3) = an73
*     7g
      call Sluggo(7,4,varmax,varupp,varned,ansats,org,lock(7,4),       
     :                             low,start(7,4),stopp(7,4),virtu)
      do 74 an74 = start(7,4),stopp(7,4),-1
         antel(7,4) = an74 + antel(7,3)
         if (antel(7,4).GT.antal) goto 74
         ansats(7,4) = an74
*     7h
      call Sluggo(7,5,varmax,varupp,varned,ansats,org,lock(7,5),       
     :                             low,start(7,5),stopp(7,5),virtu)
      do 75 an75 = start(7,5),stopp(7,5),-1
         antel(7,5) = an75 + antel(7,4)
         if (antel(7,5).GT.antal) goto 75
         ansats(7,5) = an75
*     7i
      call Sluggo(7,6,varmax,varupp,varned,ansats,org,lock(7,6),       
     :                             low,start(7,6),stopp(7,6),virtu)
      do 76 an76 = start(7,6),stopp(7,6),-1
         antel(7,6) = an76 + antel(7,5)
         if (antel(7,6).GT.antal) goto 76
         ansats(7,6) = an76
*     8s
      call Sluggo(8,0,varmax,varupp,varned,ansats,org,lock(8,0),       
     :                             low,start(8,0),stopp(8,0),virtu)
      do 80 an80 = start(8,0),stopp(8,0),-1
         antel(8,0) = an80 + antel(7,6)
         if (antel(8,0).GT.antal) goto 80
         ansats(8,0) = an80
*     8p
      call Sluggo(8,1,varmax,varupp,varned,ansats,org,lock(8,1),       
     :                             low,start(8,1),stopp(8,1),virtu)
      do 81 an81 = start(8,1),stopp(8,1),-1
         antel(8,1) = an81 + antel(8,0)
         if (antel(8,1).GT.antal) goto 81
         ansats(8,1) = an81
*     8d
      call Sluggo(8,2,varmax,varupp,varned,ansats,org,lock(8,2),       
     :                             low,start(8,2),stopp(8,2),virtu)
      do 82 an82 = start(8,2),stopp(8,2),-1
         antel(8,2) = an82 + antel(8,1)
         if (antel(8,2).GT.antal) goto 82
         ansats(8,2) = an82
*     8f
      call Sluggo(8,3,varmax,varupp,varned,ansats,org,lock(8,3),       
     :                             low,start(8,3),stopp(8,3),virtu)
      do 83 an83 = start(8,3),stopp(8,3),-1
         antel(8,3) = an83 + antel(8,2)
         if (antel(8,3).GT.antal) goto 83
         ansats(8,3) = an83
*     8g
      call Sluggo(8,4,varmax,varupp,varned,ansats,org,lock(8,4),       
     :                             low,start(8,4),stopp(8,4),virtu)
 
      do 84 an84 = start(8,4),stopp(8,4),-1
         antel(8,4) = an84 + antel(8,3)
         if (antel(8,4).GT.antal) goto 84
         ansats(8,4) = an84
*     8h
      call Sluggo(8,5,varmax,varupp,varned,ansats,org,lock(8,5),       
     :                             low,start(8,5),stopp(8,5),virtu)
      do 85 an85 = start(8,5),stopp(8,5),-1
         antel(8,5) = an85 + antel(8,4)
         if (antel(8,5).GT.antal) goto 85
         ansats(8,5) = an85
*     8i
      call Sluggo(8,6,varmax,varupp,varned,ansats,org,lock(8,6),       
     :                             low,start(8,6),stopp(8,6),virtu)
      do 86 an86 = start(8,6),stopp(8,6),-1
         antel(8,6) = an86 + antel(8,5)
         if (antel(8,6).GT.antal) goto 86
         ansats(8,6) = an86
*     8k
      call Sluggo(8,7,varmax,varupp,varned,ansats,org,lock(8,7),       
     :                             low,start(8,7),stopp(8,7),virtu)
      do 87 an87 = start(8,7),stopp(8,7),-1
         antel(8,7) = an87 + antel(8,6)
         if (antel(8,7).GT.antal) goto 87
         ansats(8,7) = an87
*     9s
      call Sluggo(9,0,varmax,varupp,varned,ansats,org,lock(9,0),       
     :                             low,start(9,0),stopp(9,0),virtu)
      do 90 an90 = start(9,0),stopp(9,0),-1
         antel(9,0) = an90 + antel(8,7)
         if (antel(9,0).GT.antal) goto 90
         ansats(9,0) = an90
*     9p
      call Sluggo(9,1,varmax,varupp,varned,ansats,org,lock(9,1),       
     :                             low,start(9,1),stopp(9,1),virtu)
      do 91 an91 = start(9,1),stopp(9,1),-1
         antel(9,1) = an91 + antel(9,0)
         if (antel(9,1).GT.antal) goto 91
         ansats(9,1) = an91
*     9d
      call Sluggo(9,2,varmax,varupp,varned,ansats,org,lock(9,2),       
     :                             low,start(9,2),stopp(9,2),virtu)
      do 92 an92 = start(9,2),stopp(9,2),-1
         antel(9,2) = an92 + antel(9,1)
         if (antel(9,2).GT.antal) goto 92
         ansats(9,2) = an92
*     9f
      call Sluggo(9,3,varmax,varupp,varned,ansats,org,lock(9,3),       
     :                             low,start(9,3),stopp(9,3),virtu)
      do 93 an93 = start(9,3),stopp(9,3),-1
         antel(9,3) = an93 + antel(9,2)
         if (antel(9,3).GT.antal) goto 93
         ansats(9,3) = an93
*     9g
      call Sluggo(9,4,varmax,varupp,varned,ansats,org,lock(9,4),       
     :                             low,start(9,4),stopp(9,4),virtu)
      do 94 an94 = start(9,4),stopp(9,4),-1
         antel(9,4) = an94 + antel(9,3)
         if (antel(9,4).GT.antal) goto 94
         ansats(9,4) = an94
*     9h
      call Sluggo(9,5,varmax,varupp,varned,ansats,org,lock(9,5),       
     :                             low,start(9,5),stopp(9,5),virtu)
      do 95 an95 = start(9,5),stopp(9,5),-1
         antel(9,5) = an95 + antel(9,4)
         if (antel(9,5).GT.antal) goto 95
         ansats(9,5) = an95
*     9i
      call Sluggo(9,6,varmax,varupp,varned,ansats,org,lock(9,6),       
     :                             low,start(9,6),stopp(9,6),virtu)
      do 96 an96 = start(9,6),stopp(9,6),-1
         antel(9,6) = an96 + antel(9,5)
         if (antel(9,6).GT.antal) goto 96
         ansats(9,6) = an96
*     9k
      call Sluggo(9,7,varmax,varupp,varned,ansats,org,lock(9,7),       
     :                             low,start(9,7),stopp(9,7),virtu)
      do 97 an97 = start(9,7),stopp(9,7),-1
         antel(9,7) = an97 + antel(9,6)
         if (antel(9,7).GT.antal) goto 97
         ansats(9,7) = an97
*     9l
      call Sluggo(9,8,varmax,varupp,varned,ansats,org,lock(9,8),       
     :                             low,start(9,8),stopp(9,8),virtu)
      do 98 an98 = start(9,8),stopp(9,8),-1
         antel(9,8) = an98 + antel(9,7)
         if (antel(9,8).GT.antal) goto 98
         ansats(9,8) = an98
*     10s
      call Sluggo(10,0,varmax,varupp,varned,ansats,org,lock(10,0),     
     :                             low,start(10,0),stopp(10,0),virtu)
      do 100 anA0 = start(10,0),stopp(10,0),-1
         antel(10,0) = anA0 + antel(9,8)
         if (antel(10,0).GT.antal) goto 100
         ansats(10,0) = anA0
*     10p
      call Sluggo(10,1,varmax,varupp,varned,ansats,org,lock(10,1),     
     :                             low,start(10,1),stopp(10,1),virtu)
      do 101 anA1 = start(10,1),stopp(10,1),-1
         antel(10,1) = anA1 + antel(10,0)
         if (antel(10,1).GT.antal) goto 101
         ansats(10,1) = anA1
*     10d
      call Sluggo(10,2,varmax,varupp,varned,ansats,org,lock(10,2),     
     :                             low,start(10,2),stopp(10,2),virtu)
      do 102 anA2 = start(10,2),stopp(10,2),-1
         antel(10,2) = anA2 + antel(10,1)
         if (antel(10,2).GT.antal) goto 102
         ansats(10,2) = anA2
*     10f
      call Sluggo(10,3,varmax,varupp,varned,ansats,org,lock(10,3),     
     :                             low,start(10,3),stopp(10,3),virtu)
      do 103 anA3 = start(10,3),stopp(10,3),-1
         antel(10,3) = anA3 + antel(10,2)
         if (antel(10,3).GT.antal) goto 103
         ansats(10,3) = anA3
*     10g
      call Sluggo(10,4,varmax,varupp,varned,ansats,org,lock(10,4),     
     :                             low,start(10,4),stopp(10,4),virtu)
      do 104 anA4 = start(10,4),stopp(10,4),-1
         antel(10,4) = anA4 + antel(10,3)
         if (antel(10,4).GT.antal) goto 104
         ansats(10,4) = anA4
*     10h
      call Sluggo(10,5,varmax,varupp,varned,ansats,org,lock(10,5),     
     :                             low,start(10,5),stopp(10,5),virtu)
      do 105 anA5 = start(10,5),stopp(10,5),-1
         antel(10,5) = anA5 + antel(10,4)
         if (antel(10,5).GT.antal) goto 105
         ansats(10,5) = anA5
*     10i
      call Sluggo(10,6,varmax,varupp,varned,ansats,org,lock(10,6),     
     :                             low,start(10,6),stopp(10,6),virtu)
      do 106 anA6 = start(10,6),stopp(10,6),-1
         antel(10,6) = anA6 + antel(10,5)
         if (antel(10,6).GT.antal) goto 106
         ansats(10,6) = anA6
*     10k
      call Sluggo(10,7,varmax,varupp,varned,ansats,org,lock(10,7),     
     :                             low,start(10,7),stopp(10,7),virtu)
      do 107 anA7 = start(10,7),stopp(10,7),-1
         antel(10,7) = anA7 + antel(10,6)
         if (antel(10,7).GT.antal) goto 107
         ansats(10,7) = anA7
*     10l
      call Sluggo(10,8,varmax,varupp,varned,ansats,org,lock(10,8),     
     :                             low,start(10,8),stopp(10,8),virtu)
      do 108 anA8 = start(10,8),stopp(10,8),-1
         antel(10,8) = anA8 + antel(10,7)
         if (antel(10,8).GT.antal) goto 108
         ansats(10,8) = anA8
*     10m
      call Sluggo(10,9,varmax,varupp,varned,ansats,org,lock(10,9),     
     :                             low,start(10,9),stopp(10,9),virtu)
      do 109 anA9 = start(10,9),stopp(10,9),-1
         antel(10,9) = anA9 + antel(10,8)
         if (antel(10,9).GT.antal) goto 109
         ansats(10,9) = anA9
*     11s
      call Sluggo(11,0,varmax,varupp,varned,ansats,org,lock(11,0),     
     :                             low,start(11,0),stopp(11,0),virtu)
      do 110 anB0 = start(11,0),stopp(11,0),-1
         antel(11,0) = anB0 + antel(10,9)
         if (antel(11,0).GT.antal) goto 110
         ansats(11,0) = anB0
*     11p
      call Sluggo(11,1,varmax,varupp,varned,ansats,org,lock(11,1),     
     :                             low,start(11,1),stopp(11,1),virtu)
      do 111 anB1 = start(11,1),stopp(11,1),-1
         antel(11,1) = anB1 + antel(11,0)
         if (antel(11,1).GT.antal) goto 111
         ansats(11,1) = anB1
*     11d
      call Sluggo(11,2,varmax,varupp,varned,ansats,org,lock(11,2),     
     :                             low,start(11,2),stopp(11,2),virtu)
      do 112 anB2 = start(11,2),stopp(11,2),-1
         antel(11,2) = anB2 + antel(11,1)
         if (antel(11,2).GT.antal) goto 112
         ansats(11,2) = anB2
*     11f
      call Sluggo(11,3,varmax,varupp,varned,ansats,org,lock(11,3),     
     :                             low,start(11,3),stopp(11,3),virtu)
      do 113 anB3 = start(11,3),stopp(11,3),-1
         antel(11,3) = anB3 + antel(11,2)
         if (antel(11,3).GT.antal) goto 113
         ansats(11,3) = anB3
*     11g
      call Sluggo(11,4,varmax,varupp,varned,ansats,org,lock(11,4),     
     :                             low,start(11,4),stopp(11,4),virtu)
      do 114 anB4 = start(11,4),stopp(11,4),-1
         antel(11,4) = anB4 + antel(11,3)
         if (antel(11,4).GT.antal) goto 114
         ansats(11,4) = anB4
*     11h
      call Sluggo(11,5,varmax,varupp,varned,ansats,org,lock(11,5),     
     :                             low,start(11,5),stopp(11,5),virtu)
      do 115 anB5 = start(11,5),stopp(11,5),-1
         antel(11,5) = anB5 + antel(11,4)
         if (antel(11,5).GT.antal) goto 115
         ansats(11,5) = anB5
*     11i
      call Sluggo(11,6,varmax,varupp,varned,ansats,org,lock(11,6),     
     :                             low,start(11,6),stopp(11,6),virtu)
      do 116 anB6 = start(11,6),stopp(11,6),-1
         antel(11,6) = anB6 + antel(11,5)
         if (antel(11,6).GT.antal) goto 116
         ansats(11,6) = anB6
*     11k
      call Sluggo(11,7,varmax,varupp,varned,ansats,org,lock(11,7),     
     :                             low,start(11,7),stopp(11,7),virtu)
      do 117 anB7 = start(11,7),stopp(11,7),-1
         antel(11,7) = anB7 + antel(11,6)
         if (antel(11,7).GT.antal) goto 117
         ansats(11,7) = anB7
*     11l
      call Sluggo(11,8,varmax,varupp,varned,ansats,org,lock(11,8),     
     :                             low,start(11,8),stopp(11,8),virtu)
      do 118 anB8 = start(11,8),stopp(11,8),-1
         antel(11,8) = anB8 + antel(11,7)
         if (antel(11,8).GT.antal) goto 118
         ansats(11,8) = anB8
*     11m
      call Sluggo(11,9,varmax,varupp,varned,ansats,org,lock(11,9),     
     :                             low,start(11,9),stopp(11,9),virtu)
      do 119 anB9 = start(11,9),stopp(11,9),-1
         antel(11,9) = anB9 + antel(11,8)
         if (antel(11,9).GT.antal) goto 119
         ansats(11,9) = anB9
*     11n
      call Sluggo(11,10,varmax,varupp,varned,ansats,org,lock(11,10),   
     :                             low,start(11,10),stopp(11,10),virtu)
      do 1110 anBA = start(11,10),stopp(11,10),-1
         antel(11,10) = anBA + antel(11,9)
         if (antel(11,10).GT.antal) goto 1110
         ansats(11,10) = anBA
*     12s
      call Sluggo(12,0,varmax,varupp,varned,ansats,org,lock(12,0),     
     :                             low,start(12,0),stopp(12,0),virtu)
      do 120 anC0 = start(12,0),stopp(12,0),-1
         antel(12,0) = anC0 + antel(11,10)
         if (antel(12,0).GT.antal) goto 120
         ansats(12,0) = anC0
*     12p
      call Sluggo(12,1,varmax,varupp,varned,ansats,org,lock(12,1),     
     :                             low,start(12,1),stopp(12,1),virtu)
      do 121 anC1 = start(12,1),stopp(12,1),-1
         antel(12,1) = anC1 + antel(12,0)
         if (antel(12,1).GT.antal) goto 121
         ansats(12,1) = anC1
*     12d
      call Sluggo(12,2,varmax,varupp,varned,ansats,org,lock(12,2),     
     :                             low,start(12,2),stopp(12,2),virtu)
      do 122 anC2 = start(12,2),stopp(12,2),-1
         antel(12,2) = anC2 + antel(12,1)
         if (antel(12,2).GT.antal) goto 122
         ansats(12,2) = anC2
*     12f
      call Sluggo(12,3,varmax,varupp,varned,ansats,org,lock(12,3),     
     :                             low,start(12,3),stopp(12,3),virtu)
      do 123 anC3 = start(12,3),stopp(12,3),-1
         antel(12,3) = anC3 + antel(12,2)
         if (antel(12,3).GT.antal) goto 123
         ansats(12,3) = anC3
*     12g
      call Sluggo(12,4,varmax,varupp,varned,ansats,org,lock(12,4),     
     :                             low,start(12,4),stopp(12,4),virtu)
      do 124 anC4 = start(12,4),stopp(12,4),-1
         antel(12,4) = anC4 + antel(12,3)
         if (antel(12,4).GT.antal) goto 124
         ansats(12,4) = anC4
*     12h
      call Sluggo(12,5,varmax,varupp,varned,ansats,org,lock(12,5),     
     :                             low,start(12,5),stopp(12,5),virtu)
      do 125 anC5 = start(12,5),stopp(12,5),-1
         antel(12,5) = anC5 + antel(12,4)
         if (antel(12,5).GT.antal) goto 125
         ansats(12,5) = anC5
*     12i
      call Sluggo(12,6,varmax,varupp,varned,ansats,org,lock(12,6),     
     :                             low,start(12,6),stopp(12,6),virtu)
      do 126 anC6 = start(12,6),stopp(12,6),-1
         antel(12,6) = anC6 + antel(12,5)
         if (antel(12,6).GT.antal) goto 126
         ansats(12,6) = anC6
*     12k
      call Sluggo(12,7,varmax,varupp,varned,ansats,org,lock(12,7),     
     :                             low,start(12,7),stopp(12,7),virtu)
      do 127 anC7 = start(12,7),stopp(12,7),-1
         antel(12,7) = anC7 + antel(12,6)
         if (antel(12,7).GT.antal) goto 127
         ansats(12,7) = anC7
*     12l
      call Sluggo(12,8,varmax,varupp,varned,ansats,org,lock(12,8),     
     :                             low,start(12,8),stopp(12,8),virtu)
      do 128 anC8 = start(12,8),stopp(12,8),-1
         antel(12,8) = anC8 + antel(12,7)
         if (antel(12,8).GT.antal) goto 128
         ansats(12,8) = anC8
*     12m
      call Sluggo(12,9,varmax,varupp,varned,ansats,org,lock(12,9),     
     :                             low,start(12,9),stopp(12,9),virtu)
      do 129 anC9 = start(12,9),stopp(12,9),-1
         antel(12,9) = anC9 + antel(12,8)
         if (antel(12,9).GT.antal) goto 129
         ansats(12,9) = anC9
*     12n
      call Sluggo(12,10,varmax,varupp,varned,ansats,org,lock(12,10),   
     :                             low,start(12,10),stopp(12,10),virtu)
      do 1210 anCA = start(12,10),stopp(12,10),-1
         antel(12,10) = anCA + antel(12,9)
         if (antel(12,10).GT.antal) goto 1210
         ansats(12,10) = anCA
*     13s
      call Sluggo(13,0,varmax,varupp,varned,ansats,org,lock(13,0),     
     :                             low,start(13,0),stopp(13,0),virtu)
      do 130 anD0 = start(13,0),stopp(13,0),-1
         antel(13,0) = anD0 + antel(12,10)
         if (antel(13,0).GT.antal) goto 130
         ansats(13,0) = anD0
*     13p
      call Sluggo(13,1,varmax,varupp,varned,ansats,org,lock(13,1),     
     :                             low,start(13,1),stopp(13,1),virtu)
      do 131 anD1 = start(13,1),stopp(13,1),-1
         antel(13,1) = anD1 + antel(13,0)
         if (antel(13,1).GT.antal) goto 131
         ansats(13,1) = anD1
*     13d
      call Sluggo(13,2,varmax,varupp,varned,ansats,org,lock(13,2),     
     :                             low,start(13,2),stopp(13,2),virtu)
      do 132 anD2 = start(13,2),stopp(13,2),-1
         antel(13,2) = anD2 + antel(13,1)
         if (antel(13,2).GT.antal) goto 132
         ansats(13,2) = anD2
*     13f
      call Sluggo(13,3,varmax,varupp,varned,ansats,org,lock(13,3),     
     :                             low,start(13,3),stopp(13,3),virtu)
      do 133 anD3 = start(13,3),stopp(13,3),-1
         antel(13,3) = anD3 + antel(13,2)
         if (antel(13,3).GT.antal) goto 133
         ansats(13,3) = anD3
*     13g
      call Sluggo(13,4,varmax,varupp,varned,ansats,org,lock(13,4),     
     :                             low,start(13,4),stopp(13,4),virtu)
      do 134 anD4 = start(13,4),stopp(13,4),-1
         antel(13,4) = anD4 + antel(13,3)
         if (antel(13,4).GT.antal) goto 134
         ansats(13,4) = anD4
*     13h
      call Sluggo(13,5,varmax,varupp,varned,ansats,org,lock(13,5),     
     :                             low,start(13,5),stopp(13,5),virtu)
      do 135 anD5 = start(13,5),stopp(13,5),-1
         antel(13,5) = anD5 + antel(13,4)
         if (antel(13,5).GT.antal) goto 135
         ansats(13,5) = anD5
*     13i
      call Sluggo(13,6,varmax,varupp,varned,ansats,org,lock(13,6),     
     :                             low,start(13,6),stopp(13,6),virtu)
      do 136 anD6 = start(13,6),stopp(13,6),-1
         antel(13,6) = anD6 + antel(13,5)
         if (antel(13,6).GT.antal) goto 136
         ansats(13,6) = anD6
*     13k
      call Sluggo(13,7,varmax,varupp,varned,ansats,org,lock(13,7),     
     :                             low,start(13,7),stopp(13,7),virtu)
      do 137 anD7 = start(13,7),stopp(13,7),-1
         antel(13,7) = anD7 + antel(13,6)
         if (antel(13,7).GT.antal) goto 137
         ansats(13,7) = anD7
*     13l
      call Sluggo(13,8,varmax,varupp,varned,ansats,org,lock(13,8),     
     :                             low,start(13,8),stopp(13,8),virtu)
      do 138 anD8 = start(13,8),stopp(13,8),-1
         antel(13,8) = anD8 + antel(13,7)
         if (antel(13,8).GT.antal) goto 138
         ansats(13,8) = anD8
*     13m
      call Sluggo(13,9,varmax,varupp,varned,ansats,org,lock(13,9),     
     :                             low,start(13,9),stopp(13,9),virtu)
      do 139 anD9 = start(13,9),stopp(13,9),-1
         antel(13,9) = anD9 + antel(13,8)
         if (antel(13,9).GT.antal) goto 139
         ansats(13,9) = anD9
*     13n
      call Sluggo(13,10,varmax,varupp,varned,ansats,org,lock(13,10),   
     :                             low,start(13,10),stopp(13,10),virtu)
      do 1310 anDA = start(13,10),stopp(13,10),-1
         antel(13,10) = anDA + antel(13,9)
         if (antel(13,10).GT.antal) goto 1310
         ansats(13,10) = anDA
*     14s
      call Sluggo(14,0,varmax,varupp,varned,ansats,org,lock(14,0),     
     :                             low,start(14,0),stopp(14,0),virtu)
      do 140 anE0 = start(14,0),stopp(14,0),-1
         antel(14,0) = anE0 + antel(13,10)
         if (antel(14,0).GT.antal) goto 140
         ansats(14,0) = anE0
*     14p
      call Sluggo(14,1,varmax,varupp,varned,ansats,org,lock(14,1),     
     :                             low,start(14,1),stopp(14,1),virtu)
      do 141 anE1 = start(14,1),stopp(14,1),-1
         antel(14,1) = anE1 + antel(14,0)
         if (antel(14,1).GT.antal) goto 141
         ansats(14,1) = anE1
*     14d
      call Sluggo(14,2,varmax,varupp,varned,ansats,org,lock(14,2),     
     :                             low,start(14,2),stopp(14,2),virtu)
      do 142 anE2 = start(14,2),stopp(14,2),-1
         antel(14,2) = anE2 + antel(14,1)
         if (antel(14,2).GT.antal) goto 142
         ansats(14,2) = anE2
*     14f
      call Sluggo(14,3,varmax,varupp,varned,ansats,org,lock(14,3),     
     :                             low,start(14,3),stopp(14,3),virtu)
      do 143 anE3 = start(14,3),stopp(14,3),-1
         antel(14,3) = anE3 + antel(14,2)
         if (antel(14,3).GT.antal) goto 143
         ansats(14,3) = anE3
*     14g
      call Sluggo(14,4,varmax,varupp,varned,ansats,org,lock(14,4),     
     :                             low,start(14,4),stopp(14,4),virtu)
      do 144 anE4 = start(14,4),stopp(14,4),-1
         antel(14,4) = anE4 + antel(14,3)
         if (antel(14,4).GT.antal) goto 144
         ansats(14,4) = anE4
*     14h
      call Sluggo(14,5,varmax,varupp,varned,ansats,org,lock(14,5),     
     :                             low,start(14,5),stopp(14,5),virtu)
      do 145 anE5 = start(14,5),stopp(14,5),-1
         antel(14,5) = anE5 + antel(14,4)
         if (antel(14,5).GT.antal) goto 145
         ansats(14,5) = anE5
*     14i
      call Sluggo(14,6,varmax,varupp,varned,ansats,org,lock(14,6),     
     :                             low,start(14,6),stopp(14,6),virtu)
      do 146 anE6 = start(14,6),stopp(14,6),-1
         antel(14,6) = anE6 + antel(14,5)
         if (antel(14,6).GT.antal) goto 146
         ansats(14,6) = anE6
*     14k
      call Sluggo(14,7,varmax,varupp,varned,ansats,org,lock(14,7),     
     :                             low,start(14,7),stopp(14,7),virtu)
      do 147 anE7 = start(14,7),stopp(14,7),-1
         antel(14,7) = anE7 + antel(14,6)
         if (antel(14,7).GT.antal) goto 147
         ansats(14,7) = anE7
*     14l
      call Sluggo(14,8,varmax,varupp,varned,ansats,org,lock(14,8),     
     :                             low,start(14,8),stopp(14,8),virtu)
      do 148 anE8 = start(14,8),stopp(14,8),-1
         antel(14,8) = anE8 + antel(14,7)
         if (antel(14,8).GT.antal) goto 148
         ansats(14,8) = anE8
*     14m
      call Sluggo(14,9,varmax,varupp,varned,ansats,org,lock(14,9),     
     :                             low,start(14,9),stopp(14,9),virtu)
      do 149 anE9 = start(14,9),stopp(14,9),-1
         antel(14,9) = anE9 + antel(14,8)
         if (antel(14,9).GT.antal) goto 149
         ansats(14,9) = anE9
*     14n
      call Sluggo(14,10,varmax,varupp,varned,ansats,org,lock(14,10),   
     :                             low,start(14,10),stopp(14,10),virtu)
      do 1410 anEA = start(14,10),stopp(14,10),-1
         antel(14,10) = anEA + antel(14,9)
         if (antel(14,10).GT.antal) goto 1410
         ansats(14,10) = anEA
*     15s
      call Sluggo(15,0,varmax,varupp,varned,ansats,org,lock(15,0),     
     :                             low,start(15,0),stopp(15,0),virtu)
      do 150 anF0 = start(15,0),stopp(15,0),-1
         antel(15,0) = anF0 + antel(14,10)
         if (antel(15,0).GT.antal) goto 150
         ansats(15,0) = anF0
*     15p
      call Sluggo(15,1,varmax,varupp,varned,ansats,org,lock(15,1),     
     :                             low,start(15,1),stopp(15,1),virtu)
      do 151 anF1 = start(15,1),stopp(15,1),-1
         antel(15,1) = anF1 + antel(15,0)
         if (antel(15,1).GT.antal) goto 151
         ansats(15,1) = anF1
*     15d
      call Sluggo(15,2,varmax,varupp,varned,ansats,org,lock(15,2),     
     :                             low,start(15,2),stopp(15,2),virtu)
      do 152 anF2 = start(15,2),stopp(15,2),-1
         antel(15,2) = anF2 + antel(15,1)
         if (antel(15,2).GT.antal) goto 152
         ansats(15,2) = anF2
*     15f
      call Sluggo(15,3,varmax,varupp,varned,ansats,org,lock(15,3),     
     :                             low,start(15,3),stopp(15,3),virtu)
      do 153 anF3 = start(15,3),stopp(15,3),-1
         antel(15,3) = anF3 + antel(15,2)
         if (antel(15,3).GT.antal) goto 153
         ansats(15,3) = anF3
*     15g
      call Sluggo(15,4,varmax,varupp,varned,ansats,org,lock(15,4),     
     :                             low,start(15,4),stopp(15,4),virtu)
      do 154 anF4 = start(15,4),stopp(15,4),-1
         antel(15,4) = anF4 + antel(15,3)
         if (antel(15,4).GT.antal) goto 154
         ansats(15,4) = anF4
*     15h
      call Sluggo(15,5,varmax,varupp,varned,ansats,org,lock(15,5),     
     :                             low,start(15,5),stopp(15,5),virtu)
      do 155 anF5 = start(15,5),stopp(15,5),-1
         antel(15,5) = anF5 + antel(15,4)
         if (antel(15,5).GT.antal) goto 155
         ansats(15,5) = anF5
*     15i
      call Sluggo(15,6,varmax,varupp,varned,ansats,org,lock(15,6),     
     :                             low,start(15,6),stopp(15,6),virtu)
      do 156 anF6 = start(15,6),stopp(15,6),-1
         antel(15,6) = anF6 + antel(15,5)
         if (antel(15,6).GT.antal) goto 156
         ansats(15,6) = anF6
*     15k
      call Sluggo(15,7,varmax,varupp,varned,ansats,org,lock(15,7),     
     :                             low,start(15,7),stopp(15,7),virtu)
      do 157 anF7 = start(15,7),stopp(15,7),-1
         antel(15,7) = anF7 + antel(15,6)
         if (antel(15,7).GT.antal) goto 157
         ansats(15,7) = anF7
*     15l
      call Sluggo(15,8,varmax,varupp,varned,ansats,org,lock(15,8),     
     :                             low,start(15,8),stopp(15,8),virtu)
      do 158 anF8 = start(15,8),stopp(15,8),-1
         antel(15,8) = anF8 + antel(15,7)
         if (antel(15,8).GT.antal) goto 158
         ansats(15,8) = anF8
*     15m
      call Sluggo(15,9,varmax,varupp,varned,ansats,org,lock(15,9),     
     :                             low,start(15,9),stopp(15,9),virtu)
      do 159 anF9 = start(15,9),stopp(15,9),-1
         antel(15,9) = anF9 + antel(15,8)
         if (antel(15,9).GT.antal) goto 159
         ansats(15,9) = anF9
*     15n
      call Sluggo(15,10,varmax,varupp,varned,ansats,org,lock(15,10),   
     :                             low,start(15,10),stopp(15,10),virtu)
      do 1510 anFA = start(15,10),stopp(15,10),-1
         antel(15,10) = anFA + antel(15,9)
         if (antel(15,10).NE.antal) goto 1510
         ansats(15,10) = anFA
         par = 0
         do 6 i=1,nmax
            do 6 j=0,min(10,i-1)
    6          par = mod(par+j*ansats(i,j),2)
         if (par.EQ.par0) then
            cf = cf+1
            do 1600 i=1,nmax
C              write(*,5000) (ansats(i,j),j=0,min(10,i-1))
 1600          write(fil,5000) (ansats(i,j),j=0,min(10,i-1))
         endif
 1510 continue
  159 continue
  158 continue
  157 continue
  156 continue
  155 continue
  154 continue
  153 continue
  152 continue
  151 continue
  150 continue
 1410 continue
  149 continue
  148 continue
  147 continue
  146 continue
  145 continue
  144 continue
  143 continue
  142 continue
  141 continue
  140 continue
 1310 continue
  139 continue
  138 continue
  137 continue
  136 continue
  135 continue
  134 continue
  133 continue
  132 continue
  131 continue
  130 continue
 1210 continue
  129 continue
  128 continue
  127 continue
  126 continue
  125 continue
  124 continue
  123 continue
  122 continue
  121 continue
  120 continue
 1110 continue
  119 continue
  118 continue
  117 continue
  116 continue
  115 continue
  114 continue
  113 continue
  112 continue
  111 continue
  110 continue
  109 continue
  108 continue
  107 continue
  106 continue
  105 continue
  104 continue
  103 continue
  102 continue
  101 continue
  100 continue
   98 continue
   97 continue
   96 continue
   95 continue
   94 continue
   93 continue
   92 continue
   91 continue
   90 continue
   87 continue
   86 continue
   85 continue
   84 continue
   83 continue
   82 continue
   81 continue
   80 continue
   76 continue
   75 continue
   74 continue
   73 continue
   72 continue
   71 continue
   70 continue
   65 continue
   64 continue
   63 continue
   62 continue
   61 continue
   60 continue
   54 continue
   53 continue
   52 continue
   51 continue
   50 continue
   43 continue
   42 continue
   41 continue
   40 continue
   32 continue
   31 continue
   30 continue
   21 continue
   20 continue
   10 continue
 5000 format(11I2)
C     write(*,*) 'cf =',cf
      return
      end
*     last edited April 21, 1993
      logical function Block(rank,pos,nr,maskal,ny,lika)
      integer noofco
Ctc   parameter (noofco=200000)
      parameter (noofco=2100000)
      integer rank(0:noofco,0:10,15),pos(0:noofco),maskal
      integer i,j,ytter,nasta,denna,ant
      logical ny,lika(8)

      denna = pos(nr)
      nasta = pos(nr+1)
      Block = .FALSE.
      ytter = 0
      do 10 i=maskal,(maskal+3)/2,-1
         do 10 j=1,2
   10       if (rank(denna,j,i) .NE. rank(nasta,j,i)) return
      do 20 i=(maskal+1)/2,1,-1
         if (ytter.EQ.0 .AND. rank(denna,3,i).GE.0) ytter=i
         do 20 j=1,3
   20       if (rank(denna,j,i) .NE. rank(nasta,j,i)) return
C     if (ny) then
C        ant = 0
C        do 30 i=1,8
C           if (i.LE.ytter) then
C              lika(i) = rank(denna,4,i) .EQ. rank(nasta,4,i)
C           else
C              lika(i) = .TRUE.
C           endif
C  30       if (.NOT.lika(i)) ant = ant+1
C        Block = ant .EQ. 1
C        return
C     else
C	 do 40 i=1,ytter
C  40	    if (lika(i) .NEQV.
C    f          rank(denna,4,i) .EQ. rank(nasta,4,i)) return
C     endif
      Block = .TRUE.
      return
      end
*     last edited April 1, 1993
      subroutine Bubbla(rank,pos,antal,maskal,nmax,lmax,opt)
      integer    noofco
Ctc   parameter (noofco=200000)
      parameter (noofco=2100000)
      integer    rank(0:noofco,0:10,15),antal,pos(0:noofco),maskal,
     :           nmax,lmax
      logical    opt
      integer    nr,starta,low,sluta,high
      logical    Change,dirty,lika(noofco),Kolla

      pos(0) = 0
      write(*,*) 'number of items =',antal
      do 10 nr=1,antal
         lika(nr) = opt .AND. Kolla(rank,nr,maskal) 
*        if (lika(nr)) write(*,*) nr,' e lik!'
   10    pos(nr) = nr

      if (opt) then
         low = 1
      else
         low = 0
      endif
      high = antal-2

   15 starta = low
      sluta  = high
      write(*,*) sluta-starta
      dirty  = .False.
      do 20 nr=starta,sluta
         if (lika(pos(nr)) .EQV. lika(pos(nr+1))) then
            if (Change(rank,pos,nr,maskal,nmax,lmax,opt)) then
               high  = nr-1
               dirty = .True.
            endif
         elseif (.NOT.lika(pos(nr)) .AND. lika(pos(nr+1))) then
            call Skifta(pos,nr)
            high  = nr-1
            dirty = .True.
         endif
   20 continue
      if (.not.dirty .OR. low.GT.high) return
      starta = high
      sluta  = low
      write(*,*) starta-sluta
      dirty  = .False.
      do 30 nr=starta,sluta,-1
         if (lika(pos(nr)) .EQV. lika(pos(nr+1))) then
            if (Change(rank,pos,nr,maskal,nmax,lmax,opt)) then
               low   = nr+1
               dirty = .True.
            endif
         elseif (.NOT.lika(pos(nr)) .AND. lika(pos(nr+1))) then
            call Skifta(pos,nr)
            low   = nr+1
            dirty = .True.
         endif
   30 continue
      if (dirty .AND. low.LT.high) goto 15

      return
      end
*     last edited April 1, 1993
      logical function Change(rank,pos,nr,maskal,nmax,lmax,opt)
      integer noofco
Ctc   parameter (noofco=200000)
      parameter (noofco=2100000)
      integer rank(0:noofco,0:10,15),antal,pos(0:noofco),maskal,nmax
      integer lmax
      logical opt
      integer i,j,temp

      if (opt) then
         do 10 i=1,maskal
            do 10 j=1,2
               if (rank(pos(nr),j,i).LT.rank(pos(nr+1),j,i)) then
                  Change = .False.
                  return
               elseif (rank(pos(nr),j,i).GT.rank(pos(nr+1),j,i)) then
                  goto 40
               endif
   10    continue
         do 20 i=1,(maskal+1)/2
            if (rank(pos(nr),3,i).LT.rank(pos(nr+1),3,i)) then
               Change = .False.
               return
            elseif (rank(pos(nr),3,i).GT.rank(pos(nr+1),3,i)) then
               goto 40
            endif
   20    continue
      else
         do 30 i=1,nmax
            do 30 j=0,min(i-1,lmax)
               if (rank(pos(nr),j,i).GT.rank(pos(nr+1),j,i)) then
                  Change = .False.
                  return
               elseif (rank(pos(nr),j,i).LT.rank(pos(nr+1),j,i)) then
                  goto 40
               endif
   30    continue
      endif
      Change    = .False.
      return
   40 Change    = .True.
      call Skifta(pos,nr)
      return
      end
*     last edited February 21, 1993
      subroutine Gen(ansats,posn,posl,resS,resL,skal,
     :                   cf,dyn,first,breit,J2min,J2max,minS,minL,check)
      integer ansats(1:15,0:10),koppl(0:10,0:5),SKVANT(0:10,0:6,21)
      integer EXTRA(0:10,0:6,21),LKVANT(0:10,0:6,21),antmax(0:10)
      integer resS,resL,pos,i,S(8),L(8),SK(8),LK(8),orbit(8),antel(8)
      integer i1,i2,i3,i4,i5,i6,i7,i8,E(8),SK1,SK2,SK3,SK4,SK5,SK6,SK7
      integer LK1,LK2,LK3,LK4,LK5,LK6,LK7,stopp,fil,antko(10),skal,cf
      integer posn(110),posl(110),J2min,J2max,minS,minL
      character rad1*80,rad2*72,L1(0:10)
      logical dyn,first,breit,check
      data (L1(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/
*     The value of antmax(l-number) is the maximum number of electrons
*     in the orbital
      data (antmax(i),i=0,10)/2,6,10,14,18,22,26,30,34,38,42/
*     The value of koppl(l-number,number of electrons) is the number of 
*     possible couplings for a certain orbital. If the orbital is
*     populated with more than half of the maximal number of electrons
*     the index "number of electrons" should be substituted with
*     "antmax(l-number) - number of electrons".
      data (koppl(0,i),i=0,1)/1,1/                                   
*     l=0
      data (koppl(1,i),i=0,3)/1,1,3,3/                               
*     l=1
      data (koppl(2,i),i=0,5)/1,1,5,8,16,16/                         
*     l=2
      data (koppl(3,i),i=0,2)/1,1,7/                                 
*     l=3
*     data (koppl(4,i),i=0,2)/1,1,9/                                 
      data (koppl(4,i),i=0,2)/1,1,7/                                 
*     l=4 Limit; doubled occoupied high angular
*     data (koppl(5,i),i=0,2)/1,1,11/                                
      data (koppl(5,i),i=0,2)/1,1,7/                                 
*     l=5 Limit  momentum orbits may not coupple
*     data (koppl(6,i),i=0,2)/1,1,13/                                
      data (koppl(6,i),i=0,2)/1,1,7/                                 
*     l=6 Limit  higher than "I".
*     data (koppl(7,i),i=0,2)/1,1,15/                                
      data (koppl(7,i),i=0,2)/1,1,7/                                 
*     l=7 Limit
*     data (koppl(8,i),i=0,2)/1,1,17/                                
      data (koppl(8,i),i=0,2)/1,1,7/                                 
*     l=8 Limit
*     data (koppl(9,i),i=0,2)/1,1,19/                                
      data (koppl(9,i),i=0,2)/1,1,7/                                 
*     l=9 Limit
*     data (koppl(10,i),i=0,2)/1,1,21/                               
      data (koppl(10,i),i=0,2)/1,1,7/                                
*     l=10 Limit
*     The couplings and there order to appear in the list.
*  SKVANT(l-number, number of electrons, coupling number) is the spin
*  LKVANT(l-number, number of electrons, coupling number) is L-number
*  EXTRA(l-number, number of electrons, coupling number) is the
*     seniority number
*     If orbitals with l-number higher than 3 should be allowed to be
*     populated with more than 2 electrons, changes in koppl(,) have to
*     be done and a lot of new SKVANT(,,), LKVANT(,,) and EXTRA(,,) have
*     to be added to the lists below.
*      no:    1   2   3   4   5   6   7   8   9  10  11  12  13  14  15 
*            16  17  18
*     s( 0): 1S0
*     s( 1): 2S1
*     p( 0): 1S0
*     p( 1): 2P1
*     p( 2): 1S0,1D2,3P2
*     p( 3): 2P1,2D3,4S3
*     d( 0): 1S0
*     d( 1): 2D1
*     d( 2): 1S0,1D2,1G2,3P2,3F2
*     d( 3): 2P3,2D1,2D3,2F3,2G3,2H3,4P3,4F3
*     d( 4): 1S0,1S4,1D2,1D4,1F4,1G2,1G4,1I4,3P2,3P4,3D4,3F2,3F4,3G4,3H4
*           ,5D4
*     d( 5): 2S5,2P3,2D1,2D3,2D5,2F3,2F5,2G3,2G5,2H3,2I5,4P3,4D5,4F3,4G5
*           ,6S5
*     :( 0): 2S0
*     :( 1): 2F1
*     :( 2): 1S0,1D2,1G2,1I2,3P2,3F2,3H2
*     g( 0): 1S0
*     g( 1): 2G1
*     g( 2): 1S0,1D2,1G2,1I2,3P2,3F2,3H2,(1L2,3K2)
*     h( 0): 1S0
*     h( 1): 2H1
*     h( 2): 1S0,1D2,1G2,1I2,3P2,3F2,3H2,(1L2,1N2,3K2,3M2)
*     i( 0): 1S0
*     i( 1): 2I1
*     i( 2): 1S0,1D2,1G2,1I2,3P2,3F2,3H2,(1L2,1N2,1Q2,3K2,3M2,3O2)
*     k( 0): 1S0
*     k( 1): 2K1
*     k( 2): 1S0,1D2,1G2,1I2,3P2,3F2,3H2,
*           (1L2,1N2,1Q2,1T2,3K2,3M2,3O2,3R2)
*     l( 0): 1S0
*     l( 1): 2L1
*     l( 2): 1S0,1D2,1G2,1I2,3P2,3F2,3H2,
*           (1L2,1N2,1Q2,1T2,1V2,3K2,3M2,3O2,3R2,3U2)
*     m( 0): 1S0
*     m( 1): 2M1
*     m( 2): 1S0,1D2,1G2,1I2,3P2,3F2,3H2,
*           (1L2,1N2,1Q2,1T2,1V2,1X2,3K2,3M2,3O2,3R2,3U2,3W2)
*     n( 0): 1S0
*     n( 1): 2N1
*     n( 2): 1S0,1D2,1G2,1I2,3P2,3F2,3H2,
*           (1L2,1N2,1Q2,1T2,1V2,1X2,1Z2,3K2,3M2,3O2,3R2,3U2,3W2,3Y2)
      data  SKVANT(0,0,1)        /1/                                 
*     l=0 #=0
      data  SKVANT(0,1,1)        /2/                                 
*     l=0 #=1
      data  SKVANT(1,0,1)        /1/                                 
*     l=1 #=0
      data  SKVANT(1,1,1)        /2/                                 
*     l=1 #=1
      data (SKVANT(1,2,i),i=1,3) /1,1,3/                             
*     l=1 #=2
      data (SKVANT(1,3,i),i=1,3) /2,2,4/                             
*     l=1 #=3
      data  SKVANT(2,0,1)        /1/                                 
*     l=2 #=0
      data  SKVANT(2,1,1)        /2/                                 
*     l=2 #=1
      data (SKVANT(2,2,i),i=1,5) /1,1,1,3,3/                         
*     l=2 #=2
      data (SKVANT(2,3,i),i=1,8) /2,2,2,2,2,2,4,4/                   
*     l=2 #=3
      data (SKVANT(2,4,i),i=1,16)/1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,5/   
*     l=2 #=4
      data (SKVANT(2,5,i),i=1,16)/2,2,2,2,2,2,2,2,2,2,2,4,4,4,4,6/   
*     l=2 #=5
      data  SKVANT(3,0,1)        /1/                                 
*     l=3 #=0
      data  SKVANT(3,1,1)        /2/                                 
*     l=3 #=1
      data (SKVANT(3,2,i),i=1,7) /1,1,1,1,3,3,3/                     
*     l=3 #=2
      data  SKVANT(4,0,1)        /1/                                 
*     l=4 #=0
      data  SKVANT(4,1,1)        /2/                                 
*     l=4 #=1
*     data (SKVANT(4,2,i),i=1,9) /1,1,1,1,1,3,3,3,3/                 
      data (SKVANT(4,2,i),i=1,7) /1,1,1,1,3,3,3/                     
*     l=4 #=2 Limit
      data  SKVANT(5,0,1)        /1/                                 
*     l=5 #=0
      data  SKVANT(5,1,1)        /2/                                 
*     l=5 #=1
*     data (SKVANT(5,2,i),i=1,11)/1,1,1,1,1,1,3,3,3,3,3/             
      data (SKVANT(5,2,i),i=1,7) /1,1,1,1,3,3,3/                     
*     l=5 #=2 Limit
      data  SKVANT(6,0,1)        /1/                                 
*     l=6 #=0
      data  SKVANT(6,1,1)        /2/                                 
*     l=6 #=1
*     data (SKVANT(6,2,i),i=1,13)/1,1,1,1,1,1,1,3,3,3,3,3,3/         
      data (SKVANT(6,2,i),i=1,7) /1,1,1,1,3,3,3/                     
*     l=6 #=2 Limit
      data  SKVANT(7,0,1)        /1/                                 
*     l=7 #=0
      data  SKVANT(7,1,1)        /2/                                 
*     l=7 #=1
*     data (SKVANT(7,2,i),i=1,15)/1,1,1,1,1,1,1,1,3,3,3,3,3,3,3/     
      data (SKVANT(7,2,i),i=1,7) /1,1,1,1,3,3,3/                     
*     l=7 #=2 Limit
      data  SKVANT(8,0,1)        /1/                                 
*     l=8 #=0
      data  SKVANT(8,1,1)        /2/                                 
*     l=8 #=1
*     data (SKVANT(8,2,i),i=1,17)/1,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3/ 
      data (SKVANT(8,2,i),i=1,7) /1,1,1,1,3,3,3/                     
*     l=8 #=2 Limit
      data  SKVANT(9,0,1)        /1/                                 
*     l=9 #=0
      data  SKVANT(9,1,1)        /2/                                 
*     l=9 #=1
*     data (SKVANT(9,2,i),i=1,19)/1,1,1,1,1,1,1,1,1,1,3,3,3,3,3,3,3,3,3/
      data (SKVANT(9,2,i),i=1,7) /1,1,1,1,3,3,3/                     
*     l=9 #=2 Limit
      data  SKVANT(10,0,1)        /1/                                
*     l=10 #=0
      data  SKVANT(10,1,1)        /2/                                
*     l=10 #=1
*     data (SKVANT(10,2,i),i=1,21)/1,1,1,1,1,1,1,1,1,1,1,
*    f                                            3,3,3,3,3,3,3,3,3,3/  
      data (SKVANT(10,2,i),i=1,7) /1,1,1,1,3,3,3/                    
*     l=10 #=2 Limit
      data  LKVANT(0,0,1)        /0/                                 
*     l=0 #=0
      data  LKVANT(0,1,1)        /0/                                 
*     l=0 #=1
      data  LKVANT(1,0,1)        /0/                                 
*     l=1 #=0
      data  LKVANT(1,1,1)        /1/                                 
*     l=1 #=1
      data (LKVANT(1,2,i),i=1,3) /0,2,1/                             
*     l=1 #=2
      data (LKVANT(1,3,i),i=1,3) /1,2,0/                             
*     l=1 #=3
      data  LKVANT(2,0,1)        /0/                                 
*     l=2 #=0
      data  LKVANT(2,1,1)        /2/                                 
*     l=2 #=1
      data (LKVANT(2,2,i),i=1,5) /0,2,4,1,3/                         
*     l=2 #=2
      data (LKVANT(2,3,i),i=1,8) /1,2,2,3,4,5,1,3/                   
*     l=2 #=3
      data (LKVANT(2,4,i),i=1,16)/0,0,2,2,3,4,4,6,1,1,2,3,3,4,5,2/   
*     l=2 #=4
      data (LKVANT(2,5,i),i=1,16)/0,1,2,2,2,3,3,4,4,5,6,1,2,3,4,0/   
*     l=2 #=5
      data  LKVANT(3,0,1)        /0/                                 
*     l=3 #=0
      data  LKVANT(3,1,1)        /3/                                 
*     l=3 #=1
      data (LKVANT(3,2,i),i=1,7) /0,2,4,6,1,3,5/                     
*     l=3 #=2
      data  LKVANT(4,0,1)        /0/                                 
*     l=4 #=0
      data  LKVANT(4,1,1)        /4/                                 
*     l=4 #=1
*     data (LKVANT(4,2,i),i=1,9) /0,2,4,6,8,1,3,5,7/                 
*     l=4 #=2
      data (LKVANT(4,2,i),i=1,7) /0,2,4,6,1,3,5/                     
*     l=4 #=2 Limit
      data  LKVANT(5,0,1)        /0/                                 
*     l=5 #=0
      data  LKVANT(5,1,1)        /5/                                 
*     l=5 #=1
*     data (LKVANT(5,2,i),i=1,11)/0,2,4,6,8,10,1,3,5,7,9/            
      data (LKVANT(5,2,i),i=1,7) /0,2,4,6,1,3,5/                     
*     l=5 #=2 Limit
      data  LKVANT(6,0,1)        /0/                                 
*     l=6 #=0
      data  LKVANT(6,1,1)        /6/                                 
*     l=6 #=1
*     data (LKVANT(6,2,i),i=1,13)/0,2,4,6,8,10,12,1,3,5,7,9,11/      
      data (LKVANT(6,2,i),i=1,7) /0,2,4,6,1,3,5/                     
*     l=6 #=2 Limit
      data  LKVANT(7,0,1)        /0/                                 
*     l=7 #=0
      data  LKVANT(7,1,1)        /7/                                 
*     l=7 #=1
*     data (LKVANT(7,2,i),i=1,15)/0,2,4,6,8,10,12,14,1,3,5,7,9,11,13/
      data (LKVANT(7,2,i),i=1,7) /0,2,4,6,1,3,5/                     
*     l=7 #=2 Limit
      data  LKVANT(8,0,1)        /0/                                 
*     l=8 #=0
      data  LKVANT(8,1,1)        /8/                                 
*     l=8 #=1
*     data (LKVANT(8,2,i),i=1,17)/0,2,4,6,8,10,12,14,16,
*    :                                         1,3,5,7,9,11,13,15/   
      data (LKVANT(8,2,i),i=1,7) /0,2,4,6,1,3,5/                     
*     l=8 #=2 Limit
      data  LKVANT(9,0,1)        /0/                                 
*     l=9 #=0
      data  LKVANT(9,1,1)        /9/                                 
*     l=9 #=1
*     data (LKVANT(9,2,i),i=1,19)/0,2,4,6,9,10,12,14,16,18,
*    :                                         1,3,5,7,9,11,13,15,17/
      data (LKVANT(9,2,i),i=1,7) /0,2,4,6,1,3,5/                     
*     l=9 #=2 Limit
      data  LKVANT(10,0,1)        /0/                                
*     l=10 #=0
      data  LKVANT(10,1,1)        /10/                               
*     l=10 #=1
*     data (LKVANT(10,2,i),i=1,21)/0,2,4,6,9,10,12,14,16,18,20
*    :                                         1,3,5,7,9,11,13,15,17,19/
      data (LKVANT(10,2,i),i=1,7) /0,2,4,6,1,3,5/                    
*     l=10 #=2 Limit
      data  EXTRA(0,0,1)        /0/
      data  EXTRA(0,1,1)        /1/
      data  EXTRA(1,0,1)        /0/
      data  EXTRA(1,1,1)        /1/
      data (EXTRA(1,2,i),i=1,3) /0,2,2/
      data (EXTRA(1,3,i),i=1,3) /1,3,3/
      data  EXTRA(2,0,1)        /0/
      data  EXTRA(2,1,1)        /1/
      data (EXTRA(2,2,i),i=1,5) /0,2,2,2,2/
      data (EXTRA(2,3,i),i=1,8) /3,1,3,3,3,3,3,3/
      data (EXTRA(2,4,i),i=1,16)/0,4,2,4,4,2,4,4,2,4,4,2,4,4,4,4/
      data (EXTRA(2,5,i),i=1,16)/5,3,1,3,5,3,5,3,5,3,5,3,5,3,5,5/
      data  EXTRA(3,0,1)        /0/
      data  EXTRA(3,1,1)        /1/
      data (EXTRA(3,2,i),i=1,7) /0,2,2,2,2,2,2/
      data  EXTRA(4,0,1)        /0/
      data  EXTRA(4,1,1)        /1/
      data (EXTRA(4,2,i),i=1,9) /0,2,2,2,2,2,2,2,2/
      data  EXTRA(5,0,1)        /0/
      data  EXTRA(5,1,1)        /1/
      data (EXTRA(5,2,i),i=1,11)/0,2,2,2,2,2,2,2,2,2,2/
      data  EXTRA(6,0,1)        /0/
      data  EXTRA(6,1,1)        /1/
      data (EXTRA(6,2,i),i=1,13)/0,2,2,2,2,2,2,2,2,2,2,2,2/
      data  EXTRA(7,0,1)        /0/
      data  EXTRA(7,1,1)        /1/
      data (EXTRA(7,2,i),i=1,15)/0,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
      data  EXTRA(8,0,1)        /0/
      data  EXTRA(8,1,1)        /1/
      data (EXTRA(8,2,i),i=1,17)/0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
      data  EXTRA(9,0,1)        /0/
      data  EXTRA(9,1,1)        /1/
      data (EXTRA(9,2,i),i=1,19)/0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
      data  EXTRA(10,0,1)        /0/
      data  EXTRA(10,1,1)        /1/
      data (EXTRA(10,2,i),i=1,21)
     :                       /0,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2/
      if (first) then
         fil = 7
      else
         fil = 8
      endif
      do 12 i=1,8
   12    antko(i) = 1
      pos = 0
      do 20 i=1,110
         if (ansats(posn(i),posl(i)).NE.0) then
            rad1(pos*8+1:pos*8+2)  = '  '
            rad1(pos*8+3:pos*8+3) = char(48+posn(i))
            rad1(pos*8+4:pos*8+4) = L1(posl(i))
            rad1(pos*8+5:pos*8+8) = '(  )'
            if (ansats(posn(i),posl(i)).GE.10) then
               rad1(pos*8+6:pos*8+6) =
     :                               char(ansats(posn(i),posl(i))/10+48)
            else
               rad1(pos*8+6:pos*8+6) = ' '
            endif
            rad1(pos*8+7:pos*8+7) =
     :                          char(mod(ansats(posn(i),posl(i)),10)+48)
            pos = pos + 1
            if (pos.GT.skal) return
            orbit(pos) = posl(i)
            antel(pos) = min(ansats(posn(i),posl(i)),
     :                          antmax(posl(i))-ansats(posn(i),posl(i)))
            antko(pos) = koppl(posl(i),antel(pos))
         endif
   20 continue
      if (pos.EQ.0) return
      if (dyn) then
         stopp =  8*pos - 4
      else
         stopp = 16*pos - 8
      endif
      do 90 i1=1,antko(1)
       do 90 i2=1,antko(2)
        do 90 i3=1,antko(3)
         do 90 i4=1,antko(4)
          do 90 i5=1,antko(5)
           do 90 i6=1,antko(6)
            do 90 i7=1,antko(7)
             do 90 i8=1,antko(8)
              S(1) = SKVANT(orbit(1),antel(1),i1)
              L(1) = LKVANT(orbit(1),antel(1),i1)
              E(1) = EXTRA(orbit(1),antel(1),i1)
              if (breit) then
               if (pos.EQ.1) then
                if (S(1).GE.minS .AND. S(1).LE.resS .AND.
     :              L(1).GE.minL .AND. L(1).LE.resL .AND.
     :              S(1)+2*L(1)-1.GE.J2min .AND.
     :              abs(S(1)-2*L(1)-1).LE.J2max) then
                 if (check) then
                  write(*,999) ' ',rad1(1:8)
                  return
                 endif
                 call KOPP1(pos,rad2,S,L,E,dyn)
                 write(fil,999) rad1(1:8)
                 write(fil,999) rad2(1:stopp)
                 cf = cf + 1
                endif
               else
                S(2) = SKVANT(orbit(2),antel(2),i2)
                L(2) = LKVANT(orbit(2),antel(2),i2)
                E(2) = EXTRA(orbit(2),antel(2),i2)
                do 25 SK1=abs(S(1)-S(2))+1,S(1)+S(2)-1,2
                 do 25 LK1=abs(L(1)-L(2)),L(1)+L(2)
*Limit; keeps the intermediate couplings less than or equal to "N".
                  if (LK1.GT.10) goto 25
                  SK(1) = SK1
                  LK(1) = LK1
                  if (pos.EQ.2) then
                   if (SK(1).GE.minS .AND. SK(1).LE.resS .AND.
     :                 LK(1).GE.minL .AND. LK(1).LE.resL .AND.
     :                 SK(1)+2*LK(1)-1.GE.J2min .AND.
     :                 abs(SK(1)-2*LK(1)-1).LE.J2max) then
                    if (check) then
                     write(*,999) ' ',rad1(1:16)
                     return
                    endif
                    call KOPP1(pos,rad2,S,L,E,dyn)
                    call KOPP2(pos,rad2,SK,LK,dyn)
                    write(fil,999) rad1(1:16)
                    write(fil,999) rad2(1:stopp)
                    cf = cf + 1
                   endif
                  else
                   S(3) = SKVANT(orbit(3),antel(3),i3)
                   L(3) = LKVANT(orbit(3),antel(3),i3)
                   E(3) = EXTRA(orbit(3),antel(3),i3)
                   do 35 SK2=abs(SK1-S(3))+1,SK1+S(3)-1,2
                    do 35 LK2=abs(LK1-L(3)),LK1+L(3)
*Limit; keeps the intermediate couplings less than or equal to "N".
                     if (LK2.GT.10) goto 35
                     SK(2) = SK2
                     LK(2) = LK2
                     if (pos.EQ.3) then
                      if (SK(2).GE.minS .AND. SK(2).LE.resS .AND.
     :                    LK(2).GE.minL .AND. LK(2).LE.resL .AND.
     :                    SK(2)+2*LK(2)-1.GE.J2min .AND.
     :                    abs(SK(2)-2*LK(2)-1).LE.J2max) then
                       if (check) then
                        write(*,999) ' ',rad1(1:24)
                        return
                       endif
                       call KOPP1(pos,rad2,S,L,E,dyn)
                       call KOPP2(pos,rad2,SK,LK,dyn)
                       write(fil,999) rad1(1:24)
                       write(fil,999) rad2(1:stopp)
                       cf = cf + 1
                      endif
                     else
                      S(4) = SKVANT(orbit(4),antel(4),i4)
                      L(4) = LKVANT(orbit(4),antel(4),i4)
                      E(4) = EXTRA(orbit(4),antel(4),i4)
                      do 45 SK3=abs(SK2-S(4))+1,SK2+S(4)-1,2
                       do 45 LK3=abs(LK2-L(4)),LK2+L(4)
*Limit; keeps the intermediate couplings less than or equal to "N".
                        if (LK3.GT.10) goto 45
                        SK(3) = SK3
                        LK(3) = LK3
                        if (pos.EQ.4) then
                         if (SK(3).GE.minS .AND. SK(3).LE.resS .AND.
     :                       LK(3).GE.minL .AND. LK(3).LE.resL .AND.
     :                       SK(3)+2*LK(3)-1.GE.J2min .AND.
     :                       abs(SK(3)-2*LK(3)-1).LE.J2max) then
                          if (check) then
                           write(*,999) ' ',rad1(1:32)
                           return
                          endif
                          call KOPP1(pos,rad2,S,L,E,dyn)
                          call KOPP2(pos,rad2,SK,LK,dyn)
                          write(fil,999) rad1(1:32)
                          write(fil,999) rad2(1:stopp)
                          cf = cf + 1
                         endif
                        else
                         S(5) = SKVANT(orbit(5),antel(5),i5)
                         L(5) = LKVANT(orbit(5),antel(5),i5)
                         E(5) = EXTRA(orbit(5),antel(5),i5)
                         do 55 SK4=abs(SK3-S(5))+1,SK3+S(5)-1,2
                          do 55 LK4=abs(LK3-L(5)),LK3+L(5)
*Limit; keeps the intermediate couplings less than or equal to "N".
                           if (LK4.GT.10) goto 55
                           SK(4) = SK4
                           LK(4) = LK4
                           if (pos.EQ.5) then
                            if (SK(4).GE.minS .AND. SK(4).LE.resS .AND.
     :                          LK(4).GE.minL .AND. LK(4).LE.resL .AND.
     :                          SK(4)+2*LK(4)-1.GE.J2min .AND. 
     :                          abs(SK(4)-2*LK(4)-1).LE.J2max) then
                              if (check) then
                              write(*,999) ' ',rad1(1:40)
                              return
                             endif
                             call KOPP1(pos,rad2,S,L,E,dyn)
                             call KOPP2(pos,rad2,SK,LK,dyn)
                             write(fil,999) rad1(1:40)
                             write(fil,999) rad2(1:stopp)
                             cf = cf + 1
                            endif
                           else
                            S(6) = SKVANT(orbit(6),antel(6),i6)
                            L(6) = LKVANT(orbit(6),antel(6),i6)
                            E(6) = EXTRA(orbit(6),antel(6),i6)
                            do 65 SK5=abs(SK4-S(6))+1,SK4+S(6)-1,2
                             do 65 LK5=abs(LK4-L(6)),LK4+L(6)
*Limit; keeps the intermediate couplings less than or equal to "N".
                              if (LK5.GT.10) goto 65
                              SK(5) = SK5
                              LK(5) = LK5
                              if (pos.EQ.6) then
                               if (SK(5).GE.minS .AND. SK(5).LE.resS 
     :                             .AND. LK(5).GE.minL .AND.
     :                             LK(5).LE.resL .AND.
     :                             SK(5)+2*LK(5)-1.GE.J2min .AND.
     :                             abs(SK(5)-2*LK(5)-1).LE.J2max) then
                                if (check) then
                                 write(*,999) ' ',rad1(1:48)
                                 return
                                endif
                                call KOPP1(pos,rad2,S,L,E,dyn)
                                call KOPP2(pos,rad2,SK,LK,dyn)
                                write(fil,999) rad1(1:48)
                                write(fil,999) rad2(1:stopp)
                                cf = cf + 1
                               endif
                              else
                               S(7) = SKVANT(orbit(7),antel(7),i7)
                               L(7) = LKVANT(orbit(7),antel(7),i7)
                               E(7) = EXTRA(orbit(7),antel(7),i7)
                               do 75 SK6=abs(SK5-S(7))+1,SK5+S(7)-1,2
                                do 75 LK6=abs(LK5-L(7)),LK5+L(7)
*Limit; keeps the intermediate couplings less than or equal to "N".
                                 if (LK6.GT.10) goto 75
                                 SK(6) = SK6
                                 LK(6) = LK6
                                 if (pos.EQ.7) then
                                 if (SK(6).GE.minS .AND. SK(6).LE.resS 
     :                               .AND. LK(6).GE.minL .AND.
     :                               LK(6).LE.resL .AND.
     :                               SK(6)+2*LK(6)-1.GE.J2min .AND.
     :                               abs(SK(6)-2*LK(6)-1).LE.J2max) then
                                   if (check) then
                                    write(*,999) ' ',rad1(1:56)
                                    return
                                   endif
                                   call KOPP1(pos,rad2,S,L,E,dyn)
                                   call KOPP2(pos,rad2,SK,LK,dyn)
                                   write(fil,999) rad1(1:56)
                                   write(fil,999) rad2(1:stopp)
                                   cf = cf + 1
                                  endif
                                 else
                                  S(8) = SKVANT(orbit(8),antel(8),i8)
                                  L(8) = LKVANT(orbit(8),antel(8),i8)
                                  E(8) = EXTRA(orbit(8),antel(8),i8)
                                  do 85 SK7=abs(SK6-S(8))+1,SK6+S(8)-1,2
                                   do 85 LK7=abs(LK6-L(8)),LK6+L(8)
*Limit; keeps the intermediate couplings less than or equal to "N".
                                    if (LK7.GT.10) goto 85
                                    SK(7) = SK7
                                    LK(7) = LK7
                                    if (SK(7).GE.minS .AND. 
     :                                  SK(7).LE.resS .AND.
     :                                  LK(7).GE.minL .AND.
     :                                  LK(7).LE.resL .AND.
     :                                  SK(7)+2*LK(7)-1.GE.J2min .AND. 
     :                               abs(SK(7)-2*LK(7)-1).LE.J2max) then
                                     if (check) then
                                      write(*,999) ' ',rad1(1:64)
                                      return
                                     endif
                                     call KOPP1(pos,rad2,S,L,E,dyn)
                                     call KOPP2(pos,rad2,SK,LK,dyn)
                                     write(fil,999) rad1(1:64)
                                     write(fil,999) rad2(1:stopp)
                                     cf = cf + 1
                                    endif
   85                             continue
                                 endif
   75                          continue
                              endif
   65                       continue
                           endif
   55                    continue
                        endif
   45                 continue
                     endif
   35              continue
                  endif
   25           continue
               endif
              else
               if (pos.EQ.1) then
                if (S(1).EQ.resS .AND. L(1).EQ.resL) then
                 if (check) then
                  write(*,999) ' ',rad1(1:8)
                  return
                 endif
                 call KOPP1(pos,rad2,S,L,E,dyn)
                 write(fil,999) rad1(1:8)
                 write(fil,999) rad2(1:stopp)
                 cf = cf + 1
                endif
               else
                S(2) = SKVANT(orbit(2),antel(2),i2)
                L(2) = LKVANT(orbit(2),antel(2),i2)
                E(2) = EXTRA(orbit(2),antel(2),i2)
                SK(pos-1) = resS
                LK(pos-1) = resL
                if (pos.EQ.2) then
                 if (resS.GE.abs(S(1)-S(2))+1 .AND.
     :               resS.LE.S(1)+S(2)-1 .AND.
     :               resL.GE.abs(L(1)-L(2)) .AND.
     :               resL.LE.L(1)+L(2)              ) then
                  if (check) then
                   write(*,999) ' ',rad1(1:16)
                   return
                  endif
                  call KOPP1(pos,rad2,S,L,E,dyn)
                  call KOPP2(pos,rad2,SK,LK,dyn)
                  write(fil,999) rad1(1:16)
                  write(fil,999) rad2(1:stopp)
                  cf = cf + 1
                 endif
                else
                 S(3) = SKVANT(orbit(3),antel(3),i3)
                 L(3) = LKVANT(orbit(3),antel(3),i3)
                 E(3) = EXTRA(orbit(3),antel(3),i3)
                 do 30 SK1=abs(S(1)-S(2))+1,S(1)+S(2)-1,2
                  do 30 LK1=abs(L(1)-L(2)),L(1)+L(2)
*Limit; keeps the intermediate couplings less than or equal to "N".
                   if (LK1 .GT.10) goto 30
                   SK(1) = SK1
                   LK(1) = LK1
                   if (pos.EQ.3) then
                    if (resS.GE.abs(SK1-S(3))+1 .AND.
     :                  resS.LE.SK1+S(3)-1 .AND.
     :                  resL.GE.abs(LK1-L(3)) .AND.
     :                  resL.LE.LK1+L(3)              ) then
                     if (check) then
                      write(*,999) ' ',rad1(1:24)
                      return
                     endif
                     call KOPP1(pos,rad2,S,L,E,dyn)
                     call KOPP2(pos,rad2,SK,LK,dyn)
                     write(fil,999) rad1(1:24)
                     write(fil,999) rad2(1:stopp)
                     cf = cf + 1
                    endif
                   else
                    S(4) = SKVANT(orbit(4),antel(4),i4)
                    L(4) = LKVANT(orbit(4),antel(4),i4)
                    E(4) = EXTRA(orbit(4),antel(4),i4)
                    do 40 SK2=abs(SK1-S(3))+1,SK1+S(3)-1,2
                     do 40 LK2=abs(LK1-L(3)),LK1+L(3)
*Limit; keeps the intermediate couplings less than or equal to "N".
                      if (LK2 .GT.10) goto 40
                      SK(2) = SK2
                      LK(2) = LK2
                      if (pos.EQ.4) then
                       if (resS.GE.abs(SK2-S(4))+1 .AND.
     :                     resS.LE.SK2+S(4)-1 .AND.
     :                     resL.GE.abs(LK2-L(4)) .AND.
     :                     resL.LE.LK2+L(4)              ) then
                        if (check) then
                         write(*,999) ' ',rad1(1:32)
                         return
                        endif
                        call KOPP1(pos,rad2,S,L,E,dyn)
                        call KOPP2(pos,rad2,SK,LK,dyn)
                        write(fil,999) rad1(1:32)
                        write(fil,999) rad2(1:stopp)
                        cf = cf + 1
                       endif
                      else
                       S(5) = SKVANT(orbit(5),antel(5),i5)
                       L(5) = LKVANT(orbit(5),antel(5),i5)
                       E(5) = EXTRA(orbit(5),antel(5),i5)
                       do 50 SK3=abs(SK2-S(4))+1,SK2+S(4)-1,2
                        do 50 LK3=abs(LK2-L(4)),LK2+L(4)
*Limit; keeps the intermediate couplings less than or equal to "N".
                         if (LK3 .GT.10) goto 50
                         SK(3) = SK3
                         LK(3) = LK3
                         if(pos.EQ.5) then
                          if (resS.GE.abs(SK3-S(5))+1 .AND.
     :                        resS.LE.SK3+S(5)-1 .AND.
     :                        resL.GE.abs(LK3-L(5)) .AND.
     :                        resL.LE.LK3+L(5)              ) then
                           if (check) then
                            write(*,999) ' ',rad1(1:40)
                            return
                           endif
                           call KOPP1(pos,rad2,S,L,E,dyn)
                           call KOPP2(pos,rad2,SK,LK,dyn)
                           write(fil,999) rad1(1:40)
                           write(fil,999) rad2(1:stopp)
                           cf = cf + 1
                          endif
                         else
                          S(6) = SKVANT(orbit(6),antel(6),i6)
                          L(6) = LKVANT(orbit(6),antel(6),i6)
                          E(6) = EXTRA(orbit(6),antel(6),i6)
                          do 60 SK4=abs(SK3-S(5))+1,SK3+S(5)-1,2
                           do 60 LK4=abs(LK3-L(5)),LK3+L(5)
*Limit; keeps the intermediate couplings less than or equal to "N".
                            if (LK4 .GT.10) goto 60
                            SK(4) = SK4
                            LK(4) = LK4
                            if(pos.EQ.6) then
                             if (resS.GE.abs(SK4-S(6))+1 .AND.
     :                           resS.LE.SK4+S(6)-1 .AND.
     :                           resL.GE.abs(LK4-L(6)) .AND.
     :                           resL.LE.LK4+L(6)              ) then
                              if (check) then
                               write(*,999) ' ',rad1(1:48)
                               return
                              endif
                              call KOPP1(pos,rad2,S,L,E,dyn)
                              call KOPP2(pos,rad2,SK,LK,dyn)
                              write(fil,999) rad1(1:48)
                              write(fil,999) rad2(1:stopp)
                              cf = cf + 1
                             endif
                            else
                             S(7) = SKVANT(orbit(7),antel(7),i7)
                             L(7) = LKVANT(orbit(7),antel(7),i7)
                             E(7) = EXTRA(orbit(7),antel(7),i7)
                             do 70 SK5=abs(SK4-S(6))+1,SK4+S(6)-1,2
                              do 70 LK5=abs(LK4-L(6)),LK4+L(6)
*Limit; keeps the intermediate couplings less than or equal to "N".
                               if (LK5 .GT.10) goto 70
                               SK(5) = SK5
                               LK(5) = LK5
                               if(pos.EQ.7) then
                                if (resS.GE.abs(SK5-S(7))+1 .AND.
     :                              resS.LE.SK5+S(7)-1 .AND.
     :                              resL.GE.abs(LK5-L(7)) .AND.
     :                              resL.LE.LK5+L(7)              ) then
                                 if (check) then
                                  write(*,999) ' ',rad1(1:56)
                                  return
                                 endif
                                 call KOPP1(pos,rad2,S,L,E,dyn)
                                 call KOPP2(pos,rad2,SK,LK,dyn)
                                 write(fil,999) rad1(1:56)
                                 write(fil,999) rad2(1:stopp)
                                 cf = cf + 1
                                endif
                               else
                                S(8) = SKVANT(orbit(8),antel(8),i8)
                                L(8) = LKVANT(orbit(8),antel(8),i8)
                                E(8) = EXTRA(orbit(8),antel(8),i8)
                                do 80 SK6=abs(SK5-S(7))+1,SK5+S(7)-1,2
                                 do 80 LK6=abs(LK5-L(7)),LK5+L(7)
*Limit; keeps the intermediate couplings less than or equal to "N".
                                  if (LK6 .GT.10) goto 80
                                  if ( resS.GE.abs(SK6-S(8))+1 .AND.
     :                                 resS.LE.SK6+S(8)-1 .AND.
     :                                 resL.GE.abs(LK6-L(8)) .AND.
     :                                 resL.LE.LK6+L(8)           ) then
                                   SK(6) = SK6
                                   LK(6) = LK6
                                   if (check) then
                                    write(*,999) ' ',rad1(1:64)
                                    return
                                   endif
                                   call KOPP1(pos,rad2,S,L,E,dyn)
                                   call KOPP2(pos,rad2,SK,LK,dyn)
                                   write(fil,999) rad1(1:64)
                                   write(fil,999) rad2(1:stopp)
                                   cf = cf + 1
                                  endif
   80                           continue
                               endif
   70                        continue
                            endif
   60                     continue
                         endif
   50                  continue
                      endif
   40               continue
                   endif
   30            continue
                endif
               endif
              endif
   90 continue
  999 format(2A)
      check = .FALSE.
      return
      end
*     last edited July 2, 1993
      subroutine Koda(rad1,rad2,rank,antal,maskal,nmax,lmax,opt)
      integer    noofco
Ctc   parameter (noofco=200000)
      parameter (noofco=2100000)
      character  rad1(0:noofco)*64,rad2(0:noofco)*60,L(0:20),orb(0:10)
      integer    rank(0:noofco,0:10,15),antal,maskal,nmax,lmax
      logical    opt
      integer    i,j,k,k1,l1,m,tal,grans,nr

      data (L(i),i=0,20)/'S','P','D','F','G','H','I','K','L','M','N'
     :                         ,'O','Q','R','T','U','V','W','X','Y','Z'/
      data (orb(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/

      do 10 j=1,15
         if (opt) then
            do 5 i=1,4
    5          rank(antal,i,j) = 999
         else
            do 7 i=0,10
    7          rank(antal,i,j) = 0
         endif
   10 continue 
      maskal = 0
      nmax   = 1
      lmax   = 0
      do 40 i=0,antal-1
         grans = 0
         do 20 j=1,15
            k  = (j-1)*4+2
            k1 = k+1
            tal = ichar(rad2(i)(k:k)) - ichar('0')
            if (tal.GT.0 .AND. tal.LT.10) then
               grans  = j
               maskal = max(j,maskal)
               if (opt) then
                  rank(i,1,j) = tal
                  rank(i,2,j) =-2
                  do 15 m=0,20
   15                if (rad2(i)(k1:k1) .EQ. L(m)) rank(i,2,j) = m
                  if (rank(i,2,j).EQ.-2) then
                     write(*,*) 
     :                   'Wrong notation in configuration number',i
                     write(*,*) 'Stopp 1' 
                     stop
                  endif
               endif
            elseif (opt) then
               rank(i,1,j) = 0
               rank(i,2,j) =-1
            endif
   20    continue
         do 30 j=1,(grans+1)/2
            k   = (j-1)*8+3
            k1  = k+1
            tal = ichar(rad1(i)(k:k)) - ichar('0')
            if (tal.GT.0 .AND. tal.LT.16) then
               nmax = max(nmax,tal)
               if (opt) rank(i,4,j) = tal 
               l1 =-1
               do 25 m=0,10
   25             if (rad1(i)(k1:k1) .EQ. orb(m)) l1 = m 
               if (opt) then
                  rank(i,3,j) = l1
               else
                  lmax = max(lmax,m)
               endif
               if (l1.EQ.-1 .OR. l1.GT.min(tal-1,10)) then 
                  write(*,*) 'Wrong notation in configuration number',i
                  write(*,*) 'Stopp 2' 
                  stop
               endif
            else
               write(*,*) 'Wrong notation in configuration number',i
               write(*,*) 'Stopp 3' 
               stop
            endif
            if (.NOT. opt) then
               k = (j-1)*8+6
               if (rad1(i)(k:k).NE.' ') then
                  nr = 10*(ichar(rad1(i)(k:k)) - ichar('0'))
               else
                  nr = 0
               endif
               k1 = k+1
               nr = nr + ichar(rad1(i)(k1:k1)) - ichar('0')
               rank(i,l1,tal) = nr
            endif
   30    continue
         if (opt) then
            do 35 j=(grans+3)/2,15
               rank(i,3,j) =-1
   35          rank(i,4,j) = 0
         endif
   40 continue

*     write(*,*) 'max number of columns',maskal
      return
      end
*     last edited April 1,1993
      logical function Kolla(rank,nr,maskal)
      integer noofco
Ctc   parameter (noofco=200000)
      parameter (noofco=2100000)
      integer rank(0:noofco,0:10,15),maskal
      integer i,j,nr

      Kolla = .FALSE.
      do 10 i=maskal,(maskal+3)/2,-1
         do 10 j=2,1,-1
   10       if (rank(nr,j,i) .NE. rank(0,j,i)) return
      do 20 i=(maskal+1)/2,1,-1
         do 20 j=3,1,-1
   20       if (rank(nr,j,i) .NE. rank(0,j,i)) return
      Kolla = .TRUE.
      return
      end
*     last edited April 1, 1993
      subroutine Kopp1(pos,rad2,S,L,E,dyn)
      integer pos,i,j,S(8),L(8),E(8)
      character rad2*72,L2(0:20)
      logical dyn
      data (L2(i),i=0,20)/'S','P','D','F','G','H','I','K','L','M','N'
     :                         ,'O','Q','R','T','U','V','W','X','Y','Z'/
      do 10 i=1,pos
         if (dyn) then
            j = 4*i
            rad2(j-3:j-3) = ' '
         else
            j = 8*i
            rad2(j-7:j-3) = '    '
         endif
         rad2(j-2:j-2) = char(48+S(i))
         rad2(j-1:j-1) = L2(L(i))
   10    rad2(j:j)     = char(48+E(i))
      return
      end
*     last edited April 1, 1993 
      subroutine Kopp2(pos,rad2,S,L,dyn)
      integer pos,i,j,S(8),L(8)
      character rad2*72,L2(0:20)
      logical dyn
      data (L2(i),i=0,20)/'S','P','D','F','G','H','I','K','L','M','N'
     :                         ,'O','Q','R','T','U','V','W','X','Y','Z'/
      do 10 i=1,pos-1
         if (dyn) then
            j             = 4*(pos+i)
            rad2(j-3:j-3) = ' '
            rad2(j:j)     = ' '
         else
            j             = 8*(pos+i)
            rad2(j-7:j-3) = '    '
            rad2(j:j)     = ' '
         endif
         rad2(j-2:j-2) = char(48 + S(i))
   10    rad2(j-1:j-1) = L2(L(i))
      return
      end
*     last edited February 21, 1993
      subroutine Lasa(rad1,rad2,antal)
      integer noofco
Ctc   parameter (noofco=200000)
      parameter (noofco=2100000)
      character rad1(0:noofco)*64,rad2(0:noofco)*60
      integer antal

      open(unit=7,file='clist.inp',status='old')
      open(unit=13,file='clist.out',status='unknown')
      antal = 0
      read(7,1000,end=666) rad1(antal)
      write(13,1000)
      read(7,1000,end=666) rad1(antal)
      write(13,1000) rad1(antal)
    7 read(7,1000,end=666) rad1(antal)
      if (rad1(antal)(5:5).NE.'(') goto 7
      read(7,1000,end=99) rad2(antal)
      antal = antal+1
      if (antal.EQ.noofco) then
         write(*,*) 'To many configurations in list!'
         close(7)
         stop
      endif
      goto 7
   99 write(*,*)
     f 'Improper truncation of inputfile at configuration number',antal
         close(7)
      stop
  666 if (antal.EQ.0) then
         write(*,*) 'No configurations in input list!'
         close(7)
         stop
      endif
      return
 1000 format(A)
      end
*     last edited February 21, 1993 
      subroutine Lasa1(fil,rad,pop,skal,slut,dyn)
      character rad*64
      integer   fil,pop(15,0:10),skal
      logical   slut,dyn
      if (.NOT.slut) then
         read(fil,999,end=10) rad
         call Reada(rad,pop,skal,slut,dyn)
         return
      endif
   10 slut = .TRUE.
  999 format(A)
      return
      end
*     last edited February 21, 1993 
      subroutine Lasa2(fil,rad,stopp,slut)
      character rad*72
      integer   fil,stopp
      logical   slut
      if (.NOT.slut) then
         read(fil,999,end=10) rad(1:stopp)
         return
      endif
   10 slut = .TRUE.
  999 format(A)
      return
      end
*     last edited February 21, 1993 
      logical function Lika(pop0,pop1)
      integer i,j,pop0(15,0:10),pop1(15,0:10)
      logical dum
      dum = .TRUE.
      do 10 i=1,15
         do 10 j=0,min(10,i-1)
   10       dum = dum .AND. pop0(i,j).EQ.pop1(i,j)
      Lika = dum
      return
      end
*     last edited February 21, 1993
      subroutine Lockad(closed,slut,expand)
      logical closed(15,0:10),slut,expand
      character rad*80,orb(0:10)
      integer i,j,n,l
      data (orb(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/
      if (expand) then
         read(8,*,end=40)
         read(8,100,end=40) rad
      else
         read(7,*,end=40)
         read(7,100,end=40) rad
      endif
      do 10 n=1,15
         do 10 l=1,min(10,n-1)
   10       closed(n,l) = .FALSE.
      do 30 i=0,19
         j = i*4
         n = ichar(rad(j+3:j+3)) - ichar('0')
         if (n.GE.1 .AND. n.LE.15) then
            do 20 l=0,min(10,n-1)
               if (rad(j+4:j+4).EQ.orb(l)) then
                  closed(n,l) = .TRUE.
                  goto 30
               endif
   20       continue
         else
            return
         endif
   30 continue
      return
   40 slut = .TRUE.
      return
  100 format(A)
      end
*     ------------------------------------------------------------------
*     LSGEN  -- program to generate LS configurations states
*
*     Written by: Lennart Sturesson
*                 Department of Physics
*                 Lund Institute of Technology
*                 P.O. Box 118, S-221 00  Lund, Sweden
*
*     Received   March, 1992
*     1:st modified version: including minimum option
*     2:nd modified version: including customized coupling order
*     3:rd modified version: including the Breit-option
*     4:th modified version: including virtual sets
*     5:th modified version: including limit of population and double excitation
*     6:th modified version: including optimzed sorting
*     7:th modified version: including standard sorting
*     8:th modified version: including simple expansion of existing list
*     9:th modified version: including advanced expansion of existing list
*                                      
*     Last edited July 2, 1993
*     ------------------------------------------------------------------
*
      program Lsgen15
      integer org(1:15,0:10),varmax,resS,resL,skal,anel,par,
     :        low(1:15,0:10),posn(110),posl(110),nmax,cfmax,
     :        J2min,J2max,minS,minL,virmax,lim(15),nold
      logical lock(1:15,0:10),closed(1:15,0:10),second,slut,dyn,breit,
     :        virtu(1:15,0:10),vir,dubbel(1:15,0:10),extra,advexp
      character X*1,Y*1
      open(unit=31,file='clist.log',status='unknown')
   10 write(*,200) 'New list, add to existing list, expand existing ',
     :             'list, optimized sorting,',
     :             ' restored order or quit? (*/a/e/s/r/q)'
      read(*,100) X
      write(31,200) ' Option : ',X
      advexp = X.EQ.'e' .OR. X.EQ.'E'
      if (X.EQ.'h' .OR. X.EQ.'H' .OR. X.EQ.'?') then
         write(*,*) 'Sorry! No online help is implemented!'
         goto 10
      elseif (X.EQ.'q' .OR. X.EQ.'Q') then
         stop
      endif
      if (X.EQ.'s' .OR. X.EQ.'S') then
         call Sort(.true.)
         write(*,200) 'The generated file is called clist.out. '
         write(*,200) 'NB! This file is not possible to use as ',
     :                'inputfile for the adding routine.'
         write(*,200) 'The option "r" will restore to the ',
     :                'standard order, enableing further adding.'
         stop
      elseif (X.EQ.'r' .OR. X.EQ.'R') then
         call Sort(.false.)
         write(*,200) 'The generated file is called clist.out. '
         write(*,200) 'NB! This file will work with ',
     :                'the adding routine.'
         stop
      endif
*     write(*,200) 'Dynamic or static version? (d/*)'
*     read(*,100) Y
*     dyn = Y.EQ.'d' .OR. Y.EQ.'D'
      dyn = .TRUE.
      if (advexp) then
         write(*,*)
         write(*,*) 'This option is only running in the MCHF-mode!'
         write(*,*)
         breit = .FALSE.
      else
         write(*,200) 'Breit or MCHF? (B/*)'
         read(*,100) Y
         breit = Y.EQ.'b' .OR. Y.EQ.'B'
         if (breit) then
            write(31,200) ' Breit'
	 else
            write(31,200) ' MCHF'
	 endif
      endif
      call Reffa(posn,posl,dyn)
      if (X.EQ.'a' .OR. X.EQ.'A' .OR. advexp) then
         call Adder(closed,slut,resS,resL,anel,par,dyn,advexp)
         if (slut) then
            if (dyn) then
               write(*,200)
               write(*,200) 'The clist.inp-file is not readable! ',
     :                      'Run the "r"-option if the file ',
     :                      'is a copy of a cfg-file!' 
            else
               write(*,200)
               write(*,200) 'The cfg.old-file is not readable!'
            endif
            stop
         endif
         if (.NOT. advexp) then
            write(*,200)
            write(*,200) 'Use new reference set or make simple ',
     :                   'expansion of the active sets in input-file? ',
     :                   '(*/e)'
            read(*,100) Y
            extra = Y.EQ.'e' .OR. Y.EQ.'E'
            if (extra) then
	       write(31,200) 
     :                ' Expansion of the active sets in input-file.'
            else
               write(31,200) ' New reference set.'
	    endif 
            if (extra) then
               call Scan(nold)
               if (nold.LT.10) then
                  write(*,300) 'The highest identified n-number is',
     :                         nold,'.'
               else
                  write(*,400) 'The highest identified n-number is',
     :                         nold,'.' 
               endif
               if (nold.GT.14) then
                  write(*,200) 'The identified n-number is to high.', 
     :                         ' - Highest possible n-number is 14!'
                  stop
               endif
               call Expand(nold,dyn)
            endif
         else
            extra = .FALSE.
         endif
         second = .TRUE.
      else
         call Matain(org,lock,closed,resS,resL,varmax,virmax,skal,nmax,
     :               anel,par,dyn,low,breit,J2min,J2max,minS,minL,virtu,
     :               vir,lim,dubbel)
         call Twolines(org,closed,.TRUE.)
         call Blanda(org,varmax,virmax,lock,resS,resL,skal,nmax,dyn,
     :               .TRUE.,low,posn,posl,breit,J2min,J2max,minS,minL,
     :               virtu,vir,lim,dubbel)
         second = .FALSE.
         extra  = .FALSE.
      endif
      if (.NOT.extra) then
         if (advexp) then
            call Matcin(lock,closed,varmax,cfmax,nmax)
            call Blandc(varmax,cfmax,lock,resS,resL,nmax,dyn,posn,
     :                  posl)
            second = .FALSE.
         else
            call Matbin(org,lock,closed,varmax,virmax,skal,second,anel,
     :            par,dyn,low,breit,J2min,J2max,resS,resL,minS,minL,
     :            virtu,vir,nmax,lim,dubbel)
         endif
      endif
      if (second) then
         if (.NOT.extra) then
            call Twolines(org,closed,.FALSE.)
            call Blanda(org,varmax,virmax,lock,resS,resL,skal,nmax,dyn,
     :                  .FALSE.,low,posn,posl,breit,J2min,J2max,minS,
     :                  minL,virtu,vir,lim,dubbel)
         endif
         call Merge(.FALSE.,dyn)
         if (dyn) then
            write(*,200) 'The merged file is called clist.out.'
         else
            write(*,200) 'The merged file is called cfg.inp.'
         endif
      else
         call Merge(.TRUE.,dyn)
         if (dyn) then
            write(*,200) 'The generated file is called clist.out.'
         else
            write(*,200) 'The generated file is called cfg.inp.'
         endif
      endif
      close(31)
      stop
  100 format(A)
  200 format(' ',10A)
  300 format(' ',A,I2,A)
  400 format(' ',A,I3,A)
      end
*     last edited July 2, 1993
      subroutine Matain(org,lock,closed,resS,resL,varmax,virmax,skal,
     :                  nmax,anel,par,dyn,low,breit,J2min,J2max,minS,
     :                  minL,virtu,vir,lim,dubbel)
      integer org(1:15,0:10),anela,low(1:15,0:10),anelb,J2min,J2max,minS
      integer varmax,resS,resL,i,j,skal,nmax,lmax,em,anel,par,nenter
      integer minL,virmax,lim(15),block,mshell,enn
      character X,orb(0:10),L(0:20),Y*3
      logical lock(1:15,0:10),closed(1:15,0:10),dyn,all,breit,
     :        virtu(1:15,0:10),vir,dubbel(1:15,0:10),lima
      data (L(i),i=0,20)/'S','P','D','F','G','H','I','K','L','M','N'
     :                         ,'O','Q','R','T','U','V','W','X','Y','Z'/
      data (orb(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/
      if (dyn) then
*  50    write(*,200) 'Maximum number of open subshells? (1..8)'
*        read(*,*,err=50) skal
*        skal = max(skal,1)
*        skal = min(skal,8)
         skal = 8
      else
   51    write(*,200) 'Maximum number of open subshells? (1..5)'
         read(*,*,err=51) skal
         skal = max(skal,1)
         skal = min(skal,5)
      endif
   60 write(*,200) 'Highest principal quantum number, n? (1..15)'
      read(*,*,err=60) nmax
      nmax = max(nmax,1)
      nmax = min(nmax,15)
      write(31,*) nmax,' Highest principal quantum number.'
   70 write(*,300) 'Highest orbital angular momentum, l? (s..',
     :                                        orb(min(10,nmax-1)),')'
      read(*,1000) X
      lmax = -1
      do 71 i=0,min(10,nmax-1)
   71    if (X.EQ.orb(i)) lmax=i
      if (lmax.EQ.-1) goto 70
      write(31,*) lmax,' Highest orbital angular momentum.' 
      write(*,200) 'Are all these nl-subshells active? (n/*)'
      read(*,1000) X
      all   = .NOT.(X.EQ.'n' .OR. X.EQ.'N')
      write(31,*) all,' all subshells active.'
      do 72 i=1,15
   72    lim(i) = 0
      if (nmax.GE.2) then
         write(*,200) 'Limitations on population of n-subshells? (y/*)'
         read(*,1000) X
         lima = X.EQ.'y' .OR. X.EQ.'Y'
         write(31,*) lima,' limitations on population of n-subshells.'
         if (lima) then
            mshell = 0
            do 85 i=1,nmax-1
               mshell = mshell + 2*i*i
   83          continue
               if (i.EQ.1) then
                  write(*,200)
     :                    'Minimum number of electrons with n=1? (0..2)'
               elseif (i.LT.10) then
                  if (mshell.LT.100) then
                     write(*,208) 'Minimum number of electrons with n<='
     :                       ,i,'? (0..',mshell,')'
                  else
                     write(*,208) 'Minimum number of electrons with n<='
     :                       ,i,'? (0..)'
                  endif
               else
                  write(*,202) 'Minimum number of electrons with n<=',i,
     :                                         '? (0..)'
               endif
               read(*,*,err=83) lim(i)
               if (lim(i).GT.mshell) lim(i) = mshell
               write(31,*) lim(i),
     :                     ' is minimum number of electrons with n =',i
   85       continue
         endif
      endif
   90 continue
      if (nmax.LT.10) then
         write(*,200)
     :     'Highest n-number in reference configuration? (1..', nmax,')'
      else
         write(*,202)
     :     'Highest n-number in reference configuration? (1..', nmax,')'
      endif
      read(*,*,err=90) nenter
      nenter = max(nenter,1)
      nenter = min(nenter,nmax)
      write(31,*) nenter,' highest n-number.'
      anela  = 0
      anelb  = 0
      anel   = 0
      par    = 0
      block  = 0
      vir    = .FALSE.
      do 160 i=1,15
         do 150 j=0,min(10,i-1)
            low(i,j)    = 0
            virtu(i,j)  = .FALSE.
            dubbel(i,j) = .FALSE.
            if (nmax.GE.i .AND. lmax.GE.j) then
               if (nenter.GE.i) then
                  if (j.EQ.3) then                                   
*Limit; special case for f-orbitals
                     em = 2
   95                continue
                     if (i.LE.9) then
                        write(*,200) 'Number of electrons in ',i,
     :                                  'f? (0,1,2 or 14)'
                     else
                        write(*,202) 'Number of electrons in ',i,
     :                                  'f? (0,1,2 or 14)'
                     endif
                     read(*,*,err=95) org(i,j)
                     if ( org(i,j).LT.0 .OR. org(i,j).GT.14 .OR.
     :                    (org(i,j).GT.2 .AND. org(i,j).LT.14) ) goto 95
                  else
                     em = 2 + 4*j
                     if (j.GE.3) em = 2                              
*Limit; high angular momentum orbits, except f
                     if (em.LT.10) then
  100                   continue
                        if (i.LE.9) then
                           write(*,200) 'Number of electrons in ',i,
     :                                          orb(j),'? (0..',em,')'
                        else
                           write(*,202) 'Number of electrons in ',i,
     :                                          orb(j),'? (0..',em,')'
                        endif
                        read(*,*,err=100) org(i,j)
                        if (org(i,j).LT.0 .OR. org(i,j).GT.em) 
     :                        goto 100
                     else
  101                   continue
                        if (i.LT.10) then
                           write(*,201) 'Number of electrons in ',i,
     :                                           orb(j),'? (0..',em,')'
                        else
                           write(*,203) 'Number of electrons in ',i,
     :                                           orb(j),'? (0..',em,')'
                        endif
                        read(*,*,err=101) org(i,j)
                        if (org(i,j).LT.0 .OR. org(i,j).GT.em) 
     :                        goto 101
                     endif
                  endif
                  write(31,*) org(i,j),' number of electrons in',i,
     :                        orb(j)
                  anel = anel + org(i,j)
                  par  = mod(par+j*org(i,j),2)
                  if (org(i,j).EQ.14) then                           
*Limit; full f subshell
                     lock(i,j)   = .TRUE.
                     closed(i,j) = .TRUE.
                     block       = block + 14
                  elseif (all) then
                     lock(i,j)   = .FALSE.
                     closed(i,j) = .FALSE.
                  else
                     if (org(i,j).EQ.em .AND. j.LE.2) then           
*Limit; only s, p, and d shells may be chosen
                        write(*,201)
     :                'Closed, inactive, active or minimum? (c/i/*/0..',
     :                               org(i,j)-1,')'
                        read(*,1000) X
                        write(31,*) X,' closed, inactive, etc...'
                        closed(i,j) = X.EQ.'c' .OR. X.EQ.'C'
                        lock(i,j)   = X.EQ.'i' .OR. X.EQ.'I'
     :                                         .OR. closed(i,j)
                        if (closed(i,j)) block = block + em
                     else
                        if (org(i,j).GT.1) then
                           write(*,201)
     :           'Inactive, active or minimum? (i/*/0..',org(i,j)-1,')'
                        elseif (org(i,j).EQ.1) then 
                           write(*,201) 'Inactive or active? (i/*)'
                        else
                           write(*,400) 'Inactive, active, double ',
     :                            'excited or virtual? (i/*/d/v)'
                        endif
                        read(*,1000) X
                        write(31,*) X,' inactive, active, etc...'
                        if (org(i,j).EQ.0) then
                           dubbel(i,j) = X.EQ.'d' .OR. X.EQ.'D'
                           virtu(i,j)  = X.EQ.'v' .OR. X.EQ.'V'
                           vir         = vir .OR. virtu(i,j)
                        endif
                        lock(i,j)   = X.EQ.'i' .OR. X.EQ.'I' .OR.
     :                                virtu(i,j)
                        closed(i,j) = .FALSE.
                     endif
                     if (X.GE.'0' .AND. X.LE.'9') then
                        if (org(i,j).GT.0) 
     :                     low(i,j) = min(org(i,j),ICHAR(X)-ICHAR('0'))
                     endif
                  endif
                  if (.NOT. lock(i,j)) anela = anela + org(i,j)
               elseif (all) then
                  org(i,j) = 0
                  lock(i,j) = .FALSE.
                  closed(i,j) = .FALSE.
               else
                  org(i,j) = 0
                  closed(i,j) = .FALSE.
                  if (i.LT.10) then
                     write(*,204) i,orb(j),' inactive, active, ',
     :                           'doubled excited or virtual? (i/*/d/v)'
                  else
                     write(*,205) i,orb(j),' inactive, active, ',
     :                           'doubled excited or virtual? (i/*/d/v)'
                  endif
                  read(*,1000) X
                  write(31,*) X,i,orb(j),' inactive, active, etc...'
                  dubbel(i,j) = X.EQ.'d' .OR. X.EQ.'D'
                  virtu(i,j)  = X.EQ.'v' .OR. X.EQ.'V'
                  lock(i,j)   = X.EQ.'i' .OR. X.EQ.'I' .OR. virtu(i,j)
                  vir         = vir .OR. virtu(i,j)      
               endif
            else
               org(i,j)  = 0
               lock(i,j) = .TRUE.
               closed(i,j) = .FALSE.
            endif
            anelb = anelb + low(i,j)
  150    continue
         lim(i) = lim(i) - block
         if (lim(i).LT.0) lim(i) = 0
  160 continue
      if (breit) then
 1154    write(*,200) 'Maximum 2*J-value? (0..)'
         read(*,*,err=1154) J2max
         write(31,*) J2max,' maximum 2*J-value.'
         if (J2max.LT.0) J2max = 0
         if (J2max.EQ.0) then
            J2min = 0
         else
 1155       write(*,200) 'Minimum 2*J-value? (0..',J2max,')'
            read(*,*,err=1155) J2min
            write(31,*) J2min,' minimum 2*J-value.'      
         endif
 1156    write(*,200) 'Maximum (2*S+1)-value? (1..9)'
         read(*,*,err=1156) resS
         if (resS.GT.9) resS = 9
         write(31,*) resS,' maximum (2*S+1)-value.'
         if (resS.LE.1) then
            resS = 1
            minS = 1
         else
 1157       write(*,200) 'Minimum (2*S+1)-value? (1..',resS,')' 
            read(*,*,err=1157) minS
            if (minS.GT.resS) minS = resS
            if (minS.LT.1)    minS = 1
            write(31,*) minS,' minimum (2*S+1)-value.'
         endif
 1158    write(*,200) 'Maximum resulting angular momentum? (S..N/N=*)'
         read(*,1000,ERR=1158) X
         resL = 10
*        do 1161 i=0,20
         do 1161 i=0,10 
*Limit; NONH does not accept higher resultant angular momentum than "N".
 1161       if (X.EQ.L(i)) resL=i 
*        if (resL.EQ.-1) then
*           write(*,200) 'S, P, D, F, G, H, I, K, L, M or N are valid!'
*           goto 1158
*        endif
         write(31,*) L(resL),' = maximium resulting angular momentum'
         minL = 0
         if (resL.GT.0) then
 1159       write(*,300) 'Minimum resulting angular momentum? (S..',
     :                   L(resL),'/S=*)'
            read(*,1000,ERR=1159) X
            do 1162 i=0,resL
 1162          if (X.EQ.L(i)) minL=i 
         write(31,*) L(minL),' = minimium resulting angular momentum'
         endif
      else
 1100    write(*,200) 'Resulting term? (1S, 3P, etc.)'
         read(*,2000,ERR=1100) resS,X
         if (X.EQ.' ') then
            write(*,*) 'What orbital is to be inspected?',
     :                 ' (1s, 14p, etc..)'
            read(*,1000,ERR=1100) Y
            if (Y(2:2).GE.'0' .AND. Y(2:2).LE.'5') then
               enn = (ichar(Y(1:1))-ichar('0'))*10 + ichar(Y(2:2)) -
     :               ichar('0')
               X   = Y(3:3)
            else
               enn = ichar(Y(1:1))-ichar('0')
               X   = Y(2:2)
            endif
            if (enn.LT.1 .OR. enn.GT.15) goto 1100
            resL = -1
            do 1125 i=0,min(enn-1,10)
 1125          if (X.EQ.orb(i)) resL = i
            if(resL.EQ.-1) goto 1100
            if (enn.LT.10) then
               if (org(enn,resL).EQ.0) then
                  write(*,200) 'There is no electron in ',enn,X,'.'
               elseif (org(enn,resL).EQ.1) then
                  write(*,200) 'There is one electron in',enn,X,'.'
               elseif (org(enn,resL).LT.10) then
                  write(*,200) 'Number of electrons in ',enn,X,'are ',
     :                        org(enn,resL),'.'
               else
                  write(*,201) 'Number of electrons in ',enn,X,'are ',
     :                        org(enn,resL),'.'
               endif
            else
               if (org(enn,resL).EQ.0) then
                  write(*,202) 'There is no electron in ',enn,X,'.'
               elseif (org(enn,resL).EQ.1) then
                  write(*,202) 'There is one electron in',enn,X,'.'
               elseif (org(enn,resL).LT.10) then
                  write(*,202) 'Number of electrons in ',enn,X,'are ',
     :                        org(enn,resL),'.'
               else
                  write(*,203) 'Number of electrons in ',enn,X,'are ',
     :                        org(enn,resL),'.'
               endif
            endif
            if (closed(enn,resL)) then
               write(*,200) 'The shell is closed.'
            elseif (virtu(enn,resL)) then
               write(*,200) 'The orbital is virtual.'
            elseif (lock(enn,resL)) then
               write(*,200) 'The orbital is inactive.'
            elseif (dubbel(enn,resL)) then
               write(*,200) 'The excitation to the orbital is',
     :                    ' restricted to be doubble.'
            else
               write(*,200) 'The orbital is active.'
            endif
            
            goto 1100
         endif
         if (mod(resS-1,2).NE.mod(anel,2)) then
            write(*,200)
     :    'The spin is not consistent with the number of electrons!'
            goto 1100
         endif
         resL = -1
*        do 1150 i=0,20
         do 1150 i=0,10 
*Limit; NONH does not accept higher resultant angular momentum than "N".
 1150       if (X.EQ.L(i)) resL = i 
         if(resL.EQ.-1) then
            write(*,200) 'S, P, D, F, G, H, I, K, L, M or N are valid!'
            goto 1100
         endif
         write(31,*) resS,X,' is the resulting term.'
      endif
      anelb = anela - anelb
 1200 continue
      if (anelb.LT.10) then
         write(*,200) 'Number of excitations = ? (0..',anelb,')'
         read(*,*,err=1200) varmax
         if (vir) then
 1210       write(*,400) 'Number of excitations from the active set',
     :                      ' to the virtual set? (0..',anelb,')'
            read(*,*,err=1210) virmax
         endif
      else
         write(*,202) 'Number of excitations = ? (0..',anelb,')'
         read(*,*,err=1200) varmax
         if (vir) then
 1220       write(*,402) 'Number of excitations from the active set',
     :                      ' to the virtual set? (0..',anelb,')'
            read(*,*,err=1220) virmax
         endif
      endif
      write(31,*) varmax,' number of excitations.'
      if (vir) write(31,*) virmax,
     :                  ' number of excitations to virtual set.'
  200 format(' ',A,I1,A,A,I1,A)
  201 format(' ',A,I1,A,A,I2,A)
  202 format(' ',A,I2,A,A,I1,A)
  203 format(' ',A,I2,A,A,I2,A)
  204 format(' ',I1,3A)
  205 format(' ',I2,3A)
  206 format(' ',I1,A,A,I2,A)
  207 format(' ',I2,A,A,I2,A)
  208 format(' ',A,I1,A,I2,A)
  300 format(' ',3A)
  400 format(' ',2A,I1,A)
  402 format(' ',2A,I2,A)
 1000 format(3A)
 2000 format(I1,2A)
 3000 format(A,I2,2A)
      close(33)
      return
      end
*     last edited June 20, 1993 
      subroutine Matbin(org,lock,closed,varmax,virmax,skal,second,anel0,
     :                  par0,dyn,low,breit,J2min,J2max,maxS,maxL,minS,
     :                  minL,virtu,vir,nmax,lim,dubbel)
      integer org(1:15,0:10),anel,par,anel0,par0,low(1:15,0:10),J2max
      integer varmax,i,j,skal,nmax,lmax,em,nenter,anela,anelb,J2min
      integer maxS,maxL,minS,minL,virmax,lim(15),block,mshell
      character X*1,orb(0:10),L(0:20)
      logical lock(1:15,0:10),closed(1:15,0:10),second,dyn,all,breit,
     :        virtu(1:15,0:10),vir,dubbel(1:15,0:10),lima
      data (L(i),i=0,20)/'S','P','D','F','G','H','I','K','L','M','N'
     :                         ,'O','Q','R','T','U','V','W','X','Y','Z'/
      data (orb(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/
   40 if (.NOT. second) then
         write(*,200) 'Generate a second list? (y/*)'
         read(*,1000) X
         second = X.EQ.'y' .OR. X.EQ.'Y'
         write(31,*) second,' Generate a second list.'
         if (.NOT.second) return 
      endif
      anel  = 0
      anela = 0
      anelb = 0
      par   = 0
      if (dyn) then
*  50    write(*,200) 'Maximim number of open subshells? (1..8)'
*        read(*,*,err=50) skal
*        skal = max(skal,1)
*        skal = min(skal,8)
         skal = 8
      else
   51    write(*,200) 'Maximum number of open subshells? (1..5)'
         read(*,*,err=51) skal
         skal = max(skal,1)
         skal = min(skal,5)
         write(31,*) 'Maximium number of open shells =',skal
      endif
   60 write(*,200) 'Highest n-number? (1..15)'
      read(*,*,err=60) nmax
      nmax = max(nmax,1)
      nmax = min(nmax,15)
      write(31,*) nmax,' Highest principal quantum number.'
   70 write(*,400) 'Highest l-number? (s..',orb(min(10,nmax-1)),')'
      read(*,1000) X
      lmax = -1
      do 71 i=0,min(10,nmax-1)
   71    if (X.EQ.orb(i)) lmax=i
      if (lmax.EQ.-1) goto 70
      write(31,*) lmax,' Highest orbital angular momentum.'
      write(*,200) 'Are all these nl-subshells active? (n/*)'
      read(*,1000) X
      all = .NOT.(X.EQ.'n' .OR. X.EQ.'N')
      write(31,*) all,' all subshells active.'
      do 72 i=1,15
   72    lim(i) = 0
      if (nmax.GE.2) then
         write(*,200) 'Limitations on population of n-subshells? (y/*)'
         read(*,1000) X
         lima = X.EQ.'y' .OR. X.EQ.'Y'
         write(31,*) lima,' limitations on population of n-subshells.'
         if (lima) then
            mshell = 0
            do 85 i=1,nmax-1
               mshell = mshell + 2*i*i
   83          continue
               if (i.EQ.1) then
                  write(*,200)
     :                    'Minimum number of electrons with n=1? (0..2)'
               elseif (i.LT.10) then
                  if (mshell.LT.100) then
                     write(*,208) 'Minimum number of electrons with n<='
     :                            ,i,'? (0..',mshell,')'
                  else
                     write(*,208) 'Minimum number of electrons with n<='
     :                            ,i,'? (0..)'
                  endif
               else
                  write(*,202) 'Minimum number of electrons with n<=',i,
     :                         '? (0..)'
               endif
               read(*,*,err=83) lim(i)
               if (lim(i).GT.mshell) lim(i) = mshell
               write(31,*) lim(i),
     :                     ' is minimum number of electrons with n =',i
   85       continue
         endif
      endif
   95 continue
      if (nmax.LT.10) then
         write(*,200)
     :      'Highest n-number in reference configuration? (1..',nmax,')'
      else
         write(*,202)
     :      'Highest n-number in reference configuration? (1..',nmax,')'
      endif
      read(*,*,err=95) nenter
      nenter = max(nenter,1)
      nenter = min(nenter,nmax)
      write(31,*) nenter,' highest n-number.'
      block  = 0
      vir    = .FALSE.
      do 160 i=1,15
         do 150 j=0,min(10,i-1)
            low(i,j)    = 0
            virtu(i,j)  = .FALSE.
            dubbel(i,j) = .FALSE.
            if (nmax.GE.i .AND. lmax.GE.j .AND. .NOT.closed(i,j)) then
               if (nenter.GE.i) then
                  em = 2 + 4*j
                  if (j.GE.3) em = 2 
*Limit; higher angular momentum orbits than d may only be populated with
*       2 or less electrons.
                  if (em.LT.10) then 
  100                continue
                     if (i.LE.9) then
                        write(*,200) 'Number of electrons in ',i,orb(j),
     :                                  '? (0..',em,')'
                     else
                        write(*,202) 'Number of electrons in ',i,orb(j),
     :                                  '? (0..',em,')'
                     endif
                     read(*,*,err=100) org(i,j)
                     if (org(i,j).LT.0 .OR. org(i,j).GT.em) goto 100
                  else
  101                continue
                     if (i.LT.10) then
                        write(*,201) 'Number of electrons in ',i,orb(j),
     :                                  '? (0..',em,')'
                     else
                        write(*,203) 'Number of electrons in ',i,orb(j),
     :                                  '? (0..',em,')'
                     endif
                     read(*,*,err=101) org(i,j)
                     if (org(i,j).LT.0 .OR. org(i,j).GT.em) goto 101
                  endif
                  write(31,*) org(i,j),' number of electrons in',i,
     :                        orb(j)
                  if (all) then
                     lock(i,j) = .FALSE.
                  else
                     if (org(i,j).GT.1) then
                        write(*,201)
     :          'Inactive, active or minimum? (i/*/0..',org(i,j)-1,')'
                        read(*,1000) X
                     elseif (org(i,j).EQ.1) then 
                        write(*,201) 'Inactive or active? (i/*)'
                        read(*,1000) X
                     else
                        write(*,400) 'Inactive, active, doubled ',
     :                               'excited or virtual? (i/*/d/v)'
                        read(*,1000) X
                        dubbel(i,j) = X.EQ.'d' .OR. X.EQ.'D' 
                        virtu(i,j)  = X.EQ.'v' .OR. X.EQ.'V' 
                        vir         = vir .OR. virtu(i,j)
                     endif
                     lock(i,j) = X.EQ.'i' .OR. X.EQ.'I' .OR. virtu(i,j)
                     if (X.GE.'0' .AND. X.LE.'9')
     :                     low(i,j) = min(org(i,j),ICHAR(X)-ICHAR('0'))
                     write(31,1000) X,' inactive, active, etc...'
                  endif
                  if (.NOT.lock(i,j)) anela = anela + org(i,j)
                  anel = anel + org(i,j)
                  par  = mod(par+j*org(i,j),2)
               elseif (all) then
                  org(i,j) = 0
                  lock(i,j) = .FALSE.
               else
                  org(i,j) = 0
                  closed(i,j) = .FALSE.
                  if (i.LT.10) then
                     write(*,204) i,orb(j),' inactive, active, ',
     :                           'doubled excited or virtual? (i/*/d/v)'
                  else
                     write(*,205) i,orb(j),' inactive, active, ',
     :                           'doubled excited or virtual? (i/*/d/v)'
                  endif
                  read(*,1000) X
                  dubbel(i,j) = X.EQ.'d' .OR. X.EQ.'D' 
                  virtu(i,j)  = X.EQ.'v' .OR. X.EQ.'V'
                  lock(i,j)   = X.EQ.'i' .OR. X.EQ.'I' .OR. virtu(i,j)
                  vir         = vir .OR. virtu(i,j)
                  write(31,*) X,i,orb(j),' inactive, active, etc...'     
               endif
            else
               org(i,j)  = 0
               lock(i,j) = .TRUE.
               if (closed(i,j)) then
                  if (i.LT.10) then
                     write(*,204) i,orb(j),' is a closed shell.'
                  else
                     write(*,205) i,orb(j),' is a closed shell.' 
                  endif
                  em    = 2 + 4*j
                  anel  = anel + em
                  block = block + em
               endif
            endif
            anelb = anelb + low(i,j)
  150    continue
         lim(i) = lim(i) - block
         if (lim(i).LT.0) lim(i) = 0
  160 continue
      if (anel.NE.anel0) then
         if (anel0.LT.10) write(*,300)
     :             'Wrong number of electrons. The first list had ',
     :                 anel0,' electrons.'
         if (anel0.GE.10) write(*,301)
     :             'Wrong number of electrons. The first list had ',
     :                 anel0,' electrons.'
         if (anel.LT.10) write(*,300) 'This list has ',anel,
     :                                 ' electrons.'
         if (anel.GE.10) write(*,301) 'This list has ',anel,
     :                                 ' electrons.'
         second = .FALSE.
         goto 40
      endif
      if (par.NE.par0) then
         write(*,200) 'Wrong parity.'
         if (par0.EQ.0) write(*,*)
     :           'The first list had even parity and this list has odd.'
         if (par0.EQ.1) write(*,*)
     :           'The first list had odd parity and this list has even.'
         second = .FALSE.
         goto 40
      endif
      if (breit) then
 1154    write(*,200) 'Maximum 2*J-value? (0..)'
         read(*,*,err=1154) J2max
         if (J2max.LT.0) J2max = 0
         write(31,*) J2max,' maximum 2*J-value.'
         if (J2max.EQ.0) then
            J2min = 0
         else
 1155       write(*,200) 'Minimum 2*J-value? (0..',J2max,')'
            read(*,*,err=1155) J2min
         write(31,*) J2min,' minimum 2*J-value.'
         endif
 1156    write(*,200) 'Maximum (2*S+1)-value? (1..9)'
         read(*,*,err=1156) maxS
         write(31,*) maxS,' maximum (2*S+1)-value.'
         if (maxS.GT.9) maxS = 9
         if (maxS.LE.1) then
            maxS = 1
            minS = 1
         else
 1157       write(*,200) 'Minimum (2*S+1)-value? (1..',maxS,')' 
            read(*,*,err=1157) minS
            if (minS.GT.maxS) minS = maxS
            if (minS.LT.1)    minS = 1
         endif
         write(31,*) minS,' minimum (2*S+1)-value.'
 1158    write(*,200) 'Maximum resulting angular momentum? (S..N/N=*)'
         read(*,1000,ERR=1158) X
         maxL = 10
*        do 1161 i=0,20
         do 1161 i=0,10 
*Limit; NONH does not accept higher resultant angular momentum than "N".
 1161       if (X.EQ.L(i)) maxL=i
         write(31,*) 'Maximum resulting angular momentum = ',L(maxL) 
*        if (maxL.EQ.-1) then
*           write(*,200) 'S, P, D, F, G, H, I, K, L, M or N are valid!'
*           goto 1158
*        endif
         write(31,*) L(maxL),' = maximium resulting angular momentum'
         minL = 0
         if (maxL.GT.0) then
 1159       write(*,400) 'Minimum resulting angular momentum? (S..',
     :                   L(maxL),'/S=*)'
            read(*,1000,ERR=1159) X
            do 1162 i=0,maxL
 1162          if (X.EQ.L(i)) minL = i 
            write(31,*) L(minL),' = minimium resulting angular momentum'
         endif
      endif
      anelb = anela - anelb
 1200 continue
      if (anelb.LT.10) then
         write(*,200) 'Number of excitations = ? (0..',anelb,')'
         read(*,*,err=1200) varmax
         if (vir) then
 1210       write(*,401) 'Number of excitations from the active set',
     :                      ' to the virtual set? (0..',anelb,')'
            read(*,*,err=1210) virmax
         endif
      else
         write(*,202) 'Number of excitations = ? (0..',anelb,')'
         read(*,*,err=1200) varmax
         if (vir) then
 1220       write(*,402) 'Number of excitations from the active set',
     :                      ' to the virtual set? (0..',anelb,')'
            read(*,*,err=1220) virmax
         endif
      endif
      write(31,*) varmax,' number of excitations.'
      if (vir) write(31,*) virmax,
     :                  ' number of excitations to virtual set.'
  200 format(' ',A,I1,A,A,I1,A)
  201 format(' ',A,I1,A,A,I2,A)
  202 format(' ',A,I2,A,A,I1,A)
  203 format(' ',A,I2,A,A,I2,A)
  204 format(' ',I1,3A)
  205 format(' ',I2,3A)
  208 format(' ',A,I1,A,I2,A)
  300 format(' ',A,I1,A)
  301 format(' ',A,I2,A)
  400 format(' ',3A)
  401 format(' ',2A,I1,A)
  402 format(' ',2A,I2,A)
 1000 format(A,A,A)
 2000 format(I1,A)
      return
      end
*     last edited February 21, 1993
      subroutine Merge(single,dyn)
      character rad11*64,rad12*64,rad21*72,rad22*72
      logical p1,p2,slut1,slut2,Lika,single,dyn
      integer pop1(15,0:10),pop2(15,0:10),skal1,skal2,popo(15,0:10)
      integer i,j,cf,stopp1,stopp2
      if (dyn) then
         open(unit=9,file='clist.out',status='unknown')
         open(unit=13,file='clist.new',status='unknown')
      else
         open(unit=9,file='cfg.inp',status='unknown')
      endif
      slut1 = .FALSE.
      slut2 = single
      cf    = 0
      call Twofirst(slut1,slut2)
      call Lasa1(7,rad11,pop1,skal1,slut1,dyn)
      call Lasa1(8,rad12,pop2,skal2,slut2,dyn)
   10 if (.NOT.slut1 .AND. .NOT.slut2) then
         call Test(p1,p2,pop1,pop2,15)
         if (p1) then
            do 20 i=1,15
               do 20 j=0,min(10,i-1)
   20             popo(i,j) = pop1(i,j)
            stopp1 = 8*skal1
            if (dyn) then
               stopp2 =  8*skal1 - 4
            else
               stopp2 = 16*skal1 - 8
            endif
   30       call Lasa2(7,rad21,stopp2,slut1)
            if (.NOT.slut1) then
               write(9,999) rad11(1:stopp1)
               write(9,999) rad21(1:stopp2)
               cf = cf + 1
            endif
            call Lasa1(7,rad11,pop1,skal1,slut1,dyn)
            if (.NOT.slut1) then
               if (Lika(popo,pop1)) goto 30
            endif
            if (p2) then
   40          call Lasa2(8,rad22,stopp2,slut2)
               call Lasa1(8,rad12,pop2,skal2,slut2,dyn)
               if (.NOT.slut2) then
                  if (Lika(popo,pop2)) goto 40
               endif
            endif
            goto  10
         elseif (p2) then
            do 50 i=1,15
               do 50 j=0,min(10,i-1)
   50             popo(i,j) = pop2(i,j)
            stopp1 = 8*skal2
            if (dyn) then
               stopp2 =  8*skal2 - 4
            else
               stopp2 = 16*skal2 - 8
            endif
   60       call Lasa2(8,rad22,stopp2,slut2)
            if (.NOT.slut2) then
               write(9,999) rad12(1:stopp1)
               write(9,999) rad22(1:stopp2)
               if (dyn) then
                  write(13,999) rad12(1:stopp1)
                  write(13,999) rad22(1:stopp2)
               endif
               cf = cf + 1
            endif
            call Lasa1(8,rad12,pop2,skal2,slut2,dyn)
            if (.NOT.slut2) then
               if (Lika(popo,pop2)) goto 60
            endif
            goto 10
         else
            write(*,*) 'fatal error'
            stop
         endif
      elseif (.NOT.slut1 .AND. slut2) then
   70    stopp1 = 8*skal1
         if (dyn) then
            stopp2 =  8*skal1 - 4
         else
            stopp2 = 16*skal1 - 8
         endif
         call Lasa2(7,rad21,stopp2,slut1)
         if (.NOT.slut1) then
            write(9,999) rad11(1:stopp1)
            write(9,999) rad21(1:stopp2)
            cf = cf + 1
         endif
         call Lasa1(7,rad11,pop1,skal1,slut1,dyn)
         if (.NOT.slut1) goto 70
      elseif (slut1 .AND. .NOT.slut2) then
   80    stopp1 = 8*skal2
         if (dyn) then
            stopp2 =  8*skal2 - 4
         else
            stopp2 = 16*skal2 - 8
         endif
         call Lasa2(8,rad22,stopp2,slut2)
         if (.NOT.slut2) then
            write(9,999) rad12(1:stopp1)
            write(9,999) rad22(1:stopp2)
            if (dyn) then
               write(13,999) rad12(1:stopp1)
               write(13,999) rad22(1:stopp2)
            endif
            cf = cf + 1
         endif
         call Lasa1(8,rad12,pop2,skal2,slut2,dyn)
         if (.NOT.slut2) goto 80
      endif
      close(7)
      close(8)
      write(9,999) '*'
      close(9)
      if (dyn) close(13)
      if (cf.EQ.0) then
         write(*,105) 'No configuration state in the final list.'
      elseif (cf.EQ.1) then
         write(*,105) 'One configuration state in the final list.'
      elseif (cf.LT.10) then
         write(*,101) cf,' configuration states in the final list.'
      elseif (cf.LT.100) then
         write(*,102) cf,' configuration states in the final list.'
      elseif (cf.LT.1000) then
         write(*,103) cf,' configuration states in the final list.'
      elseif (cf.LT.10000) then
         write(*,104) cf,' configuration states in the final list.'
      elseif (cf.LT.100000) then
         write(*,106) cf,' configuration states in the final list.'
      else
         write(*,*) cf,' configuration states in the final list.'
      endif
      return
  101 format(' ',I1,A)
  102 format(' ',I2,A)
  103 format(' ',I3,A)
  104 format(' ',I4,A)
  105 format(' ',A)
  106 format(' ',I5,A)
  999 format(A)
      end
*     last edited April 26, 1993
      subroutine Mergeb(nmax,antal)
      logical p1,p2,slut1,slut2
      integer pop1(15,0:10),pop2(15,0:10),popo(15,0:10),nmax
      integer i,j,antal
      slut1 = .FALSE.
      slut2 = .FALSE.
      antal = 0
      open(unit=22,status='scratch')
      do 1 i=1,nmax
    1    read(20,5000,end=2) (pop1(i,j),j=0,min(10,i-1))
      goto 3
    2 slut1 = .TRUE.
    3 do 4 i=1,nmax
    4    read(21,5000,end=5) (pop2(i,j),j=0,min(10,i-1))
      goto 6
    5 slut2 = .TRUE.
    6 continue
   10 if (.NOT.slut1 .AND. .NOT.slut2) then
         call Test(p1,p2,pop1,pop2,nmax)
         if (p1) then
            do 20 i=1,nmax
               do 20 j=0,min(10,i-1)
   20             popo(i,j) = pop1(i,j)
            do 120 i=1,nmax
  120          write(22,5000) (pop1(i,j),j=0,min(10,i-1))
            do 121 i=1,nmax
  121          read(20,5000,end=21) (pop1(i,j),j=0,min(10,i-1))
            goto 22
   21       slut1 = .TRUE.
   22       continue
            if (p2) then
            do 122 i=1,nmax
  122          read(21,5000,end=23) (pop2(i,j),j=0,min(10,i-1))
               goto 10
   23          slut2 = .TRUE.
            endif
         elseif (p2) then
            do 50 i=1,nmax
               do 50 j=0,min(10,i-1)
   50             popo(i,j) = pop2(i,j)
            if (.NOT.slut2) then
            do 51 i=1,nmax
   51          write(22,5000) (pop2(i,j),j=0,min(10,i-1))
            endif
            do 52 i=1,nmax
   52          read(21,5000,end=53) (pop2(i,j),j=0,min(10,i-1))
            goto 10
   53       slut2 = .TRUE.
         else
            write(*,*) 'fatal error'
            stop
         endif
         goto 10
      elseif (.NOT.slut1 .AND. slut2) then
   70    continue
         do 170 i=1,nmax
  170      write(22,5000) (pop1(i,j),j=0,min(10,i-1))
         do 171 i=1,nmax
  171       read(20,5000,end=71) (pop1(i,j),j=0,min(10,i-1))
         goto 70
   71    slut1 = .TRUE.
      elseif (slut1 .AND. .NOT.slut2) then
   80    continue
         do 180 i=1,nmax
  180       write(22,5000) (pop2(i,j),j=0,min(10,i-1)) 
         do 181 i=1,nmax
  181       read(21,5000,end=81) (pop2(i,j),j=0,min(10,i-1))
         goto 80
   81    slut2 = .TRUE.
      endif
      rewind(22)
      close(20)
      close(21)
      open(unit=20,status='scratch')
  580 continue
      do 581 i=1,nmax
  581    read(22,5000,end=999) (pop2(i,j),j=0,min(10,i-1))
      do 582 i=1,nmax
  582    write(20,5000) (pop2(i,j),j=0,min(10,i-1))
      antal = antal + 1 
      goto 580
  999 close(22)
      rewind(20)
      return
 5000 format(11I2)
      end
*     last edited April 1, 1993 
      subroutine Reada(rad1,pop,skal,slut,dyn)
      character rad1*64,orb(0:10)
      integer   skal,i,j,n,l,antal,pop(15,0:10),stopp
      logical   slut,dyn
      data (orb(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/
      slut = .false.
      if (rad1(1:1).EQ.'*' .OR. rad1(3:3).EQ.' ') then
         slut = .TRUE.
         return
      endif
      do 5 n=1,15
         do 5 l=0,min(10,n-1)
    5       pop(n,l) = 0
      if (dyn) then
         stopp = 7
      else
         stopp = 4
      endif
      do 10 i=0,stopp
         j = 8*i
         if (rad1(j+3:j+3).EQ.' ') return
         skal = i + 1
         slut = .TRUE.
         n = ichar(rad1(j+3:j+3)) - ichar('0')
         if (n.LE.15 .AND. n.GE.1) then
            do 15 l=0,min(10,n-1)
               if (rad1(j+4:j+4).EQ.orb(l)) then
                  slut = .FALSE.
                  if (rad1(j+6:j+6).EQ.' ' .OR. rad1(j+6:j+6).EQ.'0')
     :                                                             then
                     antal = 0
                  else
                     antal = 10*(ichar(rad1(j+6:j+6)) - ichar('0'))
                  endif
                  antal = antal + ichar(rad1(j+7:j+7)) - ichar('0')
                  if (antal.GT.4*l+2) then
                     slut = .TRUE.
                     return
                  endif
                  pop(n,l) = antal
                  goto 10
               endif
   15       continue
         else
            slut = .TRUE.
            return
         endif
   10 continue
      return
      end
*     last edited June 20, 1993
      subroutine Reffa(posn,posl,dyn)
      integer posn(110),posl(110),stat(15,0:10),i,n,l,num
      character orb(0:10)*1,M*1,X*1
      logical dyn,OK
      data (orb(i),i=0,10) /'s','p','d','f','g','h','i','k','l','m','n'/
      do 10 n=1,15
         do 10 l=0,min(n-1,10)
   10       stat(n,l) = 0
      write(*,200) 
     :       'Default, symmetry or user specified ordering? (*/s/u)'
      read(*,1000) X
      if (X.EQ.'u' .OR. X.EQ.'U') then
         write(31,*) 'User specified ordering.'
         if (dyn) then
            inquire(file='clist.ref',exist=OK)
            if (OK) open(unit=18,status='old',file='clist.ref')
         else
            inquire(file='cfg.ref',exist=OK)
            if (OK) open(unit=18,status='old',file='cfg.ref')
         endif
         l   = -1
         num = 1
         if (.NOT.OK) then
            write(*,200) 'No reference file present! ',
     :              'The couplings will appear in standard order.'
         else
            write(*,200) 'Reference file present!'
   20       read(18,1000,end=40) M,X
               n = ichar(M) - ichar('0')
               do 30 i=0,10
                  if (orb(i).EQ.X) l=i
   30          continue
               if (l.EQ.-1 .OR. n.LT.0 .OR. n.GT.15 .OR.
     :                                      n.LE.l .OR. l.GT.10) goto 40
               if (stat(n,l).NE.0) then
                  write(*,200) 
     :                 'The same orbital appeared more than once!'
                  l = -1
                  goto 20
               endif
               posn(num) = n
               posl(num) = l
               stat(n,l) = num
               num       = num + 1
               l         = -1
               goto 20
   40       continue
            if (num.EQ.1) then
               write(*,200) 'The program failed reading the order of ',
     :                 'the customized coupling scheme.'
            else
               write(*,200) 'The couplings will ',
     :                      'be made in the following customized order:'
               if (num.EQ.2) then
                  write(*,100) posn(1),orb(posl(1))
               else
                  write(*,100) posn(1),orb(posl(1)),
     :                      (',',posn(i),orb(posl(i)),i=2,num-1)
               endif
            endif
         endif
         do 50 n=1,15
            do 50 l=0,min(n-1,10)
               if (stat(n,l).EQ.0) then
                  posn(num) = n
                  posl(num) = l
                  num       = num + 1
               endif
   50    continue
         close(18)
         write(*,200)
      elseif (X.EQ.'s' .OR. X.EQ.'S') then
         write(31,*) 'Symetri ordering.'
         num = 1
         do 60 l=0,10
            do 60 n=l+1,15
               posn(num) = n
               posl(num) = l
   60          num       = num + 1
      else
         write(31,*) 'Standard ordering.'
         num = 1
         do 70 n=1,15
            do 70 l=0,min(n-1,10)
               posn(num) = n
               posl(num) = l
   70          num       = num + 1
      endif
      return
  100 format(' ',110(I2,A,A))
  200 format(' ',2A)
 1000 format(2A)
      end
*     last edited April 1,1993
      subroutine Skifta(pos,nr)
      integer    noofco
Ctc   parameter (noofco=200000)
      parameter (noofco=2100000)
      integer    pos(0:noofco),nr,temp

      temp       = pos(nr)
      pos(nr)    = pos(nr+1)
      pos(nr+1)  = temp

      return
      end
*     last edited February 21, 1993
      subroutine Skriva(rad1,rad2,rank,pos,antal,maskal,opt)
      integer    noofco
Ctc   parameter (noofco=200000)
      parameter (noofco=2100000)
      character  rad1(0:noofco)*64,rad2(0:noofco)*60,ref*64,
     f           reffa(noofco)*64
      integer    rank(0:noofco,4,15),antal,pos(0:noofco),i,j,k,no,
     f           ant(noofco),numb
      logical    Block,ny,lika(8),opt
      real       summa,summa2,mean,stdv

      if (opt) open(unit=17,file='cfg.lst',status='unknown')
      ny     = .TRUE.
      no     = 0
      numb   = 0
      summa  = 0.
      summa2 = 0.
      do 60 i=0,antal-1
         do 10 j=64,1,-1
   10       if (rad1(pos(i))(j:j).NE.' ') goto 20
   20    if (opt) then
            if (ny) ref = rad1(pos(i))
            no = no + 1
            if (.NOT.Block(rank,pos,i,maskal,ny,lika)) then 
               numb        = numb + 1
               reffa(numb) = ref
               ant(numb)   = no
               summa       = summa  + no
               summa2      = summa2 + no*no 
               ny          = .TRUE.
               no          = 0
            else
               ny          = .FALSE.
            endif
         endif
         write(13,1000) rad1(pos(i))(1:j)
         do 40  j=60,0,-1
   40       if (rad2(pos(i))(j:j).NE.' ') goto 50
   50    write(13,1000) rad2(pos(i))(1:j)

   60 continue

      write(13,1000) '*'
      if (opt) then
         mean = summa/numb
*        write(*,*) 'medelvarde =',mean
         stdv = sqrt((summa2 - mean*summa)/(numb-1))
*        write(*,*) 'standardavvikelse =',stdv
         write(17,*) numb,mean,stdv
         do 70 i=1,numb
   70       write(17,2000) reffa(i),ant(i)

         close(17)
      endif
      return
 1000 format(A)
 2000 format(A,I3)
      end
*     last edited February 21, 1993
      subroutine Slug(i,j,varmax,varupp,varned,ansats,org,lock,dubbel,
     :                                                 low,start,stopp)
      integer i,j,varmax,varupp(1:15,0:10),varned(1:15,0:10),minmax
      integer ansats(1:15,0:10),org(1:15,0:10),start,stopp,iold,jold
      integer low(1:15,0:10)
      logical lock,dubbel(1:15,0:10) 
      if (i.EQ.1) then
         varupp(1,0) = 0
         varned(1,0) = 0
      else
         if (j.EQ.0) then
            iold = i - 1
            jold = min(10,iold-1)
         else
            iold = i
            jold = j - 1
         endif
         varupp(i,j) = varupp(iold,jold) +
     :                 max(0,ansats(iold,jold)-org(iold,jold))
         varned(i,j) = varned(iold,jold) +
     :                 max(0,org(iold,jold)-ansats(iold,jold))
      endif
      if (lock) then
         start = org(i,j)
         stopp = org(i,j)
         return
      endif
      if (j.GE.3) then
*Limit; high angular momentum orbits may only have at the most 2
*       electrons.
         minmax = 2
      else
         minmax = 4*j + 2
      endif
      start = min(minmax,org(i,j)+(varmax-varupp(i,j)))
      if (dubbel(i,j)) start = 2*(start/2)
      stopp = max(low(i,j),org(i,j)-(varmax-varned(i,j)))
      return
      end
*     last edited February 21, 1993
      subroutine Sluggo(
     :   i,j,varmax,varupp,varned,ansats,org,lock,low,start,stopp,virtu)
      integer i,j,varmax,varupp(1:15,0:10),varned(1:15,0:10),minmax
      integer ansats(1:15,0:10),org(1:15,0:10),start,stopp,iold,jold
      integer low(1:15,0:10)
      logical lock,virtu(1:15,0:10)
      if (i.EQ.1) then
         varupp(1,0) = 0
         varned(1,0) = 0
      else
         if (j.EQ.0) then
            iold = i - 1
            jold = min(10,iold-1)
         else
            iold = i
            jold = j - 1
         endif
         varupp(i,j) = varupp(iold,jold) +
     :                 max(0,ansats(iold,jold)-org(iold,jold))
         varned(i,j) = varned(iold,jold) +
     :                 max(0,org(iold,jold)-ansats(iold,jold))
      endif
      if (lock .AND. .NOT.virtu(i,j)) then
         start = org(i,j)
         stopp = org(i,j)
         return
      endif
      if (j.GE.3) then
*Limit; high angular momentum orbits may only have at the most 2
*       electrons.
         minmax = 2
      else
         minmax = 4*j + 2
      endif
      if (virtu(i,j)) then
         if (j.GE.3) then
*Limit; high angular momentum orbits may only have at the most 2
*       electrons.
            minmax = 2
         else
            minmax = 4*j + 2
         endif
         start = min(minmax,org(i,j)+(varmax-varupp(i,j)))
      else
         start = org(i,j)
      endif
      stopp = max(low(i,j),org(i,j)-(varmax-varned(i,j)))
      return
      end
*     last edited April 1, 1993
      subroutine Sort(opt)
      integer    noofco
Ctc   parameter (noofco=200000)
      parameter (noofco=2100000)
      logical opt
      character  rad1(0:noofco)*64,rad2(0:noofco)*60
      integer    rank(0:noofco,0:10,15),antal,pos(0:noofco),maskal,
     :           nmax,lmax

      call Lasa(rad1,rad2,antal)
      call Koda(rad1,rad2,rank,antal,maskal,nmax,lmax,opt)
      call Bubbla(rank,pos,antal,maskal,nmax,lmax,opt)
      call Skriva(rad1,rad2,rank,pos,antal,maskal,opt)

      return
      end
      subroutine Test(p1,p2,pop1,pop2,nmax)
      logical p1,p2
      integer pop1(15,0:10),pop2(15,0:10),n,l,nmax
      p1 = .TRUE.
      p2 = .TRUE.
      do 10 n=1,nmax
         do 10 l=0,min(10,n-1)
            if (pop1(n,l) .LT. pop2(n,l)) then
               p1 = .false.
               return
            elseif (pop1(n,l) .GT. pop2(n,l)) then
               p2 = .false.
               return
            endif
   10 continue
      return
      end
*     last edited February 21, 1993
      subroutine Twofirst(slut1,slut2)
      character rad1*64,rad2*64
      logical slut1,slut2
      read(7,999,end=100) rad1
      if (rad1(1:1).EQ.'*') goto 100
      read(7,999,end=100) rad2
      write(9,999) rad1
      write(9,999) rad2
      if (slut2) return
      read(8,999,end=50) rad1
      if (rad1(1:1).EQ.'*') goto 50
      read(8,999,end=50) rad2
      return
   50 slut2 = .TRUE.
      return
  100 slut1 = .TRUE.
      if (slut2) return
      read(8,999,end=50) rad1
      if (rad1(1:1).EQ.'*') goto 50
      read(8,999,end=50) rad2
      write(9,999) rad1
      write(9,999) rad2
      return
  999 format(A)
      end
*     last edited February 21, 1993
      subroutine Twolines(org,closed,first)
      integer org(1:15,0:10)
      logical closed(1:15,0:10),first
      integer i,j,start,stopp
      character rad*80,orb(0:10)
      data (orb(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/
      if (first) then
         open(unit=7,status='scratch')
         write(7,999)
      else
         open(unit=8,status='scratch')
         write(8,999)
         write(8,999)
      endif
      do 10 i=1,80
   10    rad(i:i) = ' '
      start =-3
      stopp = 0
      do 20 i=1,15
         do 20 j=0,min(10,i-1)
            if (closed(i,j)) then
               start                = start + 4
               stopp                = stopp + 4
               rad(start+2:start+2) = CHAR(ICHAR('0')+i)
               rad(stopp:stopp)     = orb(j)
               org(i,j)             = 0
            endif
   20 continue
      if (.NOT.first) return
      if (stopp .EQ. 0) then
         write(7,999)
      else
         write(7,999) rad(1:stopp)
      endif
  999 format(A)
      return
      end
*     last edited April 1, 1993
      subroutine Expand(nold,dyn)
      implicit character (a-z)
      logical dyn
      integer nold

      integer cf,i,k,skal,pop(1:15,0:10),stopp1,stopp2
      logical finns,slut,closed(1:15,0:10)
      character rad1*64,rad2*72

      cf    = 0
      finns = .FALSE.
      slut  = .FALSE.
      call Lockad(closed,slut,.FALSE.)
      call Twolines(pop,closed,.FALSE.)
   10 if (.NOT.slut) call Lasa1(7,rad1,pop,skal,slut,dyn)
         if (.NOT.slut) then
            do 20 i=1,skal
               k = 8*i - 5
               if (ichar(rad1(k:k))-ichar('0') .EQ. nold) then
                  rad1(k:k) = char(nold + ichar('1'))
                  finns = .TRUE.
               endif
   20       continue
            stopp1 = 8*skal
            if (dyn) then
               stopp2 = 8*skal - 4
            else
               stopp2 = 16*skal - 8
            endif
            read(7,100,end=99) rad2(1:stopp2)
            if (finns) then
               cf = cf+1
               write(8,100) rad1(1:stopp1)
               write(8,100) rad2(1:stopp2)
               finns = .FALSE.
            endif
            goto 10
         endif
   99 write(8,100) '*'
      rewind(7)
      rewind(8)

      if (cf.EQ.0) then
         write(*,1000) 'No configuration state has been generated.'
      elseif (cf.EQ.1) then
         write(*,1000) 'One configuration state has been generated.'
      elseif (cf.LT.10) then
         write(*,1001) cf,' configuration states have been generated.'
      elseif (cf.LT.100) then
         write(*,1002) cf,' configuration states have been generated.'
      elseif (cf.LT.1000) then
         write(*,1003) cf,' configuration states have been generated.'
      elseif (cf.LT.10000) then
         write(*,1004) cf,' configuration states have been generated.'
      elseif (cf.LT.100000) then
         write(*,1005) cf,' configuration states have been generated.'
      else
         write(*,*) cf,' configuration states have been generated.'
      endif

 100  format(A)
 1000 format(' ',A)
 1001 format(' ',I1,A)
 1002 format(' ',I2,A)
 1003 format(' ',I3,A)
 1004 format(' ',I4,A)
 1005 format(' ',I5,A)

      end
*     last edited April 1, 1993
      subroutine Scan(nold)
      implicit character (a-z)
      integer pop(1:15,0:10),skal,i,j,nold
      logical closed(1:15,0:10),slut
      character rad1*64

      slut = .FALSE.
      nold = 0
      call Lockad(closed,slut,.FALSE.)
    1 if (.NOT.slut) call Lasa1(7,rad1,pop,skal,slut,.TRUE.)
      if (.NOT.slut) then
         do 10 i=nold+1,15
            do 10 j=0,min(i-1,10)
   10          if (pop(i,j).NE.0) nold = i
         read(7,*,end=99)
         if (nold.GT.14) goto 99
         goto 1 
      endif
   99 rewind(7)
      return
      end
*     last edited September 26, 1993
      subroutine Blandc(varmax,cfmax,lock,resS,resL,nmax,dyn,
     :                  posn,posl)
      implicit character (a-z)
      integer varmax,cfmax,resS,resL,nmax,posn(110),posl(110)
      integer cf,ansats(1:15,0:10),i,j,n,l,l1,antal,tal,antalc
      integer low(1:15,0:10),tot,lista(15000,15,0:10),k
      logical lock(1:15,0:10),dyn,finns,napp,fool(1:15,0:10)
      logical lik(15000)
      character rad*64,orb(0:10)
      data (orb(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/
      cf     = 0
      antalc = 0
      tot    = 0
      finns  = .FALSE.
      do 1 i=1,15
         do 1 j=0,min(10,i-1)
            ansats(i,j) = 0
            fool(i,j)   = .NOT. lock(i,j)
    1       low(i,j)    = 0
      open(unit=7,status='scratch')
      read(8,1000) rad
      write(7,1000) rad
      read(8,1000) rad
      write(7,1000) rad
    2 read(8,1000,end=200) rad
      read(8,1000,end=200)
C     write(*,*) 'raden =',rad
      tot = tot+1
C     do 10 i=1,nmax
      do 10 i=1,15
         do 10 j=0,min(10,i-1)
   10       lista(tot,i,j) = 0
      do 20 i=1,8
         n = i*8-5
         l = n+1
         tal = ichar(rad(n:n)) - ichar('0')
C        if (tal.GT.nmax) then
C           write(*,*) 'Too high n-value in input-file!'
C           stop
C        endif
         if (tal.GE.1 .AND. tal.LE.15) then
            l1 = -1
            do 15 j=0,tal-1
   15          if (orb(j).EQ.rad(l:l)) l1=j
            if (l1 .EQ. -1) goto 30
         else
            goto 30
         endif
         antal = ichar(rad(l+2:l+2)) - ichar('0')
         if (antal.GE.0 .AND. antal.LE.9) then
            antal = antal*10
         else
            antal = 0
         endif
C        lista(tot,tal,l1) = antal + ichar(rad(l+3:l+3)) - ichar('0')
   20    lista(tot,tal,l1) = antal + ichar(rad(l+3:l+3)) - ichar('0')
C  20    ansats(tal,l1) = antal + ichar(rad(l+3:l+3)) - ichar('0')
C        ansats(tal,l1) = antal + ichar(rad(l+3:l+3)) - ichar('0')
C  20    write(*,*) tal,orb(l1),lista(tot,tal,l1)
   30 continue    
      goto 2
  200 if (tot.EQ.0) then
         write(*,*) 'Nothing in inputfile!'
         stop
      endif
C     write(*,*) 'tot =',tot
      lik(1) = .FALSE.
      if (tot.GE.2) then
         do 310 i=1,tot-1
            do 310 j=i+1,tot
               lik(j) = .TRUE.
               do 300 k=1,nmax
                  do 300 l=0,min(10,k-1)
                     lik(j) = lik(j) .AND. lista(i,k,l).EQ.lista(j,k,l)
  300                if (.NOT.lik(j)) goto 310
  310    continue
      endif

      do 320 i=1,tot
C        write(*,*) 'slinga, inconf nummer',i
         if (.NOT.lik(i)) then
            do 315 k=1,nmax
C           do 315 k=1,15
               do 315 l=0,min(k-1,10)
  315             ansats(k,l) = lista(i,k,l)
            if (finns) then
C              write(*,*) 'position 1'
               open(unit=21,status='scratch')
               call Blandb(ansats,nmax,varmax,lock,21,low,fool)
C              call Blandb(ansats,15,varmax,lock,21,low,fool)
               rewind(21)
               call Mergeb(nmax,antalc)
C              call Mergeb(15,antalc)
               write(*,*) 'number of uncoupled csf:s =',antalc
            else
C              write(*,*) 'position 2'
               open(unit=20,status='scratch')
               call Blandb(ansats,nmax,varmax,lock,20,low,fool)
C              call Blandb(ansats,15,varmax,lock,20,low,fool)
               rewind(20)
               finns = .TRUE.
C              write(*,*) 'number of uncoupled csf:s =',antalc
               antalc = 0
            endif
            if (antalc .GE. cfmax) goto 350
         else
C           write(*,*) 'lik!'
         endif
  320 continue

  350 if (nmax.LT.15) then
         do 391 i=nmax+1,15
            do 391 j=0,min(10,i-1)
  391          ansats(i,j) = 0
      endif
      cf   = 0
      napp = .FALSE.
  490 do 491 i=1,nmax
  491    read(20,5000,end=492) (ansats(i,j),j=0,min(10,i-1))
      call Gen(ansats,posn,posl,resS,resL,
     :             8,cf,dyn,.TRUE.,.FALSE.,0,0,0,0,napp)
      goto 490
  492 continue
      write(7,1000) '*'
      rewind(7)
      if (cf.EQ.0) then
         write(*,1005) 'No configuration state has been generated.'
      elseif (cf.EQ.1) then
         write(*,1005) 'One configuration state has been generated.'
      elseif (cf.LT.10) then
         write(*,1001) cf,' configuration states have been generated.'
      elseif (cf.LT.100) then
         write(*,1002) cf,' configuration states have been generated.'
      elseif (cf.LT.1000) then
         write(*,1003) cf,' configuration states have been generated.'
      elseif (cf.LT.10000) then
         write(*,1004) cf,' configuration states have been generated.'
      elseif (cf.LT.100000) then
         write(*,1006) cf,' configuration states have been generated.'
      else
         write(*,*) cf,' configuration states have been generated.'
      endif
 1000 format(A)
 1001 format(' ',I1,A)
 1002 format(' ',I2,A)
 1003 format(' ',I3,A)
 1004 format(' ',I4,A)
 1005 format(' ',A)
 1006 format(' ',I5,A)
 5000 format(11I2)
      return
      end
*     last edited June 20, 1993 
      subroutine Matcin(lock,closed,varmax,cfmax,nmax)
      integer org(1:15,0:10),anel,par,anel0,par0,low(1:15,0:10),J2max
      integer varmax,i,j,skal,nmax,lmax,em,nenter,anela,anelb,J2min
      integer maxS,maxL,minS,minL,virmax,lim(15),block,mshell,cfmax
      character X*1,orb(0:10),L(0:20)
      logical lock(1:15,0:10),closed(1:15,0:10),second,dyn,all,breit,
     :        virtu(1:15,0:10),vir,dubbel(1:15,0:10)
      data (L(i),i=0,20)/'S','P','D','F','G','H','I','K','L','M','N'
     :                         ,'O','Q','R','T','U','V','W','X','Y','Z'/
      data (orb(i),i=0,10)/'s','p','d','f','g','h','i','k','l','m','n'/
   60 write(*,200) 'Highest n-number? (1..15)'
      read(*,*,err=60) nmax
      nmax = max(nmax,1)
      nmax = min(nmax,15)
      write(31,*) nmax,' Highest principal quantum number.'
   70 write(*,400) 'Highest l-number? (s..',orb(min(10,nmax-1)),')'
      read(*,1000) X
      lmax = -1
      do 71 i=0,min(10,nmax-1)
   71    if (X.EQ.orb(i)) lmax=i
      if (lmax.EQ.-1) goto 70
      write(31,*) lmax,' Highest orbital angular momentum.'
      write(*,200) 'Are all these nl-subshells active? (n/*)'
      read(*,1000) X
      all = .NOT.(X.EQ.'n' .OR. X.EQ.'N')
      write(31,*) all,' all subshells active.'
      do 150 i=1,15
         do 150 j=0,min(10,i-1)
            if (nmax.GE.i .AND. lmax.GE.j .AND. .NOT.closed(i,j)) then
               if (all) then
                  lock(i,j) = .FALSE.
               else
                  if (i.LT.10) then
                     write(*,204) i,orb(j),' inactive or active? ',
     :                           '(i/*)'
                  else
                     write(*,205) i,orb(j),' inactive or active? ',
     :                           '(i/*)'
                  endif
                  read(*,1000) X
                  write(31,*) X,i,orb(j),' inactive, active, etc...'
                  lock(i,j)   = X.EQ.'i' .OR. X.EQ.'I'
               endif
            else
               lock(i,j) = .TRUE.
               if (closed(i,j)) then
                  if (i.LT.10) then
                     write(*,204) i,orb(j),' is a closed shell.'
                  else
                     write(*,205) i,orb(j),' is a closed shell.' 
                  endif
               endif
            endif
  150 continue
  160 write(*,200) 'Number of excitations = ? (0..)'
      read(*,*,err=160) varmax
      write(31,*) varmax,' number of excitations.'
  170 write(*,400) 'Maximum number of uncoupled configuration',
     :             ' states? (0..)'
      read(*,*,err=170) cfmax
      write(31,*) cfmax,' maximum number '
  200 format(' ',A,I1,A,A,I1,A)
  201 format(' ',A,I1,A,A,I2,A)
  202 format(' ',A,I2,A,A,I1,A)
  203 format(' ',A,I2,A,A,I2,A)
  204 format(' ',I1,3A)
  205 format(' ',I2,3A)
  208 format(' ',A,I1,A,I2,A)
  300 format(' ',A,I1,A)
  301 format(' ',A,I2,A)
  400 format(' ',3A)
  401 format(' ',2A,I1,A)
  402 format(' ',2A,I2,A)
 1000 format(A,A,A)
 2000 format(I1,A)
      return
      end
