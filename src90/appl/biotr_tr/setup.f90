!
!     ------------------------------------------------------------------
!     S E T U P
!     ------------------------------------------------------------------
!
!     In order to keep NC in COMMON /RED/ as in other
!   routines, the role of NCC and NC has been reversed from
!   the earlier version of this code.
!
      subroutine setup(ja, jb, let) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  double1=>double 
      use medefn_C
      use ndims_C
      use non30_C
      use red_C
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  13:37:41  11/20/01  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setupm_I 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ja 
      integer  :: jb 
      integer , intent(out) :: let 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ia, ib, ipa, ipb, j1, j2, icase, ish, i, k, ii 
      real(double1) :: diff 
      logical :: finisha, finishb 
!-----------------------------------------------
!
!      COMMON/MEDEFN/IHSH,NJ(16),LJ(16),NOSH(16,2),J1QN(31,3,2),IJFUL(16)
!      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
!     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
!      POINTER(QIAJCMP,IAJCMP(1)),(QLJCOMP,LJCOMP(1)),
!     :       (QNJCOMP,NJCOMP(1)),(QIAJCLD,IAJCLD(1)),
!     :       (QLJCLSD,LJCLSD(1))
!      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
!      COMMON /NON30/ QIAJCMP,QNJCOMP,QLJCOMP,QIAJCLD,QLJCLSD,MAXORB
!      COMMON/RED/INDL(20),INDR(20),IORST(20),NCI1,KPL,KPR,NC,LIS(16),
!     : JIST,JFST,NCONTR
!
!     NOTICE THE DIFFERENT NAMES IN THE COMMON BLOCK MEDEFN  -  WE
!      STORE NOSH1(I=1,10) AS NOSH((I=1,10),1) AND NOSH2(I=1,10) AS
!     NOSH((I=1,10),2)   AND USE THE FACT THAT NOSH1 AND NOSH2 WILL THEN
!     BE EQUIVALENT TO THE SINGLE 2-DIMENSIONAL ARRAY NOSH.  SIMILARLY
!     FOR J1QN
!
! === GENERATES THE ARRAYS  NJ,LJ - DEFINING THE QUANTUM NUMBERS OF THE
!     SHELLS,   NOSH - DEFINING THE OCCUPATION OF THE SHELLS,  J1QN -
!     DEFINING THE COUPLING OF THE SHELLS,   FOR EACH OF THE TWO
!     CONFIGURATIONS CONSIDERED.    ONLY THOSE SHELLS OCCURRING IN AT
!     LEAST ONE CONFIGURATION ARE INCLUDED.
!                   AT LEAST TWO SHELLS MUST BE CONSIDERED OCCUPIED.
!     THUS (1S)**2    HELIUM  MUST BE TREATED AS ,E.G., (1S)**2(2S)**0
!     THE SIZE OF THE ARRAYS HERE CALCULATED IS ARRANGED TO BE NO
!     GREATER THAN IS NECESSARY TO INCLUDE ALL ORBITALS WHICH ARE
!     DEEMED TO BE OCCUPIED IN EITHER OR BOTH OF THE CONFIGURATIONS
!     JA,JB
!
! --- INITIALIZE BASIC QUANTITIES - (I1+1) RUNS OVER 1,MAXORB,  IHSH IS
!     THE CURRENT VALUE OF THE HIGHEST OCCUPIED SHELL YET CONSIDERED,
!     WHILE I2HSH=2*IHSH-1
!
      let = 1 
      diff = 0 
      ia = noccsh(ja) 
      ib = noccsh(jb) 
      if (abs(ia - ib) > 1) then 
         let = 0 
         return  
      endif 
      ipa = 1 
      ipb = 1 
      j1 = nocorb(ipa,ja) 
      j2 = nocorb(ipb,jb) 
      finisha = .FALSE. 
      finishb = .FALSE. 
 
      if (ja /= jb) then 
         ihsh = min(ia,ib) + 1 
!       ..determine which shells to advance
!      icase=1 left: =2 both; =3 right
!
   10    continue 
         if (finisha) then 
            icase = 3 
         else if (finishb) then 
            icase = 1 
         else 
            if (j1 < j2) then 
               icase = 1 
            else if (j1 == j2) then 
               icase = 2 
            else 
               icase = 3 
            endif 
         endif 
         select case (icase)  
         case (1)  
!         .. shell occupied only in left
            diff = diff + nelcsh(ipa,ja) 
            if (ipa < ia) then 
               ipa = ipa + 1 
               j1 = nocorb(ipa,ja) 
            else 
               finisha = .TRUE. 
            endif 
         case (2)  
!         .. shell occupied in both
            diff = diff + iabs(nelcsh(ipa,ja)-nelcsh(ipb,jb)) 
            if (ipa < ia) then 
               ipa = ipa + 1 
               j1 = nocorb(ipa,ja) 
            else 
               finisha = .TRUE. 
            endif 
            if (ipb < ib) then 
               ipb = ipb + 1 
               j2 = nocorb(ipb,jb) 
            else 
               finishb = .TRUE. 
            endif 
         case (3)  
!         .. shells occupied only in right
            diff = diff + nelcsh(ipb,jb) 
            if (ipb < ib) then 
               ipb = ipb + 1 
               j2 = nocorb(ipb,jb) 
            else 
               finishb = .TRUE. 
            endif 
         end select 
         if (diff > 2) then 
            let = 0 
            return  
         endif 
         if (.not.(finisha .and. finishb)) go to 10 
      else 
         ihsh = ia 
      endif 
!
!     We have two configurations differing by one electron
!     same angular configuration.
!
      ipa = 1 
      ipb = 1 
      j1 = nocorb(ipa,ja) 
      j2 = nocorb(ipb,jb) 
      finisha = .FALSE. 
      finishb = .FALSE. 
      ish = 0 
      if (finisha) then 
         icase = 3 
      else if (finishb) then 
         icase = 1 
      else 
         if (j1 < j2) then 
            icase = 1 
         else if (j1 == j2) then 
            icase = 2 
         else 
            icase = 3 
         endif 
      endif 
      ish = ish + 1 
      if (icase == 1) then 
!       .. insert shell on left; dummy for right
!       .. in the call to setupm, 1=.not.occupied; 2=occupied
         call setupm (ish, ipa, ipb, ja, jb, 2, 1) 
         if (ipa < ia) then 
            ipa = ipa + 1 
            j1 = nocorb(ipa,ja) 
         else 
            finisha = .TRUE. 
         endif 
      else if (icase == 2) then 
!       .. insert same shell for left and right
         call setupm (ish, ipa, ipb, ja, jb, 2, 2) 
         if (ipa < ia) then 
            ipa = ipa + 1 
            j1 = nocorb(ipa,ja) 
         else 
            finisha = .TRUE. 
         endif 
         if (ipb < ib) then 
            ipb = ipb + 1 
            j2 = nocorb(ipb,jb) 
         else 
            finishb = .TRUE. 
         endif 
      else 
!       .. insert dummy shell on left, shell on right
         call setupm (ish, ipa, ipb, ja, jb, 1, 2) 
         if (ipb < ib) then 
            ipb = ipb + 1 
            j2 = nocorb(ipb,jb) 
         else 
            finishb = .TRUE. 
         endif 
      endif 
      do while(.not.(finisha .and. finishb)) 
         if (finisha) then 
            icase = 3 
         else if (finishb) then 
            icase = 1 
         else 
            if (j1 < j2) then 
               icase = 1 
            else if (j1 == j2) then 
               icase = 2 
            else 
               icase = 3 
            endif 
         endif 
         ish = ish + 1 
         if (icase == 1) then 
!       .. insert shell on left; dummy for right
!       .. in the call to setupm, 1=.not.occupied; 2=occupied
            call setupm (ish, ipa, ipb, ja, jb, 2, 1) 
            if (ipa < ia) then 
               ipa = ipa + 1 
               j1 = nocorb(ipa,ja) 
            else 
               finisha = .TRUE. 
            endif 
         else if (icase == 2) then 
!       .. insert same shell for left and right
            call setupm (ish, ipa, ipb, ja, jb, 2, 2) 
            if (ipa < ia) then 
               ipa = ipa + 1 
               j1 = nocorb(ipa,ja) 
            else 
               finisha = .TRUE. 
            endif 
            if (ipb < ib) then 
               ipb = ipb + 1 
               j2 = nocorb(ipb,jb) 
            else 
               finishb = .TRUE. 
            endif 
         else 
!       .. insert dummy shell on left, shell on right
            call setupm (ish, ipa, ipb, ja, jb, 1, 2) 
            if (ipb < ib) then 
               ipb = ipb + 1 
               j2 = nocorb(ipb,jb) 
            else 
               finishb = .TRUE. 
            endif 
         endif 
      end do 
!
!     Shift the resultant coupling to be contiguous
!     to the shell coupling, if necessary
!
      if (ihsh /= ish) then 
!       if (ihsh .lt. ish) then
!       stop 'Error in SETUP'
!     else
         j1qn(ish+1:2*ish-1,:,:) = j1qn(ihsh+1:ish-1+ihsh,:,:) 
!     end if
         ihsh = ish 
      endif 
      return  
      end subroutine setup 
