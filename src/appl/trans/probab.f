*
*     -------------------------------------------------------------
*      P R O B A B 
*     -------------------------------------------------------------
*                                                                  *
*     CALCULATE THE TRASITION PROBABILITIES                        *
*                                                                  *
*     DA > 0       1->2     ABSORTION                              *
*     DA < 0       1->2     EMISSION                               *
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville           September 1997   *
*
      SUBROUTINE 
     :     PROBAB(ICI,IULS,IULSJ,NPAIR,CONFIGI,CONFIGF,LESS,IPRINT,IM)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL LESS,REL,VOK,IPRINT
      CHARACTER*1 IM
      CHARACTER*64 configi,configf
      COMMON /DEBUG/IBUG1,IBUG2,IBUG3,NBUG6,NBUG7,IFULL
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2)
      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      COMMON /NTRM/NTERMS
      COMMON /STATE/QET1,QLBL1,QWT1,QJV1,QET2,QLBL2,QWT2,QJV2,NJV(2),
     : NVC(2),LGTH(2),NCF(2),QCFG1,QCFG2
      POINTER(QET1,ET1(1)),(QLBL1,LBL1(1)),(QWT1,WT1(1)),(QJV1,JV1(1)),
     :       (QET2,ET2(1)),(QLBL2,LBL2(1)),(QWT2,WT2(1)),(QJV2,JV2(1)),
     :       (QCFG1,CFG1(8,1)),(QCFG2,CFG2(8,1))
      COMMON /MULT/QSL,QSV,QIL,QIR,QJVL,QJVR
      POINTER(QSL,SL(1)),(QSV,SV(1)),(QIL,IL(1)),(QIR,IR(1)),
     :       (QJVL,JVL(1)),(QJVR,JVR(1))
      COMMON /RYDBERG/ RYDBRG,ZMU,ENERGY
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,LCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /LSJ/LL1,LL2,IS1,IS2
      COMMON/CONSTS/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
*
    2 FORMAT( ///32X,'NUMBER OF TERMS IN THE ABOVE SUMMATION =',I10)
    3 FORMAT(/'------------------------------------------------------'/
     :'Pair number ',I3,//)
    4 FORMAT(/,' Initial CSF : ',A,' J = ',F3.1)
    5 FORMAT(' Final   CSF : ',A,' J = ',F3.1/)
    6 FORMAT(' 2*j = ',i5,' lbl = ',i5,' total energy = ',f16.7)
    7 FORMAT(/10X,'ENERGY DIFFERENCE OF THE STATES:',9X,1PD15.7,' CM-1',
     :   /51X,1PD15.7,' ANGSTROMS'/51X,1PD15.7,' A. U.')
    8 FORMAT(/,1H ,' SL =',1PE15.7,' TRPT =',1PE15.7,' AKI =',1PE15.7)
    9 FORMAT(/,1H ,' SV =',1PE15.7,' TRPT =',1PE15.7,' AKI =',1PE15.7)
   10 FORMAT(1H ,' I = ',I3,' GL = ',1PD15.7,' J = ',I3, ' D1 = ',F3.0)
*
  237 FORMAT(//,10x,'LENGTH   FORMALISM: ',/10x,'-------- ----------',/)
  233 FORMAT(/10X,'SL                                       = ',1PD15.7)
  234 FORMAT(10X,'FINAL OSCILLATOR STRENGTH (GF)           = ',1PD15.7)
  134 FORMAT(10X,'TRANSITION PROBABILITY IN EMISSION (Aki) =',1PD16.7//)
  235 FORMAT(//,10x,'VELOCITY FORMALISM: ',/10x,'-------- ----------',/)
  236 FORMAT(/10X,'SV                                       = ',1PD15.7)
  238 FORMAT(10X,'FINAL OSCILLATOR STRENGTH (GF)           = ',1PD15.7)
  239 FORMAT(10X,'TRANSITION PROBABILITY IN EMISSION (Aki) =',1PD16.7//)
*
   42 format(10x,' LS RESULTS (length   formalism):  S   = ',D14.6,/
     :   10x,' ----------                        Aki = ',D14.6,/
     :   10x,'                                   gf  = ',D14.6)
   43 format(10x,'            (velocity formalism):  S   = ',D14.6,
     :    /44x,' Aki = ',D14.6,
     :    /44x,' gf  = ',D14.6)
*
*    FORMAT FOR OUTPUT PN tr.lsj FILE
*
   38 FORMAT(F11.2,' CM-1',2X,F11.2,' ANGS(VAC)',2X,F11.2,' ANGS(AIR)'/
     :    1X,A1,I1,2X,'S = ',1PD12.5,3X,'GF = ',D12.5,3X,'AKI = ',D12.5)
*
*     FORMATS FOR  OUTPUT ON tr.ls FILE
*
   41 FORMAT(F11.2,' CM-1',2X,F11.2,' ANGS(VAC)',2X,F11.2,' ANGS(AIR)'/
     :1X,A1,I1,2X,'length:   ','S = ',1PD12.5,3X,'GF = ',D12.5,3X,
     :'AKI = ',D12.5)
   40 FORMAT(4x,'velocity: ',' S = ',1PD12.5,3X,'GF = ',D12.5,3X,
     : 'AKI = ',D12.5)
*
      WRITE(IWRITE,2) NTERMS
      DO 1 K = 1,NPAIR
*
* --- print heading for the pair (j,j')
*
        WRITE(IWRITE,3) K
        IF(ICI.NE.0) THEN
          WRITE(configi,'(8A8)') (cfg1(ii,il(k)),ii=1,8)
          WRITE(configf,'(8A8)') (cfg2(ii,ir(k)),ii=1,8)
        ENDIF 
        AJL=DBLE(JVL(K))/TWO
        AJR=DBLE(JVR(K))/TWO
        WRITE(iwrite,4) configi(1:50),ajl
        WRITE(iwrite,5) configf(1:50),ajr
        WRITE(iwrite,6) jvl(k),lbl1(il(k)),et1(il(k))
        WRITE(iwrite,6) jvr(k),lbl2(ir(k)),et2(ir(k))
*
        da = et2(ir(k)) - et1(il(k))
CGG        ryrat = rydbrg/109737.31534
        sl(k) = sl(k)**2 *(jvl(k)+1)*(jvr(k)+1)
        if(vok) then 
          if(lam.eq.1)  sv(k) = (sv(k)/da)**2 *(jvl(k)+1)*(jvr(k)+1)
          if(lam.eq.2)  sv(k) = (TWO*sv(k)/da)**2 *(jvl(k)+1)*(jvr(k)+1)
        end if
        if(im .eq. 'M') sl(k) = sl(k)*d4*lam*(2*lam-1)
      if(ibug1.ne.0)print*,' k = ',k,' sl(k) = ',sl(k),' sv(k) = ',sv(k)
        if (da .eq. d0) go to 1
        d = dabs(da)
        dd = d*d2*rydbrg
        angs = d10**8/dd
        angsa = angs
        if (angs .gt. 2000.d0) then
          sigma = (1.d8/angs)**2
          angsa = angs/(d1 +8342.13d-8 +206030./(130.d+8 - sigma)
     :                                 +15997./(38.9d+8 - sigma))
        end if
        if (Iprint) write(iwrite,7) dd,angs,d
*
* --- compute the transition probabilities (sec-1) in emission
*     and the gf(1->2) value
*
        i = 2
        jvi = jvr(k)
        if (da .lt. d0) then
          i = 1
          jvi = jvl(k)
        end if
        trpt = trp(im,lam,dd)
        aki = sl(k)*trpt/(jvi+1)
        if(vok) akiv = sv(k)*trpt/(jvi+1)
        if(ibug1.ne.0) write(iwrite,8) sl(k),trpt,aki
        if(ibug1.ne.0.and.vok) write(iwrite,9) sv(k),trpt,akiv
        fe = (angs**2)*1.499193d-16
        if (ibug1.ne.0) write(iwrite,10) i,fe,jvi,d1
        if(vok) gv = fe*akiv*(jvi+1)*(-d1)**i
        gl = fe*aki*(jvi+1)*(-d1)**i
        if (dabs(gl) .le. tol) go to 1
*
        write(iwrite,237)
        write(iwrite,233) sl(k)
        write(iwrite,234) gl
        write(iwrite,134) aki
*
      if(vok) then
        write(iwrite,235)
        write(iwrite,236) sv(k)
        write(iwrite,238) gv
        write(iwrite,239) akiv
      end if
*
*  *****  OUTPUT ON tr.lsj FILE
*
      if(rel) then
      if (less .and. et1(il(k)) .gt. et2(ir(k))) go to 1
*
      WRITE(iulsj,'(/)')
      WRITE(iulsj,'(I4,F14.8,2X,A)')jvl(k),et1(il(k)),configi(1:50)
      WRITE(iulsj,'(I4,F14.8,2X,A)')jvr(k),et2(ir(k)),configf(1:50)
      WRITE(iulsj,38) DD,ANGS,ANGSA,im,lam,SL(k),GL,AKI
      go to 1
*
* --- calculate the transition properties in LS
*
      else
*
      if(da.gt.0) then
        CALL SIXJ(LL2,JVR(K),IS2,JVL(K),LL1,LAM+LAM,1,W)
        PHAS=ONE
        flsj = (jvl(k)+1)*(ll2+1)*w*W
        g = (ll2+1)*(is2+1)
      else
        CALL SIXJ(LL1,JVL(K),IS1,JVR(K),LL2,LAM+LAM,1,W)
        PHAS=-ONE
        flsj = (jvr(k)+1)*(ll1+1)*w*W
        g = (ll1+1)*(is1+1)
      end if
      slsl = sl(k)*(is1+1)/((jvr(k)+1)*(jvl(k)+1)*w**2)
      akilsl = aki/flsj
      glsl = fe*akilsl*g*PHAS
      write(iwrite,42) slsl,akilsl,glsl
*
      if(vok)then
        slsv = sv(k)*(is1+1)/((jvr(k)+1)*(jvl(k)+1)*w*W)
        akilsv = akiv/flsj
        glsv = fe*akilsv*g*phas
        write(iwrite,43) slsv,akilsv,glsv
      end if
*
*  *****  OUTPUT ON tr.ls FILE
*
      WRITE(iuls,'(/)')
      WRITE(iuls,'(I4,F14.8,2X,A)')jvl(k),et1(il(k)),configi(1:50)
      WRITE(iuls,'(I4,F14.8,2X,A)')jvr(k),et2(ir(k)),configf(1:50)
      WRITE(iuls,41) DD,ANGS,ANGSA,im,lam,slsl,Glsl,AKIlsl
      if(vok) WRITE(iuls,40) slsv,Glsv,AKIlsv
      end if
    1 CONTINUE
      RETURN
      END
