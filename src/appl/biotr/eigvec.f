*
*     -------------------------------------------------------------
*      E I G V E C 
*     -------------------------------------------------------------
*                                                                  *
*     Written by G. Gaigalas,                                      *
*     Vanderbilt University,  Nashville           September 1997   *
*
      SUBROUTINE EIGVEC(ICI,CONFIGI,CONFIGF,LSJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL rel,vok
      CHARACTER ligne*80
      character line*70,configi*64,configf*64
      character elc(8)*3,couple(15)*3
      integer q(8),lsj(2)
      COMMON /DBG  /IBUGM
      COMMON/INOUT/IREAD,IWRITE,ISCW,
     :         iuc(2),iuw(2),iul(2),iuj(2),iut(2)
      POINTER(QNOC,NOCCSH(1)),(QNELCSH,NELCSH(8,1)),
     :       (QNOCORB,NOCORB(8,1)),(QJ1,J1QNRD(15,1))
      COMMON /NDIMS/ QNOC,QNELCSH,QNOCORB,QJ1,NCFG
      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      COMMON /NOR/NCOM,NCLOSI,NCLOSF,NORBI,NORBF,IWAR
      COMMON /STATE/QET1,QLBL1,QWT1,QJV1,QET2,QLBL2,QWT2,QJV2,NJV(2),
     : NVC(2),LGTH(2),NCF(2),QCFG1,QCFG2
      POINTER(QET1,ET1(1)),(QLBL1,LBL1(1)),(QWT1,WT1(1)),(QJV1,JV1(1)),
     :       (QET2,ET2(1)),(QLBL2,LBL2(1)),(QWT2,WT2(1)),(QJV2,JV2(1)),
     :       (QCFG1,CFG1(8,1)),(QCFG2,CFG2(8,1))
      COMMON /MULT/QSL,QSV,QIL,QIR,QJVL,QJVR
      POINTER(QSL,SL(1)),(QSV,SV(1)),(QIL,IL(1)),(QIR,IR(1)),
     :       (QJVL,JVL(1)),(QJVR,JVR(1))
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,LCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
      COMMON /FO/IWF(2),NCLOS(2),NOPEN(2),MAXNFO
      double precision MAX_wt1, MAX_wt2
*
    1 format(15x,f14.7,1x,A)
   14 FORMAT(//8X,I4,10X,I4)
   15 FORMAT(/I6,F16.9,8A8/(7F11.8))
*
* --- determine the number of different J-values (njv)
*               the total number of eigenvectors (nvc)
*               the total length of the wt vector to allocate (lgth)
*
      if (ici .eq. 0) then
        nvc(1) = 1
        nvc(2) = 1
        lgth(1) = ncf(1)
        lgth(2) = ncf(2)
      else
        call rdegvc
      end if
*
* --- allocate the /STATE/
*
      if(ibugm.ne.0) print*,' nvc(1) = ',nvc(1),' nvc(2) = ',nvc(2)
      if(ibugm.ne.0) print*,' lgth(1)= ',lgth(1),' lgth(2)= ',lgth(2)
      if(ibugm.ne.0) print*,' qet1    allocation: nvc(1)   = ',nvc(1)
      call alloc(qet1,nvc(1),8)
      if(ibugm.ne.0) print*,' qlbl1   allocation: nvc(1)   = ',nvc(1)
      call alloc(qlbl1,nvc(1),4)
      if(ibugm.ne.0) print*,' qwt1    allocation: lgth(1)  = ',lgth(1)
      call alloc(qwt1,lgth(1),8)
      if(ibugm.ne.0) print*,' qjv1    allocation: nvc(1)   = ',nvc(1)
      call alloc(qjv1,nvc(1),4)
      if(ibugm.ne.0) print*,' qcfg1   allocation: 8*nvc(1) = ',8*nvc(1)
      call alloc(qcfg1,8*nvc(1),8)
      if(ibugm.ne.0) print*,' qet2    allocation: nvc(2)   = ',nvc(2)
      call alloc(qet2,nvc(2),8)
      if(ibugm.ne.0) print*,' qlbl2   allocation: nvc(2)   = ',nvc(2)
      call alloc(qlbl2,nvc(2),4)
      if(ibugm.ne.0) print*,' qwt2    allocation: lgth(2)  = ',lgth(2)
      call alloc(qwt2,lgth(2),8)
      if(ibugm.ne.0) print*,' qjv2    allocation: nvc(2)   = ',nvc(2)
      call alloc(qjv2,nvc(2),4)
      if(ibugm.ne.0) print*,' qcfg2   allocation: 8*nvc(2) = ',8*nvc(2)
      call alloc(qcfg2,8*nvc(2),8)
*
* --- read in all the eigenvectors of the initial state
*
C ja = seniority; jb = 2*L+1 ; jc = 2*S+1
      if (ici .eq. 0) then
        read(iuc(1), 1) et1(1),ligne(1:36)
         read(iuc(1),'(A)') configi
*           .. closed shells
         lwf = iwf(1) - nclos(1)
        if(ibugm.ne.0) print*,' lwf = ',lwf
    2   read (iuc(1),'(A)') ligne
*           .. orbitals (could have 20)
         lwf = lwf -20
         if (lwf .gt. 0) go to 2

        MAX_wt1 = 0.0D0
        do 3 k = 1,ncf(1)
*          if(k.eq.1) then
           read(iuc(1),'(8(1x,a3,1x,i2,1x),f15.12)')
     :         (elc(j),q(j),j=1,8),wt1(k)
           read(iuc(1),'(15(1x,a3))') (couple(j),j=1,15)

           if (abs(wt1(k)).gt.MAX_wt1) then
              MAX_wt1 = abs(wt1(k)) 
              nocc=0
    4         if(elc(nocc+1).ne.'   ') then
                 nocc=nocc+1
                 if(nocc.lt.8) go to 4
              end if
              call pack(nocc,elc,q,couple,configi)
           end if

*          else
*            read(iuc(1),'(t65,f15.12)') wt1(k)
*            read(iuc(1),'(A)')
*          end if

    3   continue
        jd = j1qnrd(2*noccsh(1)-1,1)
        ja = mod(jd,64)
        jd = jd/64
        jb = mod(jd,64)
        jc = jd/64
        jv1(1) = jb + jc - 2
        print*,' seniority = ',ja,' (2L+1) = ',jb,' (2S+1) = ',jc
     :  ,' 2*J = ',jv1(1)
        read(iuc(2), 1) et2(1),ligne(1:36)
         read(iuc(2),'(A)') configf
*           .. closed shells
         lwf = iwf(2) - nclos(2)
        if(ibugm.ne.0) print*,' lwf = ',lwf
        print*,' lwf = ',lwf
    5   read (iuc(2),'(A)') ligne
*           .. orbitals (could have 20)
         lwf = lwf -20
         if (lwf .gt. 0) go to 5

        MAX_wt2 = 0.0D0
        do 6 k = 1,ncf(2)
*          if(k.eq.1) then

           read(iuc(2),'(8(1x,a3,1x,i2,1x),f15.12)')
     :          (elc(j),q(j),j=1,8),wt2(k)
           read(iuc(2),'(15(1x,a3))') (couple(j),j=1,15)
           if (abs(wt2(k)).gt.MAX_wt2) then  
              MAX_wt2 = abs(wt2(k))
              nocc=0
    7         if(elc(nocc+1).ne.'   ') then
                 nocc=nocc+1
                 if(nocc.lt.8) go to 7
              end if
              call pack(nocc,elc,q,couple,configf)
           end if

*          else
*            read(iuc(2),'(t65,f15.12)') wt2(k)
*            read(iuc(2),'(A)')
*          end if
    6   continuE
CGG        mc = mcfg+1
        mc = NCF(1)+1
        jd = j1qnrd(2*noccsh(mc)-1,mc)
        ja = mod(jd,64)
        jd = jd/64
        jb = mod(jd,64)
        jc = jd/64
        jv2(1) = jb + jc - 2
        print*,' seniority = ',ja,' (2L+1) = ',jb,' (2S+1) = ',jc
     :  ,' 2*J = ',jv2(1)
      else
        if (rel) then
          iu=iuj(1)
        else
          iu=iul(1)
        end if
        k1=1
        k2=ncf(1)
        k3=0
        rewind(iu)
        read(iu,'(a)') line
    8   read(iu,14,end=10) jv,nb
        if(ibugm.ne.0) print*,'     jv = ',jv,' nb = ',nb
        if(.not.rel) jv = lsj(1)
        do 9 m = 1,nb
          k3=k3+1
          jv1(k3)=jv
          read(iu,'(/i6,f16.9,2x,8A8)') lbl1(k3),et1(k3),
     :       (cfg1(k,k3),k=1,8)
          read(iu,'(7F11.8)') (wt1(k),k=k1,k2)
          k1=k1+ncf(1)
          k2=k2+ncf(1)
    9   continue
        go to 8
*
* --- read in all the eigenvectors of the final state
*
   10   close(iu)
        if(ibugm.ne.0) print*,' jv1 = ',(jv1(m),m=1,nvc(1))
        if (rel) then
          iu=iuj(2)
        else
          iu=iul(2)
        end if
        k1=1
        k2=ncf(2)
        k3=0
        read(iu,'(a)') line
   11   read(iu,14,end=13) jv,nb
        if(ibugm.ne.0) print*,'     jv = ',jv,' nb = ',nb
        if (.not.rel) jv = lsj(2)
        do 12 m = 1,nb
          k3=k3+1
          jv2(k3)=jv
          read(iu,'(/i6,f16.9,2x,8A8)') lbl2(k3),et2(k3),
     :      (cfg2(k,k3),k=1,8)
          read(iu,'(7F11.8)') (wt2(k),k=k1,k2)
          k1=k1+ncf(2)
          k2=k2+ncf(2)
   12   continue
        go to 11
   13   close(iu)
      end if
      if(ibugm.ne.0) print*,' jv2 = ',(jv2(m),m=1,nvc(2))
      RETURN
      END
