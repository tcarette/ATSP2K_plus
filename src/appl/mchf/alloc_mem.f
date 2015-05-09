*
*     ------------------------------------------------------------------
*    3-3       allocate memory for the largest block!
*     ------------------------------------------------------------------
*
*       Data concerning the number of configurations (NCFG), the number
*   and type of electrons in each  configuration,  as  well  as  data
*   associated with the energy expression are read and stored.
*
*
      SUBROUTINE alloc_mem(nwf)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
      PARAMETER (NWD=94,NWC=20,NOD=220,NOFFD=800,MTERM=20,MEIG=20)
      DIMENSION IVAR(NWD)
*
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :     (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :     (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :      QMETH,QIEPTR
*
      POINTER (pkval,kval(1)),(pvalue,value(1)),
     :        (pih,ih(1)),(pjh,jh(1)),
     :        (pcoeff,coeff(1)),(pnijptr,nijptr(1)),
     :        (pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PVALUE,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
*
      COMMON ZZ(NWD),IND(NWD),IELI(8)

      POINTER (qjptr,jptr(1))
      COMMON/COLS/qjptr
      POINTER (QWT,WT(1)),(QWP,WP(1)),  (IQWPTR,W(1))
      COMMON /CFGS/ETOTAL,QWT,QWP,IQWPTR

      POINTER (qhmx,hmx(1)),(qtm,tm(1)),(qtp,tp(1)),
     :        (qdiag,hii(1)),(qiwork,iwork(1))
      common/spd/qhmx,qtm,qtp,qdiag,qiwork

      logical leigen(meig,mterm), lguess
      integer nume(mterm),iws(mterm),iiws(mterm), 
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm),
     :        nze_max(mterm)
      character term_bl(mterm)*3
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess, nze_max
*
      character line*72, t(60)*3, str*72, tmp*5, tmpa*5, tmpb*5

      i = 0 
      j =0
      do nb = 1,nblock
         i = i + nume(nb)*ncfg_bl(nb) 
         j = j + ncfg_bl(nb)
      end do
      if (i.eq.1) i=2*nblock
      call alloc(qeigvec,i,8)
      call alloc(qjptr,j,4)
      call alloc(pen,(nblock*mterm),4) 
      call alloc(iqp,nod*nwf,8)
      call alloc(iqn,nwf,4)
      call alloc(iql,nwf,4)
      call alloc(iqaz,nwf,8)
      call alloc(iqmax,nwf,4)
      call alloc(qvard,nwf,4)
      call alloc(qsum,nwf,8)
      call alloc(qs,nwf,8)
      call alloc(qdpm,nwf,8)
      call alloc(qacc,nwf,8)
      call alloc(qmeth,nwf,4)
      call alloc(qieptr,nwf,4)
         
*
      do nb = 1,nblock
!         lim = min0(nume(nb)+20,ncfg_bl(nb))
         lim = min0(2*nume(nb)+40,ncfg_bl(nb))
         iws(nb) = (2*ncfg_bl(nb)+lim+nume(nb)+10)*lim + nume(nb)
         iiws(nb) = lim
      end do 

      nze = maxval(nze_bl(1:nblock))
      lim = maxval(iiws(1:nblock))
      iiwsz = 6*lim + maxval(nume(1:nblock))
      iworksz = maxval(iws(1:nblock))
      ncfg = maxval(ncfg_bl(1:nblock))
      call alloc(qtm,ncfg,8)
      call alloc(qtp,ncfg,8)
      call alloc(qdiag,ncfg,8)
      call alloc(qiwork,iiwsz,4)
      call alloc(qwt,iworksz,8)

      return 
      END
   
