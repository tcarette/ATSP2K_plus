*
*     ------------------------------------------------------------------
*      P R O P E R T Y
*     -----------------------------------------------------------------
*
      SUBROUTINE PRPRTY(ieigval,iblock,coef_tmp,SS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER (NWD=70,NWC=20,NOD=220,NOFFD=800)
      parameter(MEIG=20,MTERM=20)
      COMMON/PARAM/  H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,
     :                NCFG,IB,IC,ID,
     :                D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,
     :                NSCF,NCLOSD,RMASS
      POINTER (QWT,WT(1)),(QWP,WP(1)),  (IQWPTR,W(1))
      COMMON /CFGS/ETOTAL,QWT,QWP,IQWPTR
      POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :     (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :     (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
      POINTER (pkval,kval(1)),(pcoef,coef(1)),(pih,ih(1)),(pjh,jh(1)),
     :        (pcoeff,coeff(1)),(pnijptr,nijptr(1)),(pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PCOEF,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
      INTEGER       IN,OUT,PRI,OUC,OUF,OUH,OUD
      COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,iud,OUC,OUF,OUD,OUH,ISCW
      COMMON/LABEL/ EL(NWD),ATOM,TERM
      COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :           QMETH,QIEPTR
      COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR

      integer size, iblock
      dimension coef_tmp(idim)
      character string*64
      dimension eigst_weight(meig,mterm)
      integer cf_tot(nblock)
      logical leigen(meig,mterm), lguess
      integer nume(mterm),iws(mterm),iiws(mterm),
     :        nze_bl(mterm), ncfg_bl(mterm), niv_bl(mterm),
     :        nze_max(mterm)
      character term_bl(mterm)*3
      pointer   (qeigvec,eigvec(1)),(pen,en(1))
      common/st/leigen,nblock,nume,iws,iiws,ncfg_bl,nze_bl,
     :          term_bl, qeigvec, pen, lguess,nze_max
      double precision SS

*     Compute some wave function properties
*    

      WRITE(PRI,'(//15X,A/)')  'Some WaveFunction Properties'
      rmean1 = 0.d0
      rmean2 = 0.d0
      SS = 0.d0
      R_R = 0.0

      Do I = 1,NCLOSD
	sumi = 4*L(i)+2
	rmean1 = rmean1 + sumi*quadr(i,i,1)
	rmean2 = rmean2 + sumi*quadr(i,i,2)
      END DO
*
*     a) From Common core interactions: 2-body  operators
      Do I = 2,NCLOSD
        sumi = 4*L(i)+2
        Do J =  1,I-1
          sumj = 4*L(j)+2
          DO K=IABS(L(I)-L(J)),L(I)+L(J),2 
            IF (K.EQ.1) THEN            
              C= -CB(L(I),L(J),K)           
*             .. this minus is from exchange term
              SS = SS - C*SUMI*SUMJ*GRAD(I,J)**2
	      R_R = R_R + c*SUMI*SUMJ*quadr(i,j,1)**2
*	      print *, i,j, -c*sumi*sumj,-grad(i,j)**2, ss
            ENDIF
          END DO
        END DO
      END DO

*
*     b) From outer electrons
*       i) From G-integrals
      IBEGIN = INTPTR(1)+1
      IEND = INTPTR(2)
      Do I = IBEGIN, IEND
	call unpacki(2,i,kv,iel1,iel2,iel3,iel4)
	IF (kv .eq. 1) then
*         .. this minus is from exchange form
	  SS = SS - coef_tmp(i)*grad(iel1,iel2)**2
	  R_R = R_R +coef_tmp(i)*quadr(iel1,iel2,1)**2
*	  print *, iel1, iel2, coef_tmp(i),-grad(iel1,iel2)**2, ss
	END IF
      END DO

*      ii) From R-integrals
      IBEGIN = INTPTR(4)+1
      IEND = INTPTR(5)
      Do I = IBEGIN, IEND
	IF (coef_tmp(i) .ne. 0.d0 ) THEN
	  call unpacki(5,i,kv,iel1,iel2,iel3,iel4)
	  IF (kv .eq. 1) then
	    SS = SS + coef_tmp(i)*grad(iel1,iel3)*grad(iel2,iel4)
	    R_R = R_R +coef_tmp(i)*quadr(iel1,iel3,1)*quadr(iel2,iel4,1)
*	    print *, iel1, iel3, iel2, iel4, coef_tmp(i),
*    : 	    grad(iel1,iel3)*grad(iel2,iel4), ss
	  END IF
	END IF
      END DO

*     c) Contributions from outer electrons and common core

      IBEGIN = intptr(5)+1
      IEND = intptr(6)
      DO i = ibegin, iend
        call unpacki(6,i,kv,iel1,iel2,iel3,iel4)
        SUMI= -D2*coef_tmp(i)
	rmean1 = rmean1 + sumi*quadr(iel1, iel2, 1)
	rmean2 = rmean2 + sumi*quadr(iel1, iel2, 2)
*	print *, i,iel1,iel2,sumi
        DO  J=1,NCLOSD
          SUMJ=4*L(J)+2
          RR1 = 0.d0
          GR1 = 0.d0
          DO K=IABS(L(iel1)-L(J)),L(iel1)+L(J),2 
              IF (K.EQ.1) THEN        
                C= -CB(L(IEL1),L(J),K)
                IF (IEL1 .EQ. IEL2 ) THEN
*                    .. we have a diagonal case
                     RR1=C*QUADR(IEL1,J,1)*QUADR(J,IEL2,1)
*                    .. this is from exchange form
                     GR1=-C*GRAD(IEL1,J)**2
*                    Print *, iel1,j,j,iel2,-c*sumi*sumj,
*    :                     GRAD(IEL1,J)**2
                  ELSE
*                     .. we have an off-diagonal case
                     GR1=  C*GRAD(IEL1,J)*GRAD(J,IEL2)
		     RR1=  C*Quadr(IEL1,J,1)*Quadr(J,IEL2,1)
*		     Print *, iel1,j,j,iel2,c*sumi*sumj,
*    :		     GRAD(IEL1,J)*GRAD(J,IEL2)
                  ENDIF
               ENDIF
         END DO          
            SS = SS + SUMI*SUMJ*GR1
            R_R = R_R + SUMI*SUMJ*RR1
        END DO
      END DO

*     Change sign of SS to adhere to our definition
      SS = - SS
      write (pri,'(A8,A3)') ' Term ', term_bl(iblock)     
      WRITE(PRI,'(A,F15.8)') ' Mean radius             =', rmean1
      WRITE(PRI,'(A,F15.8)') ' Mean square radius      =', rmean2
      WRITE(PRI,'(A,F15.8)') ' Mean R.R parameter      =', 
     : rmean2 + 2*r_r 
      WRITE(PRI,'(A,F15.8)') ' Isotope Shift parameter =', SS
      RETURN
      END
