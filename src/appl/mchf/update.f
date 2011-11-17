        SUBROUTINE UPDATE
        IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*
        PARAMETER (NWD=70,NOD=220,NOFFD=800)
*
        CHARACTER EL*3,ATOM*6,TERM*6
        INTEGER       IN,OUT,PRI,OUC,OUF,OUH,OUD
        COMMON/INOUT/ IN,OUT,PRI,IUC,IUF,IUD,OUC,OUF,OUD,OUH,ISCW
        COMMON/LABEL/ EL(NWD),ATOM,TERM
        LOGICAL       FAIL,OMIT,EZERO,REL,ALL,TRACE,VARIED
        COMMON/TEST/  FAIL,OMIT,EZERO,REL,ALL,TRACE
       COMMON/WAVE/EC,ED,AZD,PDE(NOD),IJE(noffd),EIJ(noffd),IPR
       COMMON/PARAM/  H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,
     :                NCFG,IB,IC,ID,
     :                D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,
     :                NSCF,NCLOSD,RMASS
      COMMON/RADIAL/R(NOD),RR(NOD),R2(NOD),YK(NOD),YR(NOD),X(NOD)

      POINTER (QWT,WT(1)),(QWP,WP(1)),  (IQWPTR,W(1))
      COMMON /CFGS/ETOTAL,QWT,QWP,IQWPTR
*
        POINTER(IQP,P(NOD,1)),(IQN,N(1)),(IQL,L(1)),(IQAZ,AZ(1)),
     :     (IQMAX,MAX(1)),(QVARD,VARIED(1)),(QSUM,SUM(1)),(QS,S(1)),
     :     (QDPM,DPM(1)),(QACC,ACC(1)),(QMETH,METH(1)),(QIEPTR,IEPTR(1))
        COMMON/NEL/IQP,IQN,IQL,IQAZ,IQMAX,QVARD,QSUM,QS,QDPM,QACC,
     :             QMETH,QIEPTR
*
      POINTER (pkval,kval(1)),(pvalue,value(1)),
     :        (pih,ih(1)),(pjh,jh(1)),
     :        (pcoeff,coeff(1)),(pnijptr,nijptr(1)),
     :        (pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PVALUE,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR

*     .. Logical array passed from nonh to give the used integrals
      POINTER (plused,lused(1))
      COMMON/USED/plused
      LOGICAL lused

      LOGICAL change
*
        len = intptr(6)
        call dfill (len, 0.d0, value, 1)
        IBEGIN =  1
        IEND = INTPTR(3)
        DO 1 I = IBEGIN,IEND
          if (lused(i)) then
           call unpacki(1,i,kv,iel1,iel2,iel3,iel4)
              IF (I .LE. INTPTR(1)) THEN
                 VALUE(I) = FK(IEL1,IEL2,KV,REL)
              ELSE IF (I .LE. INTPTR(2)) THEN
                 VALUE(I) = GK(IEL1,IEL2,KV,REL)
              ELSE
                 VALUE(I) = QUADR(IEL1,IEL2,0)**KV
              END IF
          end if
  1     CONTINUE
*
        IBEGIN = IEND + 1
        IEND = INTPTR(4)
        DO 30 I = IBEGIN,IEND
          if (lused(i)) then
           call unpacki(4,i,kv,iel1,iel2,iel3,iel4)
              K1 = KVAL(I)/16
              K2 = KVAL(I) - 16*K1
              VALUE(I) = QUADR(IEL1,IEL2,0)**K1
     :                  *QUADR(IEL3,IEL4,0)**K2
          end if
 30     CONTINUE
        IBEGIN = IEND + 1
        IEND = INTPTR(5)
        DO 10 I = IBEGIN,IEND
          if (lused(i)) then
           call unpacki(5,i,kv,iel1,iel2,iel3,iel4)
           VALUE(I) = RK(IEL1,IEL2,IEL3,IEL4,KV,REL)
          end if
 10     CONTINUE
*
        IBEGIN = IEND + 1
        IEND = INTPTR(6)
        DO 20 I = IBEGIN,IEND
          if (lused(i)) then
           call unpacki(6,i,kv,iel1,iel2,iel3,iel4)
           IF (VARIED(IEL1) .OR. VARIED(IEL2))
     :              VALUE(I) = HLC(EL,IEL1,IEL2,REL)
          end if
 20     CONTINUE
*       print *, 'Values:'
*       print '(3(I6,F18.15))', (i,value(i),i=1,intptr(6))

*      ... Test if any of the core functions have changed
*
         CHANGE = .FALSE.
         DO 35 I = 1,NCLOSD
           CHANGE = CHANGE .OR. VARIED(I)
  35     CONTINUE
         IF (CHANGE .OR.  EC.EQ.D0) CALL ECORE(EL,EC,REL)
*
        RETURN
        END




