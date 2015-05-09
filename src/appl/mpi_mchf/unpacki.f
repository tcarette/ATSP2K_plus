*-----------------------------------------------------------------------
*               U N P A C K I
*-----------------------------------------------------------------------
*
*     Given the type of integral, the index of the integral,
*     determine -- kv, iel1, iel2, iel3, iel4
*
      SUBROUTINE unpacki(int,i,kv,iel1,iel2,iel3,iel4)
        IMPLICIT double precision(A-H,O-Z)
        PARAMETER (NWD=94,NOD=220,NOFFD=800)
*
      COMMON /PARAM/H,H1,H3,CH,EH,RHO,Z,TOL,NO,ND,NWF,MASS,NCFG,IB,IC,ID
     :   ,D0,D1,D2,D3,D4,D5,D6,D8,D10,D12,D16,D30,FINE,NSCF,NCLOSD,RMASS
*
      POINTER (pkval,kval(1)),(pvalue,value(1)),
     :        (pih,ih(1)),(pjh,jh(1)),
     :        (pcoeff,coeff(1)),(pnijptr,nijptr(1)),
     :        (pinptr,inptr(1))
      COMMON/STATE/INTPTR(6),IDIM,PKVAL,PVALUE,NZE,PIH,PJH,
     :             NCODIM,PCOEFF,PNIJPTR,PINPTR
*
*      ... unpack the electron data
*
ctc  NWD limited to 63 (pack root =64). Change the pack root to 95
          kv = kval(i)
          if (int .le. 3 .or. int .eq. 6) then
              iel2 = mod(kv,95) + nclosd
              kv = kv/95
              iel1 = mod(kv,95) + nclosd
              kv = kv/95
          else if (int .eq. 4 .or. int .eq. 5) then
              iel4 = mod(kv,95) + nclosd
              kv = kv/95
              iel3 = mod(kv,95) + nclosd
              kv = kv/95
              iel2 = mod(kv,95) + nclosd
              kv = kv/95
              iel1 = mod(kv,95) + nclosd
              kv = kv/95
          end if
          end
