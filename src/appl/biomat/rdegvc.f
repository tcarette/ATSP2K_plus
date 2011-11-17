*     ------------------------------------------------------------
*     R D E G V C
*     ------------------------------------------------------------
*
* --- determine the number of different J-values (njv)
*               the total number of eigenvectors (nvc)
*               the total length of the wt vector to allocate (lgth)
*
      SUBROUTINE rdegvc
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL rel,vok
      COMMON /EMS/IEM(4),IFL,JI,JF,LAM,REL,VOK
      COMMON /INOUT/IREAD,IWRITE,ISCW,iuc(2),iuw(2),iul(2),iuj(2),iut(2)
      COMMON /STATE/QET1,QLBL1,QWT1,QJV1,QET2,QLBL2,QWT2,QJV2,NJV(2),
     : NVC(2),LGTH(2),NCF(2),QCFG1,QCFG2
      POINTER(QET1,ET1(1)),(QLBL1,LBL1(1)),(QWT1,WT1(1)),(QJV1,JV1(1)),
     :       (QET2,ET2(1)),(QLBL2,LBL2(1)),(QWT2,WT2(1)),(QJV2,JV2(1)),
     :       (QCFG1,CFG1(8,1)),(QCFG2,CFG2(8,1))
*
   10 FORMAT(//8X,I4,10X,I4)
*
Cdbg  print*,' ncf(1) = ',ncf(1),' ncf(2) = ',ncf(2)
      do 1 ii = 1,2
        if (rel) then
          iu=iuj(ii)
        else
          iu=iul(ii)
        end if
        nvc(ii) = 0
        njv(ii) = 0
        nl=(ncf(ii)-1)/7+1
*
        read(iu,'(A)') 
    2   read(iu,10,end=5) jv,nb
        njv(ii) = njv(ii)  + 1
        nvc(ii) = nvc(ii) + nb
        do 3 m = 1,nb
           read(iu,'(A)') 
           read(iu,'(A)') 
           do 4 k = 1,nl
             read(iu,'(A)') 
    4      continue
    3   continue
        go to 2
    5   lgth(ii) = nvc(ii)*ncf(ii)
Cdbg  print*,' ii = ',ii,' nvc = ',nvc(ii),'ncf = ',ncf(ii)
Cdbg  print*,'             lgth = ',lgth(ii),' njv = ',njv(ii)
        rewind(unit=iu)
    1 continue
      return
      end
