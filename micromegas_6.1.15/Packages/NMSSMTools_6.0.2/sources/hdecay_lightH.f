      SUBROUTINE HDECAY_lightH(PAR)

      IMPLICIT NONE

      INTEGER K,L,N0,VFLAG,INDH,NF

      DOUBLE PRECISION PAR(*),MH,EPS,PI,aux,HF,H1,H2,SINB,COSB,GHCC
      DOUBLE PRECISION BETA,X,SP,QCD0,HQCDM,HQCD,QCDH,TQCDH,RATCOUP
      DOUBLE PRECISION FPI,THETAE,MPI,MPIC,MPI8,MPI9,META,METAP,MKc,MK0
      DOUBLE PRECISION CRUL,CRDL,CRUR,CRDR,CRLL,CRLR,CRULR,CRDLR,CRLLR
      DOUBLE PRECISION GHCHACHA(2,2),MS2,GamS2,MS3,GamS3,MSIG,GamSIG
      DOUBLE PRECISION AHee,AHmumu,AHtata,MG0,GamG0
      DOUBLE PRECISION RMS,ASH,ASH2,AS3,AS4,HIGTOP,ALPHAS,RUNM
      DOUBLE PRECISION AH2Pi,AHPiE,GamH2Pi,GamH2PiC,GamHPiE,GamHPiEP
      DOUBLE PRECISION AH2K,GamH2KC,GamH2K0,AH2E,GamH2E,GamHEEP,GamH2EP
      DOUBLE PRECISION AHuu,AHdd,AHss,GamHuu,GamHdd,GamHss,GamHjj
      DOUBLE PRECISION MD,DCC,RMC,MBm,DBB,RMB,runmb,HTWW,HTZZ
      DOUBLE PRECISION MCHIC,DMC,CMIX,MCHIB1,MCHIB2,DMB1,DMB2
      DOUBLE PRECISION MMIX(3,3),VALP(3),VECP(3,3),OMIX(3,3)
      DOUBLE PRECISION ALEM0
      DOUBLE PRECISION SMASS(3),S(3,3),PMASS(2),P2(2,2),CMASS
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),N(5,5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5),CL(5)
      DOUBLE PRECISION LAMBDA,KAPPA,ALAMBDA,AKAPPA,MUEFF,NUEFF
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MCC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION MPI0,MEL,MSTRANGE
      DOUBLE PRECISION XLAMBDA,MC0,MB0,MT0
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      DOUBLE PRECISION HUQ,HDQ,MTQ,MBQ
      DOUBLE PRECISION XIF,XIS,MUP,MSP,M3H
      DOUBLE PRECISION QSTSB
      DOUBLE PRECISION GamHGAGA,GamHee,GamHmumu,GamHtata,GamHhadr,
     . GamHcc,GamHbb,GamHinv(3),GamHWW,GamHZZ,GamHAA
      DOUBLE PRECISION DDCOS,DDSIN
* New May 2019:
      DOUBLE PRECISION ZETA2,ZETA3,ASG,HGGQCD2,CIH
* New June 2022:
      DOUBLE PRECISION AbsThetaPiPiFit,phiPiPiFit,AbsThetaKKFit
      DOUBLE PRECISION phiKKFit,RFqqPiPiFit,IFqqPiPiFit,RFssPiPiFit
      DOUBLE PRECISION IFssPiPiFit,RFqqPiEtaFit,IFqqPiEtaFit
      DOUBLE PRECISION RFqqEtaEtaFit,IFqqEtaEtaFit,RFssEtaEtaFit
      DOUBLE PRECISION IFssEtaEtaFit,RFqqKKFit,IFqqKKFit
      DOUBLE PRECISION RFudKKFit,IFudKKFit,RFssKKFit,IFssKKFit
* End New

      DOUBLE COMPLEX XC,TC,BC,CC,LC,MC,EC,CH1C,CH2C,WC,HC,PropI,PropII
      DOUBLE COMPLEX ULC,URC,DLC,DRC,T1C,T2C,B1C,B2C,LLC,LRC,L1C,L2C
      DOUBLE COMPLEX CJH,CGH,F0,FF,FS,FV,CM,CJHF

      COMMON/ALEM0/ALEM0
      COMMON/HIGGSPEC/SMASS,S,PMASS,P2,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,N
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/CB/CB,CL
      COMMON/QPAR/LAMBDA,KAPPA,ALAMBDA,AKAPPA,MUEFF,NUEFF
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MCC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/SMEXT/MPI0,MEL,MSTRANGE
      COMMON/ALS/XLAMBDA,MC0,MB0,MT0,N0
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     .      MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     .      CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/QQUARK/HUQ,HDQ,MTQ,MBQ
      COMMON/QEXT/XIF,XIS,MUP,MSP,M3H
      COMMON/STSBSCALE/QSTSB
      COMMON/VFLAG/VFLAG
      COMMON/LIGHTHDECAYS/GamHGAGA,GamHee,GamHmumu,GamHtata,GamHhadr,
     . GamHcc,GamHbb,GamHinv,GamHWW,GamHZZ,GamHAA

      CM(X)= DCMPLX(MIN(1d3,X)**2,-EPS/4d0)
      F0(XC)=-CDLOG((CDSQRT(1d0-4d0*XC)-1d0)
     .               /(CDSQRT(1d0-4d0*XC)+1d0))**2/4d0
      FF(XC)=8d0*XC*(1d0+(1d0-4d0*XC)*F0(XC))
      FS(XC)=2d0*XC*(4d0*XC*F0(XC)-1d0)
      FV(XC)=-(2d0+12d0*XC+24d0*XC*(1d0-2d0*XC)*F0(XC))
      CRUL(HF)= dsqrt(2d0)*(HF**2*H1Q*S(1,1)
     . + (g1/12d0-g2/4d0)*(H1Q*S(1,1)-H2Q*S(1,2)))
      CRDL(HF)= dsqrt(2d0)*(HF**2*H2Q*S(1,2)
     . + (g1/12d0+g2/4d0)*(H1Q*S(1,1)-H2Q*S(1,2)))
      CRUR(HF)= dsqrt(2d0)*(HF**2*H1Q*S(1,1)
     . - g1/3d0*(H1Q*S(1,1)-H2Q*S(1,2)))
      CRDR(HF)= dsqrt(2d0)*(HF**2*H2Q*S(1,2)
     . + g1/6d0*(H1Q*S(1,1)-H2Q*S(1,2)))
      CRLL(HF)= dsqrt(2d0)*(HF**2*H2Q*S(1,2)
     .      + (-g1/4d0+g2/4d0)*(H1Q*S(1,1)-H2Q*S(1,2)))
      CRLR(HF)= dsqrt(2d0)*(HF**2*H2Q*S(1,2)
     .                + g1/2d0*(H1Q*S(1,1)-H2Q*S(1,2)))
      BETA(X)= DSQRT(1d0-4d0*X)
      QCd0(X)= (1d0+X**2)*(4d0*SP((1d0-X)/(1d0+X))
     . +2d0*SP((X-1d0)/(X+1d0))
     . - 3d0*DLOG((1d0+X)/(1d0-X))*DLOG(2d0/(1d0+X))
     . - 2d0*DLOG((1d0+X)/(1d0-X))*DLOG(X))
     . - 3d0*X*DLOG(4d0/(1d0-X**2))-4d0*X*DLOG(X)
      HQCDM(X)= QCd0(X)/X+(3d0+34d0*X**2-13d0*X**4)/16d0/X**3
     . * DLOG((1d0+X)/(1d0-X))+3d0/8d0/X**2*(7d0*X**2-1d0)
* July 2010:
c      HQCD(X)=(4d0/3d0*HQCDM(BETA(X))
c     .   +2d0*(4d0/3d0-DLOG(X))*(1d0-10d0*X)/(1d0-4d0*X))*ASH/PI
c     .   + (29.14671d0 + RATCOUP*(1.570d0 - 2d0*DLOG(HIGTOP)/3d0
c     .   + DLOG(X)**2/9d0))*(ASH/PI)**2
c     .   +(164.14d0 - 25.77d0*5 + 0.259d0*5**2)*(ASH/PI)**3
* New May 2019:
      HQCD(X)=(4d0/3d0*HQCDM(BETA(X))
     .       + 2d0*(4d0/3d0-DLOG(X))*(1d0-10d0*X)/(1d0-4d0*X))*ASH/PI
     .       + (29.14671d0
     .         + RATCOUP*(1.570d0 - 2d0*DLOG(HIGTOP)/3d0
     .            + DLOG(X)**2/9d0))*(ASH/PI)**2
     .       + (164.14d0 - 25.77d0*5d0 + 0.259d0*5d0**2)*(ASH/PI)**3
     .       + (39.34d0-220.9d0*5d0+9.685d0*5d0**2
     .         - 0.0205d0*5d0**3)*(ASH/PI)**4
* End New
      QCDH(X)= 1d0+HQCD(X)
      TQCDH(X)= 1d0+4d0/3d0*HQCDM(BETA(X))*ASH/PI
* New May 2019:
      HGGQCD2(ASG,NF,MH,MT)= 1d0+ASG/PI*(95d0/4d0-NF*7d0/6d0)
     . +(ASG/PI)**2*(149533d0/288d0-363d0/8d0*ZETA2-495d0/8d0*ZETA3
     .              +19d0/8d0*DLOG(MH**2/MT**2)
     . +NF*(-4157d0/72d0+11d0/2d0*ZETA2+5d0/4d0*ZETA3
     . +2d0/3d0*DLOG(MH**2/MT**2))
     . +NF**2*(127d0/108d0-1d0/6d0*ZETA2))+(ASG/PI)**3
     . *(467.683620788d0+122.440972222d0*DLOG(MH**2/MT**2)
     .              +10.9409722222d0*DLOG(MH**2/MT**2)**2)
* End New

      EPS=1d-8
      PI=4d0*DATAN(1d0)
* New May 2019:
      ZETA2=PI**2/6d0
      ZETA3=1.202056903159594d0
* End New

      MH=SMASS(1)
      COSB=1d0/DSQRT(1d0+PAR(3)**2)
      SINB=PAR(3)*COSB
      H1=SINB/DSQRT(2d0*dsqrt(2d0)*GF)
      H2=COSB/DSQRT(2d0*dsqrt(2d0)*GF)

      FPI=0.093d0
      THETAE=-13d0*Pi/180d0
      MPI=0.135d0
      MPIC=0.1396d0
      META=0.548d0
      METAP=0.958d0
      MKc=0.494d0
      MK0=0.498d0
      MPI8=dsqrt(META**2*DDCOS(THETAE)**2+METAP**2*DDSIN(THETAE)**2)
      MPI9=dsqrt(META**2*DDSIN(THETAE)**2+METAP**2*DDCOS(THETAE)**2)

C       Coupling to photons / gluons

      TC=CM(MT/MH)
      BC=CM(MBP/MH)
      CC=CM(MCC/MH)
      LC=CM(MTAU/MH)
      MC=CM(MMUON/MH)
      EC=CM(MEL/MH)
      CH1C=CM(MCHA(1)/MH)
      CH2C=CM(MCHA(2)/MH)
      WC=CM(MW/MH)
      HC=CM(CMASS/MH)
      ULC=CM(MUL/MH)
      URC=CM(MUR/MH)
      DLC=CM(MDL/MH)
      DRC=CM(MDR/MH)
      T1C=CM(MST1/MH)
      T2C=CM(MST2/MH)
      B1C=CM(MSB1/MH)
      B2C=CM(MSB2/MH)
      LLC=CM(MLL/MH)
      LRC=CM(MLR/MH)
      L1C=CM(MSL1/MH)
      L2C=CM(MSL2/MH)

      CRULR=HUQ/dsqrt(2d0)*(PAR(12)*S(1,1)-MUEFF*S(1,2)
     .                            -LAMBDA*H2Q*S(1,3))
      CRDLR=HDQ/dsqrt(2d0)*(-MUEFF*S(1,1)+PAR(13)*S(1,2)
     .                            -LAMBDA*H1Q*S(1,3))
      CRLLR= MTAU/H2Q/dsqrt(2d0)*(-MUEFF*S(1,1)+PAR(14)*S(1,2)
     .                                      - LAMBDA*H1Q*S(1,3))

      GHCC=LAMBDA*dsqrt(2d0)*(MUEFF*S(1,3)
     .                 -LAMBDA*SINB*COSB*(H1*S(1,2)+H2*S(1,1)))
     . +MUEFF*KAPPA*2d0*dsqrt(2d0)*S(1,3)*SINB*COSB
     . +LAMBDA*ALAMBDA*dsqrt(2d0)*S(1,3)*SINB*COSB
     . +g1/(2d0*dsqrt(2d0))*(H1*S(1,1)*(COSB**2-SINB**2)
     .                      +H2*S(1,2)*(SINB**2-COSB**2))
     . +g2/(2d0*dsqrt(2d0))*(H1*(S(1,1)+2d0*S(1,2)*SINB*COSB)
     .                      +H2*(S(1,2)+2d0*S(1,1)*SINB*COSB))
     . +LAMBDA*MUP*S(1,3)*SINB*COSB*dsqrt(2d0)
     . +6d0*dsqrt(2d0)*DLOG(MAX(QSTSB,MH**2)/(MAX(MT,MH)**2))
     .     *(MT**4/H1**3*S(1,1)*COSB**2+MB**4/H2**3*S(1,2)*SINB**2
     .      +MT**2*MB**2/(H1**2*H2**2)
     .        *(H1*S(1,1)*SINB**2+H2*S(1,1)*SINB*COSB
     .         +H2*S(1,2)*COSB**2+H1*S(1,2)*SINB*COSB))/(16d0*PI**2)

      DO K=1,2
      DO L=1,2
       GHCHACHA(K,L)=LAMBDA/dsqrt(2d0)*S(1,3)*U(K,2)*V(L,2)
     .    +dsqrt(g2/2d0)*(S(1,1)*U(K,1)*V(L,2)+S(1,2)*U(K,2)*V(L,1))
      ENDDO
      ENDDO

      CJHF=DSQRT(dsqrt(2d0)*GF)/4d0*(CU(1)*(FF(TC)+FF(CC))+CB(1)*FF(BC))

      CJH=CJHF
     .  +(CRUL(0d0)*FS(ULC)/MUL**2+CRUR(0d0)*FS(URC)/MUR**2
     .   +CRDL(0d0)*FS(DLC)/MDL**2+CRDR(0d0)*FS(DRC)/MDR**2)/2d0
     .  +(CST**2*CRUL(HUQ)+(1d0-CST**2)*CRUR(HUQ)+2d0*CST
     .                  *DSQRT(1d0-CST**2)*CRULR)*FS(T1C)/(4d0*MST1**2)
     .  +((1d0-CST**2)*CRUL(HUQ)+CST**2*CRUR(HUQ)-2d0*CST
     .                  *DSQRT(1d0-CST**2)*CRULR)*FS(T2C)/(4d0*MST2**2)
     .  +(CSB**2*CRDL(HDQ)+(1d0-CSB**2)*CRDR(HDQ)+2d0*CSB
     .                  *DSQRT(1d0-CSB**2)*CRDLR)*FS(B1C)/(4d0*MSB1**2)
     .  +((1d0-CSB**2)*CRDL(HDQ)+CSB**2*CRDR(HDQ)-2d0*CSB
     .                  *DSQRT(1d0-CSB**2)*CRDLR)*FS(B2C)/(4d0*MSB2**2)

* New May 2019:
      CIH=DREAL(DCONJG(CJH)*(CJH-CJHF))
* End New

      CGH=DSQRT(dsqrt(2d0)*GF)/2d0*(4d0/3d0*CU(1)*(FF(TC)+FF(CC))
     .    +CB(1)*FF(BC)/3d0+CD(1)*(FF(LC)+FF(MC)+FF(EC))+CV(1)*FV(WC))
     .  +GHCC/CMASS**2*FS(HC)/2d0
     .  +GHCHACHA(1,1)/MCHA(1)*FF(CH1C)/2d0
     .  +GHCHACHA(2,2)/MCHA(2)*FF(CH2C)/2d0
     .  +2d0/3d0*((CST**2*CRUL(HUQ)+(1d0-CST**2)*CRUR(HUQ)+2d0*CST
     .                     *DSQRT(1d0-CST**2)*CRULR)*FS(T1C)/MST1**2
     .           +((1d0-CST**2)*CRUL(HUQ)+CST**2*CRUR(HUQ)-2d0*CST
     .                     *DSQRT(1d0-CST**2)*CRULR)*FS(T2C)/MST2**2
     .    +2d0*CRUL(0d0)*FS(ULC)/MUL**2+2d0*CRUR(0d0)*FS(URC)/MUR**2)
     .  +1d0/6d0*((CSB**2*CRDL(HDQ)+(1d0-CSB**2)*CRDR(HDQ)+2d0*CSB
     .                     *DSQRT(1d0-CSB**2)*CRDLR)*FS(B1C)/MSB1**2
     .           +((1d0-CSB**2)*CRDL(HDQ)+CSB**2*CRDR(HDQ)-2d0*CSB
     .                     *DSQRT(1d0-CSB**2)*CRDLR)*FS(B2C)/MSB2**2
     .    +2d0*CRDL(0d0)*FS(DLC)/MDL**2+2d0*CRDR(0d0)*FS(DRC)/MDR**2)
     .  +CRLL(0d0)*FS(LLC)/MLL**2+CRLR(0d0)*FS(LRC)/MLR**2
     .  +(CSL**2*CRLL(MTAU/H2Q)+(1d0-CSL**2)*CRLR(MTAU/H2Q)+2d0*CSL
     .                   *DSQRT(1d0-CSL**2)*CRLLR)*FS(L1C)/MSL1**2/2d0
     .  +((1d0-CSL**2)*CRLL(MTAU/H2Q)+CSL**2*CRLR(MTAU/H2Q)-2d0*CSL
     .                   *DSQRT(1d0-CSL**2)*CRLLR)*FS(L2C)/MSL2**2/2d0

C       Leptonic decays

C  * H -> ee
      AHee=MEL*CD(1)*DSQRT(dsqrt(2d0)*GF)
      GamHee=AHee**2/(8d0*Pi)*MH
     .               *dsqrt(max(0d0,1d0-4d0*MEL**2/MH**2))**3

C  * H -> mumu
      AHmumu=MMUON*CD(1)*DSQRT(dsqrt(2d0)*GF)
      GamHmumu=AHmumu**2/(8d0*Pi)*MH
     .               *dsqrt(max(0d0,1d0-4d0*MMUON**2/MH**2))**3

C  * H -> tautau
      AHtata=MTAU*CL(1)*DSQRT(dsqrt(2d0)*GF)
      GamHtata=AHtata**2/(8d0*Pi)*MH
     .               *dsqrt(max(0d0,1d0-4d0*MTAU**2/MH**2))**3

c       Decay to light neutralinos

      IF(MH.LE.2d0*DABS(MNEU(1)))THEN
       GamHinv(1)=0d0
      ELSE
       aux=dsqrt(2d0)*LAMBDA*(S(1,1)*N(1,4)*N(1,5)
     .                       +S(1,2)*N(1,3)*N(1,5)
     .                       +S(1,3)*N(1,3)*N(1,4))
     .    -dsqrt(2d0)*KAPPA*S(1,3)*N(1,5)*N(1,5)
     .    +dsqrt(g1)*(-S(1,1)*N(1,1)*N(1,3)+S(1,2)*N(1,1)*N(1,4))
     .    +dsqrt(g2)*(S(1,1)*N(1,2)*N(1,3)-S(1,2)*N(1,2)*N(1,4))

       GamHinv(1)=aux**2/(16d0*PI)*MH*dsqrt(1d0-4d0*(MNEU(1)/MH)**2)**3
      ENDIF

      IF(MH.LE.DABS(MNEU(1))+DABS(MNEU(2)))THEN
       GamHinv(2)=0d0
      ELSE
       aux=LAMBDA*(S(1,1)*(N(1,4)*N(2,5)+N(2,4)*N(1,5))
     .            +S(1,2)*(N(1,3)*N(2,5)+N(2,3)*N(1,5))
     .            +S(1,3)*(N(1,3)*N(2,4)+N(2,3)*N(1,4)))/dsqrt(2d0)
     .    -dsqrt(2d0)*KAPPA*S(1,3)*N(1,5)*N(2,5)
     .    +dsqrt(g1)/2d0*(-S(1,1)*(N(1,1)*N(2,3)+N(2,1)*N(1,3))
     .                    +S(1,2)*(N(1,1)*N(2,4)+N(2,1)*N(1,4)))
     .    +dsqrt(g2)/2d0*(S(1,1)*(N(1,2)*N(2,3)+N(2,2)*N(1,3))
     .                   -S(1,2)*(N(1,2)*N(2,4)+N(2,2)*N(1,4)))

       GamHinv(2)=aux**2/(8d0*PI)*MH
     .               *dsqrt(1d0-((MNEU(1)+MNEU(2))/MH)**2)**3
     .               *dsqrt(1d0-((MNEU(1)-MNEU(2))/MH)**2)
      ENDIF

      IF(MH.LE.2d0*DABS(MNEU(2)))THEN
       GamHinv(3)=0d0
      ELSE
       aux=dsqrt(2d0)*LAMBDA*(S(1,1)*N(2,4)*N(2,5)
     .                       +S(1,2)*N(2,3)*N(2,5)
     .                       +S(1,3)*N(2,3)*N(2,4))
     .    -dsqrt(2d0)*KAPPA*S(1,3)*N(2,5)*N(2,5)
     .    +dsqrt(g1)*(-S(1,1)*N(2,1)*N(2,3)+S(1,2)*N(2,1)*N(2,4))
     .    +dsqrt(g2)*(S(1,1)*N(2,2)*N(2,3)-S(1,2)*N(2,2)*N(2,4))

       GamHinv(3)=aux**2/(16d0*PI)*MH*dsqrt(1d0-4d0*(MNEU(2)/MH)**2)**3
      ENDIF

c       Decay to light pseudoscalars

       IF(MH.LE.2d0*PMASS(1))THEN
        GamHAA=0d0
       ELSE
      aux=(g1+g2)/(2d0*dsqrt(2d0))
     . *P2(1,1)**2*(H1*S(1,1)*(COSB**2-SINB**2)
     .             +H2*S(1,2)*(SINB**2-COSB**2))
     . +LAMBDA*ALAMBDA*dsqrt(2d0)*(S(1,1)*P2(1,1)*P2(1,2)*SINB
     .   +S(1,2)*P2(1,1)*P2(1,2)*COSB+S(1,3)*P2(1,1)**2*SINB*COSB)
     . -KAPPA*AKAPPA*dsqrt(2d0)*S(1,3)*P2(1,2)**2
     . +LAMBDA**2*dsqrt(2d0)*(H1*S(1,1)*(P2(1,1)**2*SINB**2+P2(1,2)**2)
     .    +H2*S(1,2)*(P2(1,1)**2*COSB**2+P2(1,2)**2)
     .    +MUEFF/LAMBDA*S(1,3)*P2(1,1)**2)
     . +KAPPA**2*2d0*dsqrt(2d0)*MUEFF/LAMBDA*S(1,3)*P2(1,2)**2
     . +LAMBDA*KAPPA*dsqrt(2d0)
     .   *(H1*(S(1,2)*P2(1,2)**2-2d0*S(1,3)*P2(1,1)*P2(1,2)*SINB)
     .    +H2*(S(1,1)*P2(1,2)**2-2d0*S(1,3)*P2(1,1)*P2(1,2)*COSB)
     .    +2d0*MUEFF/LAMBDA*(S(1,3)*P2(1,1)**2*SINB*COSB
     .     -S(1,1)*P2(1,1)*P2(1,2)*SINB-S(1,2)*P2(1,1)*P2(1,2)*COSB))
     . +LAMBDA*MUP*dsqrt(2d0)*(-S(1,1)*P2(1,1)*P2(1,2)*SINB
     .      -S(1,2)*P2(1,1)*P2(1,2)*COSB+S(1,3)*P2(1,1)**2*SINB*COSB)
     . +KAPPA*MUP*S(1,3)*P2(1,2)**2*dsqrt(2d0)
     . +6d0*dsqrt(2d0)*DLOG(MAX(QSTSB,MH**2)/(MAX(MT,MH)**2))
     .    *(MT**4/H1**3*S(1,1)*P2(1,1)**2*COSB**2
     .     +MB**4/H2**3*S(1,2)*P2(1,1)**2*SINB**2)/(16d0*PI**2)
        GamHAA=aux**2/(32d0*PI*MH)*dsqrt(1d0-4d0*(PMASS(1)/MH)**2)
       ENDIF

c       Decay to W*W*/Z*Z*

      IF(VFLAG.ne.0)then

        CALL HTOVV(MW,2.08856d0,MH,HTWW)
        GamHWW= 3d0/2d0*GF*MW**4/DSQRT(2d0)/PI/MH**3*HTWW*CV(1)**2

        CALL HTOVV(MZ,2.49581d0,MH,HTZZ)
        GamHZZ= 3d0/4d0*GF*MZ**4/DSQRT(2d0)/PI/MH**3*HTZZ*CV(1)**2

      ELSE

       GamHWW=0d0
       GamHZZ=0d0

      ENDIF

C       Hadronic decays

      GamHhadr=0d0

C   Initializing the strong coupling constant and running masses
      IF(MH.ge.1.5d0)then
       HIGTOP=(MH/MT)**2
       MT0=3d8
* New May 2019:
       ASH=ALPHAS(MH,3)
       MC0=1d8
       MB0=2d8
       AS3=ALPHAS(MH,3)
       ASH2=AS3*dsqrt(1d0+AS3/Pi*(95d0/4d0-3d0*7d0/6d0))
       MC0=MCC
       AS4=ALPHAS(MH,3)
* End New
       MB0=MBP
       MT0=MT
       RMS=RUNM(MH,3)
      ELSE
       IF(MH.lt.0.8d0)ASH2=Pi/3d0
       IF(MH.ge.0.8d0.and.MH.lt.2d0)ASH2=
     . (Pi/3d0-0.458d0)*dexp(-((MH-0.8d0)/0.732d0)**2)+0.458d0
       IF(MH.gt.1d0)RMS=RUNM(MH,3)
      ENDIF

      IF(MH.lt.4d0)then

      MS2=1d0
      GamS2=0.18d0

      MS3=1.5d0
      GamS3=0.1d0

      MSIG=1d0
      GamSIG=0.7d0

      MG0=1.6d0
      GamG0=0.06d0

C  * H -> 2Pi0
      AH2Pi=CDABS(-DSQRT(dsqrt(2d0)*GF)/2d0*PropII(MH,MS2,GamS2)
     .   *((MPI**2+MKC**2-MK0**2)*CU(1)+(MPI**2+MK0**2-MKC**2)*CD(1))
     .  -2d0*ASH2/PI/3d0*CJH*(MH**2+MPI**2)*PropII(MH,MSIG,GamSIG))**2

      GamH2Pi=AH2Pi/(32d0*PI*MH)*dsqrt(max(0d0,1d0-4d0*(MPI/MH)**2))

C  * H -> Pi+Pi-
      AH2Pi=CDABS(-DSQRT(dsqrt(2d0)*GF)/2d0*PropII(MH,MS2,GamS2)
     .   *((MPI**2+MKC**2-MK0**2)*CU(1)+(MPI**2+MK0**2-MKC**2)*CD(1))
     . -2d0*ASH2/PI/3d0*CJH*(MH**2+MPIC**2)*PropII(MH,MSIG,GamSIG))**2

      GamH2PiC=AH2Pi/(16d0*PI*MH)*dsqrt(max(0d0,1d0-4d0*(MPIC/MH)**2))
      
* New June 2022
      AH2Pi=CDABS(MPI**2/2d0*(2d0/9d0*CJH-CU(1)*dsqrt(dsqrt(2d0)*GF))
     c       *DCMPLX(-RFqqPiPiFit(MH),-IFqqPiPiFit(MH))/dsqrt(3d0) 
     c     +MPI**2/2d0*(2d0/9d0*CJH-CD(1)*dsqrt(dsqrt(2d0)*GF))
     c           *DCMPLX(-RFqqPiPiFit(MH),-IFqqPiPiFit(MH))/dsqrt(3d0) 
     c     +(MK0**2+MKc**2-MPI**2)/2d0*(2d0/9d0*CJH-
     c                   CD(1)*dsqrt(dsqrt(2d0)*GF))*dsqrt(2d0/3d0)
     c           *DCMPLX(-RFssPiPiFit(MH),-IFssPiPiFit(MH))
     c     -2d0/9d0*CJH*AbsThetaPiPiFit(MH)
     c          *DCMPLX(dcos(phiPiPiFit(MH)),dsin(phiPiPiFit(MH))))**2

      GamH2Pi=AH2Pi/(32d0*PI*MH)*dsqrt(max(0d0,1d0-4d0*(MPI/MH)**2))
      GamH2PiC=AH2Pi/(16d0*PI*MH)*dsqrt(max(0d0,1d0-4d0*(MPIC/MH)**2))
* End New

C  * H -> EtaPi0
      AHPiE=CDABS(-DSQRT(dsqrt(2d0)*GF/3d0)/2d0*PropII(MH,MS2,GamS2)
     .     *((MPI**2+MKC**2-MK0**2)*CU(1)-(MPI**2+MK0**2-MKC**2)*CD(1))
     .        -2d0*ASH2/PI/dsqrt(3d0)*CJH*PropII(MH,MSIG,GamSIG)
     .                   *(MKC**2-MK0**2-MPIC**2+MPI**2))**2

      GamHPiE=AHPiE/(16d0*PI*MH)
     .       *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))**2
     .       *dsqrt(max(0d0,1d0-(MPI+META)**2/MH**2))
     .       *dsqrt(max(0d0,1d0-(MPI-META)**2/MH**2))
      
* New June 2022
      AHPiE=CDABS(MPI**2/2d0*(2d0/9d0*CJH-CU(1)*dsqrt(dsqrt(2d0)*GF))
     c       *DCMPLX(RFqqPiEtaFit(MH),IFqqPiEtaFit(MH))/dsqrt(2d0) 
     c     -MPI**2/2d0*(2d0/9d0*CJH-CD(1)*dsqrt(dsqrt(2d0)*GF))
     c       *DCMPLX(RFqqPiEtaFit(MH),IFqqPiEtaFit(MH))/dsqrt(2d0))**2

      GamHPiE=AHPiE/(16d0*PI*MH)
     .       *dsqrt(max(0d0,1d0-(MPI+META)**2/MH**2))
     .       *dsqrt(max(0d0,1d0-(MPI-META)**2/MH**2))
* End New

C  * H -> Eta'Pi0
      AHPiE=CDABS(-DSQRT(dsqrt(2d0)*GF/3d0)/2d0*PropII(MH,MS2,GamS2)
     .     *((MPI**2+MKC**2-MK0**2)*CU(1)-(MPI**2+MK0**2-MKC**2)*CD(1))
     .        -2d0*ASH2/PI/dsqrt(3d0)*CJH*PropII(MH,MSIG,GamSIG)
     .                   *(MKC**2-MK0**2-MPIC**2+MPI**2))**2
     
* New June 2022
      AHPiE=0d0
* End New

      GamHPiEP=AHPiE/(16d0*PI*MH)
     .       *(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))**2
     .       *dsqrt(max(0d0,1d0-(MPI+METAP)**2/MH**2))
     .       *dsqrt(max(0d0,1d0-(MPI-METAP)**2/MH**2))

C  * H -> K+K-
      AH2K=CDABS(-DSQRT(dsqrt(2d0)*GF)/2d0
     .     *((MPI**2+MKC**2-MK0**2)*CU(1)*PropII(MH,MS2,GamS2)
     .      +(MKC**2+MK0**2-MPI**2)*CD(1)*PropII(MH,MS3,GamS3))
     . -2d0*ASH2/PI/3d0*CJH*(MH**2+MKC**2)*PropII(MH,MSIG,GamSIG))**2

* New June 2022
      AH2K=CDABS(MPI**2/2d0*(2d0/9d0*CJH-CU(1)*dsqrt(dsqrt(2d0)*GF))
     c           *DCMPLX(-RFqqKKFit(MH)/2d0+RFudKKFit(MH)/dsqrt(2d0),
     c                   -IFqqKKFit(MH)/2d0+IFudKKFit(MH)/dsqrt(2d0))
     c     +MPI**2/2d0*(2d0/9d0*CJH-CD(1)*dsqrt(dsqrt(2d0)*GF))
     c           *DCMPLX(-RFqqKKFit(MH)/2d0-RFudKKFit(MH)/dsqrt(2d0),
     c                   -IFqqKKFit(MH)/2d0-IFudKKFit(MH)/dsqrt(2d0)) 
     c     +(MK0**2+MKc**2-MPI**2)/2d0*(2d0/9d0*CJH-
     c                   CD(1)*dsqrt(dsqrt(2d0)*GF))/dsqrt(2d0)
     c           *DCMPLX(-RFssKKFit(MH),-IFssKKFit(MH))
     c     -2d0/9d0*CJH*AbsThetaKKFit(MH)
     c           *DCMPLX(dcos(phiKKFit(MH)),dsin(phiKKFit(MH))))**2
* End New

      GamH2KC=AH2K/(16d0*PI*MH)
     .       *dsqrt(max(0d0,1d0-4d0*MKC**2/MH**2))

C  * H -> K0K0b
      AH2K=CDABS(-DSQRT(dsqrt(2d0)*GF)/2d0
     .     *((MPI**2+MK0**2-MKC**2)*CD(1)*PropII(MH,MS2,GamS2)
     .      +(MKC**2+MK0**2-MPI**2)*CD(1)*PropII(MH,MS3,GamS3))
     . -2d0*ASH2/PI/3d0*CJH*(MH**2+MK0**2)*PropII(MH,MSIG,GamSIG))**2

* New June 2022
      AH2K=CDABS(MPI**2/2d0*(2d0/9d0*CJH-CU(1)*dsqrt(dsqrt(2d0)*GF))
     c           *DCMPLX(-RFqqKKFit(MH)/2d0-RFudKKFit(MH)/dsqrt(2d0),
     c                   -IFqqKKFit(MH)/2d0-IFudKKFit(MH)/dsqrt(2d0))
     c     +MPI**2/2d0*(2d0/9d0*CJH-CD(1)*dsqrt(dsqrt(2d0)*GF))
     c           *DCMPLX(-RFqqKKFit(MH)/2d0+RFudKKFit(MH)/dsqrt(2d0),
     c                   -IFqqKKFit(MH)/2d0+IFudKKFit(MH)/dsqrt(2d0)) 
     c     +(MK0**2+MKc**2-MPI**2)/2d0*(2d0/9d0*CJH-
     c                   CD(1)*dsqrt(dsqrt(2d0)*GF))/dsqrt(2d0)
     c           *DCMPLX(-RFssKKFit(MH),-IFssKKFit(MH))
     c     -2d0/9d0*CJH*AbsThetaKKFit(MH)
     c           *DCMPLX(dcos(phiKKFit(MH)),dsin(phiKKFit(MH))))**2
* End New

      GamH2K0=AH2K/(16d0*PI*MH)
     .       *dsqrt(max(0d0,1d0-4d0*MK0**2/MH**2))

C  * H -> 2Eta
      AH2E=CDABS(-DSQRT(dsqrt(2d0)*GF)/6d0*
     . (((MPI**2+MKC**2-MK0**2)*CU(1)+(MPI**2+MK0**2-MKC**2)*CD(1))
     . *PropII(MH,MS2,GamS2)*(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))**2
     .   +4d0*(MKC**2+MK0**2-MPI**2)*CD(1)*PropII(MH,MS3,GamS3)
     .              *(DDCOS(THETAE)+DDSIN(THETAE)/dsqrt(2d0))**2)
     .        -2d0*ASH2/PI/3d0*CJH*(PropII(MH,MSIG,GamSIG)
     .     *(MH**2-2d0*META**2+MPI**2*(DDCOS(THETAE)
     .       -dsqrt(2d0)*DDSIN(THETAE))**2+2d0*(MKC**2+MK0**2-MPI**2)
     .          *(DDCOS(THETAE)+DDSIN(THETAE)/dsqrt(2d0))**2)
     .     +PropI(MH,MG0,GamG0)*DDSIN(THETAE)**2))**2

* New June 2022
      AH2E=CDABS(MPI**2/2d0*(2d0/9d0*CJH-CU(1)*dsqrt(dsqrt(2d0)*GF))
     c           *DCMPLX(RFqqEtaEtaFit(MH),IFqqEtaEtaFit(MH))
     c     +MPI**2/2d0*(2d0/9d0*CJH-CD(1)*dsqrt(dsqrt(2d0)*GF))
     c           *DCMPLX(RFqqEtaEtaFit(MH),IFqqEtaEtaFit(MH)) 
     c     +(MK0**2+MKc**2-MPI**2)/2d0*(2d0/9d0*CJH-
     c                   CD(1)*dsqrt(dsqrt(2d0)*GF))*dsqrt(2d0)
     c           *DCMPLX(RFssEtaEtaFit(MH),IFssEtaEtaFit(MH))
     c     -2d0/9d0*CJH*(MH**2+2d0*META**2))**2
* End New

      GamH2E=AH2E/(32d0*PI*MH)*dsqrt(max(0d0,1d0-4d0*(META/MH)**2))

C  * H -> 2Eta'
      AH2E=CDABS(-DSQRT(dsqrt(2d0)*GF)/6d0*
     . (((MPI**2+MKC**2-MK0**2)*CU(1)+(MPI**2+MK0**2-MKC**2)*CD(1))
     . *PropII(MH,MS2,GamS2)*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))**2
     .   +4d0*(MKC**2+MK0**2-MPI**2)*CD(1)*PropII(MH,MS3,GamS3)
     .              *(DDSIN(THETAE)-DDCOS(THETAE)/dsqrt(2d0))**2)
     .        -2d0*ASH2/PI/3d0*CJH*(PropII(MH,MSIG,GamSIG)
     .    *(MH**2-2d0*METAP**2+MPI**2*(DDSIN(THETAE)
     .       +dsqrt(2d0)*DDCOS(THETAE))**2+2d0*(MKC**2+MK0**2-MPI**2)
     .          *(DDSIN(THETAE)-DDCOS(THETAE)/dsqrt(2d0))**2)
     .     +PropI(MH,MG0,GamG0)*DDCOS(THETAE)**2))**2
     
* New June 2022
      AH2E=0d0
* End New

      GamH2EP=AH2E/(32d0*PI*MH)*dsqrt(max(0d0,1d0-4d0*(METAP/MH)**2))

C  * H -> EtaEta'
      AH2E=CDABS(-DSQRT(dsqrt(2d0)*GF)/6d0*
     . (((MPI**2+MKC**2-MK0**2)*CU(1)+(MPI**2+MK0**2-MKC**2)*CD(1))
     .  *PropII(MH,MS2,GamS2)*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))
     .              *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
     .   +4d0*(MKC**2+MK0**2-MPI**2)*CD(1)*PropII(MH,MS3,GamS3)
     .              *(DDSIN(THETAE)-DDCOS(THETAE)/dsqrt(2d0))
     .              *(DDCOS(THETAE)+DDSIN(THETAE)/dsqrt(2d0)))
     .        -2d0*ASH2/PI/3d0*CJH*(PropII(MH,MSIG,GamSIG)*(
     .       MPI**2*(DDSIN(THETAE)+dsqrt(2d0)*DDCOS(THETAE))
     .             *(DDCOS(THETAE)-dsqrt(2d0)*DDSIN(THETAE))
     .       +2d0*(MKC**2+MK0**2-MPI**2)
     .          *(DDSIN(THETAE)-DDCOS(THETAE)/dsqrt(2d0))
     .          *(DDCOS(THETAE)+DDSIN(THETAE)/dsqrt(2d0)))
     .     +PropI(MH,MG0,GamG0)*DDSIN(THETAE)*DDCOS(THETAE)))**2

* New June 2022
      AH2E=0d0
* End New

      GamHEEP=AH2E/(16d0*PI*MH)
     .                  *dsqrt(max(0d0,1d0-((META+METAP)/MH)**2))
     .                  *dsqrt(max(0d0,1d0-((META-METAP)/MH)**2))

      GamHhadr=GamH2Pi+GamH2PiC+GamHPiE+GamHPiEP+GamH2KC+GamH2K0
     .         +GamH2E+GamH2EP+GamHEEP

C       Diphoton decays

      GamHGAGA=CDABS(ALEM0*CGH/Sqrt(2d0)/Pi
     .         +DSQRT(dsqrt(2d0)*GF)/2d0*PropI(MH,MS2,GamS2)
     .   *((MPI**2+MKC**2-MK0**2)*CU(1)+(MPI**2+MK0**2-MKC**2)*CD(1))
     .   *dsqrt(16d0*Pi*6.7d-7/MS2**3)/2.5d0
     .         +DSQRT(dsqrt(2d0)*GF)/2d0*PropI(MH,MS3,GamS3)
     .   *(MKC**2+MK0**2-MPI**2)*CD(1)*dsqrt(16d0*Pi*4d-7/MS3**3)/2.7d0
     . +2d0*ASH2/PI/3d0*CJH*
     .  (MSIG**2*PropI(MH,MSIG,GamSIG)*dsqrt(16d0*Pi*2d-7/MSIG**3)/5d0
     .   +MG0**2*PropI(MH,MG0,GamG0)*dsqrt(16d0*Pi*1d-6/MG0**3)/2d0)
     .          )**2*MH**3/(32d0*Pi)

      ENDIF

      aux=8d0*CDABS(CGH)**2*MH**3
     .          /(32d0*Pi)*(ALEM0/4d0/Pi)**2
      IF(MH.ge.3d0)then
       IF(MH.lt.4d0)THEN
        GamHGAGA=GamHGAGA*(1d0-((MH-3d0)/1d0)**2)+((MH-3d0)/1d0)**2*aux
       ELSE
        GamHGAGA=aux
       ENDIF
      ENDIF

C       Light-quark and gluon decays

      AHuu=2d-3*CU(1)*DSQRT(dsqrt(2d0)*GF)
      AHdd=4d-3*CD(1)*DSQRT(dsqrt(2d0)*GF)
      AHss=MS*CD(1)*DSQRT(dsqrt(2d0)*GF)

      IF(MH.LT.1.5d0)THEN
       GamHuu=0d0
       GamHdd=0d0
       GamHss=0d0
       GamHjj=0d0
      ELSE

       RATCOUP=1d0
       GamHuu=3d0*dabs(AHuu)**2/(8d0*Pi)*MH
     .        *(4d0*(2d-3/MH)**2*TQCDH((2d-3/MH)**2)
     .   +(1d0-4d0*(2d-3/MH)**2)*max(0d0,QCDH((2d-3/MH)**2)))
     .               *dsqrt(max(0d0,1d0-4d0*(2d-3/MH)**2))**3

       RATCOUP=0d0
       IF(CD(1).NE.0d0)RATCOUP=CU(1)/CD(1)
       GamHdd=3d0*dabs(AHdd)**2/(8d0*Pi)*MH
     .         *(4d0*(4d-3/MH)**2*TQCDH((4d-3/MH)**2)
     .   +(1d0-4d0*(4d-3/MH)**2)*max(0d0,QCDH((4d-3/MH)**2)))
     .               *dsqrt(max(0d0,1d0-4d0*(4d-3/MH)**2))**3
       GamHss=3d0*dabs(AHss)**2/(8d0*Pi)*MH
     .         *(4d0*(MS/MH)**2*TQCDH((MS/MH)**2)
     .   +(1d0-4d0*(MS/MH)**2)*max(0d0,(RMS/MS)**2*QCDH((RMS/MH)**2)))
     .               *dsqrt(max(0d0,1d0-4d0*(MS/MH)**2))**3
* New May 2019:
       GamHjj=AS3**2/(8d0*PI**3)*MH**3
     .  *Max(0d0,CDABS(CJH)**2*(1d0+AS3/Pi*(95d0/4d0-3d0*7d0/6d0))
     .           +CIH*AS3/PI*7d0/2d0)
* End New
      ENDIF
      
* New June 2022: (boundaries chiral region)
      aux=GamHhadr
      IF(MH.lt.1.5d0)then
       GamHhadr=aux
      ELSE
       IF(MH.lt.2d0)THEN
        GamHhadr=aux*(1d0-((MH-1.5d0)/0.5d0))
     .        +((MH-1.5d0)/0.5d0)*(0d0*GamHuu+0d0*GamHdd+GamHss+GamHjj)
       ELSE
        GamHhadr=0d0*GamHuu+0d0*GamHdd+GamHss+GamHjj
       ENDIF
      ENDIF
* End New

C       Mixing with the chi_c(1P)  - hep-ph/9503356

      MCHIC=3.41475d0
      DMC=dsqrt(27d0*dsqrt(2d0)/Pi*GF*MCHIC*0.1d0)*CU(1)    !0.007d0*CU(1)

c      IF(MH.gt.2d0.and.MH.lt.2d0*MTAU)THEN
       CMIX=0d0
       aux=Max(0d0,min(1d0,MH-1d0))
       IF(DMC.ne.0d0)THEN
        CMIX=datan(2d0*DMC
     .        /(MH**2-MCHIC**2+dsqrt((MH**2-MCHIC**2)**2+4d0*DMC**2)))
       ENDIF
       IF(DDCOS(CMIX)**2.gt.DDSIN(CMIX)**2)then
        GamHhadr=GamHhadr*DDCOS(CMIX)**2+10.5d-3*aux*DDSIN(CMIX)**2
        GamHee=GamHee*DDCOS(CMIX)**2
        GamHmumu=GamHmumu*DDCOS(CMIX)**2
        GamHtata=GamHtata*DDCOS(CMIX)**2
        GamHgaga=GamHgaga*DDCOS(CMIX)**2+10.5d-3*2.23d-4*aux
     .                                                 *DDSIN(CMIX)**2
        GamHinv=GamHinv*DDCOS(CMIX)**2
       ELSE
        GamHhadr=GamHhadr*DDSIN(CMIX)**2+10.5d-3*aux*DDCOS(CMIX)**2
        GamHee=GamHee*DDSIN(CMIX)**2
        GamHmumu=GamHmumu*DDSIN(CMIX)**2
        GamHtata=GamHtata*DDSIN(CMIX)**2
        GamHgaga=GamHgaga*DDSIN(CMIX)**2+10.5d-3*2.23d-4*aux
     .                                                 *DDCOS(CMIX)**2
        GamHinv=GamHinv*DDSIN(CMIX)**2
       ENDIF
c      ENDIF

C       H -> cc decays

      MD=1.865d0

      IF(MH.LE.2d0*MD)THEN
       GamHcc= 0d0
      ELSE
       RMC=RUNM(MH,4)
       RATCOUP= 1d0

* New July 2019:
       DCC=MH**3/(8d0*PI**3)*(AS4**2*Max(0d0,CDABS(CJH)**2*
     .     (1d0+AS4/Pi*(95d0/4d0-4d0*7d0/6d0))+CIH*AS4/PI*7d0/2d0)
     .   -AS3**2*Max(0d0,CDABS(CJH)**2*
     .    (1d0+AS3/Pi*(95d0/4d0-3d0*7d0/6d0))+CIH*AS3/PI*7d0/2d0))

       GamHcc=4d0*(MCC/MH)**2*
     .        3d0*GF*(MCC*CU(1))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .        *TQCDH((MCC/MH)**2)*dsqrt(1d0-4d0*MCC**2/MH**2)**3
     .         +(1d0-4d0*(MCC/MH)**2)*max(0d0,
     .         3d0*GF*(RMC*CU(1))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .         *QCDH((RMC/MH)**2)*dsqrt(1d0-4d0*RMC**2/MH**2)**3!+DCC
     .         )
* End new

       IF(MH.le.50d0)then
        GamHcc=GamHcc*((1d0-(MH-2d0*MD)/(50d0-2d0*MD))*
     .  (dsqrt(1d0-4d0*MD**2/MH**2)/dsqrt(1d0-4d0*MCC**2/MH**2))**3
     .  +(MH-2d0*MD)/(50d0-2d0*MD))
* New July 2019:
        GamHhadr=GamHhadr+DCC*
     .  ((1d0-4d0*(MCC/MH)**2)*(1d0-(MH-2d0*MD)/(50d0-2d0*MD))*
     .  (dsqrt(1d0-4d0*MD**2/MH**2)/dsqrt(1d0-4d0*MCC**2/MH**2))**3
     .  +(MH-2d0*MD)/(50d0-2d0*MD))
* End new
       ENDIF
      ENDIF

C       Mixing with the chi_b(1,2P)   - Drees, Hikasa (1990)

      MCHIB1=9.85944d0
      MCHIB2=10.2325d0

      DMB1=dsqrt(27d0*dsqrt(2d0)/Pi*GF*MCHIB1*1.7d0)*CB(1) ! 0.049d0*CB(1)
      DMB2=dsqrt(27d0*dsqrt(2d0)/Pi*GF*MCHIB2*2.0d0)*CB(1) ! 0.054d0*CB(1)

c      IF(MH.gt.4d0.and.MH.lt.2d0*5.2795d0)THEN

       MMIX(1,1)=MCHIB1**2
       MMIX(1,2)=0d0
       MMIX(1,3)=DMB1
       MMIX(2,1)=0d0
       MMIX(2,2)=MCHIB2**2
       MMIX(2,3)=DMB2
       MMIX(3,1)=DMB1
       MMIX(3,2)=DMB2
       MMIX(3,3)=MH**2

       CALL DIAGN(3,MMIX,VALP,VECP,EPS)
       CALL SORTN(3,VALP,VECP)
       DO K=1,3
        DO L=1,3
         OMIX(K,L)=VECP(L,K)
        ENDDO
       ENDDO

       INDH=1
       DO K=1,3
        IF(dabs(OMIX(INDH,3))**2.lt.dabs(OMIX(K,3))**2)INDH=K
       ENDDO
c chi_b widths from 1601.05093
       aux=Max(0d0,min(1d0,MH/2d0-1.5d0))
       GamHhadr=GamHhadr*OMIX(INDH,3)**2+aux*(2.03d-3*OMIX(INDH,1)**2
     .                                       +2.39d-3*OMIX(INDH,2)**2)
       GamHee=GamHee*OMIX(INDH,3)**2
       GamHmumu=GamHmumu*OMIX(INDH,3)**2
       GamHtata=GamHtata*OMIX(INDH,3)**2
       GamHgaga=GamHgaga*OMIX(INDH,3)**2+aux*(0.12d-3*OMIX(INDH,1)**2
     .                                       +0.14d-3*OMIX(INDH,2)**2)
       GamHcc=GamHcc*OMIX(INDH,3)**2
       GamHinv=GamHinv*OMIX(INDH,3)**2
! GamAhadr also contains contributions to the cc final state
c      ENDIF

C       H -> bb decays

      MBm=5.2795d0

      IF(MH.LE.2d0*MBm)THEN
       GamHbb= 0d0
      ELSE
       RMB=RUNMB(MH)
       IF(CB(1).ne.0d0)THEN
        RATCOUP= CU(1)/CB(1)
       ELSE
        RATCOUP=0d0
       ENDIF

* New July 2019:
       DBB=MH**3/(8d0*PI**3)*(ASH**2*Max(0d0,CDABS(CJH)**2*
     .    (1d0+ASH/Pi*(95d0/4d0-5d0*7d0/6d0))+CIH*ASH/PI*7d0/2d0)
     .   -AS4**2*Max(0d0,CDABS(CJH)**2*
     .    (1d0+AS4/Pi*(95d0/4d0-4d0*7d0/6d0))+CIH*AS4/PI*7d0/2d0))

       DBB=DBB+MH**3/(8d0*PI**3)*CDABS(CJHF)**2*ASH**2*(
     .      HGGQCD2(ASH,5,MH,MT)
     .   -(1d0+ASH/Pi*(95d0/4d0-5d0*7d0/6d0)))

       GamHbb=4d0*(MBP/MH)**2*
     .        3d0*GF*(MBP*CB(1))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .        *TQCDH((MBP/MH)**2)*dsqrt(1d0-4d0*MBP**2/MH**2)**3
     .         +(1d0-4d0*(MBP/MH)**2)*max(0d0,
     .         3d0*GF*(RMB*CB(1))**2*MH/(4d0*dsqrt(2d0)*Pi)
     .      *QCDH((RMB/MH)**2)*dsqrt(1d0-4d0*RMB**2/MH**2)**3!+DBB
     .       )
* End new
       IF(MH.le.50d0)THEN
        GamHbb=GamHbb*((1d0-(MH-2d0*MBm)/(50d0-2d0*MBm))*
     .  (dsqrt(1d0-4d0*MBm**2/MH**2)/dsqrt(1d0-4d0*MBP**2/MH**2))**3
     .  +(MH-2d0*MBm)/(50d0-2d0*MBm))
* New July 2019:
        GamHhadr=GamHhadr+DBB*
     .  ((1d0-4d0*(MBP/MH)**2)*(1d0-(MH-2d0*MBm)/(50d0-2d0*MBm))*
     .  (dsqrt(1d0-4d0*MBm**2/MH**2)/dsqrt(1d0-4d0*MBP**2/MH**2))**3
     .  +(MH-2d0*MBm)/(50d0-2d0*MBm))
* End new
       ENDIF
      ENDIF

      RETURN
      END

*********************************************************************

      DOUBLE COMPLEX FUNCTION PropI(MH,M,G)

      DOUBLE PRECISION MH,M,G

       PropI=M**2/DCMPLX(MH**2-M**2,M*G)  

      END

*********************************************************************

      DOUBLE COMPLEX FUNCTION PropII(MH,M,G)

      DOUBLE PRECISION MH,M,G

      IF(MH.le.M)THEN
       PropII=M**2/DCMPLX(MH**2-M**2,M*G)    
      ELSE
       PropII=((MH**2-M**2)**3+M**6)**(1d0/3d0)/DCMPLX(MH**2-M**2,M*G)
      ENDIF 

      END

*********************************************************************
c   Scalar Energy-momentum Tensor Form factors from arXiv:1809.01876
*********************************************************************

      DOUBLE PRECISION FUNCTION AbsThetaPiPiFit(MH)

      DOUBLE PRECISION MH

       AbsThetaPiPiFit=0d0  
      IF(MH.ge.0d0.and.MH.le.0.9790628397170731d0)THEN
       AbsThetaPiPiFit=0.035372435008801224d0+0.7190233449441249d0*MH
     c         -23.996794796757335d0*MH**2+346.86954362729125d0*MH**3
     c         -2500.8738105055777d0*MH**4+10466.913355718469d0*MH**5
     c         -26455.81261827014d0*MH**6+40884.961796922325d0*MH**7
     c         -37747.344496205194d0*MH**8+19100.2380699125d0*MH**9
     c         -4068.5398648001346d0*MH**10
      ENDIF
      IF(MH.gt.0.9790628397170731d0.and.MH.le.2d0)THEN
       AbsThetaPiPiFit=792321.8205523676d0-5629423.6094466485d0*MH
     c        +17887717.97870951d0*MH**2-33467993.693830714d0*MH**3
     c        +40825802.191426344d0*MH**4-33923329.70547809d0*MH**5
     c        +19444167.638917923d0*MH**6-7591105.801080975d0*MH**7
     c        +1931877.272082622d0*MH**8-289415.3018455889d0*MH**9
     c        +19382.9939802244d0*MH**10
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION phiPiPiFit(MH)

      DOUBLE PRECISION MH

       phiPiPiFit=0d0  
      IF(MH.ge.0.2851978214598723d0.and.
     c             MH.le.0.9921591335920872d0)THEN
       phiPiPiFit=-70.10415342002676d0+1224.7163623370416d0*MH
     c       -9432.460040975142d0*MH**2+41819.51719525536d0*MH**3
     c       -117294.49900829424d0*MH**4+215758.11607247705d0*MH**5
     c       -260186.2280515026d0*MH**6+198291.93152492034d0*MH**7
     c       -86658.09071993448d0*MH**8+16550.656540620956d0*MH**9
      ENDIF
      IF(MH.gt.0.9921591335920872d0.and.MH.le.2d0)THEN
       phiPiPiFit=189445.41331880973d0-1243868.766630274d0*MH
     c        +3591607.8335345024d0*MH**2-5986874.173883149d0*MH**3
     c        +6350488.95697874d0*MH**4-4446522.389896632d0*MH**5
     c        +2055730.2364366285d0*MH**6-605315.8369538989d0*MH**7
     c        +103040.1166864158d0*MH**8-7728.179370567146d0*MH**9
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION AbsThetaKKFit(MH)

      DOUBLE PRECISION MH

       AbsThetaKKFit=0d0  
      IF(MH.ge.0d0.and.MH.le.0.9824963645861275d0)THEN
       AbsThetaKKFit=0.47026841255227014d0+1.4131724216350454d0*MH
     c      -44.20292180337171d0*MH**2+562.3433940303258d0*MH**3
     c      -3653.5431709007953d0*MH**4+13929.40400462669d0*MH**5
     c      -32630.682023639954d0*MH**6+47387.910605680685d0*MH**7
     c      -41478.46530401207d0*MH**8+19990.869624107727d0*MH**9
     c      -4060.587912385784d0*MH**10
      ENDIF
      IF(MH.gt.0.9824963645861275d0.and.MH.le.2d0)THEN
       AbsThetaKKFit=1555921.0933963354d0-10843917.847846407d0*MH
     c       +33804108.094041914d0*MH**2-62068795.067201085d0*MH**3
     c       +74338187.13856113d0*MH**4-60683172.04295773d0*MH**5
     c       +34194210.92848184d0*MH**6-13133879.260693451d0*MH**7
     c       +3291121.107592999d0*MH**8-485875.8577210487d0*MH**9
     c       +32094.566810581477d0*MH**10
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION phiKKFit(MH)

      DOUBLE PRECISION MH

       phiKKFit=0d0  
      IF(MH.ge.0.2851978214598723d0.and.
     c             MH.le.0.9927068987104044d0)THEN
       phiKKFit=303.5934420038648d0-5934.828580634104d0*MH
     c      +51123.25880773344d0*MH**2-255582.79276919903d0*MH**3
     c      +821418.3392307195d0*MH**4-1774113.8536080194d0*MH**5
     c      +2609680.6292278985d0*MH**6-2583884.396660344d0*MH**7
     c      +1649644.2369654302d0*MH**8-613858.6437072959d0*MH**9
     c      +101206.94232097174d0*MH**10
      ENDIF
      IF(MH.gt.0.9927068987104044d0.and.MH.le.2d0)THEN
       phiKKFit=98854.13071112418d0-645072.4089060655d0*MH
     c      +1864697.6939710448d0*MH**2-3139039.2003097474d0*MH**3
     c      +3400256.1983084083d0*MH**4-2468726.5188846593d0*MH**5
     c      +1211160.7621937892d0*MH**6-393693.3390883526d0*MH**7
     c      +80205.06406669847d0*MH**8-9049.988232542877d0*MH**9
     c      +409.68922827712294d0*MH**10
      ENDIF 

      END
      
*********************************************************************
c   Scalar qqbar Form factors from arXiv:2011.00921
*********************************************************************

      DOUBLE PRECISION FUNCTION RFqqPiPiFit(MH)

      DOUBLE PRECISION MH

       RFqqPiPiFit=0d0  
      IF(MH.ge.0d0.and.MH.le.0.37093093093093094d0)THEN
       RFqqPiPiFit=-2.224007116483618d0+8.741872962623434d0*MH
     c        -401.76633423935357d0*MH**2+8354.744171574033d0*MH**3
     c        -93201.7788157771d0*MH**4+581066.1627191864d0*MH**5
     c        -2.023255387986879d6*MH**6+3.6546034899301054d6*MH**7
     c        -2.657283053455844d6*MH**8
      ENDIF
      IF(MH.gt.0.37093093093093094d0
     c            .and.MH.le.0.9415115115115116d0)THEN
       RFqqPiPiFit=410.8572358050331d0-5741.097696312264d0*MH
     c        +34756.13331818974d0*MH**2-119547.90885613284d0*MH**3
     c        +254159.8762799177d0*MH**4-340302.21202045644d0*MH**5
     c        +279456.9901845174d0*MH**6-128573.08755906425d0*MH**7
     c        +25381.027807103652d0*MH**8
      ENDIF 
      IF(MH.gt.0.9415115115115116d0
     c            .and.MH.le.0.9832032032032032d0)THEN
       RFqqPiPiFit=1.5535763789365967d10-9.68893629085741d10*MH
     c      +2.5175556264937296d11*MH**2-3.488602107854444d11*MH**3
     c      +2.719048351467953d11*MH**4-1.1301906190409894d11*MH**5
     c      +1.9572474042003315d10*MH**6
      ENDIF 
      IF(MH.gt.0.9832032032032032d0
     c            .and.MH.le.1.0225125125125125d0)THEN
       RFqqPiPiFit=5.976341528837345d9-2.5339735904760395d10*MH
     c      +3.2750041926757286d10*MH**2+1.3500995745012934d9*MH**3
     c      -2.94768099754552d10*MH**4+1.620654321475292d9*MH**5
     c      +3.2044172046400337d10*MH**6-2.472639731838865d10*MH**7
     c      +5.801633800763721d9*MH**8
      ENDIF 
      IF(MH.gt.1.0225125125125125d0
     c            .and.MH.le.1.2d0)THEN
       RFqqPiPiFit=-2.3365821791540998d8+1.4664609969400995d9*MH
     c      -3.7616054991692863d9*MH**2+4.750521610974811d9*MH**3
     c      -2.2182971369648995d9*MH**4-1.7247135227080576d9*MH**5
     c      +3.232833599506887d9*MH**6-2.0805781843905506d9*MH**7
     c      +6.523399750053192d8*MH**8-8.330362272475152d7*MH**9
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       RFqqPiPiFit=-0.6447884841997111d0*(1.2d0**2-0.99d0**2)
     c                /(MH**2-0.99d0**2)
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION IFqqPiPiFit(MH)

      DOUBLE PRECISION MH

       IFqqPiPiFit=0d0  
      IF(MH.ge.0.26253253253253256d0
     c            .and.MH.le.0.9617617617617618d0)THEN
       IFqqPiPiFit=-143.30743540227138d0+3042.869905061221d0*MH
     c     -27732.272942080806d0*MH**2+143652.114864149d0*MH**3
     c     -468911.3591168181d0*MH**4+1.0076676542158348d6*MH**5
     c     -1.446826952105172d6*MH**6+1.376011951606933d6*MH**7
     c     -833589.8475282837d0*MH**8+291925.4588915927d0*MH**9
     c     -45098.65897193111d0*MH**10
      ENDIF
      IF(MH.gt.0.9617617617617618d0
     c            .and.MH.le.0.9855855855855856d0)THEN
       IFqqPiPiFit=-1.7852436307092912d9+7.371672968361654d9*MH
     c     -9.525192384865444d9*MH**2+5.372509570446173d7*MH**3
     c     +1.0052993176636566d10*MH**4-8.313311552491243d9*MH**5
     c     +2.1453564196820405d9*MH**6
      ENDIF 
      IF(MH.gt.0.9855855855855856d0
     c            .and.MH.le.1.0248948948948948d0)THEN
       IFqqPiPiFit=-7.319048156065587d8+2.900524003105904d9*MH
     c     -3.589692556508995d9*MH**2-8.933256728283413d6*MH**3
     c     +3.538298008410879d9*MH**4-2.80221697587655d9*MH**5
     c     +6.939255927284353d8*MH**6
      ENDIF 
      IF(MH.gt.1.0248948948948948d0
     c            .and.MH.le.1.2d0)THEN
       IFqqPiPiFit=1.1525221990277072d7-6.249972027452044d7*MH
     c     +1.3209413804378164d8*MH**2-1.2090896957947361d8*MH**3
     c     +4.508848837771294d6*MH**4+9.130180983130822d7*MH**5
     c     -8.349062955888985d7*MH**6+3.2316770171020553d7*MH**7
     c     -4.84746969144681d6*MH**8
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       IFqqPiPiFit=-0.6052538681839782d0*(1.2d0**2-0.99d0**2)
     c                /(MH**2-0.99d0**2)
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION RFssPiPiFit(MH)

      DOUBLE PRECISION MH

       RFssPiPiFit=0d0  
      IF(MH.ge.0d0.and.MH.le.0.6937437437437438d0)THEN
       RFssPiPiFit=-0.015418800638912789d0-0.6227922786068617d0*MH
     c     +22.907843102297647d0*MH**2-398.2047787267778d0*MH**3
     c     +3787.6769021284626d0*MH**4-21408.689736986256d0*MH**5
     c     +74721.83596193549d0*MH**6-162334.53820974397d0*MH**7
     c     +213408.89446457734d0*MH**8-155049.41340382805d0*MH**9
     c     +47676.505372644926d0*MH**10
      ENDIF
      IF(MH.gt.0.6937437437437438d0
     c            .and.MH.le.0.9415115115115116d0)THEN
       RFssPiPiFit=-4522.521235909181d0+33974.9039261146d0*MH
     c     -106101.7365892397d0*MH**2+176360.8491950196d0*MH**3
     c     -164594.34673513463d0*MH**4+81790.65235551102d0*MH**5
     c     -16908.748654522824d0*MH**6
      ENDIF 
      IF(MH.gt.0.9415115115115116d0
     c            .and.MH.le.0.9832032032032032d0)THEN
       RFssPiPiFit=4.026744936475702d10-2.5112277041158286d11*MH
     c      +6.524952654577241d11*MH**2-9.041444869453922d11*MH**3
     c      +7.046788797646075d11*MH**4-2.9289621797961835d11*MH**5
     c      +5.072188082531002d10*MH**6
      ENDIF 
      IF(MH.gt.0.9832032032032032d0
     c            .and.MH.le.1.0248948948948948d0)THEN
       RFssPiPiFit=-5.04099619346091d7+1.894473533053543d8*MH
     c       -2.184228613785049d8*MH**2-1.3112708179371137d7*MH**3
     c       +2.1156567724739408d8*MH**4-1.5487842703262043d8*MH**5
     c       +3.581092981454866d7*MH**6
      ENDIF 
      IF(MH.gt.1.0248948948948948d0
     c            .and.MH.le.1.1047047047047047d0)THEN
       RFssPiPiFit=-2.6013665769498634d9+1.4585872737804817d10*MH
     c      -3.1805024635742714d10*MH**2+2.9732718204535027d10*MH**3
     c      -4.6334958139721334d7*MH**4-2.5873265666451794d10*MH**5
     c      +2.4180830252540047d10*MH**6-9.68014525387902d9*MH**7
     c      +1.5067158840276325d9*MH**8
      ENDIF 
      IF(MH.gt.1.1047047047047047d0
     c            .and.MH.le.1.2d0)THEN
       RFssPiPiFit=1.1409545567619689d6-6.273774441483194d6*MH
     c     +1.4320067484254692d7*MH**2-1.737393103709432d7*MH**3
     c     +1.182106776803527d7*MH**4-4.277802382258467d6*MH**5
     c     +643407.78265745d0*MH**6
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       RFssPiPiFit=-0.18465217607297763d0*(1.2d0**2-0.99d0**2)
     c                /(MH**2-0.99d0**2)
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION IFssPiPiFit(MH)

      DOUBLE PRECISION MH

       IFssPiPiFit=0d0  
      IF(MH.ge.0.26253253253253256d0
     c            .and.MH.le.0.7235235235235236d0)THEN
       IFssPiPiFit=-363.0164905072285d0+8092.397929999013d0*MH
     c   -79677.83712937417d0*MH**2+455952.0761752348d0*MH**3
     c   -1.6778176491626576d6*MH**4+4.1440278377091843d6*MH**5
     c   -6.948363983913389d6*MH**6+7.79699738071957d6*MH**7
     c   -5.592338268225534d6*MH**8+2.308924562003486d6*MH**9
     c   -415193.147133364d0*MH**10
      ENDIF
      IF(MH.ge.0.7235235235235236d0
     c            .and.MH.le.0.9617617617617618d0)THEN
       IFssPiPiFit=789105.194297382d0-7.457095817393999d6*MH
     c   +3.074329976019017d7*MH**2-7.220982201349224d7*MH**3
     c   +1.0566978610997672d8*MH**4-9.863429740243419d7*MH**5
     c   +5.733595365385374d7*MH**6-1.897192323338045d7*MH**7
     c   +2.734989435519957d6*MH**8
      ENDIF
      IF(MH.gt.0.9617617617617618d0
     c            .and.MH.le.0.9855855855855856d0)THEN
       IFssPiPiFit=-6.678540131899431d9+2.7547031320665337d10*MH
     c   -3.554615516392238d10*MH**2+1.6192328714839274d8*MH**3
     c   +3.750262105187552d10*MH**4-3.0970360034640415d10*MH**5
     c   +7.983479957558158d9*MH**6
      ENDIF 
      IF(MH.gt.0.9855855855855856d0
     c            .and.MH.le.1.0332332332332332d0)THEN
       IFssPiPiFit=-5.740903791653013d9+2.3805585961846123d10*MH
     c   -2.9931127053433678d10*MH**2-1.7782783982727406d9*MH**3
     c   +2.6399914857173866d10*MH**4-9.010878830483193d8*MH**5
     c   -2.8118526925628757d10*MH**6+2.1106941350062878d10*MH**7
     c   -4.842518114675784d9*MH**8
      ENDIF 
      IF(MH.gt.1.0332332332332332d0
     c            .and.MH.le.1.2d0)THEN
       IFssPiPiFit=1.06599457350483d7-5.650350609258485d7*MH
     c    +1.1651584752359103d8*MH**2-1.033263336491175d8*MH**3
     c    +1.3930772400478045d6*MH**4+7.814833691182758d7*MH**5
     c    -6.917639328467607d7*MH**6+2.6117227402440857d7*MH**7
     c    -3.828198445841053d6*MH**8
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       IFssPiPiFit=0.6261087496644424d0*(1.2d0**2-0.99d0**2)
     c                /(MH**2-0.99d0**2)
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION RFqqPiEtaFit(MH)

      DOUBLE PRECISION MH

       RFqqPiEtaFit=0d0  
      IF(MH.ge.0d0.and.MH.le.0.9736736736736737d0)THEN
       RFqqPiEtaFit=-1.2460458971331339d0+1.2939405229614886d0*MH
     c       -29.97118389853018d0*MH**2+296.76582376216356d0*MH**3
     c       -1570.9690582968399d0*MH**4+4701.597179965452d0*MH**5
     c       -8115.856071636114d0*MH**6+7738.790804573481d0*MH**7
     c       -3378.2067766537975d0*MH**8+59.37972289714424d0*MH**9
     c       +295.43994822426134d0*MH**10
      ENDIF
      IF(MH.gt.0.9736736736736737d0
     c            .and.MH.le.1.2d0)THEN
       RFqqPiEtaFit=1.4464626562276587d8-9.30914614397502d8*MH
     c     +2.4511989919471726d9*MH**2-3.189445043873542d9*MH**3
     c     +1.5743221885307086d9*MH**4+1.117154717541024d9*MH**5
     c     -2.2339639907826934d9*MH**6+1.4805815397057052d9*MH**7
     c     -4.757231307023046d8*MH**8+6.2143074289532006d7*MH**9
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       RFqqPiEtaFit=-1.2675468933451968d0*(1.2d0**2-0.99d0**2)
     c                /(MH**2-0.99d0**2)
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION IFqqPiEtaFit(MH)

      DOUBLE PRECISION MH

       IFqqPiEtaFit=0d0  
      IF(MH.ge.0.6865965965965967d0
     c               .and.MH.le.0.9939239239239239d0)THEN
       IFqqPiEtaFit=-1.9260786531745765d8+2.3559691681323714d9*MH
     c   -1.2945186776546362d10*MH**2+4.2075987374903465d10*MH**3
     c   -8.959041861223022d10*MH**4+1.3057657602776733d11*MH**5
     c   -1.3192938327430016d11*MH**6+9.124272497614377d10*MH**7
     c   -4.133940430146686d10*MH**8+1.1079737202586342d10*MH**9
     c   -1.3339939230837686d9*MH**10
      ENDIF
      IF(MH.gt.0.9939239239239239d0
     c            .and.MH.le.1.2d0)THEN
       IFqqPiEtaFit=-2.7523881616984032d7+1.7437083587144455d8*MH
     c   -4.5171771311724484d8*MH**2+5.769915510886595d8*MH**3
     c   -2.752037433869808d8*MH**4-2.0696047735397726d8*MH**5
     c   +3.9738834285749173d8*MH**6-2.5869065218764776d8*MH**7
     c   +8.189978664424552d7*MH**8-1.0554051146483568d7*MH**9
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       IFqqPiEtaFit=-1.079788588424158d0*(1.2d0**2-0.99d0**2)
     c                /(MH**2-0.99d0**2)
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION RFqqEtaEtaFit(MH)

      DOUBLE PRECISION MH

       RFqqEtaEtaFit=0d0  
      IF(MH.ge.0d0.and.MH.le.0.962952952952953d0)THEN
       RFqqEtaEtaFit=0.937825508055591d0-4.042793992788487d0*MH
     c     +100.72621383978884d0*MH**2-1130.4922777562551d0*MH**3
     c     +6905.7313514004845d0*MH**4-25193.773966449975d0*MH**5
     c     +57836.51005703213d0*MH**6-84807.80864081168d0*MH**7
     c     +77341.8927241942d0*MH**8-39995.59693049815d0*MH**9
     c     +8951.931052580356d0*MH**10
      ENDIF
      IF(MH.gt.0.962952952952953d0
     c            .and.MH.le.1.0213213213213213d0)THEN
       RFqqEtaEtaFit=-2.1684139833210205d10+1.0069099121822414d11*MH
     c   -1.462737321182392d11*MH**2+1.0601844633036467d10*MH**3
     c   +1.2473824071718013d11*MH**4-2.177991674223826d9*MH**5
     c   -1.2549605699188197d11*MH**6-6.682863349007509d9*MH**7
     c   +1.476676185488853d11*MH**8-1.0424248303279448d11*MH**9
     c   +2.2858571880574493d10*MH**10
      ENDIF 
      IF(MH.gt.1.0213213213213213d0
     c            .and.MH.le.1.2d0)THEN
       RFqqEtaEtaFit=9.4671348767898d6-5.127563132234187d7*MH
     c   +1.0823206104556014d8*MH**2-9.891207350320348d7*MH**3
     c   +3.584858240047766d6*MH**4+7.467331064957827d7*MH**5
     c   -6.817875650752455d7*MH**6+2.63585251652871d7*MH**7
     c   -3.9494286682399227d6*MH**8
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       RFqqEtaEtaFit=1.4350918189015411d0*(1.2d0**2-0.99d0**2)
     c                /(MH**2-0.99d0**2)
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION IFqqEtaEtaFit(MH)

      DOUBLE PRECISION MH

       IFqqEtaEtaFit=0d0  
      IF(MH.ge.0.26253253253253256d0
     c               .and.MH.le.0.6996996996996997d0)THEN
       IFqqEtaEtaFit=-691.2812531165536d0+16393.398461572633d0*MH
     c   -172395.6181056707d0*MH**2+1.0602796305458865d6*MH**3
     c   -4.228933542140724d6*MH**4+1.1441138795489782d7*MH**5
     c   -2.127790613793412d7*MH**6+2.687143742790305d7*MH**7
     c   -2.2057508258240454d7*MH**8+1.0627053504449926d7*MH**9
     c   -2.281757267080227d6*MH**10
      ENDIF
      IF(MH.gt.0.6996996996996997d0
     c            .and.MH.le.0.9617617617617618d0)THEN
       IFqqEtaEtaFit=-3.022480575805993d6+3.0237569792189866d7*MH
     c   -1.2830267012253772d8*MH**2+2.910098217678973d8*MH**3
     c   -3.442824304463291d8*MH**4+8.92764846327866d7*MH**5
     c   +3.388899025445695d8*MH**6-5.403836812246559d8*MH**7
     c   +3.880155783579015d8*MH**8-1.4368851153642368d8*MH**9
     c   +2.225042290805304d7*MH**10
      ENDIF 
      IF(MH.gt.0.9617617617617618d0
     c            .and.MH.le.0.9832032032032032d0)THEN
       IFqqEtaEtaFit=3.395117735924447d9-1.3976150555584154d10*MH
     c   +1.798404129380557d10*MH**2-2.181060283010773d7*MH**3
     c   -1.900556003142175d10*MH**4+1.5650451224260662d10*MH**5
     c   -4.0260891342788906d9*MH**6
      ENDIF 
      IF(MH.gt.0.9832032032032032d0
     c            .and.MH.le.1.0213213213213213d0)THEN
       IFqqEtaEtaFit=-1.1022525204747696d9+4.393928724097582d9*MH
     c   -5.47353190756923d9*MH**2+913021.8715836589d0*MH**3
     c   +5.43297845654411d9*MH**4-4.33069667278019d9*MH**5
     c   +1.0786608989405327d9*MH**6
      ENDIF 
      IF(MH.gt.1.0213213213213213d0
     c            .and.MH.le.1.2d0)THEN
       IFqqEtaEtaFit=-8.130605221903458d8+5.868379199409519d9*MH
     c   -1.8521580104802837d10*MH**2+3.3387705669767277d10*MH**3
     c   -3.7597774774395546d10*MH**4+2.7083515756144394d10*MH**5
     c   -1.218753117921555d10*MH**6+3.132396092152576d9*MH**7
     c   -3.520501377098154d8*MH**8
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       IFqqEtaEtaFit=-0.4399285146855076d0*(1.2d0**2-0.99d0**2)
     c                /(MH**2-0.99d0**2)
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION RFssEtaEtaFit(MH)

      DOUBLE PRECISION MH

       RFssEtaEtaFit=0d0  
      IF(MH.ge.0d0.and.MH.le.0.9677177177177178d0)THEN
       RFssEtaEtaFit=1.2585845135495943d0-13.355390598312976d0*MH
     c   +350.7715233557971d0*MH**2-4225.905486242332d0*MH**3
     c   +27855.189511616765d0*MH**4-109462.76464183061d0*MH**5
     c   +267488.2695309214d0*MH**6-409787.970822952d0*MH**7
     c   +382484.8921902139d0*MH**8-198711.91379287394d0*MH**9
     c   +44038.04435510852d0*MH**10
      ENDIF
      IF(MH.gt.0.9677177177177178d0
     c            .and.MH.le.0.9974974974974975d0)THEN
       RFssEtaEtaFit=-1.6739678984056358d9+6.846570029435989d9*MH
     c   -8.762046753334974d9*MH**2+4.7553206404884025d7*MH**3
     c   +9.073664658422716d9*MH**4-7.430845427975229d9*MH**5
     c   +1.8990721863849933d9*MH**6
      ENDIF 
      IF(MH.gt.0.9974974974974975d0
     c            .and.MH.le.1.0451451451451452d0)THEN
       RFssEtaEtaFit=-1.2280123972646408d10+7.221078414656856d10*MH
     c   -1.7690264548347595d11*MH**2+2.3110540145808344d11*MH**3
     c   -1.6980627541946695d11*MH**4+6.6533751982396385d10*MH**5
     c   -1.0860892712744898d10*MH**6
      ENDIF 
      IF(MH.gt.1.0451451451451452d0
     c            .and.MH.le.1.2d0)THEN
       RFssEtaEtaFit=-6.94360941574731d7+3.6816176823091626d8*MH
     c   -7.596924435456283d8*MH**2+6.747499101896144d8*MH**3
     c   -1.0689171008747408d7*MH**4-5.0888463415085316d8*MH**5
     c   +4.514446251744869d8*MH**6-1.7071864538997525d8*MH**7
     c   +2.50646644506233d7*MH**8
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       RFssEtaEtaFit=0.5636016545512804d0*(1.2d0**2-0.99d0**2)
     c                /(MH**2-0.99d0**2)
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION IFssEtaEtaFit(MH)

      DOUBLE PRECISION MH

       IFssEtaEtaFit=0d0  
      IF(MH.ge.0.26253253253253256d0
     c               .and.MH.le.0.6996996996996997d0)THEN
       IFssEtaEtaFit=-81.5686401827932d0+1956.6377822290074d0*MH
     c   -20858.30607649679d0*MH**2+130131.1847796296d0*MH**3
     c   -526182.5975258689d0*MH**4+1.4408481831866994d6*MH**5
     c   -2.705927219564996d6*MH**6+3.4413590586865153d6*MH**7
     c   -2.8365079589824835d6*MH**8+1.3682759532134603d6*MH**9
     c   -293347.11442348995d0*MH**10
      ENDIF
      IF(MH.gt.0.6996996996996997d0
     c            .and.MH.le.0.937937937937938d0)THEN
       IFssEtaEtaFit=-854464.578101375d0+8.653280578156931d6*MH
     c   -3.715097392543537d7*MH**2+8.514058575957462d7*MH**3
     c   -1.0122462484307261d8*MH**4+2.41164569474938d7*MH**5
     c   +1.0773814976250134d8*MH**6-1.714582924990717d8*MH**7
     c   +1.2436129358517678d8*MH**8-4.664619149249967d7*MH**9
     c   +7.324785096648426d6*MH**10
      ENDIF 
      IF(MH.gt.0.937937937937938d0
     c            .and.MH.le.0.9843943943943945d0)THEN
       IFssEtaEtaFit=-3.8659297919644875d10+2.4125167298350262d11*MH
     c   -6.272415022094594d11*MH**2+8.696747622233992d11*MH**3
     c   -6.782033667257452d11*MH**4+2.82045883216843d11*MH**5
     c   -4.886815152520404d10*MH**6
      ENDIF 
      IF(MH.gt.0.9843943943943945d0
     c            .and.MH.le.1.0094094094094095d0)THEN
       IFssEtaEtaFit=-4.856820249411857d9+1.336899087654061d10*MH
     c   -4.851555462448209d9*MH**2-1.0983799409656363d10*MH**3
     c   -1.7884604367058456d7*MH**4+1.101465465163478d10*MH**5
     c   +4.916858352313243d9*MH**6-1.3506186328953474d10*MH**7
     c   +4.915742178383084d9*MH**8
      ENDIF 
      IF(MH.gt.1.0094094094094095d0
     c            .and.MH.le.1.2d0)THEN
       IFssEtaEtaFit=-2.229499422751938d9+1.614070435286528d10*MH
     c   -5.109642162790135d10*MH**2+9.238362466426796d10*MH**3
     c   -1.043407373872383d11*MH**4+7.538189895994182d10*MH**5
     c   -3.4020021045527622d10*MH**6+8.768765836170832d9*MH**7
     c   -9.883143285590179d8*MH**8
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       IFssEtaEtaFit=0.7344458918613559d0*(1.2d0**2-0.99d0**2)
     c                /(MH**2-0.99d0**2)
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION RFqqKKFit(MH)

      DOUBLE PRECISION MH

       RFqqKKFit=0d0  
      IF(MH.ge.0d0.and.MH.le.0.9605705705705706d0)THEN
       RFqqKKFit=-1.6566964400083182d0-5.810015093295681d0*MH
     c   +127.93656759698938d0*MH**2-1347.583622640478d0*MH**3
     c   +7801.095767531852d0*MH**4-28286.114629382213d0*MH**5
     c   +64548.77687206958d0*MH**6-89556.80300468474d0*MH**7
     c   +71501.25052759216d0*MH**8-29332.150454772884d0*MH**9
     c   +4543.07076259626d0*MH**10
      ENDIF
      IF(MH.gt.0.9605705705705706d0
     c            .and.MH.le.1.0046446446446446d0)THEN
       RFqqKKFit=2.491124445003588d10-1.520256727301268d11*MH
     c   +3.8654708514143d11*MH**2-5.241573990994536d11*MH**3
     c   +3.997777397034362d11*MH**4-1.626107257690385d11*MH**5
     c   +2.755772830638619d10*MH**6
      ENDIF 
      IF(MH.gt.1.0046446446446446d0
     c            .and.MH.le.1.2d0)THEN
       RFqqKKFit=-5.334741554496561d8+3.4061662904865828d9*MH
     c   -8.893874270878363d9*MH**2+1.1459747908667297d10*MH**3
     c   -5.54997596601195d9*MH**4-4.0774850220980353d9*MH**5
     c   +7.966752852749553d9*MH**6-5.228033949130905d9*MH**7
     c   +1.6661989232511418d9*MH**8-2.1602260797919503d8*MH**9
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       RFqqKKFit=-1.3430809005892712d0*(1.2d0**2-0.99d0**2)
     c                /(MH**2-0.99d0**2)
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION IFqqKKFit(MH)

      DOUBLE PRECISION MH

       IFqqKKFit=0d0  
      IF(MH.ge.0.26253253253253256d0
     c               .and.MH.le.0.8426426426426427d0)THEN
       IFqqKKFit=-107.7911479549693d0+3730.6047746975696d0*MH
     c   -49152.97486011699d0*MH**2+349875.38417657313d0*MH**3
     c   -1.5343019098017942d6*MH**4+4.394010699259236d6*MH**5
     c   -8.391565002415549d6*MH**6+1.0606278449785901d7*MH**7
     c   -8.519290035924591d6*MH**8+3.9366789672789145d6*MH**9
     c   -796345.3932072594*MH**10
      ENDIF 
      IF(MH.gt.0.8426426426426427d0
     c            .and.MH.le.0.9843943943943945d0)THEN
       IFqqKKFit=-2.2953310339852786d8+1.7080976558638318d9*MH
     c   -4.873419256002429d9*MH**2+5.7758495052681465d9*MH**3
     c   +1.949372272105773d8*MH**4-6.758225079428388d9*MH**5
     c   +2.8047976220985785d9*MH**6+7.287592824583931d9*MH**7
     c   -1.0033152431478794d10*MH**8+5.0847021703274975d9*MH**9
     c   -9.616471941701502d8*MH**10
      ENDIF 
      IF(MH.gt.0.9843943943943945d0
     c            .and.MH.le.1.0213213213213213d0)THEN
       IFqqKKFit=6.805447463906928d8-2.7212412975287776d9*MH
     c   +3.403008799573414d9*MH**2-1.1893289687419102d7*MH**3
     c   -3.378157098901343d9*MH**4+2.7028967063077655d9*MH**5
     c   -6.751585672843113d8*MH**6
      ENDIF 
      IF(MH.gt.1.0213213213213213d0
     c            .and.MH.le.1.2d0)THEN
       IFqqKKFit=1.0709686805044966d9-7.739452686565952d9*MH
     c   +2.4456667349383514d10*MH**2-4.4138949449361374d10*MH**3
     c   +4.976256921612676d10*MH**4-3.5887239961746765d10*MH**5
     c   +1.6167182239700733d10*MH**6-4.159758783612083d9*MH**7
     c   +4.6801339593398046d8*MH**8
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       IFqqKKFit=0.2754690032909191d0*(1.2d0**2-0.99d0**2)
     c                /(MH**2-0.99d0**2)
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION RFudKKFit(MH)

      DOUBLE PRECISION MH

       RFudKKFit=0d0  
      IF(MH.ge.0d0.and.MH.le.0.6901701701701702d0)THEN
       RFudKKFit=1.244826298940701d0-4.588141542998997d0*MH
     c   +157.67973890029995d0*MH**2-2520.6302836159766d0*MH**3
     c   +22349.106201821352d0*MH**4-118799.8694173221d0*MH**5
     c   +394004.28934227134d0*MH**6-820925.8003057675d0*MH**7
     c   +1.0435004270418182d6*MH**8-738947.2972337224d0*MH**9
     c   +223331.31161485874d0*MH**10
      ENDIF
      IF(MH.gt.0.6901701701701702d0
     c            .and.MH.le.0.9855855855855856d0)THEN
       RFudKKFit=-10854.102914567855d0+86141.97985970833d0*MH
     c   -281992.66447918274d0*MH**2+488034.7077752426d0*MH**3
     c   -471387.48258050624d0*MH**4+241115.99840383418d0*MH**5
     c   -51056.457623846094d0*MH**6
      ENDIF 
      IF(MH.gt.0.9855855855855856d0
     c            .and.MH.le.1.2d0)THEN
       RFudKKFit=-2.504392614189769d8+1.8408474676940396d9*MH
     c   -5.915024469445698d9*MH**2+1.085194074464175d10*MH**3
     c   -1.243349238854185d10*MH**4+9.109994695147066d9*MH**5
     c   -4.1685664261415944d9*MH**6+1.089143291526758d9*MH**7
     c   -1.2440365192826805d8*MH**8
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       RFudKKFit=0.32798320575996587d0*(1.2d0**2-1d0**2)
     c                /(MH**2-1d0**2)
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION IFudKKFit(MH)

      DOUBLE PRECISION MH

       IFudKKFit=0d0  
      IF(MH.ge.0.6865965965965967d0
     c               .and.MH.le.0.9939239239239239d0)THEN
       IFudKKFit=1.9189337411396074d8-2.3528851418511295d9*MH
     c   +1.2957544818825844d10*MH**2-4.22057596015162d10*MH**3
     c   +9.004575664672466d10*MH**4-1.3148432883090709d11*MH**5
     c   +1.3307652365433398d11*MH**6-9.218385991139536d10*MH**7
     c   +4.182774283601784d10*MH**8-1.122590523682989d10*MH**9
     c   +1.3532773974509213d9*MH**10
      ENDIF 
      IF(MH.gt.0.9939239239239239d0
     c            .and.MH.le.1.2d0)THEN
       IFudKKFit=1.1221413462642291d8-8.107772122023996d8*MH
     c   +2.5618382468070493d9*MH**2-4.623604223100679d9*MH**3
     c   +5.213211734753592d9*MH**4-3.7603209975430317d9*MH**5
     c   +1.69448859571236d9*MH**6-4.361420146811104d8*MH**7
     c   +4.9091739385064684d7*MH**8
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       IFudKKFit=1.9480835104422614d0*(1.2d0**2-1d0**2)
     c                /(MH**2-1d0**2)
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION RFssKKFit(MH)

      DOUBLE PRECISION MH

       RFssKKFit=0d0  
      IF(MH.ge.0d0.and.MH.le.0.962952952952953d0)THEN
       RFssKKFit=-2.379915908326461d0+19.94102586190353d0*MH
     c   -521.1904235608562d0*MH**2+6248.7136755730135d0*MH**3
     c   -41036.65675744621d0*MH**4+160850.25829980886d0*MH**5
     c   -392574.26250946446d0*MH**6+601452.9656752597d0*MH**7
     c   -562018.3560646309d0*MH**8+292539.35411807935d0*MH**9
     c   -64983.29262333879d0*MH**10
      ENDIF
      IF(MH.gt.0.962952952952953d0
     c            .and.MH.le.1.007027027027027d0)THEN
       RFssKKFit=8.162549177073906d10-4.9747449717145483d11*MH
     c   +1.263205579714393d12*MH**2-1.7105876738970906d12*MH**3
     c   +1.3028920994716575d12*MH**4-5.2922396127348114d11*MH**5
     c   +8.95629613867896d10*MH**6
      ENDIF 
      IF(MH.gt.1.007027027027027d0
     c            .and.MH.le.1.2d0)THEN
       RFssKKFit=-7.258641714369085d8+4.635836120203144d9*MH
     c   -1.2106669951068125d10*MH**2+1.5600025571920063d10*MH**3
     c   -7.553323707020802d9*MH**4-5.553250428637707d9*MH**5
     c   +1.0844867095072563d10*MH**6-7.114473162852861d9*MH**7
     c   +2.266590390759649d9*MH**8-2.937377518218257d8*MH**9
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       RFssKKFit=-2.0429950987937358d0*(1.2d0**2-0.99d0**2)
     c                /(MH**2-0.99d0**2)
      ENDIF 

      END

*********************************************************************

      DOUBLE PRECISION FUNCTION IFssKKFit(MH)

      DOUBLE PRECISION MH

       IFssKKFit=0d0  
      IF(MH.ge.0.26253253253253256d0
     c               .and.MH.le.0.7235235235235236d0)THEN
       IFssKKFit=-195.26115870235503+4257.646983173696*MH
     c   -40823.19135723105*MH**2+226208.3983974684*MH**3
     c   -800090.2207473927*MH**4+1.880522624789403d6*MH**5
     c   -2.9587134029338956d6*MH**6+3.0515513030726244d6*MH**7
     c   -1.94713046084997d6*MH**8+675857.0929406591*MH**9
     c   -91077.6191902315*MH**10
      ENDIF 
      IF(MH.gt.0.7235235235235236d0
     c            .and.MH.le.0.9617617617617618d0)THEN
       IFssKKFit=6.745261835772425d6-6.7038930148773186d7*MH
     c   +2.8254382352356863d8*MH**2-6.359862026620796d8*MH**3
     c   +7.437951386929188d8*MH**4-1.784600077008128d8*MH**5
     c   -7.536298736293993d8*MH**6+1.1822366080143297d9*MH**7
     c   -8.429199469784213d8*MH**8+3.106298721330664d8*MH**9
     c   -4.791575560631714d7*MH**10
      ENDIF 
      IF(MH.gt.0.9617617617617618d0
     c            .and.MH.le.0.9843943943943945d0)THEN
       IFssKKFit=-3.2615121501321373d10+1.3425660024294904d11*MH
     c   -1.7280743792553262d11*MH**2+4.399121199530934d8*MH**3
     c   +1.8220007707128268d11*MH**4-1.500825969497722d11*MH**5
     c   +3.8608567820651436d10*MH**6
      ENDIF 
      IF(MH.gt.0.9843943943943945d0
     c            .and.MH.le.1.007027027027027d0)THEN
       IFssKKFit=4.089929579314977d10-1.640008170451977d11*MH
     c   +2.0546573094605533d11*MH**2+1.4889481914333048d8*MH**3
     c   -2.0672233307355643d11*MH**4+1.6574336821003082d11*MH**5
     c   -4.153413965436318d10*MH**6
      ENDIF 
      IF(MH.gt.1.007027027027027d0
     c            .and.MH.le.1.2d0)THEN
       IFssKKFit=2.450261610013141d9-1.7737842114255066d10*MH
     c   +5.6147786008366196d10*MH**2-1.0150604375656705d11*MH**3
     c   +1.1462959069164407d11*MH**4-8.280323629879483d10*MH**5
     c   +3.736319060405493d10*MH**6-9.628730167991293d9*MH**7
     c   +1.085023421149642d9*MH**8
      ENDIF 
      IF(MH.gt.1.2d0
     c            .and.MH.le.2d0)THEN
       IFssKKFit=-0.668149569419343d0*(1.2d0**2-0.99d0**2)
     c                /(MH**2-0.99d0**2)
      ENDIF 

      END
