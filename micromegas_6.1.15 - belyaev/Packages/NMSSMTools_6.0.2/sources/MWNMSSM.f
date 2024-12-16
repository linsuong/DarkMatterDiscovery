C 	Subroutine for MW
C
C
C	Literature:
C	[1] Awramik,Czakon,Freitas,Weiglein,
C	    'Precise prediction of the W-Boson mass in the standard 
C	     model', hep-ph/0311148
C
C	[2] Arbuzov,Awramik, 'ZFITTER: A Semi-analytical program for
C	    fermion pair production in e+ e- annihilation, from 	    
C	    version 6.21 to version 6.42',
C	    Comput.Phys.Commun.174:728-758,2006
C	
C	[3] Passarino,Veltman,' One Loop Corrections for e+ e-
C           Annihilation Into mu+ mu- in the Weinberg Model',
C	    Nucl.Phys.B160:151,1979
C
C	[4]Hagiwara,Matsumoto,'A Novel Approach to confront 
C	   electroweak Data and Theory',Z.Phys.C64:559-620,1994,
C	   hep-ph/9409380 	
C
C	[5]P.Langacker,'Precision Tests of the Standard Electroweak
C	   Model', World Scientific (1995)
C
C
C
C*********************************************************************


	SUBROUTINE MWNMSSM(PAR,PROB)

	implicit none
	INTEGER i,j,It,MWFLAG
	DOUBLE PRECISION PAR(*),PROB(*)
	DOUBLE PRECISION S2TW,C2TW,g2,pi,aux,alpha,Gmu,alSMZ,dalSMZ
	DOUBLE PRECISION MZ,dMZ,mt,dmt,Dalph,dDalph,tanb,sinb,cosb
	DOUBLE PRECISION MSU(2),MSD(2),MSE(2),MST(2),MSB(2),MSL(2)
	DOUBLE PRECISION DELT(2,2),UT(2,2),UB(2,2),UL(2,2),mBL,mHSM
	DOUBLE PRECISION DrSM,dDrSM,drHiggs,drSusy
	DOUBLE PRECISION dvSusy,dbSusy,DrhoNP1L,DrhoNP2L
	DOUBLE PRECISION vertexcorr,SigmaGamPrSusy,SigmaGamZSusy
	DOUBLE PRECISION SigmaWHiggs,SigmaZHiggs
	DOUBLE PRECISION SigmaWrenSusy,SigmaWrenHiggs
	DOUBLE PRECISION SigmaWSusy,SigmaZSusy,B21,aux1
	DOUBLE PRECISION C24zdec,B1zdec,B0zdec,CC0,DD0
	DOUBLE PRECISION selfen,M1,M2,M3,M4,D27,fone
c	DOUBLE PRECISION C12,C23,SigmaGamZHiggs,SigmaGamZSusy,Li_2
c     c ,PhiGL

	DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2)
	DOUBLE PRECISION CMASS
	DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
	DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
	DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
	DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	DOUBLE PRECISION MW,dumw,dMW,DrNMSSM,MWSM,dMWSM,decztt,
     C   deltadecztt,deczee,deltadeczee,BRZTauTau,BRZTTmin,
     C   BRZTTmax,ratio,deltaratio,S2TWeffTau,deltaS2TWeffTau,
     c   S2TWeffSM
	DOUBLE PRECISION Drv
	DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0,MWEXP

	COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
	COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
	COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C		MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C		CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	COMMON/EWPO/MW,dumw,dMW,DrNMSSM,MWSM,dMWSM,decztt,
     C   deltadecztt,deczee,deltadeczee,BRZTauTau,BRZTTmin,
     C   BRZTTmax,ratio,deltaratio,S2TWeffTau,deltaS2TWeffTau,
     c   S2TWeffSM
	COMMON/VEVSHIFT/Drv
        COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
        COMMON/MWFLAG/MWFLAG

C	0 - Constants and experimental data

C  Pi
	Pi=4.d0*datan(1.d0)

C  fine structure constant
	alpha=ALEMMZ0   !1d0/137.03599911d0
C  exp. value for G_Fermi(muon decay)
	Gmu=GF0   !1.16637d-5
C  exp. value for mt
	mt=MT0     !172.58d0
	dmt=0.45d0 !sqrt((.52d0)**2+(.72d0)**2)
C  exp. value for alpha_S(MZ)
	alSMZ=ALSMZ0   !.1177d0
	dalSMZ=1d-3
C  exp. value for MZ
	MZ=MZ0     !91.1875d0
	dMZ=2.1d-3
C  Dalpha(lept.+hadr.)
	Dalph=0.0591579d0  !0.05907d0
	dDalph=1d-4   !4.d-4

C  Trig. Functions of Betab 
        TANB=PAR(3)
        sinb=tanb/dsqrt(1.d0+tanb**2)
        cosb=sinb/tanb

C  World average measured value of MW
        MWEXP=80.4133d0

C  Sfermion sector
        MSU(1)=MUL
        MSU(2)=MUR
        MSD(1)=MDL
        MSD(2)=MDR
        MSE(1)=MLL
        MSE(2)=MLR

        MST(1)=dsqrt(RMST1)
	MST(2)=dsqrt(RMST2)
        UT(1,1)=CST
        UT(1,2)=+dsqrt(1.d0-CST**2)
        UT(2,1)=-dsqrt(1.d0-CST**2)
        UT(2,2)=CST

        MSB(1)=dsqrt(RMSB1)
	MSB(2)=dsqrt(RMSB2)
        UB(1,1)=CSB
        UB(1,2)=+dsqrt(1.d0-CSB**2)
        UB(2,1)=-dsqrt(1.d0-CSB**2)
        UB(2,2)=CSB
	If(CSB**2.ge.0.5d0)THEN
	 mBL=MSB(1)
	ELSE
	 mBL=MSB(2)
	ENDIF

        MSL(1)=MSL1
	MSL(2)=MSL2
        UL(1,1)=CSL
        UL(1,2)=+dsqrt(1.d0-CSL**2)
        UL(2,1)=-dsqrt(1.d0-CSL**2)
        UL(2,2)=CSL

        DELT(1,1)=1.d0
        DELT(1,2)=0.d0
        DELT(2,1)=0.d0
        DELT(2,2)=1.d0

C	I - Standard Model
C Chosen Higgs mass
	mHSM=125.2d0
C MW in SM
	MWSM=80.3799d0-0.05429d0*dlog(mHSM/100.d0)
     C           -8.939d-3*dlog(mHSM/100.d0)**2
     C           +8.90d-5*dlog(mHSM/100.d0)**4
     C           +1.61d-4*((mHSM/100.d0)**2-1.d0)
     C           -1.070d0*(Dalph/0.05907d0-1.d0)
     C           +.5256d0*((mt/174.3)**2-1.d0)
     C           -.0678d0*((mt/174.3)**2-1.d0)**2
     C  -1.79d-3*dlog(mHSM/100.d0)*((mt/174.3)**2-1.d0)
     C  +6.59d-5*(mHSM/100.d0)**2*((mt/174.3)**2-1.d0)
     C           -7.37d-2*(alSMZ/0.119d0-1.d0)
     C           +114.9d0*(MZ/91.1875d0-1.d0)

C param. error / mt
	aux=0.5256d0*(((mt+dmt)/174.3d0)**2-(mt/174.3)**2)
     C    -0.0678d0*((((mt+dmt)/174.3)**2-1.d0)**2
     C                -((mt/174.3)**2-1.d0)**2)
     C    -1.79d-3*dlog(mHSM/100.d0)
     C            *(((mt+dmt)/174.3)**2-(mt/174.3)**2)
     C    +6.59d-5*(mHSM/100.d0)**2
     C       *(((mt+dmt)/174.3)**2-(mt/174.3)**2)
	dMWSM=dabs(aux)
C param. error / alpha_S(MZ)
	aux=7.37d-2*((alSMZ+dalSMZ)/0.119d0-1.d0)
     C       -7.37d-2*((alSMZ)/0.119d0-1.d0)
	dMWSM=dsqrt(dMWSM**2+aux**2)
C param. error / MZ
	aux=114.9d0*((MZ+dMZ)/91.1875d0-1.d0)
     C        -114.9d0*((MZ)/91.1875d0-1.d0)
	dMWSM=dsqrt(dMWSM**2+aux**2)
C parm. error / Dalpha
	aux=1.070d0*((Dalph+dDalph)/0.05907d0-1.d0)
     C      -1.070d0*(Dalph/0.05907d0-1.d0)
	dMWSM=dsqrt(dMWSM**2+aux**2)
C error higher orders
	dMWSM=4d-3+2*dMWSM
c	print*,MWSM,dMWSM
C Delta r in SM
	DrSM=dsqrt(2.d0)*Gmu/(Pi*alpha)*MWSM**2
     C         *(1.d0-MWSM**2/MZ**2)-1.d0
	dDrSM=2.d0*dsqrt(2.d0)*Gmu/(Pi*alpha)*MWSM
     C         *(1.d0-2.d0*MWSM**2/MZ**2)*dMWSM

C Iteration
	MW=MWSM

	DO It=1,10

C Weinberg angle
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW
	g2=4.d0*sqrt(2.d0)*Gmu*MW**2

C First order dependence of DrDM in MW (Fit at 10^-4 for MW=80..80.7 GeV)
	DrNMSSM=DrSM-8.212643109d-3*(MW-MWSM)
     C -1.425198118d-3*(MW-MWSM)**2-1.918891636d-4*(MW-MWSM)**3

C	II - NMSSM specific contributions

C	 1) Higgs sector (only "gauge" since Ye,Ymu~0)
	drHiggs=SigmaWrenHiggs(0.d0,tanb,MW)/MW**2

C	 2) SUSY particles
C	- "gauge" part
	drSUSY=SigmaWrenSusy(0.d0,MW)/MW**2

C	- Vertexcorrections and selfenergies
C	   (Mass-degeneracy of the first two generations)

C		->Vertexcorrections:
	aux=0.d0
	do i=1,5
	aux=aux+(NEU(i,1)*dsqrt(S2TW/C2TW)-NEU(i,2))
     C		*(NEU(i,1)*dsqrt(S2TW/C2TW)+NEU(i,2))
     C	        *C24zdec(0.d0,0.d0,MNL,MNEU(i),MLL)
	enddo

	do i=1,5
	do j=1,2
	aux=aux-U(j,1)*(Neu(i,1)*dsqrt(S2TW/C2TW)+NEU(i,2))*
     C		((U(j,1)*Neu(i,2)+dsqrt(1.d0/2.d0)*U(j,2)*NEU(i,4))
     C     *(1.d0/2.d0-2.d0*C24zdec(0.d0,0.d0,MCHA(j),MLL,MNEU(i)))
     C	       +(NEU(i,2)*V(j,1)-dsqrt(1.d0/2.d0)*NEU(i,3)*V(j,2))
     C	   *MCHA(j)*MNEU(i)*CC0(0.d0,0.d0,MCHA(j),MLL,MNEU(i)))
	enddo
	enddo

	do i=1,5
	do j=1,2
	aux=aux+V(j,1)*(NEU(i,1)*dsqrt(S2TW/C2TW)-NEU(i,2))
     C	    *((NEU(i,2)*V(j,1)-dsqrt(1.d0/2.d0)*Neu(i,3)*V(j,2))
     C 	*(1.d0/2.d0-2.d0*C24zdec(0.d0,0.d0,MNEU(i),MNL,MCHA(j)))
     C	    +(U(j,1)*Neu(i,2)+dsqrt(1.d0/2.d0)*U(j,2)*Neu(i,4))
     C	*MCHA(j)*MNEU(i)*CC0(0.d0,0.d0,MNEU(i),MNL,MCHA(j)))			
	enddo
	enddo

	vertexcorr=g2*dsqrt(g2/2.d0)/(16.d0*pi**2)*aux

C		->self-energies (both e and nue)
C 	-> Neutralino diagrams
	aux=0.d0
	do i=1,5
 	aux=aux+1.d0/2.d0*(NEU(i,1)*dsqrt(S2TW/C2TW)
     C		+NEU(i,2))**2
     C	    *(B0zdec(0.d0,MLL,MNEU(i))+B1zdec(0.d0,MLL,MNEU(i)))
	aux=aux+1.d0/2.d0*(NEU(i,1)*dsqrt(S2TW/C2TW)-NEU(i,2))**2
     C	        *(B0zdec(0.d0,MNL,MNEU(i))+B1zdec(0.d0,MNL,MNEU(i)))
	enddo

C	->chargino diagrams
	do i=1,2
	aux=aux+U(i,1)**2
     C	*(B0zdec(0.d0,MLL,MCHA(i))+B1zdec(0.d0,MLL,MCHA(i)))
	aux=aux+V(i,1)**2
     C	     *(B0zdec(0.d0,MNL,MCHA(i))+B1zdec(0.d0,MNL,MCHA(i)))
	enddo

	selfen=g2/(16.d0*pi**2)*aux

	dvSusy=dsqrt(2.d0/g2)*vertexcorr-1.d0/2.d0*(selfen)

	drSUSY=drSUSY+2.d0*dvSusy

C	- Box-corrections 
C	  (Mass-degeneracy of the first two generations)

	aux=0.d0		
	do i=1,5
	do j=1,2
	aux=aux+U(j,1)**2*(NEU(i,1)*dsqrt(S2TW/C2TW)+NEU(i,2))**2
     C	                  *D27(MSE(1),MSE(1),MCHA(j),MNEU(i))
	enddo
	enddo
	M1=g2**2/(2.d0*16.d0*pi**2)*aux

	aux=0.d0		
	do i=1,5
	do j=1,2
	aux=aux+V(j,1)**2*(NEU(i,1)*dsqrt(S2TW/C2TW)-NEU(i,2))**2
     C	                  *D27(MNL,MNL,MCHA(j),MNEU(i))
	enddo
	enddo
	M2=g2**2/(2.d0*16.d0*pi**2)*aux

	aux=0.d0		
	do i=1,5
	do j=1,2
	aux=aux+aux+MNEU(i)*MCHA(j)*V(j,1)*U(j,1)
     C	         *(NEU(i,1)*dsqrt(S2TW/C2TW)+NEU(i,2))
     C		 *(NEU(i,1)*dsqrt(S2TW/C2TW)-NEU(i,2))
     C	                  *DD0(MNL,MSE(1),MCHA(j),MNEU(i))
	enddo
	enddo
	M3=g2**2/(4.d0*16.d0*pi**2)*aux

	aux=0.d0		
	do i=1,5
	do j=1,2
	aux=aux+MNEU(i)*MCHA(j)*V(j,1)*U(j,1)
     C	         *(NEU(i,1)*dsqrt(S2TW/C2TW)+NEU(i,2))
     C		 *(NEU(i,1)*dsqrt(S2TW/C2TW)-NEU(i,2))
     C	                  *DD0(MSE(1),MNL,MCHA(j),MNEU(i))
	enddo
	enddo
	M4=g2**2/(4.d0*16.d0*pi**2)*aux

	dbSusy=-2.d0*MW**2/g2*(M1+M2+M3+M4)

	drSUSY=drSUSY+dbSusy
c	print*,'v+b',dvSusy,dbSusy
C	 3) 1-loop NP in Drho
	DrhoNP1L=(SigmaZHiggs(0.d0,tanb,MW)+SigmaZsusy(0.d0,MW))
     C  /MZ**2-(SigmaWHiggs(0.d0,tanb,MW)+SigmaWsusy(0.d0,MW))/MW**2

C	 4) NP 2-loop effects
C	- gluon/squark
	aux=0.d0
	do i=1,2
	 do j=1,2
	 aux=aux+(UT(i,1)*UT(j,1)-4.d0/3.d0*S2TW*DELT(i,j))**2/2.d0
     C                               *Fone(MST(i)**2,MST(j)**2)
     C      +(UB(i,1)*UB(j,1)-2.d0/3.d0*S2TW*DELT(i,j))**2/2.d0
     C                               *Fone(MSB(i)**2,MSB(j)**2)
     C      -(UT(i,1)*UB(j,1))**2*Fone(MST(i)**2,MSB(j)**2)
	 enddo
	enddo
	aux=aux+(1.d0-4.d0/3.d0*S2TW)**2*Fone(MUL**2,MUL**2)
     C       +(4.d0/3.d0*S2TW)**2*Fone(MUR**2,MUR**2)
     C      +(1.d0-2.d0/3.d0*S2TW)**2*Fone(MDL**2,MDL**2)
     C       +(2.d0/3.d0*S2TW)**2*Fone(MDR**2,MDR**2)
	aux=aux
     C      -2.d0*Fone(MUL**2,MDL**2)

	DrhoNP2L=-Gmu*alSMZ/(4.d0*dsqrt(2.d0)*Pi**3)*aux

C	- gluino/squarks: ~0 for heavy gluinos
c      MGL=1000.d0
c	aux=-2.d0*mt/MGL*CST*UT(1,2)*(2.d0*CST**2-1.d0)*
c     C (MST(1)**2*CST**2-MST(2)**2*UT(1,2)**2-mBL**4
c     C *(MST(1)**2-MST(2)**2)/(mBL**2-MST(1)**2)/(mBL**2-MST(2)**2)
c     C *dlog(MBL**2)+MST(1)**2*dlog(MST(1)**2)/(mBL**2-MST(1)**2)
c     C /(MST(1)**2-MST(2)**2)*(MBL**2*MST(1)**2-2*CST**2*MST(2)**2
c     C *mBL**2+(2.d0*CST**2-1.d0)*MST(1)**2*MST(2)**2)
c     C +MST(2)**2*dlog(MST(2)**2)/(mBL**2-MST(2)**2)
c     C /(MST(1)**2-MST(2)**2)*(MBL**2*MST(2)**2-2*U(1,2)**2*MST(1)**2
c     C *mBL**2-(2.d0*CST**2-1.d0)*MST(1)**2*MST(2)**2))
c	print*,aux
c	aux=aux
c     C +1.d0/(3.d0*MGL**2)*(2.d0*((mBL**2-MST(1)**2)**2
c     C ))
c	print*,aux,mBL,MST(1),MGL
c	aux=aux+1.d0/(3.d0*MGL**2)*(2.d0*(
c     C -U(1,2)**2
c     C *((MST(1)**2-MST(2)**2)**2*CST**2+(mBL**2-MST(1)**2)**2
c     C -(mBL**2-MST(2)**2)**2)
c     C ))
c	print*,aux
c	aux=aux+1.d0/(3.d0*MGL**2)*(2.d0*(
c     C +3.d0*mt**2*(-MST(1)**2*(3-2.d0*CST**2)
c     C +2*mBL**2-4.d0*U(1,2)**2*MST(2)**2+3.d0*U(1,2)**4
c     C *(MST(1)**2+MST(2)**2)))
c     C )
c	print*,aux
c	aux=aux+1.d0/(3.d0*MGL**2)*(
c     C -6.d0*mBL**4*mt**2*dlog(mBL**2)
c     C *(1.d0/(mBL**2-MST(1)**2)-U(1,2)**2*(MST(1)**2-MST(2)**2)
c     C /((mBL**2-MST(1)**2)*(mBL**2-MST(2)**2)))
c     C )
c	print*,aux
c	aux=aux+1.d0/(3.d0*MGL**2)*(
c     C +6.d0*MST(1)**2
c     C *mt**2*CST**2*dlog(MST(1)**2)*((2*mBL**2-MST(1)**2)
c     C /(mBL**2-MST(1)**2)+6.d0*U(1,2)**2*MST(2)**2
c     C /(MST(1)**2-MST(2)**2))
c     C )
c	print*,aux
c	aux=aux+1.d0/(3.d0*MGL**2)*(
c     C +6.d0*MST(2)**2*mt**2*U(1,2)**2
c     C *dlog(MST(2)**2)*((MST(2)**4-4.d0*mBL**2*MST(1)**2-2.d0*mBL**2
c     C *MST(2)**2+5.d0*MST(1)**2*MST(2)**2)/((mBL**2-MST(2)**2)
c     C *(MST(1)**2-MST(2)**2))+6.d0*U(1,2)**2*MST(1)**2
c     C /(MST(1)**2-MST(2)**2))
c     C )
c	print*,aux
c	aux=aux+1.d0/(3.d0*MGL**2)*(
c     C +3.d0*mt**2*dlog(mt**2)*(mBL**2+MST(1)**2-2.d0*mt**2-U(1,2)**2
c     C *(5.d0*MST(1)**2+3.d0*MST(2)**2)+4.d0*U(1,2)**4*(MST(1)**2
c     C +MST(2)**2))
c     C )
c	print*,aux
c	aux=aux+1.d0/(3.d0*MGL**2)*(
c     C -6.d0*mBL**2*mt**2*dlog(mBL**2)*dlog(mt**2/MGL**2)
c     C *(MST(1)**2/(mBL**2-MST(1)**2)-U(1,2)**2*mBL**2
c     C *(MST(1)**2-MST(2)**2)/((mBL**2-MST(1)**2)
c     C *(mBL**2-MST(2)**2)))
c     C )
c	print*,aux
c	aux=aux+1.d0/(3.d0*MGL**2)*(
c     C +6.d0*MST(1)**2*mt**2*CST**2
c     C *dlog(MST(1)**2)*dlog(mt**2/MGL**2)*(mBL**2
c     C /(mBL**2-MST(1)**2)+4.d0*U(1,2)**2*MST(2)**2
c     C /(MST(1)**2-MST(2)**2))
c     C )
c	print*,aux
c	aux=aux+1.d0/(3.d0*MGL**2)*(
c     C +6.d0*mt**2*U(1,2)**2*MST(2)**2
c     C *dlog(MST(2)**2)*dlog(mt**2/MGL**2)*(mBL**2/(mBL**2-MST(2)**2)
c     C -4.d0*CST**2*MST(1)**2/(MST(1)**2-MST(2)**2))
c     C )
c	print*,aux
c	aux=aux+1.d0/(3.d0*MGL**2)*(
c     C +3.d0*mt**2*dlog(MGL**2)*(mBL**2-3.d0*MST(1)**2+2.d0*mt**2
c     C +4.d0*U(1,2)**2*CST**2*(MST(1)**2+MST(2)**2)+3.d0*U(1,2)**2
c     C *(MST(1)**2-MST(2)**2)))
c	print*,aux
c      aux=aux
c     C +mt*U(1,2)*CST/(6.d0*MGL**3)*(12.d0*mt**2*Pi**2
c     C *(MST(1)**2-MST(2)**2)+36.d0*mt**2*dlog(MGL**2)**2
c     C *(mST(1)**2-MST(2)**2)-(4.d0*mBL**2*MST(1)**2+4.d0*MST(1)**4
c     C -4.d0*mBL**2*MST(2)**2+4.d0*MST(1)**2*MST(2)**2+57.d0*MST(1)**2
c     C *mt**2-21.d0*MST(2)**2*mt**2-20.d0*MST(1)**4*U(1,2)**2
c     C -16.d0*MST(1)**2*MST(2)**2*U(1,2)**2+4.d0*MST(2)**4*U(1,2)**2
c     C -108.d0*MST(1)**2*mt**2*U(1,2)**2-36.d0*MST(2)**2*mt**2
c     C *U(1,2)**2+16.d0*MST(1)**4*U(1,2)**4+16.d0*MST(1)**2*MST(2)**2
c     C *U(1,2)**4+72.d0*MST(1)**2*mt**2*U(1,2)**4+72.d0*MST(2)**2
c     C *mt**2*U(1,2)**4)-24.d0*mt**2*dlog(MGL**2)*(3.d0*MST(1)**2
c     C -4.d0*MST(2)**2+3.d0*MST(1)**2*U(1,2)**2+MST(2)**2*U(1,2)**2
c     C -2.d0*MST(1)**2*U(1,2)**4-2.d0*MST(2)**2*U(1,2)**4)
c     C +4.d0*(MST(1)**2-MST(2)**2)*mBL**4*dlog(mBL**2)
c     C /((mBL**2-MST(1)**2)*(mBL**2-MST(2)**2))*(2.d0*MST(1)**2
c     C +15.d0*mt**2-4.d0*MST(1)**2*U(1,2)**2-18.d0*mt**2*U(1,2)**2)
c     C -4.d0*MST(1)**2*dlog(MST(1)**2)/(mBL**2-MST(1)**2)
c     C /(MST(1)**2-MST(2)**2)*(2.d0*mBL**2*MST(1)**4-4.d0*mBL**2
c     C *MST(1)**2*MST(2)**2+2.d0*MST(1)**4*MST(2)**2+3.d0*mBL**2
c     C *MST(1)**2*mt**2+12.d0*MST(1)**4*mt**2-12.d0*mBL**2*MST(2)**2
c     C *mt**2-3.d0*MST(1)**2*MST(2)**2*mt**2-4.d0*mBL**2*MST(1)**4
c     C *U(1,2)**2+12.d0*mBL**2*MST(1)**2*MST(2)**2*U(1,2)**2
c     C -8.d0*MST(1)**4*MST(2)**2*U(1,2)**2-18.d0*mBL**2*MST(1)**2
c     C *mt**2*U(1,2)**2+54.d0*mBL**2*MST(2)**2*mt**2*U(1,2)**2
c     C -36.d0*MST(1)**2*MST(2)**2*mt**2*U(1,2)**2-8.d0*mBL**2
c     C *MST(1)**2*MST(2)**2+8.d0*MST(1)**4*MST(2)**2*U(1,2)**4
c     C -36.d0*mBL**2*MST(2)**2*mt**2*U(1,2)**4+36.d0*MST(1)**2
c     C *MST(2)**2*mt**2*U(1,2)**4)+4.d0*MST(2)**2*dlog(MST(2)**2)
c     C /(mBL**2-MST(2)**2)/(MST(1)**2-MST(2)**2)*(-2.d0*mBL**2
c     C *MST(1)**2*MST(2)**2+2.d0*MST(1)**4*MST(2)**2-6.d0*mBL**2
c     C *MST(1)**2*mt**2-3.d0*mBL**2*MST(2)**2*mt**2+21.d0*mt**2
c     C *MST(1)**2*MST(2)**2-12.d0*mt**2*MST(2)**4+4.d0*mBL**2
c     C *MST(1)**2*U(1,2)**2+4.d0*mBL**2*MST(1)**2*MST(2)**2
c     C *U(1,2)**2-8.d0*MST(1)**4*MST(2)**2*U(1,2)**2+18.d0*mBL**2
c     C *MST(1)**2*mt**2*U(1,2)**2-36.d0*MST(1)**2*MST(2)**2*mt**2
c     C *U(1,2)**2-8.d0*mBL**2*MST(1)**4*U(1,2)**4+8.d0*MST(1)**4
c     C *MST(2)**2*U(1,2)**2+18.d0*mBL**2*MST(1)**2*mt**2
c     C *U(1,2)**2-36.d0*MST(1)**2*MST(2)**2*mt**2*U(1,2)**2
c     C -8.d0*mBL**2*MST(1)**4*U(1,2)**4+8.d0*MST(1)**4
c     C *MST(2)**2*U(1,2)**4-36.d0*mBL**2*MST(1)**2*mt**2
c     C *U(1,2)**4+36.d0*MST(1)**2*MST(2)**2*mt**2*U(1,2)**4)
c     C -24.d0*mt**2*dlog(mt**2)*(-MST(1)**2+2.d0*MST(2)**2
c     C -3.d0*MST(1)**2*U(1,2)**2-MST(2)**2*U(1,2)**2
c     C +2.d0*MST(1)**2*U(1,2)**4+2.d0*MST(2)**2*U(1,2)**4)
c     C -36.d0*mt**2*(MST(1)**2-MST(2)**2)*mt**2*dlog(MGL**2)
c     C *dlog(mt**2)+12.d0*(MST(1)**2-MST(2)**2)*mBL**4
c     C *dlog(mBL**2)*dlog(mt**2/MGL**2)*(4.d0*CST**2-1.d0)
c     C /(mBL**2-MST(1)**2)/(mBL**2-MST(2)**2)+12.d0*mt**2
c     C *MST(1)**2*dlog(MST(1)**2)*dlog(mt**2/MGL**2)
c     C /(mBL**2-MST(1)**2)/(MST(1)**2-MST(2)**2)*(-3.d0
c     C *MST(1)**4+2.d0*mBL**2*MST(2)**2+MST(1)**2*MST(2)**2
c     C +4.d0*mBL**2*MST(1)**2*U(1,2)**2-12.d0*mBL**2*MST(2)**2
c     C *U(1,2)**2+8.d0*MST(1)**2*MST(2)**2*U(1,2)**2+8.d0
c     C *mBL**2*MST(2)**2*U(1,2)**4-8.d0*MST(1)**2*MST(2)**2
c     C *U(1,2)**4)-12.d0*mt**2*MST(2)**2*dlog(MST(2)**2)
c     C *dlog(mt**2/MGL**2)/(mBL**2-MST(2)**2)/(MST(1)**2
c     C -MST(2)**2)*(2.d0*mBL**2*MST(1)**2-5.d0*MST(1)**2
c     C *MST(2)**2+3.d0*MST(2)**4-4.d0*mBL**2*MST(1)**2
c     C *U(1,2)**2-4.d0*mBL**2*MST(2)**2*U(1,2)**2+8.d0
c     C *MST(1)**2*MST(2)**2*U(1,2)**2+8.d0*mBL**2
c     C *MST(1)**2*U(1,2)**4-8.d0*MST(1)**2*MST(2)**2
c     C *U(1,2)**4))
c      print*,aux

c	aux=13.d0*MGL**2-4.d0*mBL**2-9.d0/2.d0*MST(1)**2
c     C -MST(2)**2/2.d0-mt**2*(MST(1)**2-5.d0*mBL**2)
c     C /(mBL**2-MST(1)**2)
c      print*,aux
c	DO i=1,2
c	 aux=aux+(MGL**2-MST(i)**2+mt**2)*B0zdec(mt,MGL,MST(i))
c     C                                /2.d0
c	ENDDO
c      print*,aux
c	aux=aux-(MGL**2-MST(1)**2+mt**2)/(mBL**2-MST(1)**2)**2
c     C *B0zdec(MST(1),MGL,mt)*(3.d0*mBL**4-4.d0*mBL**2*MST(1)**2
c     C +MST(1)**4-2.d0*mBL**4*dlog(mBL**2/MST(1)**2))
c     C -(MGL**2-mBL**2)/(mBL**2-MST(1)**2)**2
c     C *B0zdec(MbL,MGL,0.d0)*(mBL**4-4.d0*mBL**2*MST(1)**2
c     C +3.d0*MST(1)**4-2.d0*MST(1)**4*dlog(MST(1)**2/mBL**2))
c     C -(MGL**2-mBL**2)/mt**2/(MST(1)**2-mBL**2)**2
c     C *Li_2(1.d0-MGL**2/mBL**2)*((MGL**2-mBL**2)
c     C *(MST(1)**2-mBL**2)*(mBL**2-2.d0*MGL**2+MST(1)**2)
c     C +2.d0*mt**2*MST(1)**2*(2.d0*MST(1)**2-mBL**2-MGL**2))
c     C-2.d0*(MGL**2-MST(1)**2)**2/mt**2/(MST(1)**2-mBL**2)**2
c     C *Li_2(1.d0-MGL**2/MST(1)**2)*(MGL**2*(MST(1)**2-mBL**2)
c     C +MST(1)**2*(mBL**2-MST(1)**2+mt**2))
c     C +(MST(2)**2-MGL**2)**2/mt**2*Li_2(1.d0-MGL**2/MST(2)**2)
c      print*,aux
c	DO i=1,2
c	 aux=aux+dlog(MGL**2)*(3.d0/2.d0*MGL**2+(MGL**6-2.d0
c     C *MGL**4*MST(i)**2+MGL**2*MST(i)**4-MGL**4*mt**2-MGL**2
c     C *MST(i)**2*mt**2)/((MST(i)**2-MGL**2-mt**2)**2-4.d0*mt**2
c     C *MGL**2))
c	ENDDO
c      print*,aux
c	aux=aux-2.d0/(mBL**2-MST(1)**2)**2*dlog(mBL**2)*(2.d0
c     C *mt**2*mBL**4+(mBL**2-MST(1)**2)*(MGL**2*(3.d0*mBL**2
c     C +MST(1)**2)-2.d0*mBL**2*MST(1)**2))+dlog(MST(1)**2)
c     C *((4.d0*MGL**2*(mBL**2+3.d0*MST(1)**2)-11.d0*mBL**2
c     C *MST(1)**2+3.d0*MST(1)**4)/2.d0/(mBL**2-MST(1)**2)+(2.d0*
c     C mBL**2*mt**2*(mBL**2+MST(1)**2))/(mBL**2-MST(1)**2)**2
c     C -(2.d0*MGL**2*MST(1)**2*mt**2)/((MST(1)**2-MGL**2
c     C -mt**2)**2-4.d0*MGL**2*mt**2))-dlog(MST(2)**2)*(3.d0/2.d0
c     C *MST(2)**2+2.d0*MGL**2*MST(2)**2*mt**2/((MST(2)**2-MGL**2
c     C -mt**2)**2-4.d0*MGL**2*mt**2))+(MGL**2-MST(2)**2)**2/2.d0
c     C /mt**2*dlog(MST(2)**2)*dlog(MST(2)**2/MGL**2/mt**2)
c      print*,aux
c	DO i=1,2
c	 aux=aux+(-2.d0*MGL**2+mBL**2+2.d0*MST(1)**2+MST(2)**2
c     C +2.d0*MGL**2*(mt**2*(3.d0*MST(i)**2+MGL**2)-(MGL**2
c     C -MST(i)**2)**2)/((MST(i)**2-mt**2-MGL**2)**2-4.d0*mt**2
c     C *MGL**2))/2.d0*dlog(mt**2)
c	ENDDO
c      print*,aux
c	aux=aux+(MGL**2-MST(1)**2)**2/mt**2
c     C /(mBL**2-MST(1)**2)**2*dlog(MST(1)**2)**2
c     C *(MGL**2*(mBL**2-MST(1)**2)-MST(1)**2*(mBL**2-MST(1)**2
c     C +mt**2))
c      print*,aux
c	aux=aux-(MGL**2-mBL**2)*dlog(mBL**2)**2/2.d0/mt**2
c     C /(mBL**2-MST(1)**2)**2*((mBL**2-MST(1)**2)
c     C *(mBL**2-MGL**2)*(mBL**2+MST(1)**2-2.d0*MGL**2)+2.d0
c     C *mt**2*MST(1)**2*(2.d0*MST(1)**2-mBL**2-MGL**2))
c      print*,aux
c	aux=aux
c     C -dlog(MST(1)**2)*dlog(MGL**2)*(mt**2/2.d0+1.d0/mt**2
c     C /(mBL**2-MST(1)**2)**2*(mt**2*(MGL**4*(3.d0*mBL**2
c     C -4.d0*MST(1)**2)-MST(1)**2*(mBL**2-MST(1)**2)**2
c     C +MST(1)**4*(2.d0*MGL**2-mBL**2))+(MGL**2-MST(1)**2)**3
c     C *(mBL**2-MST(1)**2)+mt**4*(mBL**2*(MST(1)**2+mBL**2)
c     C -MGL**2*(mBL**2-MST(1)**2))-mt**6*mBL**2))
c      print*,aux
c	aux=aux
c     C -1.d0/2.d0/mt**2/(mBL**2-MST(1)**2)**2*
c     C ((MGL**2-mBL**2)**2*(2.d0*MGL**2-MST(1)**2-mBL**2)
c     C *(MST(1)**2-mBL**2)+2.d0*mt**2*(MGL**4*(4.d0*MST(1)**2
c     C -3.d0*mBL**2)+MST(1)**4*(mBL**2-2.d0*MGL**2))+mt**4
c     C *((mBL**2+MGL**2)**2-(MST(1)**2+MGL**2)**2-4.d0*mBL**4
c     C )+2.d0*mt**6*mBL**2)*dlog(MGL**2)*dlog(mBL**2)
c      print*,aux
c	aux=aux
c     C +(5.d0*(MGL**2-MST(1)**2)**2+(MGL**2-MST(2))**2
c     C -4.d0*(MGL**2-MST(1)**2)*(mBL**2-MST(1)**2)
c     C +(mBL**2-MST(1)**2)**2-2.d0*mt**2*MST(1)**2)
c     C *dlog(MGL**2)*dlog(mt**2)/2.d0/mt**2
c      print*,aux
c	aux=aux
c     C +((MGL**2-mBL**2)**2*(2.d0*MGL**2-MST(1)**2-mBL**2)
c     C *(mBL**2-MST(1)**2)-2.d0*mt**2*(MGL**2*mBL**2
c     C *(3.d0*MGL**2-2.d0*mBL**2)-2.d0*MGL**4*MST(1)**2
c     C +mBL**2*MST(1)**4)+mt**4*(mBL**2-MST(1)**2)
c     C *(2.d0*MGL**2+MST(1)**2+MBL**2)+2.d0*mt**6*mBL**2)
c     C *dlog(mBL**2)*dlog(mt**2)/2.d0/mt**2
c     C /(mBL**2-MST(1)**2)**2
c      print*,aux
c	aux=aux+dlog(mt**2)*dlog(MST(1)**2)
c     C *(mt**2/2.d0-MST(1)**2-((MGL**2-MST(1)**2)**3
c     C *(mBL**2-MST(1)**2)+mt**2*(MGL**2-MST(1)**2)
c     C *(2.d0*(MST(1)**2-mBL**2)**2+mBL**2*MST(1)**2+MGL**2
c     C *(2.d0*MST(1)**2-3.d0*mBL**2))+mt**4*(MGL**2+mBL**2)
c     C *(mBL**2-MST(1)**2)+mt**6*mBL**2)/mt**2
c     C /(mBL**2-MST(1)**2)**2)
c      print*,aux
c	aux=aux-((mBL**2-mt**2-MGL**2)**2
c     C -4.d0*MGL**2*mt**2)/2.d0/mt**2/mBL**2
c     C /(mBL**2-MST(1)**2)**2*(mt**2*(MST(1)**2
c     C *(MST(1)**2+2.d0*MGL**2)+mBL**2*(mBL**2-2.d0*mt**2
c     C -4.d0*MGL**2))+(MGL**2-mBL**2)*(mBL**2-MST(1)**2)
c     C *(2.d0*MGL**2-mBL**2-MST(1)**2))
c     C *PhiGL(MGL**2,mt**2,mBL**2)
c      print*,aux
c	aux=aux
c     C +(((MGL**2-MST(2)**2)**3+(MGL**2-MST(2)**2)**2*mt**2
c     C -2.d0*MGL**4*mt**2)/2.d0/mt**2/MST(2)**2+MGL**2
c     C *((MGL**2-MST(2)**2)**2+(MGL**2-MST(2)**2)*mt**2
c     C -4.d0*MGL**2*mt**2)/((MST(2)**2-MGL**2-mt**2)**2
c     C -4.d0*mt**2*MGL**2))*PhiGL(MGL**2,mt**2,mST(2)**2)
c      print*,aux,PhiGL(MGL**2,mt**2,mST(2)**2)
c	aux=aux
c     C +((MGL**2*(MGL**2-MST(1)**2)**2+MGL**2*mt**2
c     C *(MGL**2-MST(1)**2)-4.d0*MGL**4*mt**2)
c     C /((MST(1)**2-mt**2-MGL**2)**2-4.d0*mt**2*MGL**2)
c     C +(2.d0*(MGL**2-MST(1)**2)**4*(mt**2+mBL**2-MST(1)**2)
c     C +2.d0*(MGL**2-MST(1)**2)**3*mt**2*(3.d0*mt**2-MGL**2)
c     C +2.d0*mt**4*(mBL**2-MST(1)**2)*(4.d0*MGL**2
c     C *(MGL**2-MST(1)**2)+4.d0*MGL**4-4.d0*mt**2
c     C *(MGL**2-MST(1)**2)+4.d0*MGL**2*mt**2-mt**4)
c     C +2.d0*mt**2*(MGL**2-MST(1)**2)*((mBL**2-MST(1)**2)**2
c     C -4.d0*MGL**2*(mBL**2-MST(1)**2)-4.d0*mt**2
c     C *(mBL**2-MST(1)**2)-7.d0*mt**2*MGL**2+3.d0*mt**4)
c     C +mt**2*(mBL**2-MST(1)**2)**2*(4.d0*MGL**2
c     C *(MGL**2-MST(1)**2)**2-2.d0*MGL**4+5.d0*mt**2
c     C *(MGL**2-MST(1)**2)-4.d0*MGL**2*mt**2+3.d0*mt**4)
c     C +2.d0*mt**4*(4.d0*MGL**4*(MGL**2-MST(1)**2)
c     C -7.d0*mt**2*MGL**2*(MGL**2-MST(1)**2)+4.d0*mt**2
c     C *MGL**4+mt**4*(MGL**2-MST(1)**2)-MGL**2*mt**4))
c     C /2.d0/mt**2/MST(1)**2/(mBL**2-MST(1)**2)**2)
c     C *PhiGL(MGL**2,mt**2,MST(1)**2)

c	print*,DrhoNP2L,aux,PhiGL(MGL**2,mt**2,MST(1)**2),MGL,MST(2)
c	DrhoNP2L=DrhoNP2L+Gmu*alSMZ/(4.d0*dsqrt(2.d0)*Pi**3)*aux
c	print*,DrhoNP2L

C	- O(Yt^2,Yb^2,Yt.Yb): neglected

C	 5) DrNMSSM "2L"
C	- 1st order
	DrNMSSM=DrNMSSM+drHiggs+drSUSY

C	- 2nd order
	DrNMSSM=DrNMSSM+C2TW/S2TW*Dalph*DrhoNP1L
     C                 -C2TW/S2TW*DrhoNP2L
     C                 +(drHiggs+drSUSY)**2

C	- MW
	MW=dsqrt(MZ**2/2.d0*(1.d0+dsqrt(1.d0-4*Pi*alpha
     C        /(dsqrt(2.d0)*Gmu*MZ**2)*(1.d0+DrNMSSM))))
	ENDDO

C	- error estimate
	aux=min(MST(1),MST(2),MSB(1),MSB(2),MUL,MUR,MDL,MDR)
	dMW=dMWSM+Pi*alpha/(dsqrt(2.d0)*Gmu)
     C         *dabs((1.d0-MW/MWSM)*(-8.212643109d-3*(MW-MWSM)
     C -1.425198118d-3*(MW-MWSM)**2-1.918891636d-4*(MW-MWSM)**3))
     C +10.d-3/Pi*(3.d0*Pi/2.d0-datan((aux-500.d0)/100.d0))
c	dMW=dsqrt(MZ**2/2.d0*(1.d0+dsqrt(1.d0-4*Pi*alpha
c     C        /(dsqrt(2.d0)*Gmu*MZ**2)*(1.d0+DrSUSY))))
c     c -dsqrt(MZ**2/2.d0*(1.d0+dsqrt(1.d0-4*Pi*alpha
c     C        /(dsqrt(2.d0)*Gmu*MZ**2)*(1.d0))))
        dumw=dsqrt(0.008d0**2+(0.1d0*(MW-MWSM))**2+0.00624d0**2)
c	 print*,MW,dumw,MWSM,dMW,dMWSM
        Drv=DrNMSSM-Dalph
c        print*,Drv,DrNMSSM,Dalph
        IF(MWFLAG.GT.0)
     .   PROB(87)=DDIM((MW-MWEXP)/dumw,2d0)-DDIM(-2d0,(MW-MWEXP)/dumw)

	return

	end
C**********************************************************************


	SUBROUTINE BRZLepLep(PAR,PROB)

	implicit none
	INTEGER i,j,k
	DOUBLE PRECISION PAR(*),PROB(*)
	DOUBLE PRECISION C2TW,S2TW,g2,pi,aux,mHSM,ALEMMZ,DrSM
	DOUBLE PRECISION auxiVH,auxiAH,auxiVS,auxiAS,alsmz,dalsmz
	DOUBLE PRECISION rhoSM,MZ,dMZ,alpha,Dalph
	DOUBLE PRECISION deltaS2TWeffsm,deltarhosm,dDalph
	DOUBLE PRECISION deltagVeffSM,deltagAeffSM,Gmu
	DOUBLE PRECISION gVeffSM,gAeffSM,mt,dmt
	DOUBLE PRECISION gVefftauNMSSM,gAefftauNMSSM
	DOUBLE PRECISION GamZllMZV,GamZllMZA,g1
	DOUBLE PRECISION dgVtauHiggs,dgAtauHiggs
	DOUBLE PRECISION dgVtauSusy,dgAtauSusy
	DOUBLE PRECISION gVZleplep,gAZleplep
	DOUBLE PRECISION YukTau,YukMuon,YukE,B0Yuk
	DOUBLE PRECISION SigmatauL,SigmatauR
	DOUBLE PRECISION TOTDWZmin,TOTDWZmax,TOTDWZ
	DOUBLE PRECISION BRZTTexpmin,BRZTTexpmax
	DOUBLE PRECISION SigmaGamZHiggs
	DOUBLE PRECISION SigmaWHiggs,SigmaZHiggs
	DOUBLE PRECISION SigmaZPrHiggs,SigmaGamPrHiggs
	DOUBLE PRECISION SigmaGamZSusy
	DOUBLE PRECISION SigmaWSusy,SigmaZSusy
	DOUBLE PRECISION SigmaZPrSusy,SigmaGamPrSusy
	DOUBLE PRECISION C12,C23,C24zdec,B1zdec,B0zdec,CC0,ME
	DOUBLE PRECISION dgVEHiggs,dgAEHiggs,dgVESusy,dgAESusy
	DOUBLE PRECISION SigmaEL,SigmaER
	DOUBLE PRECISION gVeffENMSSM,gAeffENMSSM
	DOUBLE PRECISION dgVtauHiggsgauge,dgAtauHiggsgauge
	DOUBLE PRECISION dgVtauSUSYgauge,dgAtauSUSYgauge
	DOUBLE PRECISION dgVEHiggsgauge,dgAEHiggsgauge
	DOUBLE PRECISION dgVESUSYgauge,dgAESUSYgauge
	DOUBLE PRECISION deltadgVtauHiggs,deltadgAtauHiggs
	DOUBLE PRECISION deltadgVtauSUSY,deltadgAtauSUSY
	DOUBLE PRECISION deltadgVEHiggs,deltadgAEHiggs
	DOUBLE PRECISION deltadgVESUSY,deltadgAESUSY
	DOUBLE PRECISION deltagVefftauNMSSM,deltagAefftauNMSSM
	DOUBLE PRECISION deltagVeffENMSSM,deltagAeffENMSSM
c	DOUBLE PRECISION vertexcorr,selfen,M1,M2,M3,M4,D27,DD0
c     C SigmaWrenHiggs,SigmaWrenSusy

	DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2)
	DOUBLE PRECISION CMASS
	DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
	DOUBLE PRECISION MSU(2),MSD(2),MSE(2),MST(2),MSB(2),MSL(2)
	DOUBLE PRECISION UT(2,2),UB(2,2),UL(2,2),DELT(2,2)
	DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
	DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
	DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	DOUBLE PRECISION TANB,COSB,SINB
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	DOUBLE PRECISION MW,dumw,dMW,DrNMSSM,MWSM,dMWSM,decztt,
     C   deltadecztt,deczee,deltadeczee,BRZTauTau,BRZTTmin,
     C   BRZTTmax,ratio,deltaratio,S2TWeffTau,deltaS2TWeffTau,
     c   S2TWeffSM
	DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0

	COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
	COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
	COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C		MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C		CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	COMMON/EWPO/MW,dumw,dMW,DrNMSSM,MWSM,dMWSM,decztt,
     C   deltadecztt,deczee,deltadeczee,BRZTauTau,BRZTTmin,
     C   BRZTTmax,ratio,deltaratio,S2TWeffTau,deltaS2TWeffTau,
     c   S2TWeffSM
        COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0

C	0 - Constants and experimental data

C  Pi
	Pi=4.d0*datan(1.d0)

C  fine structure constant
	alpha=1/137.03599911d0
	ALEMMZ=ALEMMZ0 !1.d0/127.92d0
C  exp. value for G_Fermi(muon decay)
	Gmu=GF0  !1.16637d-5
C  exp. value for mt
	mt=MT0  !172.d0
	dmt=0.45d0
C  exp. value for alpha_S(MZ)
	alSMZ=ALSMZ0 !.1184d0
	dalSMZ=1.d-4
C  exp. value for MZ
	MZ=MZ0 !91.1876d0
	dMZ=2.1d-3
C  Dalpha(lept.+hadr.)
	Dalph=0.0591579  !0.05907d0
	dDalph=1d-4   !4.d-4
C Weinberg angle
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW
	g2=4.d0*sqrt(2.d0)*Gmu*MW**2
	g1=g2*S2TW/C2TW

C	experimental data
	ME=0.5109989d-3
	TOTDWZ=2.4952d0
	TOTDWZmin=2.4929d0
	TOTDWZmax=2.4975d0

	BRZTTexpmin=0.03354d0
	BRZTTexpmax=0.03386d0

C  Trig. Functions of Betab 
        TANB=PAR(3)
        sinb=tanb/dsqrt(1.d0+tanb**2)
        cosb=sinb/tanb

C  Sfermion sector
        MSU(1)=MUL
        MSU(2)=MUR
        MSD(1)=MDL
        MSD(2)=MDR
        MSE(1)=MLL
        MSE(2)=MLR

        MST(1)=dsqrt(RMST1)
	MST(2)=dsqrt(RMST2)
        UT(1,1)=CST
        UT(1,2)=+dsqrt(1.d0-CST**2)
        UT(2,1)=-dsqrt(1.d0-CST**2)
        UT(2,2)=CST

        MSB(1)=dsqrt(RMSB1)
	MSB(2)=dsqrt(RMSB2)
        UB(1,1)=CSB
        UB(1,2)=+dsqrt(1.d0-CSB**2)
        UB(2,1)=-dsqrt(1.d0-CSB**2)
        UB(2,2)=CSB

        MSL(1)=MSL1
	MSL(2)=MSL2
        UL(1,1)=CSL
        UL(1,2)=+dsqrt(1.d0-CSL**2)
        UL(2,1)=-dsqrt(1.d0-CSL**2)
        UL(2,2)=CSL

        DELT(1,1)=1.d0
        DELT(1,2)=0.d0
        DELT(2,1)=0.d0
        DELT(2,2)=1.d0
	
C	tree-level Zll-coupling (without g2/(2*CW),included in GF) 
        gVZleplep=(-1.d0/2.d0+2.d0*S2TW)
        gAZleplep=(-1.d0/2.d0)

C*********************************************************************
C*********************************************************************
C	I.) SM-Input-Paramters to calculate the BR(Z->tau antitatu)
C		(the same for Z -> e- e+)
C*********************************************************************
C*********************************************************************

C Chosen SM Higgs mass
	mHSM=125.5d0

C		 from: Arbuzov,Awramik,.. 
C			Comput.Phys.Commun.174:728-758,2006 
	rhosm=1.005165d0
	deltarhosm=2.d-4

C		calculated from: Arbuzov,Awramik,.. 
C			Comput.Phys.Commun.174:728-758,2006 
C	rhosm=1.005286905d0

C		from:Arbuzov,Awramik,.. 
C			Comput.Phys.Commun.174:728-758,2006 
	S2TWeffsm=0.2312527d0+4.729d-4*dlog(mHSM/100.d0)
     C +2.07d-5*dlog(mHSM/100.d0)**2+3.85d-6*dlog(mHSM/100.d0)**4
     C -1.85d-6*((mHSM/100.d0)**2-1.d0)
     C +2.07d-2*(Dalph/0.05907d0-1.d0)
     C -2.851d-3*((mt/178.d0)**2-1.d0)
     C +1.82d-4*((mt/178.d0)**2-1.d0)**2
     C -9.74d-6*((mt/178.d0)**2-1.d0)*(mHSM/100.d0-1.d0)
     C +3.98d-4*(alSMZ/0.117d0-1.d0)
     C +6.55d-1*(MZ/91.1876d0-1.d0)

C param. error / mt
	aux=-2.851d-3*(((mt+dmt)/178.d0)**2-(mt/178.d0)**2)
     C +1.82d-4*((((mt+dmt)/178.d0)**2-1.d0)**2
     C               -((mt/178.d0)**2-1.d0)**2)
     C -9.74d-6*(((mt+dmt)/178.d0)**2-(mt/178.d0)**2)
     C              *(mHSM/100.d0-1.d0)
	deltaS2TWeffsm=dabs(aux)
C param. error / alpha_S(MZ)
	aux=3.98d-4*((alSMZ+dalSMZ)/0.117d0-alSMZ/0.117d0)
	deltaS2TWeffsm=dsqrt(deltaS2TWeffsm**2+aux**2)
C param. error / MZ
	aux=6.55d-1*((MZ+dMZ)/91.1876d0-MZ/91.1876d0)
	deltaS2TWeffsm=dsqrt(deltaS2TWeffsm**2+aux**2)
C parm. error / Dalpha
	aux=2.07d-2*((Dalph+dDalph)/0.05907d0-Dalph/0.05907d0)
	deltaS2TWeffsm=dsqrt(deltaS2TWeffsm**2+aux**2)
C error higher orders
	deltaS2TWeffsm=4.9d-5+2.d0*deltaS2TWeffsm

	DrSM=dsqrt(2.d0)*Gmu/(Pi*alpha)*MWSM**2
     C         *(1.d0-MWSM**2/MZ**2)-1.d0
	gVeffSM=(2.d0*S2TWeffsm-1.d0/2.d0
     C                        +2.d0*(MWSM**2-MW**2)/MZ**2)
     C                         *dsqrt(rhosm/(1.d0-DrSM))
	deltagVeffSM=dabs((2.d0*S2TWeffsm-1.d0/2.d0)/2.d0
     C                        /dsqrt(rhosm))*deltarhosm
     C +2.d0*dsqrt(rhosm)*(deltaS2TWeffsm+2.d0*MW*dMW/MZ**2)
c	deltagVeffSM=deltagVeffSM/dsqrt(1.d0-DrNMSSM)
c     C   +dabs(gVeffSM/(1.d0-DrNMSSM))*dsqrt(2.d0)*Gmu
c     C        *MW*dMW/Pi/alpha

	gAeffSM=(-1.d0/2.d0)*dsqrt(rhosm/(1.d0-DrSM))
	deltagAeffSM=dabs(1.d0/4.d0/dsqrt(rhosm))*deltarhosm
     c +dabs(gAeffSM/(1.d0-DrNMSSM))*dsqrt(2.d0)*Gmu/(Pi*alpha)*
     c dabs(MWSM*(1.d0-2.d0*MWSM**2/MZ**2)*dMWSM
     C         -MW*(1.d0-2.d0*MW**2/MZ**2)*dMW)
c     C                   /dsqrt(1.d0-DrNMSSM)
c     C +dabs(gAeffSM/(1.d0-DrNMSSM))*dsqrt(2.d0)*Gmu
c     C        *MW*dMW/Pi/alpha

C*********************************************************************
C*********************************************************************
C	II.) New Contribution (-old) in NMSSM Higgs Sector
C*********************************************************************
C*********************************************************************
C	1.) Modified Yukawa Coupling
C*********************************************************************

	YukE=dsqrt(g2/2.d0)*ME/(MW*cosb)
	YukTau=dsqrt(g2/2.d0)*mtau/(MW*cosb)
	YukMuon=dsqrt(g2/2.d0)*MMUON/(MW*cosb)

	aux=0.d0
	DO i=1,2
	 aux=aux+dsqrt(2.d0)*V(i,1)*U(i,2)*mcha(i)
     C                            *B0Yuk(mcha(i),MNL)
	ENDDO
	DO i=1,5
	 aux=aux+mneu(i)*neu(i,4)*(
     C 2.d0*dsqrt(g1/g2)*neu(i,1)*B0Yuk(mneu(i),MLR)
     C -(neu(i,1)*dsqrt(g1/g2)+neu(i,2))*B0Yuk(mneu(i),MLL))
	ENDDO
	aux=aux*g2/(32.d0*pi**2*cosb*MW)
	YukE=YukE/(1.d0+aux)
	YukMuon=YukMuon/(1.d0+aux)

	aux=0.d0
	DO i=1,2
	 aux=aux+dsqrt(2.d0)*V(i,1)*U(i,2)*mcha(i)
     C                            *B0Yuk(mcha(i),MSNT)
	ENDDO
	DO i=1,5
	DO j=1,2
	 aux=aux+dsqrt(2.d0/g2)/YukTau*mneu(i)*
     C (Yuktau*UL(j,1)*neu(i,4)+dsqrt(2.d0*g1)*UL(j,2)
     C  *neu(i,1))*(YukTau*neu(i,4)*UL(j,2)-dsqrt(g1/2.d0)
     C  *neu(i,1)*UL(j,1)-dsqrt(g2/2.d0)*neu(i,2)*UL(j,1))
     C *B0Yuk(mneu(i),MSL(j))
	ENDDO
	ENDDO
	aux=aux*g2/(32.d0*pi**2*cosb*MW)
	YukTau=Yuktau/(1.d0+aux)
 
C*********************************************************************
C	Higgs-Contribution 	
C*********************************************************************
C	2.) Calculation of effective axial and vector couplings NP 
C	(-SM-Contributions + NMSSM-contributions) for Z->Tau- Tau-
C*********************************************************************	

C	-> gauge part

	dgVtauHiggsgauge=((-1.d0/2.d0+2.d0*S2TW)*(-1.d0/2.d0
     C	*SigmaZPrHiggs(MZ,tanb,MW)+1.d0/2.d0*SigmaGamPrHiggs(0.d0,MW)
     C	   +(C2TW-S2TW)/dsqrt(S2TW*C2TW)*SigmaGamZHiggs(0.d0,MW)
     C	         /MZ**2-1.d0/2.d0*(C2TW-S2TW)/S2TW
     C	         *(SigmaZHiggs(MZ,tanb,MW)/MZ**2
     C		  -SigmaWHiggs(MW,tanb,MW)/MW**2))
     C		 +2.d0*C2TW*(SigmaZHiggs(MZ,tanb,MW)/MZ**2
     C		            -SigmaWHiggs(MW,tanb,MW)/MW**2
     C		 -dsqrt(S2TW/C2TW)
     C	     *(SigmaGamZHiggs(MZ,MW)+SigmaGamZHiggs(0.d0,MW))/MZ**2)
     C	+1.d0/2.d0*dsqrt(C2TW/S2TW)*SigmaGamZHiggs(0.d0,MW)/MZ**2)
	dgAtauHiggsgauge=(-1.d0/2.d0*(-SigmaZPrHiggs(MZ,tanb,MW)/2.d0
     C	         +SigmaGamPrHiggs(0.d0,MW)/2.d0
     C	         +(C2TW-S2TW)/dsqrt(S2TW*C2TW)*SigmaGamZHiggs(0.d0,MW)
     C	         /MZ**2-1.d0/2.d0*(C2TW-S2TW)/S2TW
     C	         *(SigmaZHiggs(MZ,tanb,MW)/MZ**2
     C		  -SigmaWHiggs(MW,tanb,MW)/MW**2))
     C	+1.d0/2.d0*dsqrt(C2TW/S2TW)*SigmaGamZHiggs(0.d0,MW)/MZ**2)


C	-> "fermion part"

C	 a) Vertexcorrections
	GAMZllMZV=0.d0
	GAMZllMZA=0.d0

C		->SM-contribution (M(Higgs)=150 GeV)
	GAMZllMZV=GAMZllMZV-((-YukE**2*cosb**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0+2.d0*S2TW)*(1.d0/2.d0-2.d0
     C		*C24zdec(ME,MZ,ME,mHSM,ME)-MZ**2
     C		*(C12(ME,MZ,ME,mHSM,ME)
     C		+C23(ME,MZ,ME,mHSM,ME))))

	GAMZllMZA=GAMZllMZA-((YukE**2*cosb**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0)*(1.d0/2.d0-2.d0
     C		*C24zdec(ME,MZ,ME,mHSM,ME)-MZ**2
     C		*(C12(ME,MZ,ME,mHSM,ME)
     C		+C23(ME,MZ,ME,mHSM,ME))))

C		-> Goldstone-Contribution
	GAMZllMZA=GAMZllMZA-2.d0*((YukE**2*cosb**2)/(8.d0*pi**2)
     C		*(-1.d0/2.d0)*C24zdec(ME,MZ,mHSM,ME,MZ))

C		->NMSSM (CP-odd Contribution) 
	aux=0.d0
	do i=1,2
	aux=aux+((-(sinb*PCOMP(i,1))**2*YukTau**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0+2.d0*S2TW)*
     C		(1.d0/2.d0-2.d0*C24zdec(MTAU,MZ,MTAU,PMASS(i),MTAU)
     C		-MZ**2*(C12(MTAU,MZ,MTAU,PMASS(i),MTAU)+
     C 		C23(MTAU,MZ,MTAU,PMASS(i),MTAU))))
  	enddo
	GAMZllMZV=GAMZllMZV+aux

	aux=0.d0
	do i=1,2
	aux=aux+((YukTau**2*(sinb*PCOMP(i,1))**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0)*(1.d0/2.d0-2.d0
     C		*C24zdec(MTAU,MZ,MTAU,PMASS(i),MTAU)
     C		-MZ**2*(C12(MTAU,MZ,MTAU,PMASS(i),MTAU)
     C		+C23(MTAU,MZ,MTAU,PMASS(i),MTAU))))
	enddo
	GAMZllMZA=GAMZllMZA+aux

C		->NMSSM (CP-even Contribution)
	aux=0.d0
	do i=1,3
	aux=aux+((-YukTau**2*(SCOMP(i,2))**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0+2.d0*S2TW)*
     C		(1.d0/2.d0-2.d0*C24zdec(MTAU,MZ,MTAU,SMASS(i),MTAU)
     C		-MZ**2*(C12(MTAU,MZ,MTAU,SMASS(i),MTAU)+
     C 		C23(MTAU,MZ,MTAU,SMASS(i),MTAU))))
	enddo
	GAMZllMZV=GAMZllMZV+aux

	aux=0.d0
	do i=1,3
		aux=aux+((YukTau**2*(SCOMP(i,2))**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0)*(1.d0/2.d0-2.d0
     C		*C24zdec(MTAU,MZ,MTAU,SMASS(i),MTAU)
     C		-MZ**2*(C12(MTAU,MZ,MTAU,SMASS(i),MTAU)
     C		+C23(MTAU,MZ,MTAU,SMASS(i),MTAU))))
	enddo
	GAMZllMZA=GAMZllMZA+aux

C		->NMSSM with Zah-Vertex 
C		   (with summation over CP-even Higgs and CP-odd)
	aux=0.d0
	do j=1,2
	do i=1,3
	aux=aux+2.d0*((YukTau**2)/(8.d0*pi**2)*(-1.d0/2.d0)*sinb
     C	*PCOMP(j,1)*SCOMP(i,2)*(SCOMP(i,2)*sinb*PCOMP(j,1)
     C	-SCOMP(i,1)*cosb*PCOMP(j,1))
     C	*C24zdec(MTAU,MZ,SMASS(i),MTAU,PMASS(j))) 
	enddo
	enddo
	GAMZllMZA=GAMZllMZA+aux

C		->NMSSM (with ZGoh-Vertex)
C		   (with summation over CP-even Higgs)
	aux=0.d0
	do i=1,3
	aux=aux+2.d0*(YukTau**2)/(8.d0*pi**2)*(-1.d0/2.d0)
     C		*cosb*SCOMP(i,2)*(SCOMP(i,2)*cosb+SCOMP(i,1)*sinb)
     C		*C24zdec(MTAU,MZ,SMASS(i),MTAU,MZ) 
	enddo
	GAMZllMZA=GAMZllMZA+aux

C		-> NMSSM (with two charged Higgs in the loop)
	GAMZllMZV=GAMZllMZV+((YukTau**2*sinb**2*(S2TW-C2TW))
     C		/(16.d0*pi**2)*C24zdec(MTAU,MZ,CMASS,0.d0,CMASS))

	GAMZllMZA=GAMZllMZA+((-YukTau**2*sinb**2*(S2TW-C2TW))
     C		/(16.d0*pi**2)*C24zdec(MTAU,MZ,CMASS,0.d0,CMASS))	

C		-> NMSSM (with one charged Higgs in the loop)
	GAMZllMZV=GAMZllMZV+((-YukTau**2*sinb**2)/(32.d0*pi**2)
     C		*(1.d0/2.d0-2.d0*C24zdec(MTAU,MZ,0.d0,CMASS,0.d0)
     C		-MZ**2*(C23(MTAU,MZ,0.d0,CMASS,0.d0)
     C		+C12(MTAU,MZ,0.d0,CMASS,0.d0))))

	GAMZllMZA=GAMZllMZA+((YukTau**2*sinb**2)/(32.d0*pi**2)
     C		*(1.d0/2.d0-2.d0*C24zdec(MTAU,MZ,0.d0,CMASS,0.d0)
     C		-MZ**2*(C23(MTAU,MZ,0.d0,CMASS,0.d0)
     C		+C12(MTAU,MZ,0.d0,CMASS,0.d0))))


C*********************************************************************
C	 b) selfenergies 
C		->SM-contribution
	GAMZllMZV=GAMZllMZV-((YukTau**2*cosb**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0+2.d0*S2TW)*B1zdec(MTAU,MTAU,mHSM))

	GAMZllMZA=GAMZllMZA-((YukTau**2*cosb**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0)*B1zdec(MTAU,MTAU,mHSM))

C		->NMSSM (CP-odd contributions)
	aux=0.d0
	do i=1,2
	aux=aux+((YukTau**2*(sinb*PCOMP(i,1))**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0+2.d0*S2TW)*B1zdec(MTAU,MTAU,PMASS(i)))
	enddo
	GAMZllMZV=GAMZllMZV+aux

	aux=0.d0
	do i=1,2
	aux=aux+((YukTau**2*(sinb*PCOMP(i,1))**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0)*B1zdec(MTAU,MTAU,PMASS(i)))
	enddo
	GAMZllMZA=GAMZllMZA+aux

C		->NMSSM (CP-even contribution)
	aux=0.d0
	do i=1,3
	aux=aux+((YukTau**2*(SCOMP(i,2))**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0+2.d0*S2TW)*B1zdec(MTAU,MTAU,SMASS(i)))
	enddo
	GAMZllMZV=GAMZllMZV+aux

	aux=0.d0
	do i=1,3
	aux=aux+((YukTau**2*(SCOMP(i,2))**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0)*B1zdec(MTAU,MTAU,SMASS(i)))
	enddo
	GAMZllMZA=GAMZllMZA+aux

C		-> NMSSM (charged Higgs)
	GAMZllMZV=GAMZllMZV+((-YukTau**2*S2TW*sinb**2)/(16.d0*pi**2)
     C		*(B0zdec(MTAU,CMASS,0.d0)+B1zdec(MTAU,CMASS,0.d0)))

	GAMZllMZA=GAMZllMZA+((YukTau**2*S2TW*sinb**2)/(16.d0*pi**2)
     C		*(B0zdec(MTAU,CMASS,0.d0)+B1zdec(MTAU,CMASS,0.d0)))

C	TOGETHER:
	dgVtauHiggs=dgVtauHiggsgauge+GAMZllMZV
	dgAtauHiggs=dgAtauHiggsgauge+GAMZllMZA
	auxiVH=GAMZllMZV
	auxiAH=GAMZllMZA

C	theoretical uncertainty of the Higgs-Sector:
	deltadgVtauHiggs=0.3d0*dabs(dgVtauHiggsgauge)
     C		        +0.1d0*dabs(GAMZllMZV)
	deltadgAtauHiggs=0.3d0*dabs(dgAtauHiggsgauge)
     C			+0.1d0*dabs(GAMZllMZA)


C*********************************************************************
C*********************************************************************
C	3.) Calculation of effective axial and vector couplings NP 
C	(-SM-Contributions + NMSSM-contributions) for Z -> e- e+
C*********************************************************************	

C	-> gauge part

	dgVEHiggsgauge=((-1.d0/2.d0+2.d0*S2TW)*(-1.d0/2.d0
     C	*SigmaZPrHiggs(MZ,tanb,MW)+1.d0/2.d0*SigmaGamPrHiggs(0.d0,MW)
     C	   +(C2TW-S2TW)/dsqrt(S2TW*C2TW)*SigmaGamZHiggs(0.d0,MW)
     C	         /MZ**2-1.d0/2.d0*(C2TW-S2TW)/S2TW
     C	         *(SigmaZHiggs(MZ,tanb,MW)/MZ**2
     C		  -SigmaWHiggs(MW,tanb,MW)/MW**2))
     C		 +2.d0*C2TW*(SigmaZHiggs(MZ,tanb,MW)/MZ**2
     C		            -SigmaWHiggs(MW,tanb,MW)/MW**2
     C		 -dsqrt(S2TW/C2TW)
     C	       *(SigmaGamZHiggs(MZ,MW)+SigmaGamZHiggs(0.d0,MW))/MZ**2)
     C	+1.d0/2.d0*dsqrt(C2TW/S2TW)*SigmaGamZHiggs(0.d0,MW)/MZ**2)

	dgAEHiggsgauge=((-1.d0/2.d0)*(-SigmaZPrHiggs(MZ,tanb,MW)/2.d0
     C	         +1.d0/2.d0*SigmaGamPrHiggs(0.d0,MW)
     C	         +(C2TW-S2TW)/dsqrt(S2TW*C2TW)*SigmaGamZHiggs(0.d0,MW)
     C	         /MZ**2-1.d0/2.d0*(C2TW-S2TW)/S2TW
     C	         *(SigmaZHiggs(MZ,tanb,MW)/MZ**2
     C		  -SigmaWHiggs(MW,tanb,MW)/MW**2))
     C	+1.d0/2.d0*dsqrt(C2TW/S2TW)*SigmaGamZHiggs(0.d0,MW)/MZ**2)

C	-> "fermion part"

C		a) Vertexcorrections

	GAMZllMZV=0.d0
	GAMZllMZA=0.d0

C		->SM-contribution (M(Higgs)=150 GeV)
	GAMZllMZV=GAMZllMZV-((-YukE**2*cosb**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0+2.d0*S2TW)*(1.d0/2.d0-2.d0
     C		*C24zdec(ME,MZ,ME,mHSM,ME)-MZ**2
     C		*(C12(ME,MZ,ME,mHSM,ME)
     C		+C23(ME,MZ,ME,mHSM,ME))))

	GAMZllMZA=GAMZllMZA-((YukE**2*cosb**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0)*(1.d0/2.d0-2.d0
     C		*C24zdec(ME,MZ,ME,mHSM,ME)-MZ**2
     C		*(C12(ME,MZ,ME,mHSM,ME)
     C		+C23(ME,MZ,ME,mHSM,ME))))

C		-> Goldstone-Contribution
	GAMZllMZA=GAMZllMZA-2.d0*((YukE**2*cosb**2)/(8.d0*pi**2)
     C		*(-1.d0/2.d0)*C24zdec(ME,MZ,mHSM,ME,MZ))

C		->NMSSM (CP-odd Contribution) 
	aux=0.d0
	do i=1,2
	aux=aux+((-(sinb*PCOMP(i,1))**2*YukE**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0+2.d0*S2TW)*
     C		(1.d0/2.d0-2.d0*C24zdec(ME,MZ,ME,PMASS(i),ME)
     C		-MZ**2*(C12(ME,MZ,ME,PMASS(i),ME)+
     C 		C23(ME,MZ,ME,PMASS(i),ME))))
  	enddo
	GAMZllMZV=GAMZllMZV+aux

	aux=0.d0
	do i=1,2
	aux=aux+((YukE**2*(sinb*PCOMP(i,1))**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0)*(1.d0/2.d0-2.d0
     C		*C24zdec(ME,MZ,ME,PMASS(i),ME)
     C		-MZ**2*(C12(ME,MZ,ME,PMASS(i),ME)
     C		+C23(ME,MZ,ME,PMASS(i),ME))))

	enddo
	GAMZllMZA=GAMZllMZA+aux

C		->NMSSM (CP-even Contribution)
	aux=0.d0
	do i=1,3
	aux=aux+((-YukE**2*(SCOMP(i,2))**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0+2.d0*S2TW)*
     C		(1.d0/2.d0-2.d0*C24zdec(ME,MZ,ME,SMASS(i),ME)
     C		-MZ**2*(C12(ME,MZ,ME,SMASS(i),ME)+
     C 		C23(ME,MZ,ME,SMASS(i),ME))))
	enddo
	GAMZllMZV=GAMZllMZV+aux

	aux=0.d0
	do i=1,3
		aux=aux+((YukE**2*(SCOMP(i,2))**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0)*(1.d0/2.d0-2.d0
     C		*C24zdec(ME,MZ,ME,SMASS(i),ME)
     C		-MZ**2*(C12(ME,MZ,ME,SMASS(i),ME)
     C		+C23(ME,MZ,ME,SMASS(i),ME))))
	enddo
	GAMZllMZA=GAMZllMZA+aux

C		->NMSSM with Zah-Vertex 
C		   (with summation over CP-even Higgs and CP-odd)
	aux=0.d0
	do j=1,2
	do i=1,3
	aux=aux+2.d0*((YukE**2)/(8.d0*pi**2)*(-1.d0/2.d0)*sinb
     C	*PCOMP(j,1)*SCOMP(i,2)*(SCOMP(i,2)*sinb*PCOMP(j,1)
     C	-SCOMP(i,1)*cosb*PCOMP(j,1))
     C	*C24zdec(ME,MZ,SMASS(i),ME,PMASS(j))) 
	enddo
	enddo
	GAMZllMZA=GAMZllMZA+aux

C		->NMSSM (with ZGoh-Vertex)
C		   (with summation over CP-even Higgs)
	aux=0.d0
	do i=1,3
	aux=aux+2.d0*(YukE**2)/(8.d0*pi**2)*(-1.d0/2.d0)
     C		*cosb*SCOMP(i,2)*(SCOMP(i,2)*cosb+SCOMP(i,1)*sinb)
     C		*C24zdec(ME,MZ,SMASS(i),ME,MZ) 
	enddo
	GAMZllMZA=GAMZllMZA+aux

C		-> NMSSM (with two charged Higgs in the loop)
	GAMZllMZV=GAMZllMZV+((YukE**2*sinb**2*(S2TW-C2TW))
     C		/(16.d0*pi**2)*C24zdec(ME,MZ,CMASS,0.d0,CMASS))

	GAMZllMZA=GAMZllMZA+((-YukE**2*sinb**2*(S2TW-C2TW))
     C		/(16.d0*pi**2)*C24zdec(ME,MZ,CMASS,0.d0,CMASS))	

C		-> NMSSM (with one charged Higgs in the loop)
	GAMZllMZV=GAMZllMZV+((-YukE**2*sinb**2)/(32.d0*pi**2)
     C		*(1.d0/2.d0-2.d0*C24zdec(ME,MZ,0.d0,CMASS,0.d0)
     C		-MZ**2*(C23(ME,MZ,0.d0,CMASS,0.d0)
     C		+C12(ME,MZ,0.d0,CMASS,0.d0))))

	GAMZllMZA=GAMZllMZA+((YukE**2*sinb**2)/(32.d0*pi**2)
     C		*(1.d0/2.d0-2.d0*C24zdec(ME,MZ,0.d0,CMASS,0.d0)
     C		-MZ**2*(C23(ME,MZ,0.d0,CMASS,0.d0)
     C		+C12(ME,MZ,0.d0,CMASS,0.d0))))


C*********************************************************************

C	 b) selfenergies 

C		->SM-contribution
	GAMZllMZV=GAMZllMZV-((YukE**2*cosb**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0+2.d0*S2TW)*B1zdec(0.d0,0.d0,mHSM))

	GAMZllMZA=GAMZllMZA-((YukE**2*cosb**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0)*B1zdec(ME,ME,mHSM))

C		->NMSSM (CP-odd contributions)
	aux=0.d0
	do i=1,2
	aux=aux+((YukE**2*(sinb*PCOMP(i,1))**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0+2.d0*S2TW)*B1zdec(ME,ME,PMASS(i)))
	enddo
	GAMZllMZV=GAMZllMZV+aux

	aux=0.d0
	do i=1,2
	aux=aux+((YukE**2*(sinb*PCOMP(i,1))**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0)*B1zdec(ME,ME,PMASS(i)))
	enddo
	GAMZllMZA=GAMZllMZA+aux

C		->NMSSM (CP-even contribution)
	aux=0.d0
	do i=1,3
	aux=aux+((YukE**2*(SCOMP(i,2))**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0+2.d0*S2TW)*B1zdec(ME,ME,SMASS(i)))
	enddo
	GAMZllMZV=GAMZllMZV+aux

	aux=0.d0
	do i=1,3
	aux=aux+((YukE**2*(SCOMP(i,2))**2)/(32.d0*pi**2)
     C		*(-1.d0/2.d0)*B1zdec(ME,ME,SMASS(i)))
	enddo
	GAMZllMZA=GAMZllMZA+aux

C		-> NMSSM (charged Higgs)
	GAMZllMZV=GAMZllMZV+((-YukE**2*S2TW*sinb**2)/(16.d0*pi**2)
     C		*(B0zdec(ME,CMASS,0.d0)+B1zdec(ME,CMASS,0.d0)))

	GAMZllMZA=GAMZllMZA+((YukE**2*S2TW*sinb**2)/(16.d0*pi**2)
     C		*(B0zdec(ME,CMASS,0.d0)+B1zdec(ME,CMASS,0.d0)))

C	TOGETHER:
	dgVEHiggs=dgVEHiggsgauge+GAMZllMZV
	dgAEHiggs=dgAEHiggsgauge+GAMZllMZA
	auxiVH=0.1*dabs(GAMZllMZV-auxiVH)
	auxiAH=0.1*dabs(GAMZllMZA-auxiAH)

C	Theoretical Uncertainty of the Higgs-Sector:
	deltadgVEHiggs=0.3d0*dabs(dgVEHiggsgauge)
     C		      +0.1d0*dabs(GAMZllMZV)
	deltadgAEHiggs=0.3d0*dabs(dgAEHiggsgauge)
     C		      +0.1d0*dabs(GAMZllMZA)

C*********************************************************************
C*********************************************************************
C*********************************************************************
C	III.) SUSY Contribution 
C		(Sfermion-Chargino/Neutralino contributions
C*********************************************************************
C*********************************************************************
C	1.) Modified Yukawa Coupling (not studied yet)
C*********************************************************************

C	Yukstau=dsqrt(g2/2.d0)*MSL1/(MW*cosb)
C	Yukstau=dsqrt(g2/2.d0)*MSL2/(MW*cosb)
C	Yuksmuon=dsqrt(g2/2.d0)*MLL/(MW*cosb)
C	Yuksmuon=dsqrt(g2/2.d0)*MLR/(MW*cosb)
 
C*********************************************************************


C*********************************************************************
C	2.) Calculation of effective axial and vector couplings 
C	    SUSY-contributions for Z-> tau- tau+
C*********************************************************************	
C	-> gauge part

	dgVtauSusygauge=(gVZleplep*(-1.d0/2.d0*SigmaZPrSusy(MZ,MW)
     C	       +1.d0/2.d0*SigmaGamPrSusy(0.d0,MW)
     C	       +(C2TW-S2TW)/dsqrt(S2TW*C2TW)*SigmaGamZSusy(0.d0,MW)
     C	       /MZ**2-1.d0/2.d0*(C2TW-S2TW)/S2TW
     C	       *(SigmaZSusy(MZ,MW)/MZ**2-SigmaWSusy(MW,MW)/MW**2))
     C +2.d0*C2TW*(SigmaZSusy(MZ,MW)/MZ**2-SigmaWSusy(MW,MW)/MW**2
     C		 -dsqrt(S2TW/C2TW)
     C	      *(SigmaGamZSusy(MZ,MW)+SigmaGamZSusy(0.d0,MW))/MZ**2)
     C	      -gAZleplep*dsqrt(C2TW/S2TW)*SigmaGamZSusy(0.d0,MW)
     C		/MZ**2)

	dgAtauSusygauge=(gAZleplep*(-1.d0/2.d0*SigmaZPrSusy(MZ,MW)
     C	      +1.d0/2.d0*SigmaGamPrSusy(0.d0,MW)
     C	      +(C2TW-S2TW)/dsqrt(S2TW*C2TW)*SigmaGamZSusy(0.d0,MW)
     C	      /MZ**2-1.d0/2.d0*(C2TW-S2TW)/S2TW
     C	      *(SigmaZSusy(MZ,MW)/MZ**2-SigmaWSusy(MW,MW)/MW**2))
     C	      -gAZleplep*dsqrt(C2TW/S2TW)*SigmaGamZSusy(0.d0,MW)
     C		/MZ**2)

C	-> "fermion part"

C		a) Contributions to the Zll vertex function

C	->Vector-Part

C	Chargino corrections:
       aux=0.d0
       do i=1,2
       do j=1,2
        aux=aux+MCHA(i)*MCHA(j)*(g2*V(i,1)*V(j,1)
     C          *(U(i,1)*U(j,1)+U(i,2)*U(j,2)/2.d0-S2TW*delt(i,j))
     C          +YukTau**2*U(i,2)*U(j,2)
     C          *(V(i,1)*V(j,1)+V(i,2)*V(j,2)/2.d0-S2TW*delt(i,j)))
     C           *CC0(mtau,MZ,MCHA(j),MSNT,MCHA(i))
     C      -(g2*V(i,1)*V(j,1)
     C          *(V(i,1)*V(j,1)+V(i,2)*V(j,2)/2.d0-S2TW*delt(i,j))
     C          +YukTau**2*U(i,2)*U(j,2)
     C          *(U(i,1)*U(j,1)+U(i,2)*U(j,2)/2.d0-S2TW*delt(i,j)))
     C  *(2.d0*C24zdec(mtau,MZ,MCHA(j),MSNT,MCHA(i))
     C		-1.d0/2.d0+MZ**2
     C      *(C12(mtau,MZ,MCHA(j),MSNT,MCHA(i))
     C                       +C23(mtau,MZ,MCHA(j),MSNT,MCHA(i))))

       enddo
        aux=aux+(g2*V(i,1)**2+YukTau**2*U(i,2)**2)
     C                 *C24zdec(mtau,MZ,MSNT,MCHA(i),MSNT)
       enddo

C	Neutralino corrections:
       do i=1,5
       do j=1,5
       do k=1,2
        aux=aux-(NEU(i,4)*NEU(j,4)-NEU(i,3)*NEU(j,3))/2.d0
     C   *((NEU(i,1)*UL(k,1)*dsqrt(g1/2.d0)
     C      +NEU(i,2)*UL(k,1)*dsqrt(g2/2.d0)-YukTau*NEU(i,4)*UL(k,2))
     C     *(NEU(j,1)*UL(k,1)*dsqrt(g1/2.d0)
     C      +NEU(j,2)*UL(k,1)*dsqrt(g2/2.d0)-YukTau*NEU(j,4)*UL(k,2))
     C      -(NEU(i,1)*UL(k,2)*dsqrt(g1*2.d0)
     C                         +YukTau*NEU(i,4)*UL(k,1))
     C      *(NEU(j,1)*UL(k,2)*dsqrt(g1*2.d0)
     C                         +YukTau*NEU(j,4)*UL(k,1)))
     C   *(MNEU(i)*MNEU(j)*CC0(mtau,MZ,MNEU(j),MSL(k),MNEU(i))
     C     +2.d0*C24zdec(mtau,MZ,MNEU(j),MSL(k),MNEU(i))
     C		-1.d0/2.d0
     C     +MZ**2*(C12(mtau,MZ,MNEU(j),MSL(k),MNEU(i))
     C             +C23(mtau,MZ,MNEU(j),MSL(k),MNEU(i))))
	enddo
	enddo
	enddo

       do i=1,2
       do j=1,2
       do k=1,5
	aux=aux-2.d0*((1.d0/2.d0-S2TW)*UL(i,1)*UL(j,1)
     C                                 -S2TW*UL(i,2)*UL(j,2))
     C    *((NEU(k,1)*UL(i,1)*dsqrt(g1/2.d0)
     C      +NEU(k,2)*UL(i,1)*dsqrt(g2/2.d0)-YukTau*NEU(k,4)*UL(i,2))
     C     *(NEU(k,1)*UL(j,1)*dsqrt(g1/2.d0)
     C      +NEU(k,2)*UL(j,1)*dsqrt(g2/2.d0)-YukTau*NEU(k,4)*UL(j,2))
     C      +(NEU(k,1)*UL(i,2)*dsqrt(g1*2.d0)
     C                         +YukTau*NEU(k,4)*UL(i,1))
     C      *(NEU(k,1)*UL(j,2)*dsqrt(g1*2.d0)
     C                         +YukTau*NEU(k,4)*UL(j,1)))
     C        *C24zdec(mtau,MZ,MSL(j),MNEU(k),MSL(i))
       enddo
       enddo
       enddo	

       GamZllMZV=dsqrt(g2)/(32.d0*Pi**2*dsqrt(1.d0-S2TW))*aux

C	->Axial-Vector-Part
C	Chargino corrections:
	aux=0.d0
	do i=1,2
	do j=1,2
       aux=aux+MCHA(i)*MCHA(j)*(g2*V(i,1)*V(j,1)
     C          *(U(i,1)*U(j,1)+U(i,2)*U(j,2)/2.d0-S2TW*delt(i,j))
     C          -YukTau**2*U(i,2)*U(j,2)
     C          *(V(i,1)*V(j,1)+V(i,2)*V(j,2)/2.d0-S2TW*delt(i,j)))
     C           *CC0(mtau,MZ,MCHA(j),MSNT,MCHA(i))
     C      +(g2*V(i,1)*V(j,1)
     C          *(V(i,1)*V(j,1)+V(i,2)*V(j,2)/2.d0-S2TW*delt(i,j))
     C          -YukTau**2*U(i,2)*U(j,2)
     C          *(U(i,1)*U(j,1)+U(i,2)*U(j,2)/2.d0-S2TW*delt(i,j)))
     C  *(-2.d0*C24zdec(mtau,MZ,MCHA(j),MSNT,MCHA(i))+1.d0/2.d0-MZ**2
     C      *(C12(mtau,MZ,MCHA(j),MSNT,MCHA(i))
     C                      +C23(mtau,MZ,MCHA(j),MSNT,MCHA(i))))
	enddo
        aux=aux+(g2*V(i,1)**2-YukTau**2*U(i,2)**2)
     C                     *C24zdec(mtau,MZ,MSNT,MCHA(i),MSNT)
	enddo

C	Neutralino corrections:
	do i=1,5
	do j=1,5
	do k=1,2
	aux=aux-(NEU(i,4)*NEU(j,4)-NEU(i,3)*NEU(j,3))/2.d0
     C    *((NEU(i,1)*UL(k,1)*dsqrt(g1/2.d0)
     C      +NEU(i,2)*UL(k,1)*dsqrt(g2/2.d0)-YukTau*NEU(i,4)*UL(k,2))
     C     *(NEU(j,1)*UL(k,1)*dsqrt(g1/2.d0)
     C      +NEU(j,2)*UL(k,1)*dsqrt(g2/2.d0)-YukTau*NEU(j,4)*UL(k,2))
     C      +(NEU(i,1)*UL(k,2)*dsqrt(g1*2.d0)
     C                         +YukTau*NEU(i,4)*UL(k,1))
     C      *(NEU(j,1)*UL(k,2)*dsqrt(g1*2.d0)
     C                         +YukTau*NEU(j,4)*UL(k,1)))
     C   *(MNEU(i)*MNEU(j)*CC0(mtau,MZ,MNEU(j),MSL(k),MNEU(i))
     C     +2.d0*C24zdec(mtau,MZ,MNEU(j),MSL(k),MNEU(i))
     C		-1.d0/2.d0
     C     +MZ**2*(C12(mtau,MZ,MNEU(j),MSL(k),MNEU(i))
     C             +C23(mtau,MZ,MNEU(j),MSL(k),MNEU(i))))
	enddo	
	enddo
	enddo

	do i=1,2
	do j=1,2
	do k=1,5
	aux=aux-2.d0*((1.d0/2.d0-S2TW)*UL(i,1)*UL(j,1)
     C                                 -S2TW*UL(i,2)*UL(j,2))
     C    *((NEU(k,1)*UL(i,1)*dsqrt(g1/2.d0)
     C      +NEU(k,2)*UL(i,1)*dsqrt(g2/2.d0)-YukTau*NEU(k,4)*UL(i,2))
     C     *(NEU(k,1)*UL(j,1)*dsqrt(g1/2.d0)
     C      +NEU(k,2)*UL(j,1)*dsqrt(g2/2.d0)-YukTau*NEU(k,4)*UL(j,2))
     C      -(NEU(k,1)*UL(i,2)*dsqrt(g1*2.d0)
     C                         +YukTau*NEU(k,4)*UL(i,1))
     C      *(NEU(k,1)*UL(j,2)*dsqrt(g1*2.d0)
     C                         +YukTau*NEU(k,4)*UL(j,1)))*
     C        C24zdec(mtau,MZ,MSL(j),MNEU(k),MSL(i))
	enddo
	enddo
	enddo	

	GamZllMZA=dsqrt(g2)/(32.d0*Pi**2*dsqrt(1.d0-S2TW))*aux

C*********************************************************************

C		b) Contributions to the leptonic self-energies

C	SigmatauL  
	aux=0.d0
	do i=1,5
	do j=1,2
        aux=aux+(NEU(i,1)*UL(j,1)*dsqrt(g1/2.d0)
     C      +NEU(i,2)*UL(j,1)*dsqrt(g2/2.d0)
     C	    -YukTau*NEU(i,4)*UL(j,2))**2
     C	       *(B0zdec(mtau,MSL(j),MNEU(i))
     C		+B1zdec(mtau,MSL(j),MNEU(i)))
	enddo
	enddo

	do i=1,2
        aux=aux+g2*(V(i,1))**2
     C       *(B0zdec(mtau,MSNT,MCHA(i))
     C		+B1zdec(mtau,MSNT,MCHA(i)))
	enddo
	SigmatauL=aux/(16.d0*Pi**2)

C	SigmatauR
	aux=0.d0
	do i=1,5
	do j=1,2
	aux=aux+(-NEU(i,1)*UL(j,2)*dsqrt(g1*2.d0)
     C                         -YukTau*NEU(i,4)*UL(j,1))**2
     C  *(B0zdec(mtau,MSL(j),MNEU(i))
     C		+B1zdec(mtau,MSL(j),MNEU(i)))
	enddo
	enddo
	do i=1,2
	aux=aux+YukTau**2*(U(i,2))**2
     C  *(B0zdec(mtau,MSNT,MCHA(i))
     C		+B1zdec(mtau,MSNT,MCHA(i)))
	enddo

	SigmatauR=aux/(16.d0*Pi**2)

C		c) Correction to the vector/axial Zll couplings
	dgVtauSusy=dgVtauSusygauge+GamZllMZV*(2.d0*dsqrt(C2TW/g2))
     C	          -(gVZleplep+gAZleplep)*SigmatauL/2.d0
     C            -(gVZleplep-gAZleplep)*SigmatauR/2.d0
	dgAtauSUSY=dgAtauSUSYgauge+GamZllMZA*(2.d0*dsqrt(C2TW/g2))
     C            -(gVZleplep+gAZleplep)*SigmatauL/2.d0
     C            +(gVZleplep-gAZleplep)*SigmatauR/2.d0
	auxiVS=GamZllMZV*(2.d0*dsqrt(C2TW/g2))
     C	          -(gVZleplep+gAZleplep)*SigmatauL/2.d0
     C            -(gVZleplep-gAZleplep)*SigmatauR/2.d0
	auxiAS=GamZllMZA*(2.d0*dsqrt(C2TW/g2))
     C            -(gVZleplep+gAZleplep)*SigmatauL/2.d0
     C            +(gVZleplep-gAZleplep)*SigmatauR/2.d0

C	Theoretical Uncertainty of the Susy-Sector:
	deltadgVtauSusy=0.3d0*dabs(dgVtauSUSYgauge)
     C		   +0.1d0*dabs(GamZllMZV*(2.d0*dsqrt(C2TW/g2))
     C            -(gVZleplep+gAZleplep)*SigmatauL/2.d0
     C            -(gVZleplep-gAZleplep)*SigmatauR/2.d0)
	deltadgAtauSusy=0.3d0*dabs(dgAtauSUSYgauge)
     C		   +0.1d0*dabs(GamZllMZA*(2.d0*dsqrt(C2TW/g2))
     C            -(gVZleplep+gAZleplep)*SigmatauL/2.d0
     C            +(gVZleplep-gAZleplep)*SigmatauR/2.d0)

C*********************************************************************
C	multiplication with "2.d0*dsqrt(C2TW/g2)", because this 
C	prefactor in in the decay width formula inside GF, the 
C	selfenergies already independent because 
C	gVZleplep=(-1/2+2*S2TW)
C*********************************************************************
C
C*********************************************************************
C	3.) Calculation of effective axial and vector couplings 
C	    SUSY-contributions for Z-> e- e+
C*********************************************************************	
C	-> gauge part
	dgVESusygauge=(gVZleplep*(-1.d0/2.d0*SigmaZPrSusy(MZ,MW)
     C	    +1.d0/2.d0*SigmaGamPrSusy(0.d0,MW)
     C	    +(C2TW-S2TW)/dsqrt(S2TW*C2TW)*SigmaGamZSusy(0.d0,MW)
     C	    /MZ**2-1.d0/2.d0*(C2TW-S2TW)/S2TW
     C	    *(SigmaZSusy(MZ,MW)/MZ**2-SigmaWSusy(MW,MW)/MW**2))
     C	   +2.d0*C2TW*(SigmaZSusy(MZ,MW)/MZ**2-SigmaWSusy(MW,MW)
     C			/MW**2
     C	      -dsqrt(S2TW/C2TW)
     C	   *(SigmaGamZSusy(MZ,MW)+SigmaGamZSusy(0.d0,MW))/MZ**2)
     C	   -gAZleplep*dsqrt(C2TW/S2TW)*SigmaGamZSusy(0.d0,MW)
     C		/MZ**2)

	dgAESusygauge=(gAZleplep*(-1.d0/2.d0*SigmaZPrSusy(MZ,MW)
     C	    +1.d0/2.d0*SigmaGamPrSusy(0.d0,MW)
     C	    +(C2TW-S2TW)/dsqrt(S2TW*C2TW)*SigmaGamZSusy(0.d0,MW)
     C	         /MZ**2-1.d0/2.d0*(C2TW-S2TW)/S2TW
     C	    *(SigmaZSusy(MZ,MW)/MZ**2-SigmaWSusy(MW,MW)/MW**2))
     C	    -gAZleplep*dsqrt(C2TW/S2TW)*SigmaGamZSusy(0.d0,MW)
     C		/MZ**2)

C	-> "fermion part"

C		a) Contributions to the Zll vertex function

C	->Vector-Part

C	Chargino corrections:
       aux=0.d0
       do i=1,2
       do j=1,2
        aux=aux+MCHA(i)*MCHA(j)*(g2*V(i,1)*V(j,1)
     C          *(U(i,1)*U(j,1)+U(i,2)*U(j,2)/2.d0-S2TW*delt(i,j))
     C          +YukE**2*U(i,2)*U(j,2)
     C          *(V(i,1)*V(j,1)+V(i,2)*V(j,2)/2.d0-S2TW*delt(i,j)))
     C           *CC0(ME,MZ,MCHA(j),MNL,MCHA(i))
     C      -(g2*V(i,1)*V(j,1)
     C          *(V(i,1)*V(j,1)+V(i,2)*V(j,2)/2.d0-S2TW*delt(i,j))
     C          +YukE**2*U(i,2)*U(j,2)
     C          *(U(i,1)*U(j,1)+U(i,2)*U(j,2)/2.d0-S2TW*delt(i,j)))
     C  *(2.d0*C24zdec(ME,MZ,MCHA(j),MNL,MCHA(i))
     C		-1.d0/2.d0+MZ**2
     C      *(C12(ME,MZ,MCHA(j),MNL,MCHA(i))
     C                       +C23(ME,MZ,MCHA(j),MNL,MCHA(i))))

       enddo
        aux=aux+(g2*V(i,1)**2+YukE**2*U(i,2)**2)
     C                 *C24zdec(ME,MZ,MNL,MCHA(i),MNL)
       enddo

C	Neutralino corrections:
       do i=1,5
       do j=1,5
       do k=1,2
        aux=aux-(NEU(i,4)*NEU(j,4)-NEU(i,3)*NEU(j,3))/2.d0
     C   *((NEU(i,1)*delt(k,1)*dsqrt(g1/2.d0)
     C      +NEU(i,2)*delt(k,1)*dsqrt(g2/2.d0)-YukE*NEU(i,4)
     C	      *delt(k,2))
     C     *(NEU(j,1)*delt(k,1)*dsqrt(g1/2.d0)
     C      +NEU(j,2)*delt(k,1)*dsqrt(g2/2.d0)-YukE*NEU(j,4)
     C	      *delt(k,2))
     C      -(NEU(i,1)*delt(k,1)*dsqrt(g1*2.d0)
     C                         +YukE*NEU(i,4)*delt(k,1))
     C      *(NEU(j,1)*delt(k,1)*dsqrt(g1*2.d0)
     C                         +YukE*NEU(j,4)*delt(k,1)))
     C   *(MNEU(i)*MNEU(j)*CC0(ME,MZ,MNEU(j),MSE(k),MNEU(i))
     C     +2.d0*C24zdec(ME,MZ,MNEU(j),MSE(k),MNEU(i))
     C		-1.d0/2.d0
     C     +MZ**2*(C12(ME,MZ,MNEU(j),MSE(k),MNEU(i))
     C             +C23(ME,MZ,MNEU(j),MSE(k),MNEU(i))))
	enddo
	enddo
	enddo

       do i=1,2
       do j=1,2
       do k=1,5
	aux=aux-2.d0*((1.d0/2.d0-S2TW)*delt(i,1)*delt(j,1)
     C                                 -S2TW*delt(i,2)*delt(j,2))
     C    *((NEU(k,1)*delt(i,1)*dsqrt(g1/2.d0)
     C      +NEU(k,2)*delt(i,1)*dsqrt(g2/2.d0)-YukE*NEU(k,4)
     C	*delt(i,2))
     C     *(NEU(k,1)*delt(j,1)*dsqrt(g1/2.d0)
     C      +NEU(k,2)*delt(j,1)*dsqrt(g2/2.d0)-YukE*NEU(k,4)
     C	*delt(j,2))
     C      +(NEU(k,1)*delt(i,2)*dsqrt(g1*2.d0)
     C                         +YukE*NEU(k,4)*delt(i,1))
     C      *(NEU(k,1)*delt(j,2)*dsqrt(g1*2.d0)
     C                         +YukE*NEU(k,4)*delt(j,1)))
     C        *C24zdec(ME,MZ,MSE(j),MNEU(k),MSE(i))
       enddo
       enddo
       enddo

       GamZllMZV=dsqrt(g2)/(32.d0*Pi**2*dsqrt(1.d0-S2TW))*aux

C	->Axial-Vector-Part
C	Chargino corrections:
	aux=0.d0
	do i=1,2
	do j=1,2
        aux=aux+MCHA(i)*MCHA(j)*(g2*V(i,1)*V(j,1)
     C          *(U(i,1)*U(j,1)+U(i,2)*U(j,2)/2.d0-S2TW*delt(i,j))
     C          -YukE**2*U(i,2)*U(j,2)
     C          *(V(i,1)*V(j,1)+V(i,2)*V(j,2)/2.d0-S2TW*delt(i,j)))
     C           *CC0(ME,MZ,MCHA(j),MNL,MCHA(i))
     C      +(g2*V(i,1)*V(j,1)
     C          *(V(i,1)*V(j,1)+V(i,2)*V(j,2)/2.d0-S2TW*delt(i,j))
     C          -YukE**2*U(i,2)*U(j,2)
     C          *(U(i,1)*U(j,1)+U(i,2)*U(j,2)/2.d0-S2TW*delt(i,j)))
     C  *(-2.d0*C24zdec(ME,MZ,MCHA(j),MNL,MCHA(i))+1.d0/2.d0-MZ**2
     C      *(C12(ME,MZ,MCHA(j),MNL,MCHA(i))
     C                      +C23(ME,MZ,MCHA(j),MNL,MCHA(i))))
	enddo
        aux=aux+(g2*V(i,1)**2-YukE**2*U(i,2)**2)
     C                     *C24zdec(ME,MZ,MNL,MCHA(i),MNL)
	enddo

C	Neutralino corrections:
	do i=1,5
	do j=1,5
	do k=1,2
	aux=aux-(NEU(i,4)*NEU(j,4)-NEU(i,3)*NEU(j,3))/2.d0
     C    *((NEU(i,1)*delt(k,1)*dsqrt(g1/2.d0)
     C      +NEU(i,2)*delt(k,1)*dsqrt(g2/2.d0)-YukE*NEU(i,4)
     C	*delt(k,2))
     C     *(NEU(j,1)*delt(k,1)*dsqrt(g1/2.d0)
     C      +NEU(j,2)*delt(k,1)*dsqrt(g2/2.d0)-YukE*NEU(j,4)
     C	*delt(k,2))
     C      +(NEU(i,1)*delt(k,2)*dsqrt(g1*2.d0)
     C                         +YukE*NEU(i,4)*delt(k,1))
     C      *(NEU(j,1)*delt(k,2)*dsqrt(g1*2.d0)
     C                         +YukE*NEU(j,4)*delt(k,1)))
     C   *(MNEU(i)*MNEU(j)*CC0(ME,MZ,MNEU(j),MSE(k),MNEU(i))
     C     +2.d0*C24zdec(ME,MZ,MNEU(j),MSE(k),MNEU(i))
     C		-1.d0/2.d0
     C     +MZ**2*(C12(ME,MZ,MNEU(j),MSE(k),MNEU(i))
     C             +C23(ME,MZ,MNEU(j),MSE(k),MNEU(i))))
	enddo	
	enddo
	enddo

	do i=1,2
	do j=1,2
	do k=1,5
	aux=aux-2.d0*((1.d0/2.d0-S2TW)*delt(i,1)*delt(j,1)
     C                                 -S2TW*delt(i,2)*delt(j,2))
     C    *((NEU(k,1)*delt(i,1)*dsqrt(g1/2.d0)
     C      +NEU(k,2)*delt(i,1)*dsqrt(g2/2.d0)-YukE*NEU(k,4)
     C	*delt(i,2))
     C     *(NEU(k,1)*delt(j,1)*dsqrt(g1/2.d0)
     C      +NEU(k,2)*delt(j,1)*dsqrt(g2/2.d0)-YukE*NEU(k,4)
     C	*delt(j,2))
     C      -(NEU(k,1)*delt(i,2)*dsqrt(g1*2.d0)
     C                         +YukE*NEU(k,4)*delt(i,1))
     C      *(NEU(k,1)*delt(j,2)*dsqrt(g1*2.d0)
     C                         +YukE*NEU(k,4)*delt(j,1)))
     C        *C24zdec(ME,MZ,MSE(j),MNEU(k),MSE(i))
	enddo
	enddo
	enddo	

	GamZllMZA=dsqrt(g2)/(32.d0*Pi**2*dsqrt(1.d0-S2TW))*aux

C*********************************************************************

C		b) Contributions to the leptonic self-energies
C
C	SigmaEL  
	aux=0.d0
	do i=1,5
	do j=1,2
        aux=aux+(NEU(i,1)*delt(j,1)*dsqrt(g1/2.d0)
     C      +NEU(i,2)*delt(j,1)*dsqrt(g2/2.d0)
     C	    -YukE*NEU(i,4)*delt(j,2))**2
     C	       *(B0zdec(ME,MSE(j),MNEU(i))
     C		+B1zdec(ME,MSE(j),MNEU(i)))
	enddo
	enddo
	do i=1,2
        aux=aux+g2*(V(i,1))**2
     C       *(B0zdec(ME,MNL,MCHA(i))
     C		+B1zdec(ME,MNL,MCHA(i)))
	enddo
	SigmaEL=aux/(16.d0*Pi**2)

C	SigmaER
	aux=0.d0
	do i=1,5
	do j=1,2
	aux=aux+(-NEU(i,1)*delt(j,2)*dsqrt(g1*2.d0)
     C                         -YukE*NEU(i,4)*delt(j,1))**2
     C  *(B0zdec(ME,MSE(j),MNEU(i))
     C		+B1zdec(ME,MSE(j),MNEU(i)))
	enddo
	enddo
	do i=1,2
	aux=aux+YukE**2*(U(i,2))**2
     C  *(B0zdec(ME,MNL,MCHA(i))
     C		+B1zdec(ME,MNL,MCHA(i)))
	enddo

	SigmaER=aux/(16.d0*Pi**2)

C		c) Correction to the vector/axial Zll couplings
	dgVESusy=dgVESusygauge+GamZllMZV*(2.d0*dsqrt(C2TW/g2))
     C	          -(gVZleplep+gAZleplep)*SigmaEL/2.d0
     C            -(gVZleplep-gAZleplep)*SigmaER/2.d0
	dgAESUSY=dgAESUSYgauge+GamZllMZA*(2.d0*dsqrt(C2TW/g2))
     C            -(gVZleplep+gAZleplep)*SigmaEL/2.d0
     C            +(gVZleplep-gAZleplep)*SigmaER/2.d0
	auxiVS=0.1d0*dabs(auxiVS-GamZllMZV*(2.d0*dsqrt(C2TW/g2))
     C	          +(gVZleplep+gAZleplep)*SigmaEL/2.d0
     C            +(gVZleplep-gAZleplep)*SigmaER/2.d0)
	auxiAS=0.1d0*dabs(auxiAS-GamZllMZA*(2.d0*dsqrt(C2TW/g2))
     C            +(gVZleplep+gAZleplep)*SigmaEL/2.d0
     C            -(gVZleplep-gAZleplep)*SigmaER/2.d0)

C	Theoretical Uncertainty of the Susy-Sector:
	deltadgVESusy=0.3d0*dabs(dgVESUSYgauge)
     C		   +0.1d0*dabs(GamZllMZV*(2.d0*dsqrt(C2TW/g2))
     C            -(gVZleplep+gAZleplep)*SigmaEL/2.d0
     C            -(gVZleplep-gAZleplep)*SigmaER/2.d0)
	deltadgAESusy=0.3d0*dabs(dgAESUSYgauge)
     C		   +0.1d0*dabs(
     C		   GamZllMZA*(2.d0*dsqrt(C2TW/g2))
     C            -(gVZleplep+gAZleplep)*SigmaEL/2.d0
     C            +(gVZleplep-gAZleplep)*SigmaER/2.d0)

C*********************************************************************
C	multiplication with "2.d0*dsqrt(C2TW/g2)", because this 
C	prefactor in in the decay width formula inside GF, the 
C	selfenergies already independent because 
C	gVZleplep=(-1/2+2*S2TW)
C*********************************************************************
C*********************************************************************
C*********************************************************************
C	Now I have: SM-Couplings: gVtauSM, gAtauSM
C	NP-Contributions, devided in Higgs-,und Susy part: 
C	dgVtauHiggs,dgVtauSusy;dgAtauHiggs,dgAtauSusy
C
C	Contributions to deltar
C	SM: drSM
C	Susy: drSUSY
C	Higgs: drHiggs
C
C	everything together will then be: 
C	gVefftauNMSSM=gVtauSM+dgVtauHiggs+dgVtauSusy;
C	gAefftauNMSSM=gAtauSM+dgAtauHiggs+dgAtauSusy
C	drNMSSM=drSM+drSUSY+drHiggs
C*********************************************************************
C*********************************************************************
C	IV. a) Parameters of Z -> tau- tau+ and Z -> e- e+
C*********************************************************************

C*********************************************************************
C	Effective Couplings for Z->Tau Tau:

	gVefftauNMSSM=gVeffSM
     C		+dgVtauHiggs
     C		+dgVtauSusy
	gAefftauNMSSM=gAeffSM
     C		+dgAtauHiggs
     C		+dgAtauSusy

C	Theoretical Uncertainty on Effective Couplings*sqrt(1-Dr)
	deltagVefftauNMSSM=deltagVeffSM
     C +(deltadgVtauHiggs+deltadgVtauSusy)*dsqrt(1.d0-DrNMSSM)
     C +dabs(dgVtauHiggs+dgVtauSusy)/dsqrt(1.d0-DrNMSSM)
     C        *dsqrt(2.d0)*Gmu*MW*dMW/Pi/alpha
	deltagAefftauNMSSM=deltagAeffSM
     C +(deltadgAtauHiggs+deltadgAtauSusy)*dsqrt(1.d0-DrNMSSM)
     C +dabs(dgAtauHiggs+dgAtauSusy)/dsqrt(1.d0-DrNMSSM)
     C        *dsqrt(2.d0)*Gmu*MW*dMW/Pi/alpha

C*********************************************************************
C	Effective Couplings for Z->e e:

	gVeffENMSSM=gVeffSM
     C		+dgVEHiggs
     C		+dgVESusy
	gAeffENMSSM=gAeffSM
     C		+dgAEHiggs
     C		+dgAESusy

C	Theoretical Uncertainty on Effective Couplings*sqrt(1-Dr)
	deltagVeffENMSSM=deltagVeffSM
     C     +(deltadgVEHiggs+deltadgVESusy)*dsqrt(1.d0-DrNMSSM)
     C +dabs(dgVEHiggs+dgVESusy)/dsqrt(1.d0-DrNMSSM)*dsqrt(2.d0)*Gmu
     C        *MW*dMW/Pi/alpha
	deltagAeffENMSSM=deltagAeffSM
     C   +(deltadgAEHiggs+deltadgAESusy)*dsqrt(1.d0-DrNMSSM)
     C +dabs(dgAEHiggs+dgAESusy)/dsqrt(1.d0-DrNMSSM)*dsqrt(2.d0)*Gmu
     C        *MW*dMW/Pi/alpha

C*********************************************************************
C	 IV. b) Observables of the Z -> tau- tau+ and Z -> e- e+	
C*********************************************************************
C	Decay Width Z->Tau Tau

	decztt=1.d0/(12.d0*pi)*MZ**3*dsqrt(2.d0)*Gmu*(1.d0-drNMSSM)
     C		*(gVefftauNMSSM**2*(dsqrt(1.d0-4.d0*MTAU**2/MZ**2)
     C		*(1.d0+2.d0*MTAU**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi)
     C		+gAefftauNMSSM**2*(dsqrt(1.d0-4.d0*MTAU**2/MZ**2)
     C		*(1.d0-4.d0*MTAU**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi))

C	theoretical error of Decay Width Z->Tau Tau 
	
	deltadecztt=0.d0
Cdabs(1.d0/(12.d0*pi)*MZ**3*dsqrt(2.d0)*Gmu
C     C		*(gVefftauNMSSM**2*(dsqrt(1.d0-4.d0*MTAU**2/MZ**2)
C     C		*(1.d0+2.d0*MTAU**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi)
C     C		+gAefftauNMSSM**2*(dsqrt(1.d0-4.d0*MTAU**2/MZ**2)
C     C		*(1.d0-4.d0*MTAU**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi)))
C     C		*dsqrt(2.d0)*Gmu*MW*dMW/Pi/alpha
     C        +dabs(1.d0/(12.d0*pi)*MZ**3*dsqrt(2.d0)
     C		*Gmu*(1.d0-DrNMSSM)*2.d0*gVefftauNMSSM
     C		*(dsqrt(1.d0-4.d0*MTAU**2/MZ**2)
     C		*(1.d0+2.d0*MTAU**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi))*
     C		deltagVefftauNMSSM
     C        +dabs(1.d0/(12.d0*pi)*MZ**3
     C		*dsqrt(2.d0)*Gmu*(1.d0-drNMSSM)*2.d0*gAefftauNMSSM
     C		*(dsqrt(1.d0-4.d0*MTAU**2/MZ**2)
     C		*(1.d0-4.d0*MTAU**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi))*
     C		deltagAefftauNMSSM
	deltadecztt=dsqrt(deltadecztt**2
     C +4.d0*decztt**2*((3.d0*dMZ/MZ)**2+(1.d-10/Gmu)**2))

C		->BR(Z->Tau+Tau-)
	BRZTauTau=decztt/TOTDWZ
	BRZTTmin=(decztt-deltadecztt)/TOTDWZmax
	BRZTTmax=(decztt+deltadecztt)/TOTDWZmin

C*********************************************************************
C	Decay Width Z->e e

	deczee=1.d0/(12.d0*pi)*MZ**3*dsqrt(2.d0)*Gmu*(1.d0-drNMSSM)
     C		*(gVeffENMSSM**2*(dsqrt(1.d0-4.d0*ME**2/MZ**2)
     C		*(1.d0+2.d0*ME**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi)
     C		+gAeffENMSSM**2*(dsqrt(1.d0-4.d0*ME**2/MZ**2)
     C		*(1.d0-4.d0*ME**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi))

C	theoretical error of Decay Width Z->e e 
	
	deltadeczee=0.d0
Cdabs(1.d0/(12.d0*pi)*MZ**3*dsqrt(2.d0)*Gmu
C     C		*(gVeffENMSSM**2*(dsqrt(1.d0-4.d0*ME**2/MZ**2)
C     C		*(1.d0+2.d0*ME**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi)
C     C		+gAeffENMSSM**2*(dsqrt(1.d0-4.d0*ME**2/MZ**2)
C     C		*(1.d0-4.d0*ME**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi)))
C     C		*dsqrt(2.d0)*Gmu*MW*dMW/Pi/alpha
     C		+dabs(1.d0/(12.d0*pi)*MZ**3*dsqrt(2.d0)
     C		*Gmu*(1.d0-drNMSSM)*2.d0*gVeffENMSSM
     C		*(dsqrt(1.d0-4.d0*ME**2/MZ**2)
     C		*(1.d0+2.d0*ME**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi))*
     C		deltagVeffENMSSM+dabs(1.d0/(12.d0*pi)*(MZ+dMZ)**3
     C		*dsqrt(2.d0)*(Gmu+1.d-10)*(1.d0-drNMSSM)*2.d0
     C		*gAeffENMSSM*(dsqrt(1.d0-4.d0*ME**2/MZ**2)
     C		*(1.d0-4.d0*ME**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi))*
     C		deltagAeffENMSSM
	deltadeczee=dsqrt(deltadeczee**2
     C +4.d0*deczee**2*((3.d0*dMZ/MZ)**2+(1.d-10/Gmu)**2))

C*********************************************************************
C		->effective Tau mixing angle

	S2TWeffTau=1.d0/4.d0*(1.d0-gVefftauNMSSM/gAefftauNMSSM)

C	Theoretical Uncertainty on S2TWeffTau
	deltaS2TWeffTau=
     C		dabs(1.d0/(4.d0*gAefftauNMSSM))
     C		*deltagVefftauNMSSM/dsqrt(1.d0-DrNMSSM)
     C		+dabs(1.d0/4.d0
     C		*gVefftauNMSSM/gAefftauNMSSM**2)*deltagAefftauNMSSM
     C          /dsqrt(1.d0-DrNMSSM)
	deltaS2TWeffTau=(4.d0*deltaS2TWeffsm+(deltadgVtauHiggs+
     C deltadgVtaususy)/dabs(gAeffSM)/dsqrt(1-DrNMSSM)
     C +deltagAeffSM*dabs(dgVtauHiggs+dgVtauSUSY)/gAeffSM**2
     C /dsqrt(1-DrNMSSM))/dabs(1.d0+(dgAtauHiggs+dgAtauSUSY)
     c /gAeffsm)
	deltaS2TWeffTau=deltaS2TWeffTau+dabs(1.d0-4.d0*
     c s2tweffsm+(dgVtauHiggs+dgVtauSUSY)/gAeffSM)/(1.d0+
     c (dgAtauHiggs+dgAtauSUSY)/gAeffsm)**2*(
     c (deltadgAtauHiggs+deltadgAtaususy)/dabs(gAeffSM)
     c /dsqrt(1-DrNMSSM)+deltagAeffSM*dabs(dgAtauHiggs
     c +dgAtauSUSY)/dsqrt(1-DrNMSSM)/gAeffSM**2)
	deltaS2TWeffTau=deltaS2TWeffTau/4.d0

C*********************************************************************
C	Ratio: (Z->Tau Tau)/(Z-> e e)-1

	ratio=decztt/deczee-1.d0

C	Theoretical Uncertainty on this Ratio
	deltaratio=dabs(1.d0/((1.d0+3.d0/4.d0*ALEMMZ/pi)
     C		*(gVeffENMSSM**2+gAeffENMSSM**2)))*(dabs(2.d0
     C		*gVefftauNMSSM*(dsqrt(1.d0-4.d0*MTAU**2/MZ**2)
     C		*(1.d0+2.d0*MTAU**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi))
     C		*deltagVefftauNMSSM
     C		+dabs(2.d0*gAefftauNMSSM*(dsqrt(1.d0-4.d0*MTAU**2
     C		/MZ**2)*(1.d0-4.d0*MTAU**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi))
     C		*deltagAefftauNMSSM
     C		+dabs(2.d0*gVeffENMSSM*((gVefftauNMSSM**2
     C		*(dsqrt(1.d0-4.d0*MTAU**2/MZ**2)
     C		*(1.d0+2.d0*MTAU**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi)
     C		+gAefftauNMSSM**2*(dsqrt(1.d0-4.d0*MTAU**2/MZ**2)
     C		*(1.d0-4.d0*MTAU**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi)))
     C		/((1.d0+3.d0/4.d0*ALEMMZ/pi)*(gVeffENMSSM**2
     C		+gAeffENMSSM**2)))*deltagVeffENMSSM
     C		+dabs(2.d0*gAeffENMSSM*((gVefftauNMSSM**2
     C		*(dsqrt(1.d0-4.d0*MTAU**2/MZ**2)
     C		*(1.d0+2.d0*MTAU**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi)
     C		+gAefftauNMSSM**2*(dsqrt(1.d0-4.d0*MTAU**2/MZ**2)
     C		*(1.d0-4.d0*MTAU**2/MZ**2)+3.d0/4.d0*ALEMMZ/pi)))
     C		/((1.d0+3.d0/4.d0*ALEMMZ/pi)
     C		*(gVeffENMSSM**2+gAeffENMSSM**2)))*deltagAeffENMSSM
     C		)

	deltaratio=2.d0*(2.d0*(dabs(gVeffENMSSM)+deltagVeffENMSSM)
     C   *(auxiVH+auxiVS+0.001d0*dabs(gVeffENMSSM))
     C                  +2.d0*(dabs(gAeffENMSSM)+deltagAeffENMSSM)
     C   *(auxiAH+auxiAS+0.001d0*dabs(gAeffENMSSM))
     C                  +12.d0*MTAU**2/(MZ**2*(1+3.d0*Pi*ALEMMZ/4.d0))
     C         *(dabs(gAeffENMSSM)+deltagAeffENMSSM)*deltagAeffENMSSM)
     C              /((dabs(gVeffENMSSM)-deltagVeffENMSSM)**2
     C                +(dabs(gAeffENMSSM)-deltagAeffENMSSM)**2)

c	PROB(32)=0.d0
c	IF(BRZTTmin.GE.BRZTTexpmax)
c     C	   PROB(32)=(BRZTTmin-BRZTTexpmax)*1.D5
c	IF(BRZTTmax.LE.BRZTTexpmin)
c     C	   PROB(32)=(BRZTTmax-BRZTTexpmin)*1.D5

C	In case I want this not as a constraint
c	PROB(32)=0.d0
c      print*,MW,dMW,DrNMSSM,MWSM,BRZTauTau,BRZTTmin,BRZTTmax,
c     C   deczee,deltadeczee,ratio,deltaratio,S2TWeffTau,
c     C   deltaS2TWeffTau, deltadeczee/deczee
	dMW=dMW/2.d0
	If(SCOMP(1,3)**2.le.0.5d0)then
	mHSM=SMASS(1)
	else
	mHSM=SMASS(2)
	endif
	MWSM=80.3799d0-0.05429d0*dlog(mHSM/100.d0)
     C           -8.939d-3*dlog(mHSM/100.d0)**2
     C           +8.90d-5*dlog(mHSM/100.d0)**4
     C           +1.61d-4*((mHSM/100.d0)**2-1.d0)
     C           -1.070d0*(Dalph/0.05907d0-1.d0)
     C           +.5256d0*((mt/174.3)**2-1.d0)
     C           -.0678d0*((mt/174.3)**2-1.d0)**2
     C  -1.79d-3*dlog(mHSM/100.d0)*((mt/174.3)**2-1.d0)
     C  +6.59d-5*(mHSM/100.d0)**2*((mt/174.3)**2-1.d0)
     C           -7.37d-2*(alSMZ/0.119d0-1.d0)
     C           +114.9d0*(MZ/91.1875d0-1.d0)
	S2TWeffsm=0.2312527d0+4.729d-4*dlog(mHSM/100.d0)
     C +2.07d-5*dlog(mHSM/100.d0)**2+3.85d-6*dlog(mHSM/100.d0)**4
     C -1.85d-6*((mHSM/100.d0)**2-1.d0)
     C +2.07d-2*(Dalph/0.05907d0-1.d0)
     C -2.851d-3*((mt/178.d0)**2-1.d0)
     C +1.82d-4*((mt/178.d0)**2-1.d0)**2
     C -9.74d-6*((mt/178.d0)**2-1.d0)*(mHSM/100.d0-1.d0)
     C +3.98d-4*(alSMZ/0.117d0-1.d0)
     C +6.55d-1*(MZ/91.1876d0-1.d0)
	return
	end

	
C*********************************************************************
C*********************************************************************
C	one-,two-,three-point-functions (relations between the  
C	several functions + definition of B5:[3],conventions:[4])
C*********************************************************************
C*********************************************************************
C		->A0(m)		
	DOUBLE PRECISION function A0zdec(a)
	
	implicit none
	DOUBLE PRECISION a,aux,div

	div=0.d0
	   
	if(a.eq.0.d0) then
	aux=0.d0
	else
	aux=a**2*(-dlog(a**2)+1.d0)
	endif
	A0zdec=aux+a**2*div

	return
	end

C********************************************************************
C		->B0(p,m1,m2)
	DOUBLE PRECISION function B0zdec(x,a,b)

	implicit none
	DOUBLE PRECISION x,a,b,ME
	DOUBLE PRECISION delta,sigma,beta,L,F0,div
	div=0.d0

	ME=0.5109989d-3
	
	if(x.eq.0.d0.or.x.eq.ME) then
	
	If(dabs(a**2-b**2).le.a**2*7.d-6)then

		if(a.eq.0.d0)then
		B0zdec=div+0.d0
		else
		B0zdec=div-2.d0*dlog(dabs(a))
		endif

	elseif((a.eq.0.d0.or.a.eq.ME).and.b.ne.0.d0)then
	B0zdec=div+1.d0-2.d0*dlog(dabs(b))
	elseif((b.eq.0.d0.or.b.eq.ME).and.a.ne.0.d0)then
	B0zdec=div+1.d0-2.d0*dlog(dabs(a))
	else
	B0zdec=div+1.d0/(b**2-a**2)*(a**2*dlog(a**2)-b**2
     C		*dlog(b**2))+1.d0
	endif
	
	else
	if(a.eq.0.d0)then
	  if(b.eq.0.d0)then
	  B0zdec=div-2.d0*dlog(dabs(x))+2.d0
	  else
	  B0zdec=div-1.d0/x**2*(2.d0*dlog(dabs(b))*b**2-2.d0*x**2
     C		+(x**2-b**2)*dlog(dabs(b**2-x**2)))
          endif
	
	elseif(b.eq.0.d0)then
	  if(a.eq.0.d0)then
	  B0zdec=div-2.d0*dlog(dabs(x))+2.d0
	  else
	  B0zdec=div-1.d0/x**2*(2.d0*dlog(dabs(a))*a**2-2.d0*x**2
     C		+(x**2-a**2)*dlog(dabs(a**2-x**2)))
	  endif	
	else

C		definition of sigma and delta	
	sigma=(a**2+b**2)/x**2
	delta=(a**2-b**2)/x**2
	If(dabs(a**2-b**2).le.a**2*7.d-6)Then
	delta=0.d0
c	sigma=2.d0*a**2/x**2
	ENDIF

C		various cases for beta
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.(x**2)
     C	         .lt.((dabs(a)-dabs(b))**2)) then
	beta=dsqrt(1.d0-2.d0*sigma+delta**2)
	else
	beta=dsqrt(2.d0*sigma-delta**2-1.d0)
	endif	

C		various cases for L 
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.(x**2)
     C	         .lt.((dabs(a)-dabs(b))**2)) then
	L=1.d0/2.d0*dlog(dabs((1.d0+beta-sigma)/(1.d0-beta-sigma)))
	else
	L=datan((1.d0-delta)/dabs(beta))
     C	+datan((1.d0+delta)/dabs(beta))	
	endif

	F0=dlog(dabs(a*b))-delta*dlog(dabs(b/a))-2.d0+beta*L
	
	B0zdec=div-F0
	endif
	endif
	
	return
	end

C********************************************************************
C		B1(p,m1,m2)	
	DOUBLE PRECISION function B1zdec(x,a,b)

	implicit none
	DOUBLE PRECISION x,a,b,ME
	DOUBLE PRECISION sigma,delta,beta,L,div
	DOUBLE PRECISION F0,FA
	div=0.d0

	ME=0.5109989d-3

	if(x.eq.0.d0.or.x.eq.ME)then

	If(dabs(a**2-b**2).le.a**2*5.d-6)then
	          if(a.eq.0.d0)then
		  B1zdec=-1.d0/2.d0*div+0.d0
		  else
	  	  B1zdec=-1.d0/2.d0*div+dlog(dabs(a))
		  endif
	elseif((a.eq.0.d0.or.a.eq.ME).and.b.ne.0.d0)then
	B1zdec=-1.d0/2.d0*div+dlog(dabs(b))-1.d0/4.d0
	elseif((b.eq.0.d0.or.b.eq.ME).and.a.ne.0.d0)then
	B1zdec=-1.d0/2.d0*div+dlog(dabs(a))-3.d0/4.d0	
	else
	B1zdec=-1.d0/2.d0*div+(-3.d0*a**4+4.d0*b**2*a**2-b**4+4.d0
     C	*a**4*dlog(dabs(a))+4.d0*b**2*(b**2-2.d0*a**2)
     C	*dlog(dabs(b)))/(4.d0*(a**2-b**2)**2)
	endif
	
	else	
	if(a.eq.0.d0)then
	
	  if(b.eq.0.d0)then
	  B1zdec=-1.d0/2.d0*div+1.d0*dlog(dabs(x))-1.d0
	  else
	  B1zdec=-1.d0/2.d0*div+1.d0/(2.d0*x**4)*((b**2-x**2)**2
     C		*dlog(dabs(b**2-x**2))-(b**2-2.d0*x**2)*(2.d0*b**2
     C		*dlog(dabs(b))-x**2))
	  endif

	elseif(b.eq.0.d0)then
	  if(a.eq.0.d0)then
	  B1zdec=-1.d0/2.d0*div+1.d0*dlog(x)-1.d0
	  else
	  B1zdec=-1.d0/2.d0*div+1.d0/(2.d0*x**4)*(2.d0
     C		*dlog(dabs(a))*a**4-x**2*(a**2+2.d0*x**2)
     C		+(x**4-a**4)*dlog(dabs(a**2-x**2)))
	  endif

	elseif(a.eq.x.and.b.ge.100.d0)then
	B1zdec=-1.d0/2.d0*div+(3.d0*a**12)/b**12+(2.d0*a**10)/b**10
     C	       -(21.d0*a**8)/(4.d0*b**8)-(5.d0*a**6)/b**6
     C	       -(4.d0*dlog(a)*a**4)/b**4+(4.d0*dlog(b)*a**4)/b**4
     C	       -(17.d0*a**4)/(8.d0*b**4)+(9.d0*a**2)/(4.d0*b**2)
     C		+dlog(b)-1.d0/4.d0
	else

C		definition of sigma and delta	
	sigma=(a**2+b**2)/x**2
	delta=(a**2-b**2)/x**2
	If(dabs(a**2-b**2).le.a**2*5.d-6)Then
	delta=0.d0
c	sigma=2.d0*a**2/x**2
	ENDIF

C		various cases for beta
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.(x**2)
     C		 .lt.((dabs(a)-dabs(b))**2)) then
	beta=dsqrt(1.d0-2.d0*sigma+delta**2)
	else
	beta=dsqrt(2.d0*sigma-delta**2-1.d0)
	endif

C		various cases for L 
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.(x**2)
     C	         .lt.((dabs(a)-dabs(b))**2)) then
	L=1.d0/2.d0*dlog(dabs((1.d0+beta-sigma)/(1.d0-beta-sigma)))
	else
	L=datan((1.d0-delta)/dabs(beta))
     C	+datan((1.d0+delta)/dabs(beta))	
	endif	

	F0=dlog(dabs(a*b))-delta*dlog(dabs(b/a))-2.d0+beta*L
	FA=-(sigma-delta**2)*dlog(dabs(b/a))+delta*(1.d0-beta*L)

	B1zdec=-1.d0/2.d0*div+1.d0/2.d0*(F0-FA)
	endif
	endif

	return
	end

C********************************************************************
C		->B5(p,m1,m2)
	DOUBLE PRECISION function B5(x,a,b)

	implicit none
	DOUBLE PRECISION x,a,b,ME
	DOUBLE PRECISION sigma,delta,beta,L,div
	DOUBLE PRECISION F0,FA,F3
	div=0.d0

	ME=0.5109989d-3

	if(x.eq.0.d0.or.x.eq.ME)then

	if(dabs(a**2-b**2).le.5.d-6*a**2)then
	B5=0.d0+x**2/3.d0*div
	elseif((a.eq.0.d0.or.a.eq.ME).and.b.ne.0.d0)then
	B5=-b**2/2.d0+x**2/3.d0*div
	elseif((b.eq.0.d0.or.b.eq.ME).and.a.ne.0.d0)then
	B5=-a**2/2.d0+x**2/3.d0*div
	else
	B5=-(a**4+4.d0*a**2*b**2*dlog(dabs(b/a))-b**4)
     C		/(2.d0*(a**2-b**2))+x**2/3.d0*div
	endif
	else	
C		definition of sigma and delta	
	sigma=(a**2+b**2)/x**2
	delta=(a**2-b**2)/x**2
	if(dabs(a**2-b**2).le.5.d-6*a**2)then
	delta=0.d0
	endif

C		various cases for beta
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.(x**2)
     C	         .lt.((dabs(a)-dabs(b))**2)) then
	beta=dsqrt(1.d0-2.d0*sigma+delta**2)
	else
	beta=dsqrt(2.d0*sigma-delta**2-1.d0)
	endif

C		various cases for L 
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.(x**2)
     C		 .lt.((dabs(a)-dabs(b))**2)) then
	L=1.d0/2.d0*dlog(dabs((1.d0+beta-sigma)
     C		/(1.d0-beta-sigma)))
	else
	L=datan((1.d0-delta)/dabs(beta))
     C	+datan((1.d0+delta)/dabs(beta))	
	endif

	F0=dlog(abs(a*b))-delta*dlog(dabs(b/a))-2.d0+beta*L
	FA=-(sigma-delta**2)*dlog(dabs(b/a))+delta*(1.d0-beta*L)
	F3=1.d0/6.d0*dlog(dabs(a*b))-(3.d0*sigma-2.d0*delta**2)
     C	*delta*dlog(dabs(b/a))/6.d0-5.d0/18.d0-(sigma-delta**2)
     C	/3.d0+(1.d0+sigma-2.d0*delta**2)*beta*L/6.d0

	B5=x**2/3.d0*div-(x**2*(F0-4.d0*F3)+(a**2-b**2)
     C	*FA)
	endif

	return
	end

C*******************************************************************
C		->B21(p,m1,m2)
	DOUBLE PRECISION function B21(x,a,b)

	implicit none
	DOUBLE PRECISION x,a,b,aux,ME
	DOUBLE PRECISION sigma,delta,beta,L,div
	DOUBLE PRECISION F0,FA,F3
	div=0.d0

	ME=0.5109989d-3

	if(x.eq.0.d0.or.x.eq.ME)then

	if(dabs(a**2-b**2).le.5.d-6*a**2)then
	aux=-(2.d0*dlog(dabs(a)))/3.d0
	elseif((a.eq.0.d0.or.a.eq.ME).and.b.ne.0.d0)then
	aux=1.d0/9.d0*(1.d0-6.d0*dlog(dabs(b)))
	elseif((b.eq.0.d0.or.b.eq.ME).and.a.ne.0.d0)then
	aux=11.d0/18.d0-(2.d0*dlog(dabs(a)))/3.d0
	else
	aux=-(12.d0*dlog(dabs(a))*a**6-11.d0*a**6+18.d0*b**2*a**4
     C      -9.d0*b**4*a**2+2.d0*b**6-12.d0*(b**6-3.d0*a**2*b**4
     C      +3.d0*a**4*b**2)*dlog(dabs(b)))/(18.d0*(a**2-b**2)**3)
	endif

	else	
C		definition of sigma and delta	
	sigma=(a**2+b**2)/x**2
	delta=(a**2-b**2)/x**2
	if(dabs(a**2-b**2).le.5.d-6*a**2)then
	delta=0.d0
	endif
C		various cases for beta
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.
     C 	   (x**2).lt.((dabs(a)-dabs(b))**2)) then
	beta=dsqrt(1.d0-2.d0*sigma+delta**2)
	else
	beta=dsqrt(2.d0*sigma-delta**2-1.d0)
	endif

C		various cases for L 
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.(x**2)
     C	         .lt.((dabs(a)-dabs(b))**2)) then
	L=1.d0/2.d0*dlog(dabs((1.d0+beta-sigma)
     C  		/(1.d0-beta-sigma)))
	else
	L=datan((1.d0-delta)/dabs(beta))
     C	+datan((1.d0+delta)/dabs(beta))	
	endif

	F0=dlog(dabs(a*b))-delta*dlog(dabs(b/a))-2.d0+beta*L
	FA=-(sigma-delta**2)*dlog(dabs(b/a))+delta*(1.d0-beta*L)
	F3=1.d0/6.d0*dlog(dabs(a*b))-(3.d0*sigma-2.d0*delta**2)
     C	*delta*dlog(dabs(b/a))/6.d0-5.d0/18.d0-(sigma-delta**2)
     C	/3.d0+(1.d0+sigma-2.d0*delta**2)*beta*L/6.d0

	aux=-(1.d0/2.d0*(F0-FA)-F3)
	endif

	B21=1.d0/3.d0*div+aux

	return
	end
C********************************************************************
C		->B22(p,m1,m2)
	DOUBLE PRECISION function B22(x,a,b)

	implicit none
	DOUBLE PRECISION x,a,b
	DOUBLE PRECISION B5,A0zdec

	B22=1.d0/4.d0*(A0zdec(a)+A0zdec(b)-B5(x,a,b))

	return	
	end

C********************************************************************
C	C24 with arguments (p1,(p1+p2),m1,m2,m3)
C	in our case either
C	(MTAU,MZ,SMASS,MTAU,PMASS)or (MTAU,MZ,MTAU,PMASS,MTAU)
C	for those functions see refenrence [1]

C		->C24(p1,(p1+p2),m1,m2,m3)
       DOUBLE PRECISION function C24zdec(x,y,a,b,c)

	implicit none
	DOUBLE PRECISION x,y,a,b,c,aux1,aux2,div
	DOUBLE PRECISION B0zdec,CC0,C11,C12

	div=0.d0

	If(x.eq.0.d0.and.y.eq.0.d0)then
	aux1=(b**2-c**2)/c**2
	aux2=(a**2-b**2)/c**2
	C24zdec=div/4.d0-1.d0/4.d0*dlog(c**2)-1.d0/2.d0*((-2.d0
     C		*(1.d0+aux1)**2*dlog(1.d0+aux1))/(4.d0*aux1*aux2)
     C		+(-3*aux2*(aux1+aux2)+2.d0*(1.d0+aux1+aux2)**2
     C		*dlog(1.d0+aux1+aux2))/(4.d0*aux2*(aux1+aux2)))
	else
	C24zdec=1.d0/4.d0+1.d0/4.d0*B0zdec(x,b,c)+a**2/2.d0
     C	    *CC0(x,y,a,b,c)-(b**2-a**2-x**2)/4.d0*C11(x,y,a,b,c)
     C      -(c**2-b**2-y**2+x**2)/4.d0*C12(x,y,a,b,c)
	endif

	return	
	end

C********************************************************************
C		->C23(p1,(p1+p2),m1,m2,m3)
	DOUBLE PRECISION function C23(x,y,a,b,c)

	implicit none
	DOUBLE PRECISION x,y,a,b,c
	DOUBLE PRECISION B0zdec,B1zdec,C11,C24zdec
	
	C23=1.d0/(2.d0*((y**2/2.d0-x**2)**2-x**4))*((y**2/2.d0-x**2)
     C	    *(B1zdec(y,a,c)+B0zdec(x,b,c)+(b**2-a**2-x**2)
     C	    *C11(x,y,a,b,c)-2.d0*C24zdec(x,y,a,b,c))-x**2
     C		*(B1zdec(x,a,b)
     C	    -B1zdec(y,a,c)+(c**2-b**2-y**2+x**2)*C11(x,y,a,b,c)))

	return
	end

C********************************************************************
C		->C22(p1,(p1+p2),m1,m2,m3)
	DOUBLE PRECISION function C22(x,y,a,b,c)

	implicit none
	DOUBLE PRECISION x,y,a,b,c
	DOUBLE PRECISION B1zdec,C12,C24zdec

	C22=1.d0/(2.d0*((y**2/2.d0-x**2)**2-x**4))*((y**2/2.d0-x**2)
     C	    *(B1zdec(y,a,c)-B1zdec(x,b,c)+(b**2-a**2-x**2)
     C	    *C12(x,y,a,b,c))-x**2*(-B1zdec(y,a,c)+(c**2-b**2-y**2
     C	    +x**2)*C12(x,y,a,b,c)-2.d0*C24zdec(x,y,a,b,c)))

	return
	end

C********************************************************************
C		->C12(p1,(p1+p2),m1,m2,m3)
	DOUBLE PRECISION function C12(x,y,a,b,c)

	implicit none
	DOUBLE PRECISION x,y,a,b,c
	DOUBLE PRECISION B0zdec,CC0
	
	C12=1.d0/(2.d0*((y**2/2.d0-x**2)**2-x**4))*((y**2/2.d0-x**2)
     C	    *(B0zdec(y,a,c)-B0zdec(x,b,c)+(b**2-a**2-x**2)
     C	    *CC0(x,y,a,b,c))-x**2*(B0zdec(x,a,b)-B0zdec(y,a,c)+
     C	    (c**2-b**2-y**2+x**2)*CC0(x,y,a,b,c)))

	return
	end

C********************************************************************
C		->C11(p1,(p1+p2),m1,m2,m3)
	DOUBLE PRECISION function C11(x,y,a,b,c)

	implicit none
	DOUBLE PRECISION x,y,a,b,c
	DOUBLE PRECISION B0zdec,CC0
	
	C11=1.d0/(2.d0*((y**2/2.d0-x**2)**2-x**4))*(-x**2
     C	    *(B0zdec(y,a,c)-B0zdec(x,b,c)+(b**2-a**2-x**2)
     C	    *CC0(x,y,a,b,c))+(y**2/2.d0-x**2)*(B0zdec(x,a,b)
     C	    -B0zdec(y,a,c)+(c**2-b**2-y**2+x**2)*CC0(x,y,a,b,c))) 
		    
	return	
	end

C********************************************************************
C 	Three-Point-Function CC0(mH,mTau,mA)
C	with arguments (sqrt(p1^2),sqrt((p1+p2)^2),m1,m2,m3)

C	in our case either 
C	(MTAU,MZ,MH,MTAU,MA) or (MTAU,MZ,MTAU,MA,MTAU)(MA has to be
C		 >1Gev)

C		->CC0(p1,(p1+p2),m1,m2,m3)
	DOUBLE PRECISION function CC0(x,y,a,b,c)
	
	
	implicit none
	DOUBLE PRECISION x,y,a,b,c,aux1,aux2,Li_2,pi,ME
	DOUBLE PRECISION im,re,im1,re1,im2,re2,l,u,lsec,usec
	DOUBLE PRECISION part1,part2,part3
	DOUBLE PRECISION a1,a2,a3,a4,a5,a6,a7,a8,a9
	DOUBLE PRECISION b1,b2,b3,b4,b5,b6,b7,b8,b9
	DOUBLE PRECISION u1,u2,u3,u4,u5,u6,u7,u8,u9 
	DOUBLE PRECISION l1,l2,l3,l4,l5,l6,l7,l8,l9
	DOUBLE PRECISION u12,u22,u32,u42,u52,u62,u72,u82,u92
	DOUBLE PRECISION l12,l22,l32,l42,l52,l62,l72,l82,l92
	DOUBLE PRECISION Integral0,Integral1,Integral2,Integral3,
     C			 Integral4,Integral5,Integral6,Integral7,
     C			 Integral8,Integral
	DOUBLE PRECISION Integral1a,Integral2a,Integral3a,
     C			 Integral4a,Integral5a,Integral6a,Integral7a,
     C			 Integral8a
	DOUBLE PRECISION MS,MC,MBNP,MB,MT,MTAU,MMUON,MZ,MW
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT,MTAU,MMUON,MZ,MW

        pi=4.d0*datan(1.d0)
	ME=0.5109989d-3
C********************************************************************
C	CC0(mTau,mZ,mTau,M,mTau) 

C	valid for M>1GeV and M<4700GeV (error <3*10^-5)
C		for higher masses it couldn't be compared
C		but probably also ok, because the last thing (in if)
C		is an approximation for large M 

	if(a.eq.MTAU.and.c.eq.MTAU) then
		if(b.lt.1.d0)then
		CC0=0.d0 !dlog(0.d0)
C		Print*,"MASS too small for the expansion of CC0"
		return
		endif

	if(b.lt.40.d0)then

	CC0=2.422162848961252d-25*(dlog(b**2/y**2))**20
     C	-1.026655364615622d-7*(dlog(b**2/y**2))**4
     C 	+1.351245758336911d-7*(dlog(b**2/y**2)**3)
     C	+0.00006921786603430583d0*(dlog(b**2/y**2))**2
     C	+0.00006321053292996573d0*dlog(b**2/y**2)
     C	-0.00009458391344069682d0*(b**2/y**2)**2
     C	+0.00002760266556219259d0*(b**2/y**2)
     C	-0.00006498227061315327d0


	elseif(b.ge.40.d0.and.b.lt.200.d0)then

	CC0=-1.2708278246373165d-10*(dlog(b**2/y**2))**10
     C	+8.474331037533179d-8*(dlog(b**2/y**2))**6
     C	-5.94960306674913d-6*(dlog(b**2/y**2))**4
     C	-0.00004189754315502147d0*(dlog(b**2/y**2))**3
     C	-0.00007979947907903371d0*(dlog(b**2/y**2))**2
     C	-0.00023317749608367592d0*dlog(b**2/y**2)
     C	-1.015175150608902d-8*(b**2/y**2)**4
     c  +3.754591953002911d-7
     C	*(b**2/y**2)**3-8.440677324024552d-6
     C	*(b**2/y**2)**2
     C	+0.0002496934112170587d0*(b**2/y**2)-0.0003413246732846144d0


	else
	CC0=(1.d0/y**2)*(-2.d0*dlog(dabs(x))
     C	   *dlog(dabs(b**2-2.d0*x**2))
     C	   +2.d0*dlog(dabs(b**2-x**2))*dlog(dabs(b**2-2.d0*x**2))
     C	   -dlog(dabs(y-dsqrt(y**2-4.d0*x**2)))
     C	   *dlog(dabs(b**2-2.d0*x**2))
     C	   -dlog(dabs(y+dsqrt(y**2-4.d0*x**2)))
     C	   *dlog(dabs(b**2-2.d0*x**2))+dlog(4.d0)
     C	   *dlog(dabs(b**2-2.d0*x**2))-2.d0*dlog(dabs(x))
     C	   *dlog(dabs(b**2-2.d0*x**2+y**2))-dlog(4.d0)
     C	   *dlog(dabs(b**2-2.d0*x**2+y**2))
     C	   +dlog(dabs(b**2-2.d0*x**2+y**2))
     C	   *dlog(dabs(y-dsqrt(y**2-4.d0*x**2)))
     C	   +dlog(dabs(b**2-2.d0*x**2+y**2))
     C	   *dlog(dabs(y+dsqrt(y**2-4.d0*x**2)))
     C	   +dlog(dabs(y-dsqrt(y**2-4.d0*x**2)))
     C	   *dlog(dabs(2.d0*b**2-4.d0*x**2+y**2
     C	   -y*dsqrt(y**2-4.d0*x**2)))
     C	   -dlog(dabs(y+dsqrt(y**2-4.d0*x**2)))
     C	   *dlog(dabs(2.d0*b**2-4.d0*x**2+y**2
     C	   -y*dsqrt(y**2-4.d0*x**2)))+2.d0*dlog(dabs(x))
     C	   *dlog(dabs(b**4+(y**2-4.d0*x**2)*b**2+4.d0*x**4-x**2
     C	   *y**2))-dlog(dabs(b**2-x**2))
     C	   *dlog(dabs(b**4+(y**2-4.d0*x**2)*b**2+4.d0*x**4-x**2
     C	   *y**2))-dlog(dabs(y-dsqrt(y**2-4.d0*x**2)))
     C	   *dlog(dabs(2.d0*b**2-4.d0*x**2+y*(y
     C	   +dsqrt(y**2-4.d0*x**2))))
     C	   +dlog(dabs(y+dsqrt(y**2-4.d0*x**2)))
     C	   *dlog(dabs(2.d0*b**2-4.d0*x**2+y*(y
     C	   +dsqrt(y**2-4.d0*x**2))))
     C	   +Li_2(-((y*(y+dsqrt(y**2-4.d0*x**2)))/(2.d0*b**2-4.d0
     C	   *x**2+y**2-y*dsqrt(y**2-4.d0*x**2))))
     C	   -Li_2((x**2*y**2)/(b**4+(y**2-4.d0*x**2)*b**2
     C	   +4.d0*x**4-x**2*y**2))
     C	   +Li_2(((b**2-x**2)*y**2)/(b**4+(y**2-4.d0*x**2)*b**2
     C	   +4.d0*x**4-x**2*y**2))
     C	   -Li_2((y*(dsqrt(y**2-4.d0*x**2)-y))/(-2.d0*b**2+4.d0*x**2
     C	   +y*(dsqrt(y**2-4.d0*x**2)-y)))
     C	   +Li_2((y*(dsqrt(y**2-4.d0*x**2)-y))/(2.d0*b**2-4.d0*x**2
     C	   +y*(y+dsqrt(y**2-4.d0*x**2))))
     C	   -Li_2((y*(y+dsqrt(y**2-4.d0*x**2)))/(2.d0*b**2-4.d0*x**2
     C	   +y*(y+dsqrt(y**2-4.d0*x**2)))))
	
	endif

C********************************************************************
C	CC0(mE,mZ,mE,M,mE)
	elseif(a.eq.ME.and.c.eq.ME) then
		if(b.lt.1.d0)then
		CC0=0.d0 !dlog(0.d0)
C		Print*,"MASS too small for the expansion of CC0"
		return
		endif

	if(b.lt.40.d0)then

	CC0=0.00008574719585836061d0*(b**2/y**2)**3
     C	-0.00007547106735635636d0*(b**2/y**2)**2
     C	-0.0002947217252411331d0*(b**2/y**2)+0.001643402374029581d0
     C	*(b**2/y**2)**0.5-0.01852767384672546d0*(b**2/y**2)**0.2
     C	+0.09876471392235593d0*(b**2/y**2)**0.1-2.9797518812515396d0
     C	*1.d-6*dlog(b**2/y**2)**3-0.0001620313390967642d0
     C	*dlog(b**2/y**2)**2-0.006686775093102159d0*dlog(b**2/y**2)
     C	-0.08166603095558493d0


	elseif(b.ge.40.d0.and.b.lt.200.d0)then

	CC0=0.004617968622607519d0*dlog(b**2/y**2)**3
     C	+0.20821093696780776d0*dlog(b**2/y**2)**2
     C	+4.468152092977522d0*dlog(b**2/y**2)+6.748667132203143d-7
     C *(b**2/y**2)**2-0.00008127241562309546d0
     C	*(b**2/y**2)-0.020016553218186434d0*(b**2/y**2)**0.5
     C	+4.794776678249852d0*(b**2/y**2)**0.2-66.47576546581614d0
     C	*(b**2/y**2)**0.1+24.614214765920487d0*(b**2/y**2)**0.05
     C	+0.00015482292195980802d0/(b**2/y**2)-7.008471193915071d-7
     c  /(b**2/y**2)**2+37.0866181375125d0



	else
	CC0=(1.d0/y**2)*(-2.d0*dlog(dabs(x))
     C	   *dlog(dabs(b**2-2.d0*x**2))
     C	   +2.d0*dlog(dabs(b**2-x**2))*dlog(dabs(b**2-2.d0*x**2))
     C	   -dlog(dabs(y-dsqrt(y**2-4.d0*x**2)))
     C	   *dlog(dabs(b**2-2.d0*x**2))
     C	   -dlog(dabs(y+dsqrt(y**2-4.d0*x**2)))
     C	   *dlog(dabs(b**2-2.d0*x**2))+dlog(4.d0)
     C	   *dlog(dabs(b**2-2.d0*x**2))-2.d0*dlog(dabs(x))
     C	   *dlog(dabs(b**2-2.d0*x**2+y**2))-dlog(4.d0)
     C	   *dlog(dabs(b**2-2.d0*x**2+y**2))
     C	   +dlog(dabs(b**2-2.d0*x**2+y**2))
     C	   *dlog(dabs(y-dsqrt(y**2-4.d0*x**2)))
     C	   +dlog(dabs(b**2-2.d0*x**2+y**2))
     C	   *dlog(dabs(y+dsqrt(y**2-4.d0*x**2)))
     C	   +dlog(dabs(y-dsqrt(y**2-4.d0*x**2)))
     C	   *dlog(dabs(2.d0*b**2-4.d0*x**2+y**2
     C	   -y*dsqrt(y**2-4.d0*x**2)))
     C	   -dlog(dabs(y+dsqrt(y**2-4.d0*x**2)))
     C	   *dlog(dabs(2.d0*b**2-4.d0*x**2+y**2
     C	   -y*dsqrt(y**2-4.d0*x**2)))+2.d0*dlog(dabs(x))
     C	   *dlog(dabs(b**4+(y**2-4.d0*x**2)*b**2+4.d0*x**4-x**2
     C	   *y**2))-dlog(dabs(b**2-x**2))
     C	   *dlog(dabs(b**4+(y**2-4.d0*x**2)*b**2+4.d0*x**4-x**2
     C	   *y**2))-dlog(dabs(y-dsqrt(y**2-4.d0*x**2)))
     C	   *dlog(dabs(2.d0*b**2-4.d0*x**2+y*(y
     C	   +dsqrt(y**2-4.d0*x**2))))
     C	   +dlog(dabs(y+dsqrt(y**2-4.d0*x**2)))
     C	   *dlog(dabs(2.d0*b**2-4.d0*x**2+y*(y
     C	   +dsqrt(y**2-4.d0*x**2))))
     C	   +Li_2(-((y*(y+dsqrt(y**2-4.d0*x**2)))/(2.d0*b**2-4.d0
     C	   *x**2+y**2-y*dsqrt(y**2-4.d0*x**2))))
     C	   -Li_2((x**2*y**2)/(b**4+(y**2-4.d0*x**2)*b**2
     C	   +4.d0*x**4-x**2*y**2))
     C	   +Li_2(((b**2-x**2)*y**2)/(b**4+(y**2-4.d0*x**2)*b**2
     C	   +4.d0*x**4-x**2*y**2))
     C	   -Li_2((y*(dsqrt(y**2-4.d0*x**2)-y))/(-2.d0*b**2+4.d0*x**2
     C	   +y*(dsqrt(y**2-4.d0*x**2)-y)))
     C	   +Li_2((y*(dsqrt(y**2-4.d0*x**2)-y))/(2.d0*b**2-4.d0*x**2
     C	   +y*(y+dsqrt(y**2-4.d0*x**2))))
     C	   -Li_2((y*(y+dsqrt(y**2-4.d0*x**2)))/(2.d0*b**2-4.d0*x**2
     C	   +y*(y+dsqrt(y**2-4.d0*x**2)))))
	
	endif

C********************************************************************
C	CC0(MTAU,MZ,0,CMASS,0) 

	elseif(a.eq.0.d0.and.c.eq.0.d0.and.x.eq.mtau)then

	   if(b.lt.100.d0)then
	   CC0=-2.6581791132934073d-6*(dlog(b**2/y**2))**4
     C	   -0.00002053438809404773d0*(dlog(b**2/y**2))**3
     C	   -1.260712348861513d-6*(dlog(b**2/y**2))**2
     C	   -0.00006278129911601018d0*dlog(b**2/y**2)
     C	   +0.00006273647510478084d0*(b**2/y**2)
     C	   -0.0001616317329574659d0

  	   elseif(b.lt.200.d0.and.b.ge.100.d0)then
	   CC0=3.687391502175193d-15*(b**2/y**2)**4
     C	   -1.2908901772437909d-12*(b**2/y**2)**3
     C	   -5.769358559836516d-9*(b**2/y**2)**2
     C	   +0.00012428473502946574d0*(b**2/y**2)
     C	   -1.0388960301509819d-10*(dlog(b**2/y**2))**10
     C	   -4.151314242891012d-9*(dlog(b**2/y**2))**8
     C	   -4.223653583097421d-7*(dlog(b**2/y**2))**6
     C	   -5.850866381193399d-6*(dlog(b**2/y**2))**4
     C	   -0.000030219060475583603d0*(dlog(b**2/y**2))**3
     C	   -0.00003226204077055742d0*(dlog(b**2/y**2))**2
     C	   -0.00012426142777669476d0*(dlog(b**2/y**2))
     C	   -0.0002231795532464488d0


 	   else
	   CC0=-(1.d0/(6.d0*y**2))*(3.d0*dlog(dabs(b**2-x**2))**2
     C	      +6.d0*dlog(dabs(1.d0/(x**2-b**2)))
     C	      *dlog(dabs(b**2-x**2))+3.d0
     C	      *dlog(dabs(1/(x**2-b**2)))**2+12.d0*dlog(y)**2
     C	      -3.d0*dlog(dabs((b**2-x**2)/(b**2-x**2-y**2)))**2
     C	      -3.d0*dlog(dabs(1.d0/(-b**2+x**2+y**2)))**2
     C         -12.d0*dlog(dabs(y))*dlog(dabs(b**2-x**2+y**2))
     C	      +6.d0*dlog(dabs((b**2-x**2)/(b**2-x**2-y**2)))
     C	      *dlog(dabs(b**2-x**2+y**2))+6.d0
     C	      *dlog(dabs((b**2-x**2)/(b**2-x**2-y**2)))
     C	      *dlog(dabs(1.d0/(-b**2+x**2+y**2)))-6.d0
     C	      *dlog(dabs(b**2-x**2+y**2))
     C	      *dlog(dabs(1.d0/(-b**2+x**2+y**2)))
     C	      +6.d0*Li_2((x**2-b**2)/y**2)-dlog(4.d0)*dlog(8.d0)
     C	      +6.d0*dlog(2.d0)**2+Pi**2)

	   endif

C********************************************************************
C	CC0(ME,MZ,0,CMASS,0) 
C	is approximated by ME~0, precision heigher than 1 Permil

	elseif(a.eq.0.d0.and.c.eq.0.d0.and.x.eq.ME)then

	CC0=1.d0/y**2*(dlog(y**2/b**2)*dlog(y**2/b**2+1.d0)
     C		+Li_2(-y**2/b**2))

C********************************************************************
C	CC0(MTAU,MZ,CMASS,0,CMASS)
 
C	here was no analytic approximation for large CMASS made
	elseif(b.eq.0.d0.and.x.eq.mtau)then
	   if(a.lt.100.d0)then
	   CC0=-7.036290232116209d-6*(dlog(a**2/y**2))**10
     C	   -0.00007827711463961622d0*(dlog(a**2/y**2))**4
     C	   -0.00012468133759555692d0*(dlog(a**2/y**2))**3
     C	   -0.0005593194452141794d0*(dlog(a**2/y**2))**2
     C	   -0.0007939255576352166d0*dlog(a**2/y**2)
     C	   +0.0009397059398562709d0*(a**2/y**2)
     C	   -0.0010716021973681098d0

	   elseif(a.lt.1000.d0.and.a.ge.100.d0)then


	   CC0=8.643334009196724d-16*(dlog(a**2/y**2))**20
     C	   +9.290386413384689d-13*(dlog(a**2/y**2))**16
     C	   -1.918135637679234d-10*(dlog(a**2/y**2))**12
     C	   -1.7028329227032502d-8*(dlog(a**2/y**2))**10
     C	   -4.454240273260743d-7*(dlog(a**2/y**2))**8
     C	   -5.953508220381455d-6*(dlog(a**2/y**2))**6
     C	   -0.0000315976348080128d0*(dlog(a**2/y**2))**4
     C	   +0.000044545718648148424d0*(dlog(a**2/y**2))**3
     C	   -9.917090550526695d-8*(dlog(a**2/y**2))**2
     C	   +0.00040437219710217894d0*(dlog(a**2/y**2))
     C	   -1.0698084672457089d-27*(a**2/y**2)**10
     C	   +1.8216000420923203d-16*(a**2/y**2)**6
     C	   -2.4619719700618396d-13*(a**2/y**2)**5
     C	   +2.2210967144877502d-10*(a**2/y**2)**4
     C	   -2.1198746248746692d-7*(a**2/y**2)**3
     C	   +0.0000417536722113765d0*(a**2/y**2)**2
     C	   -0.00034179706655195944d0*(a**2/y**2)
     C	   +0.00016834917895022482d0

	   elseif(a.lt.5000.d0.and.a.ge.1000.d0)then
	   CC0=-1.1052741214174954d-18*(dlog(a**2/y**2))**14
     C	   -2.001214666456646d-14*(dlog(a**2/y**2))**10
     C	   -7.094657711644969d-8*(dlog(a**2/y**2))**4
     C	   +1.2857919021074231d-6*(dlog(a**2/y**2))**3
     C	   -0.00001009166456006905d0*(dlog(a**2/y**2))**2
     C	   +0.00003834068898729404d0*dlog(a**2/y**2)
     C	   +1.9107966415303365d-8*(a**2/y**2)
     C	   -0.00005923870062006238d0

	else
	CC0=0.d0 !dlog(0.d0)
	   endif

C********************************************************************
C	CC0(ME,MZ,CMASS,0,CMASS)
C	is approximated by ME~0, precision heigher than 1 Permil
C	for CMASS>

	elseif(b.eq.0.d0.and.x.eq.ME)then

	CC0=(2.d0*dlog(dabs(a/y))*dlog(dabs(a**2/(a**2-y**2)))
     C	-Li_2(y**2/a**2))/y**2-((dlog(a**2/y**2-1.d0/4.d0)
     C	*dlog(dabs(a**2/(a**2-y**2))))/y**2)+
     C	(64.d0*dsqrt((4.d0*a**2)/y**2-1.d0)*dlog(dabs((-2.d0*a**2
     C	+y**2+(y
     C	*dsqrt(4.d0*a**2-y**2))/dsqrt((4.d0*a**2)/y**2-1.d0))/(-2.d0
     C	*a**2+y**2-(y*dsqrt(4.d0*a**2-y**2))/dsqrt((4.d0*a**2)/y**2
     C	-1.d0))))*a**12)/(3.d0*y**7*dsqrt((4.d0*a**2-y**2)**7))+(96.d0
     C	*dsqrt((4.d0*a**2)/y**2-1.d0)*dlog(dabs((-2.d0*a**2+y**2-(y
     C	*dsqrt(4.d0*a**2-y**2))/dsqrt((4.d0*a**2)/y**2-1.d0))/(-2.d0
     C	*a**2+y**2+(y*dsqrt(4.d0*a**2-y**2))/dsqrt((4.d0*a**2)/y**2
     C	-1.d0))))*a**10)/(y**5*dsqrt((4.d0*a**2-y**2)**7))-(64.d0
     C	*a**10)/(3.d0*y**6*(y**2-4.d0*a**2)**3)+(216.d0*dsqrt((4.d0
     C	*a**2)/y**2-1.d0)*dlog(dabs((-2.d0*a**2+y**2+(y*dsqrt(4.d0*a**2
     C	-y**2))/dsqrt((4.d0*a**2)/y**2-1.d0))/(-2.d0*a**2+y**2-(y
     C	*dsqrt(4.d0*a**2-y**2))/dsqrt((4.d0*a**2)/y**2-1.d0))))*a**8)
     C	/(y**3*dsqrt((4.d0*a**2-y**2)**7))+(256.d0*a**8)/(3.d0*y**4
     C	*(y**2-4.d0*a**2)**3)+(640.d0*dsqrt((4.d0*a**2)/y**2-1.d0)
     C	*dlog(dabs((-2.d0*a**2+y**2-(y*dsqrt(4.d0*a**2-y**2))
     C	/dsqrt((4.d0
     C	*a**2)/y**2-1.d0))/(-2.d0*a**2+y**2+(y*dsqrt(4.d0*a**2-y**2))
     C	/dsqrt((4.d0*a**2)/y**2-1.d0))))*a**6)/(3.d0*y*dsqrt((4.d0
     C	*a**2-y**2)**7))+(16.d0*a**6)/(9.d0*((4.d0*a**2)/y**2-1.d0)
     C	*y**4*(y**2-4.d0*a**2)**2)-(520.d0*a**6)/(3.d0*y**2*(y**2
     C	-4.d0*a**2)**3)+(100.d0*dsqrt((4.d0*a**2)/y**2-1.d0)*y
     C	*dlog(dabs((-2.d0*a**2+y**2+(y*dsqrt(4.d0*a**2-y**2))
     C	/dsqrt((4.d0*a**2)/y**2-1.d0))/(-2.d0*a**2+y**2-(y
     C	*dsqrt(4.d0*a**2-y**2))/dsqrt((4.d0*a**2)/y**2-1.d0))))
     C	*a**4)/dsqrt((4.d0*a**2-y**2)**7)-(16.d0*a**4)/(3.d0*((4.d0
     C	*a**2)/y**2-1.d0)*y**2*(y**2-4.d0*a**2)**2)+(380.d0*a**4)
     C	/(3.d0*(y**2-4.d0*a**2)**3)+(22.d0*dsqrt((4.d0*a**2)/y**2
     C	-1.d0)*y**3*dlog(dabs((-2.d0*a**2+y**2-(y
     C	*dsqrt(4.d0*a**2-y**2))
     C	/dsqrt((4.d0*a**2)/y**2-1.d0))/(-2.d0*a**2+y**2+(y*dsqrt(4.d0
     C	*a**2-y**2))/dsqrt((4.d0*a**2)/y**2-1.d0))))*a**2)/dsqrt((4.d0
     C	*a**2-y**2)**7)-(4.d0*a**2)/(15.d0*((4.d0*a**2)/y**2-1.d0)**2
     C	*y**2*(y**2-4.d0*a**2))+(10.d0*a**2)/(3.d0*((4.d0*a**2)/y**2
     C	-1.d0)*(y**2-4.d0*a**2)**2)-(110.d0*y**2*a**2)/(3.d0*(y**2
     C	-4.d0*a**2)**3)+(11.d0*dsqrt((4.d0*a**2)/y**2-1.d0)*y**5
     C	*dlog(dabs((-2.d0*a**2+y**2+(y*dsqrt(4.d0*a**2-y**2))
     C	/dsqrt((4.d0*a**2)/y**2-1.d0))/(-2.d0*a**2+y**2-(y
     C	*dsqrt(4.d0*a**2-y**2))/dsqrt((4.d0*a**2)/y**2-1.d0)))))
     C	/(6.d0*dsqrt((4.d0*a**2-y**2)**7))+2/(15.d0*((4.d0*a**2)
     C	/y**2-1.d0)**2*(y**2-4.d0*a**2))-(5.d0*y**2)/(9.d0*((4.d0
     C	*a**2)/y**2-1.d0)*(y**2-4.d0*a**2)**2)+(11.d0*y**4)/(3.d0
     C	*(y**2-4.d0*a**2)**3)

C	If(a.lt.y) CC0=dlog(0.d0)

C********************************************************************
C	CC0(0.d0,0.d0,m1,m2,m3)
	
	elseif(x.eq.0.d0.and.y.eq.0.d0)then

	aux1=(b**2-c**2)/c**2
	aux2=(a**2-b**2)/c**2

	CC0=-1.d0/c**2*(-((1.d0+aux1)*dlog(1.d0+aux1))/(aux1*aux2)
     C	+((1.d0+aux1+aux2)*dlog(1.d0+aux1+aux2))/((aux1+aux2)*aux2))

C********************************************************************
C	CC0(mH,mTau(ME),mA)	
	elseif(b.eq.mtau.or.b.eq.ME)then

	if(a.lt.c)then
	aux1=c
	aux2=a
	else
	aux1=a
	aux2=c
	endif


C		Solution of part 1, Form: a/(b*y+c)
	part1=(dlog(b**2/(y**2))*dlog(dabs((-y**2-b**2+aux1**2)/
     C	(-2.d0*b**2+aux1**2))))/(b**2-y**2)
	


C		Definition of imaginary and real part of the roots
	im1=1.d0/(2.d0*y**2)*dsqrt(4.d0*aux1**2.d0*y**2	
     C	-(aux1**2+y**2-aux2**2)**2)
	re1=1.d0/(2.d0*y**2)*(aux1**2+y**2-aux2**2)
	im2=1.d0/(2.d0*b**2)*dsqrt(4.d0*b**4
     C		-(aux2**2-2.d0*b**2)**2)
	re2=1.d0/(2.d0*b**2)*(2.d0*b**2-aux2**2)
	
	
C		Solution of part 2
	if(((aux1**2+y**2-aux2**2)**2-4.d0*aux1**2*y**2.d0)
     C		.ge.0.d0)then
C		Evaluation for real roots
	part2=-(Integral0(1.d0/(2.d0*y**2)*(aux1**2+y**2-aux2**2
     C	+dsqrt((aux1**2+y**2-aux2**2)**2-4.d0*aux1**2*y**2)),
     C	(b**2-y**2),(aux1**2-2.d0*b**2))+Integral0(1.d0
     C	/(2.d0*y**2)*(aux1**2+y**2-aux2**2-dsqrt((aux1**2+y**2
     C	-aux2**2)**2-4.d0*aux1**2*y**2)),(b**2-y**2),		
     C	(aux1**2-2.d0*b**2)))



	else
C		calculating the upper and lower bound of the Integral
	l=(-re1)/im1
	u=(1.d0-re1)/im1


C		Coefficients in front of Integral to set them =0 if 
C		its not necessary to integrate in this intervall
	a1=1.d0
	a2=1.d0
	a3=1.d0
	a4=1.d0
	a5=1.d0
	a6=1.d0
	a7=1.d0
	a8=1.d0
	a9=1.d0
	
	u1=-3.7d0
	l2=-3.7d0
	u2=-1.53d0
	l3=-1.53d0
	u3=-0.68d0
	l4=-0.68d0
	u4=-0.23d0
	l5=-0.23d0
	u5=0.23d0
	l6=0.23d0
	u6=0.68d0
	l7=0.68d0
	u7=1.53d0
	l8=1.53d0
	u8=3.7d0
	l9=3.7d0

C	 just to define that the program can handle the formula
C	 in which they appear	
	l1=4.d0
	u9=4.d0
	
	
	
C		checking in which intervall u lies
	If(u.le.-3.7d0)then
	u1=u
	a2=0.d0
	a3=0.d0
	a4=0.d0
	a5=0.d0
	a6=0.d0
	a7=0.d0
	a8=0.d0
	a9=0.d0
	endif
	If(-3.7d0.lt.u.and.u.le.-1.53d0)then
	u2=u
	a3=0.d0
	a4=0.d0
	a5=0.d0
	a6=0.d0
	a7=0.d0
	a8=0.d0
	a9=0.d0
	endif
	If(-1.53d0.lt.u.and.u.le.-0.68d0)then
	u3=u
	a4=0.d0
	a5=0.d0
	a6=0.d0
	a7=0.d0
	a8=0.d0
	a9=0.d0
	endif
	If(-0.68d0.lt.u.and.u.le.-0.23d0)then
	u4=u
	a5=0.d0
	a6=0.d0
	a7=0.d0
	a8=0.d0
	a9=0.d0
	endif
	If(-0.23d0.lt.u.and.u.le.0.23d0)then
	u5=u
	a6=0.d0
	a7=0.d0
	a8=0.d0
	a9=0.d0
	endif
	If(0.23d0.lt.u.and.u.le.0.68d0)then
	u6=u
	a7=0.d0
	a8=0.d0
	a9=0.d0
	endif
	If(0.68d0.lt.u.and.u.le.1.53d0)then
	u7=u
	a8=0.d0
	a9=0.d0
	endif
	If(1.53d0.lt.u.and.u.le.3.7d0)then
	u8=u
	a9=0.d0
	endif
	If(3.7d0.lt.u)then
	u9=u
	endif

C		checking in which intervall l lies

	If(l.le.-3.7d0)then
	l1=l
	endif
	If(-3.7d0.lt.l.and.l.le.-1.53d0)then
	l2=l
	a1=0.d0
	endif
	If(-1.53d0.lt.l.and.l.le.-0.68d0)then
	l3=l
	a1=0.d0
	a2=0.d0
	endif
	If(-0.68d0.lt.l.and.l.le.-0.23d0)then
	l4=l
	a1=0.d0
	a2=0.d0
	a3=0.d0
	endif
	If(-0.23d0.lt.l.and.l.le.0.23d0)then
	l5=l
	a1=0.d0
	a2=0.d0
	a3=0.d0
	a4=0.d0
	endif
	
	


	part2=-(1.d0/(b**2-y**2)
     C	*(a1*(Integral1(u1,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1)-Integral1(l1,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1))
     C	+a9*(Integral1(u9,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1)-Integral1(l9,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1))
     C  +a2*(Integral2(u2,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1)-Integral2(l2,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1))
     C	+a3*(Integral3(u3,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1)-Integral3(l3,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1))
     C	+a4*(Integral4(u4,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1)-Integral4(l4,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1))
     C	+a5*(Integral5(u5,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1)-Integral5(l5,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1))
     C	+a6*(Integral6(u6,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1)-Integral6(l6,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1))
     C	+a7*(Integral7(u7,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1)-Integral7(l7,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1))
     C	+a8*(Integral8(u8,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1)-Integral8(l8,(aux1**2-2.d0*b**2)
     C	/((b**2-y**2)*im1)+re1/im1))-2.d0*dlog(dabs(im1))
     C	*dlog(dabs((2.d0*b**2-aux1**2)/(-y**2-b**2+aux1**2)))))
	endif
	
C		Solution of part 3
	if(((aux2**2-2.d0*b**2)**2-4.d0*b**4.d0).ge.(0.d0))then

C		Evaluation for real roots

	If(b.eq.mtau)then
	part3=Integral0(1.d0/(2.d0*b**2)*(2.d0*b**2-aux2**2
     C	+dsqrt((aux2**2-2.d0*b**2)**2-4.d0*b**4)),(b**2-y**2),
     C	(aux1**2-2.d0*b**2))+Integral0(1.d0/(2.d0*b**2)
     C	*(2.d0*b**2-aux2**2
     C	-dsqrt((aux2**2-2.d0*b**2)**2-4.d0*b**4)),(b**2-y**2),
     C	(aux1**2-2.d0*b**2))

	elseif(b.eq.ME)then
	part3=Integral0(-b**2/(aux2**2-2.d0*b**2),(b**2-y**2),
     C	(aux1**2-2.d0*b**2))+Integral0(1.d0/(2.d0*b**2)
     C	*(2.d0*b**2-aux2**2
     C	-dsqrt((aux2**2-2.d0*b**2)**2-4.d0*b**4)),(b**2-y**2),
     C	(aux1**2-2.d0*b**2))
	endif

	else

C		Calculating the lower and upper bound of the Integral
	lsec=(-re2)/im2
	usec=(1.d0-re2)/im2



C		Coefficients in front of Integral to set them =0 if 
C		its not necessary to integrate in this intervall*)
	b1=1.d0
	b2=1.d0
	b3=1.d0
	b4=1.d0
	b5=1.d0
	b6=1.d0
	b7=1.d0
	b8=1.d0
	b9=1.d0

	u12=-3.7d0
	l22=-3.7d0
	u22=-1.53d0
	l32=-1.53d0
	u32=-0.68d0
	l42=-0.68d0
	u42=-0.23d0
	l52=-0.23d0
	u52=0.23d0
	l62=0.23d0
	u62=0.68d0
	l72=0.68d0
	u72=1.53d0
	l82=1.53d0
	u82=3.7d0
	l92=3.7d0

	l12=3.d0
	u92=3.d0

C		checking in which intervall usec lies
	If(usec.le.-3.7d0)then
	u12=usec
	b2=0.d0
	b3=0.d0
	b4=0.d0
	b5=0.d0
	b6=0.d0
	b7=0.d0
	b8=0.d0
	b9=0.d0
	endif
	If(-3.7d0.lt.usec.and.usec.le.-1.53d0)then
	u22=usec
	b3=0.d0
	b4=0.d0
	b5=0.d0
	b6=0.d0
	b7=0.d0
	b8=0.d0
	b9=0.d0
	endif
	If(-1.53d0.lt.usec.and.usec.le.-0.68d0)then
	u32=usec
	b4=0.d0
	b5=0.d0
	b6=0.d0
	b7=0.d0
	b8=0.d0
	b9=0.d0
	endif
	If(-0.68d0.lt.usec.and.usec.le.-0.23d0)then
	u42=usec
	b5=0.d0
	b6=0.d0
	b7=0.d0
	b8=0.d0
	b9=0.d0
	endif
	If(-0.23d0.lt.usec.and.usec.le.0.23d0)then
	u52=usec
	b6=0.d0
	b7=0.d0
	b8=0.d0	
	b9=0.d0
	endif
	If(0.23d0.lt.usec.and.usec.le.0.68d0)then
	u62=usec
	b7=0.d0
	b8=0.d0
	b9=0.d0
	endif
	If(0.68d0.lt.usec.and.usec.le.1.53d0)then
	u72=usec
	b8=0.d0
	b9=0.d0
	endif
	If(1.53d0.lt.usec.and.usec.le.3.7d0)then
	u82=usec
	b9=0.d0
	endif
	If(3.7.lt.usec)then
	u92=usec
	endif

C		checking in which intervall lsec lies
	If(lsec.le.-3.7d0)then
	l12=lsec
	endif
	If(-3.7d0.lt.lsec.and.lsec.le.-1.53d0)then
	l22=lsec
	b1=0.d0
	endif
	If(-1.53d0.lt.lsec.and.lsec.le.-0.68d0)then
	l32=lsec
	b1=0.d0
	b2=0.d0
	endif
	If(-0.68d0.lt.lsec.and.lsec.le.-0.23d0)then
	l42=lsec
	b1=0.d0
	b2=0.d0
	b3=0.d0
	endif
	If(-0.23d0.lt.lsec.and.lsec.le.0.23d0)then
	l52=lsec
	b1=0.d0
	b2=0.d0
	b2=0.d0
	b4=0.d0
	endif
	If(0.23d0.lt.lsec.and.lsec.le.0.68d0)then
	l62=lsec
	b1=0.d0
	b2=0.d0
	b3=0.d0
	b4=0.d0
	b5=0.d0
	endif
	If(0.68d0.lt.lsec.and.lsec.le.1.53d0)then
	l72=lsec
	b1=0.d0
	b2=0.d0
	b3=0.d0
	b4=0.d0
	b5=0.d0
	b6=0.d0
	endif
	If(1.53d0.lt.lsec.and.lsec.le.3.7)then
	l82=lsec
	b1=0.d0
	b2=0.d0
	b3=0.d0
	b4=0.d0
	b5=0.d0
	b6=0.d0
	b7=0.d0
	endif
	If(3.7.lt.lsec)then
	l92=lsec
	b1=0.d0
	b2=0.d0
	b3=0.d0
	b4=0.d0
	b5=0.d0
	b6=0.d0
	b7=0.d0
	b8=0.d0
	endif


	part3=1.d0/(b**2-y**2)*
     C	(b1*(Integral1(u12,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)
     C		+re2/im2)
     C	-Integral1(l12,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)+re2/im2))
     C	+b9*(Integral1(u92,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)
     C		+re2/im2)
     C	-Integral1(l92,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)+re2/im2))
     C	+b2*(Integral2(u22,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)
     C 		+re2/im2)
     C	-Integral2(l22,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)+re2/im2))
     C	+b3*(Integral3(u32,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)
     C		+re2/im2)
     C	-Integral3(l32,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)+re2/im2))
     C	+b4*(Integral4(u42,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)
     C		+re2/im2)
     C	-Integral4(l42,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)+re2/im2))
     C	+b5*(Integral5(u52,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)
     C		+re2/im2)
     C	-Integral5(l52,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)+re2/im2))
     C	+b6*(Integral6(u62,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)
     C		+re2/im2)
     C	-Integral6(l62,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)+re2/im2))
     C	+b7*(Integral7(u72,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)
     C		+re2/im2)
     C	-Integral7(l72,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)+re2/im2))
     C	+b8*(Integral8(u82,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)
     C		+re2/im2)
     C	-Integral8(l82,(aux1**2-2.d0*b**2)/((b**2-y**2)*im2)+re2/im2))
     C	+2.d0*dlog(dabs(im2))
     C	*dlog(dabs((-y**2-b**2+aux1**2)/(-2.d0*b**2+aux1**2))))
	endif
	
	CC0=part1+part2+part3
C********************************************************************

C	CC0(m0,m1,m2) with m0,m1,m2>>mtau
	else
	

	aux1=a
	aux2=c

	If(dabs(MZ**2+b**2-aux1**2).le.100000.d0
     C	.and.dabs(MZ**2+b**2-aux1**2).le.dabs(MZ**2+b**2-aux2**2))
     C	then
	aux1=c
	aux2=a
	endif

	part1=-dlog(dabs(aux2**2-b**2))/y**2
     C		*dlog(dabs((y**2+b**2-aux1**2)/(b**2-aux1**2)))
     C		-Integral0(-b**2/(aux2**2-b**2),y**2,(b**2-aux1**2))

C	for real roots	
	if((aux2**2-aux1**2-y**2)**2-4.d0*aux1**2*y**2.ge.0.d0)then

	part2=dlog(y**2)/y**2
     C	      *dlog(dabs((y**2+b**2-aux1**2)/(b**2-aux1**2)))
     C	      +Integral0(1.d0/(2.d0*y**2)*(aux1**2+y**2-aux2**2
     C	      +dsqrt((aux2**2-aux1**2-y**2)**2-4.d0*aux1**2*y**2))
     C		,y**2,b**2-aux1**2)
     C	      +Integral0(1.d0/(2.d0*y**2)*(aux1**2+y**2-aux2**2
     C	      -dsqrt((aux2**2-aux1**2-y**2)**2-4.d0*aux1**2*y**2))
     C		,y**2,b**2-aux1**2)

C	for complex roots
	else
	im=1.d0/(2.d0*y**2)
     C	   *dsqrt(4.d0*aux1**2*y**2-(aux2**2-aux1**2-y**2)**2)
	re=1.d0/(2.d0*y**2)*(aux1**2+y**2-aux2**2)

C	calculating the upper and lower bound of the Integral
	l=(-re)/im
	u=(1.d0-re)/im


C		Coefficients in front of Integral to set them =0 if 
C		its not necessary to integrate in this intervall
	a1=1.d0
	a2=1.d0
	a3=1.d0
	a4=1.d0
	a5=1.d0
	a6=1.d0
	a7=1.d0
	a8=1.d0
	a9=1.d0
	
	u1=-3.7d0
	l2=-3.7d0
	u2=-1.53d0
	l3=-1.53d0
	u3=-0.68d0
	l4=-0.68d0
	u4=-0.23d0
	l5=-0.23d0
	u5=0.23d0
	l6=0.23d0
	u6=0.68d0
	l7=0.68d0
	u7=1.53d0
	l8=1.53d0
	u8=3.7d0
	l9=3.7d0

C	 just to define that the program can handle the formula
C	 in which they appear	
	l1=4.d0
	u9=4.d0
	
	
	
C		checking in which intervall u lies
	If(u.le.-3.7d0)then
	u1=u
	a2=0.d0
	a3=0.d0
	a4=0.d0
	a5=0.d0
	a6=0.d0
	a7=0.d0
	a8=0.d0
	a9=0.d0
	endif
	If(-3.7d0.lt.u.and.u.le.-1.53d0)then
	u2=u
	a3=0.d0
	a4=0.d0
	a5=0.d0
	a6=0.d0
	a7=0.d0
	a8=0.d0
	a9=0.d0
	endif
	If(-1.53d0.lt.u.and.u.le.-0.68d0)then
	u3=u
	a4=0.d0
	a5=0.d0
	a6=0.d0
	a7=0.d0
	a8=0.d0
	a9=0.d0
	endif
	If(-0.68d0.lt.u.and.u.le.-0.23d0)then
	u4=u
	a5=0.d0
	a6=0.d0
	a7=0.d0
	a8=0.d0
	a9=0.d0
	endif
	If(-0.23d0.lt.u.and.u.le.0.23d0)then
	u5=u
	a6=0.d0
	a7=0.d0
	a8=0.d0
	a9=0.d0
	endif
	If(0.23d0.lt.u.and.u.le.0.68d0)then
	u6=u
	a7=0.d0
	a8=0.d0
	a9=0.d0
	endif
	If(0.68d0.lt.u.and.u.le.1.53d0)then
	u7=u
	a8=0.d0
	a9=0.d0
	endif
	If(1.53d0.lt.u.and.u.le.3.7d0)then
	u8=u
	a9=0.d0
	endif
	If(3.7d0.lt.u)then
	u9=u
	endif

C		checking in which intervall l lies

	If(l.le.-3.7d0)then
	l1=l
	endif
	If(-3.7d0.lt.l.and.l.le.-1.53d0)then
	l2=l
	a1=0.d0
	endif
	If(-1.53d0.lt.l.and.l.le.-0.68d0)then
	l3=l
	a1=0.d0
	a2=0.d0
	endif
	If(-0.68d0.lt.l.and.l.le.-0.23d0)then
	l4=l
	a1=0.d0
	a2=0.d0
	a3=0.d0
	endif
	If(-0.23d0.lt.l.and.l.le.0.23d0)then
	l5=l
	a1=0.d0
	a2=0.d0
	a3=0.d0
	a4=0.d0
	endif
	If(0.23d0.lt.l.and.l.le.0.68d0)then
	l6=l
	a1=0.d0
	a2=0.d0
	a3=0.d0
	a4=0.d0
	a5=0.d0
	endif
	If(0.68d0.lt.l.and.l.le.1.53d0)then
	l7=l
	a1=0.d0
	a2=0.d0
	a3=0.d0
	a4=0.d0
	a5=0.d0
	a6=0.d0
	endif
	If(1.53d0.lt.l.and.l.le.3.7d0)then
	l8=l
	a1=0.d0
	a2=0.d0
	a3=0.d0
	a4=0.d0
	a5=0.d0
	a6=0.d0
	a7=0.d0
	endif
	If(3.7d0.lt.l)then
	l9=l
	a1=0.d0
	a2=0.d0
	a3=0.d0
	a4=0.d0
	a5=0.d0
	a6=0.d0
	a7=0.d0
	a8=0.d0
	endif
	


C 	discremenate between large g and small g large g= g>100GeV	
	If((b**2-aux1**2+re*y**2)/(im*y**2).le.100.d0)then
	
	part2=1.d0/y**2*
     C 		(a1*(Integral1(u1,(b**2-aux1**2+re*y**2)/(im*y**2))
     C	   	    -Integral1(l1,(b**2-aux1**2+re*y**2)/(im*y**2)))
     C		+a9*(Integral1(u9,(b**2-aux1**2+re*y**2)/(im*y**2))
     C	    	    -Integral1(l9,(b**2-aux1**2+re*y**2)/(im*y**2)))
     C		+a2*(Integral2(u2,(b**2-aux1**2+re*y**2)/(im*y**2))
     C	    	    -Integral2(l2,(b**2-aux1**2+re*y**2)/(im*y**2)))
     C		+a3*(Integral3(u3,(b**2-aux1**2+re*y**2)/(im*y**2))
     C	    	    -Integral3(l3,(b**2-aux1**2+re*y**2)/(im*y**2)))
     C		+a4*(Integral4(u4,(b**2-aux1**2+re*y**2)/(im*y**2))
     C	    	    -Integral4(l4,(b**2-aux1**2+re*y**2)/(im*y**2)))
     C		+a5*(Integral5(u5,(b**2-aux1**2+re*y**2)/(im*y**2))
     C	    	    -Integral5(l5,(b**2-aux1**2+re*y**2)/(im*y**2)))
     C		+a6*(Integral6(u6,(b**2-aux1**2+re*y**2)/(im*y**2))
     C	    	    -Integral6(l6,(b**2-aux1**2+re*y**2)/(im*y**2)))
     C		+a7*(Integral7(u7,(b**2-aux1**2+re*y**2)/(im*y**2))
     C	    	    -Integral7(l7,(b**2-aux1**2+re*y**2)/(im*y**2)))
     C		+a8*(Integral8(u8,(b**2-aux1**2+re*y**2)/(im*y**2))
     C	    	    -Integral8(l8,(b**2-aux1**2+re*y**2)/(im*y**2)))
     C		    +dlog(dabs((y**2+b**2-aux1**2)/(b**2-aux1**2)))
     C		    *(dlog(y**2)+dlog(im**2)))


	else
	part2=1.d0/y**2*
     C 		(a1*Integral1a((b**2-aux1**2+re*y**2)/(im*y**2),l1,u1) 
     C		+a9*Integral1a((b**2-aux1**2+re*y**2)/(im*y**2),l9,u9)
     C		+a2*Integral2a((b**2-aux1**2+re*y**2)/(im*y**2),l2,u2)
     C		+a3*Integral3a((b**2-aux1**2+re*y**2)/(im*y**2),l3,u3)
     C		+a4*Integral4a((b**2-aux1**2+re*y**2)/(im*y**2),l4,u4)
     C		+a5*Integral5a((b**2-aux1**2+re*y**2)/(im*y**2),l5,u5)
     C		+a6*Integral6a((b**2-aux1**2+re*y**2)/(im*y**2),l6,u6)
     C		+a7*Integral7a((b**2-aux1**2+re*y**2)/(im*y**2),l7,u7)
     C		+a8*Integral8a((b**2-aux1**2+re*y**2)/(im*y**2),l8,u8)
     C		    +dlog(dabs((y**2+b**2-aux1**2)/(b**2-aux1**2)))
     C		    *(dlog(y**2)+dlog(im**2)))
	endif



	if(dabs(((b**2-aux1**2+re*y**2)/(im*y**2))/l).lt.0.2d0
     C	.and.dabs(((b**2-aux1**2+re*y**2)/(im*y**2))/u).lt.0.2d0)
     C	then
	if((l.lt.0.d0.and.u.lt.0.d0).or.(l.gt.0.d0.and.u.gt.0.d0))
     C  then

	part2=1.d0/y**2*
     C		(Integral((b**2-aux1**2+re*y**2)/(im*y**2),l,u)
     C		    +dlog(dabs((y**2+b**2-aux1**2)/(b**2-aux1**2)))
     C		    *(dlog(y**2)+dlog(im**2)))
	endif
	endif

	endif
	CC0=part1+part2




	endif



	return
	end


C*********************************************************************
C	solutions of important Integrals for CC0 and in the various 
C	regions for the approximate Log[1+x^2]
C*********************************************************************
C Integral for real roots, dlog(y-g)/(e*y+f)	

	DOUBLE PRECISION function Integral0(g,e,f)
	implicit none
	DOUBLE PRECISION g,e,f,Li_2

	Integral0=((-dlog(dabs(-g)))*dlog(dabs(f/(g*e+f)))
     C	+dlog(dabs(1.d0-g))*dlog(dabs((e+f)/(g*e+f)))
     C	+Li_2(((g-1.d0)*e)/(g*e+f))-Li_2((g*e)/(g*e+f)))/e
	return
	end

C*********************************************************************
C Integral1 from -/+inf til Abs(3.7)

	DOUBLE PRECISION function Integral1(x,g)
	implicit none
	DOUBLE PRECISION g,x,Li_2
	external Li_2

	Integral1=-(dlog(dabs(x))/g**2)+dlog(dabs(g+x))/g**2+2.d0
     C	*(dlog(dabs(x))*dlog(dabs(x/g+1.d0))+Li_2(-(x/g)))
     C		-1.d0/(g*x)

	Integral1=(dlog(dabs(x))*(1.d0-2.d0*g**2))/(2.d0*g**4)
     C	+(1.d0-2.d0*g**2)/(2.d0*g**3*x)+((2.d0*g**2-1.d0)
     C	*dlog(dabs(g+x)))/(2.d0*g**4)+2.d0*(dlog(dabs(x))
     C	*dlog(dabs(x/g+1.d0))+Li_2(-(x/g)))-1.d0/(4.d0*g**2*x**2)
     C	+1.d0/(6.d0*g*x**3);

	return
	end

C*********************************************************************
C Integral2 from -3.7 til -1.53

	DOUBLE PRECISION function Integral2(x,g) 
	implicit none
	DOUBLE PRECISION g,x

	Integral2=(1.d0/1631461442.d0)*((149375.d0*x**4)/4.d0
     C	-(125.d0/3.d0)*(1195.d0*g+175104.d0)*x**3.d0+(25.d0/2.d0)
     C	*(5975.d0*g**2.d0+875520.d0*g-13307758.d0)*x**2-5.d0
     C	*(29875.d0*g**3+4377600.d0*g**2-66538790.d0*g+473776128.d0)
     C	*x-(-149375.d0*g**4-21888000.d0*g**3+332693950.d0*g**2
     C	-2368880640.d0*g-1631461442.d0*dlog(dabs(169.d0/25.d0))
     C	+4076532000.d0)*dlog(dabs(g+x)))
	return
	end

C*********************************************************************
C Integral3 from -1.53 til -0.68

	DOUBLE PRECISION function Integral3(x,g)
	implicit none
	DOUBLE PRECISION g,x

	Integral3=(703919.d0*(20.d0*x+21.d0)**4)
     C	/4001971303688.d0+((32655231.d0-14078380.d0*g)*(20.d0*x
     C	+21.d0)**3)/3001478477766.d0+((281567600.d0*g**2
     C	-948750600.d0*g+627762809.d0)*(20.d0*x+21.d0)**2)
     C	/2000985651844.d0-(5.d0*(1126270400.d0*g**3-4977586320.d0*g**2
     C	+6495803756.d0*g+7356427995.d0)*(20.d0*x+21.d0))
     C	/1000492825922.d0+((112627040000.d0*g**4-616017024000.d0*g**3
     C	+1172226939200.d0*g**2+53583405120.d0*g+1000492825922.d0
     C	*dlog(dabs(841.d0/400.d0))-772424939475.d0)
     C	*dlog(dabs(20.d0*g+20.d0*x)))/1000492825922.d0
	return
	end

C*********************************************************************
C Integral4 from -0.68 til -0.23

	DOUBLE PRECISION function Integral4(x,g)
	implicit none
	DOUBLE PRECISION g,x

	Integral4=(1.d0/4243686.d0)*(-((76875.d0*x**4)/4.d0)
     C	+(125.d0/3.d0)*(615.d0*g+15488.d0)*x**3-(25.d0/2.d0)*(3075.d0
     C	*g**2+77440.d0*g-201846.d0)*x**2+5.d0*(15375.d0*g**3+387200.d0
     C	*g**2-1009230.d0*g+32256.d0)*x+(-76875.d0*g**4-1936000.d0*g**3
     C	+5046150.d0*g**2-161280.d0*g+4243686.d0
     C	*dlog(dabs(29.d0/25.d0))
     C	-617000.d0)*dlog(dabs(g+x)))
	return
	end

C*********************************************************************
C Integral5 from -0.23 til 0.23

	DOUBLE PRECISION function Integral5(x,g)
	implicit none
	DOUBLE PRECISION g,x

	Integral5=-(x**4/8.d0)+(g*x**3)/6.d0-(1.d0/4.d0)
     C	*(g**2-2.d0)*x**2+(1.d0/2.d0)*g*(g**2-2.d0)*x-(1.d0/2.d0)
     C	*(g**4-2.d0*g**2)*dlog(dabs(g+x))
	return
	end

C*********************************************************************
C Integral6 from 0.23 til 0.68

	DOUBLE PRECISION function Integral6(x,g)
	implicit none
	DOUBLE PRECISION g,x

	Integral6=(1.d0/4243686.d0)*(-((76875.d0*x**4)/4.d0)
     C	+(125.d0/3.d0)*(615.d0*g-15488.d0)*x**3-(25.d0/2.d0)
     C	*(3075.d0*g**2-77440.d0*g-201846.d0)*x**2+5.d0*(15375.d0*g**3
     C	-387200.d0*g**2-1009230.d0*g-32256.d0)*x-(76875.d0*g**4
     C	-1936000.d0*g**3-5046150.d0*g**2-161280.d0*g-4243686.d0
     C	*dlog(dabs(29.d0/25.d0))+617000.d0)*dlog(dabs(g+x)))
	return
	end

C*********************************************************************
C Integral7 from 0.68 till 1.53
	
	DOUBLE PRECISION function Integral7(x,g)
	implicit none
	DOUBLE PRECISION g,x

	Integral7=(1.d0/1000492825922.d0)*(28156760000.d0*x**4
     C 	-(32000.d0/3.d0)*(3519595.d0*g+19250532.d0)*x**3+800.d0
     C	*(70391900.d0*g**2+385010640.d0*g+732641837.d0)*x**2-320.d0
     C 	*(351959500.d0*g**3+1925053200.d0*g**2+3663209185.d0*g
     C	-167448141.d0)*x+(112627040000.d0*g**4+616017024000.d0*g**3
     C	+1172226939200.d0*g**2-53583405120.d0*g+1000492825922.d0
     C	*dlog(dabs(841.d0/400.d0))-772424939475.d0)*dlog(dabs(g+x)))
	return
	end

C*********************************************************************
C Integral8 from 1.53 til 3.7
	
	DOUBLE PRECISION function Integral8(x,g) 
	implicit none
	DOUBLE PRECISION g,x

	Integral8=(1.d0/1631461442.d0)*((149375.d0*x**4)/4.d0
     C	-(125.d0/3.d0)*(1195.d0*g-175104.d0)*x**3+(25.d0/2.d0)
     C	*(5975.d0*g**2-875520.d0*g-13307758.d0)*x**2-5.d0*(29875.d0
     C	*g**3-4377600.d0*g**2-66538790.d0*g-473776128.d0)*x+(149375.d0
     C	*g**4-21888000.d0*g**3-332693950.d0*g**2-2368880640.d0*g
     C	+1631461442.d0*dlog(dabs(169.d0/25.d0))-4076532000.d0)
     C	*dlog(dabs(g+x)))


	Integral8=(318600187109375.d0*dlog(dabs(-5.d0*g-5.d0*x))
     C	*g**8)
     C	/2661666436732719364.d0-(318600187109375.d0*x*g**7)
     C	/2661666436732719364.d0+(12244912560000000.d0
     C	*dlog(dabs(-5.d0*g-5.d0*x))*g**7)/4657916264282258887.d0
     C	+(191160112265625.d0*g**7)/665416609183179841.d0
     C	+(318600187109375.d0*x**2*g**6)/5323332873465438728.d0
     C	-(12244912560000000.d0*x*g**6)/4657916264282258887.d0
     C	+(51167625589890625.d0*dlog(dabs(-5.d0*g-5.d0*x))*g**6)
     C	/1996249827549539523.d0+(27782045200968750.d0*g**6)
     C	/4657916264282258887.d0-(318600187109375.d0*x**3*g**5)
     C	/7984999310198158092.d0+(6122456280000000.d0*x**2*g**5)
     C	/4657916264282258887.d0-(51167625589890625.d0*x*g**5)
     C	/1996249827549539523.d0+(96552009123840000.d0
     C	*dlog(dabs(-5.d0*g-5.d0*x))*g**5)/665416609183179841.d0
     C	+(253842547039437500.d0*g**5)/4657916264282258887.d0
     C	+(318600187109375.d0*x**4*g**4)/10646665746930877456.d0
     C	-(4081637520000000.d0*x**3*g**4)/4657916264282258887.d0
     C	+(51167625589890625.d0*x**2*g**4)/3992499655099079046.d0
     C	-(96552009123840000.d0*x*g**4)/665416609183179841.d0
     C	+(692887912190909375.d0*dlog(dabs(-5.d0*g-5.d0*x))*g**4)
     C	/1330833218366359682.d0+(102309793919769000.d0*g**4)
     C	/358301251098635299.d0-(63720037421875.d0*x**5*g**3)
     C	/2661666436732719364.d0+(3061228140000000.d0*x**4*g**3)
     C	/4657916264282258887.d0-(51167625589890625.d0*x**3*g**3)
     C	/5988749482648618569.d0+(48276004561920000.d0*x**2*g**3)
     C	/665416609183179841.d0-(692887912190909375.d0*x*g**3)
     C	/1330833218366359682.d0+(794443426136064000.d0
     C	*dlog(dabs(-5.d0*g-5.d0*x))*g**3)/665416609183179841.d0
     C	+(4331239193308849950.d0*g**3)/4657916264282258887.d0
     C	+(318600187109375.d0*x**6*g**2)/15969998620396316184.d0
     C	-(2448982512000000.d0*x**5*g**2)/4657916264282258887.d0
     C	+(51167625589890625.d0*x**4*g**2)/7984999310198158092.d0
     C	-(32184003041280000.d0*x**3*g**2)/665416609183179841.d0
     C  +(692887912190909375.d0*x**2*g**2)/2661666436732719364.d0
     C	-(794443426136064000.d0*x*g**2)/665416609183179841.d0
     C	+(1035683341174645825.d0*dlog(dabs(-5.d0*g-5.d0*x))*g**2)
     C	/665416609183179841.d0+(8663688266665928220.d0*g**2)
     C	/4657916264282258887.d0-(318600187109375.d0*x**7*g)
     C	/18631665057129035548.d0+(2040818760000000.d0*x**6*g)
     C	/4657916264282258887.d0-(10233525117978125.d0*x**5*g)
     C	/1996249827549539523.d0+(24138002280960000.d0*x**4*g)
     C	/665416609183179841.d0-(692887912190909375.d0*x**3*g)
     C	/3992499655099079046.d0+(397221713068032000.d0*x**2*g)
     C	/665416609183179841.d0-(1035683341174645825.d0*x*g)
     C	/665416609183179841.d0+(11739945039298560.d0
     C	*dlog(dabs(-5.d0*g-5.d0*x))*g)/665416609183179841.d0
     C	+(653874464745081204.d0*g)/358301251098635299.d0
     C	+(318600187109375.d0*x**8)/21293331493861754912.d0
     C	-(12244912560000000.d0*x**7)/32605413849975812209.d0
     C	+(51167625589890625.d0*x**6)/11977498965297237138.d0
     C	-(19310401824768000.d0*x**5)/665416609183179841.d0
     C	+(692887912190909375.d0*x**4)/5323332873465438728.d0
     C	-(264814475378688000.d0*x**3)/665416609183179841.d0
     C	+(1035683341174645825.d0*x**2)/1330833218366359682.d0
     C	-(11739945039298560.d0*x)/665416609183179841.d0
     C	+dlog(169.d0/25.d0)*dlog(dabs(-5.d0*g-5.d0*x))
     C	-(45705946425790759344.d0*dlog(dabs(-5.d0*g-5.d0*x)))
     C	/23289581321411294435.d0-1309905000276164682408.d0
     C	/815135346249395305225.d0

     

	return 
	end

C********************************************************************
C	Integral for small g in CC0(m1,m2,m3) for susy-loops

	DOUBLE PRECISION function Integral(g,a,b) 
	implicit none
	DOUBLE PRECISION g,a,b,Li_2

	Integral=(1.d0/12.d0)*(-((3.d0*(b**2*a**4-b**4
     C		*dlog(dabs(1.d0+1.d0/b**2))*a**4
     C		+dlog(dabs(1.d0+1.d0/b**2))*a**4+2.d0*dlog(dabs(b))
     C		*a**4-b**4*a**2+(a**4-1.d0)*b**4
     C		*dlog(dabs(1.d0+1.d0/a**2))-2.d0*b**4
     C		*dlog(dabs(a)))*g**4)/(a**4*b**4))
     C		-4.d0*(-2.d0*datan(1.d0/a)+2.d0*datan(1.d0/b)
     C		+dlog(dabs(1.d0+1.d0/a**2))/a**3+(2.d0
     C		*dlog(dabs(a)))/a**3-dlog(dabs(1.d0+1.d0/b**2))/b**3
     C		-(2.d0*dlog(dabs(b)))/b**3+2.d0/a-2.d0/b)*g**3
     C		+(6.d0*((-((b**2+1.d0)*dlog(dabs(1.d0+1.d0/b**2))
     C		+2.d0*dlog(dabs(b))))*a**2+(a**2+1.d0)*b**2
     C		*dlog(dabs(1.d0+1.d0/a**2))+2.d0*b**2
     C		*dlog(dabs(a)))*g**2)/(a**2*b**2)-12.d0*
     C		(2.d0*datan(1.d0/a)-2.d0*datan(1.d0/b)
     C		+dlog(dabs(a**2+1.d0))/a
     C		-dlog(dabs(b**2+1.d0))/b)*g
     C		+6.d0*(Li_2(-a**2)-Li_2(-b**2)))

	return
	end


C********************************************************************
C Dilogarithm

        DOUBLE PRECISION function Li_2(x)
C       Li_2(x)
	implicit none
        DOUBLE PRECISION x,pi,fhilf,bla
        external fhilf
        pi=4.d0*datan(1.d0)
	bla=0.d0
        if(x.ge.-1.d0.and.x.le.0.5d0)then
         bla=fhilf(x)
        else if(x.gt.0.5d0.and.x.lt.1.d0)then
         bla=-fhilf(1.d0-x)+pi**2/6.d0-dlog(x)*dlog(1.d0-x)
        else if(x.lt.-1.d0)then
         bla=-fhilf(1.d0/x)-pi**2/6.d0-.5d0*(dlog(-x))**2
	else if(x.gt.1.d0.and.x.lt.2.d0)then
	 bla=fhilf((x-1)/x)-0.5d0*dlog((x-1.d0)**2/x)*dlog(x)
     C		+pi**2/6.d0
	else if(x.ge.2.d0)then
	 bla=pi**2/3.d0-0.5d0*(dlog(x))**2-fhilf(1.d0/x)
C	else if(x.gt.-500.d0.and.x.lt.-10.d0)then
C	 bla=1.50616d0+0.153765d0*x-0.0000484249d0*x**2
C     C        -2.69934d-8*x**3
C     C        -1.97807d0*dlog(dabs(x))-0.0245271d0*x*dlog(x)
        endif
	Li_2=bla
        return
        end

C********************************************************************
C	Integrals with approximation for large g (g>100GeV) in the 
C		part2 of CC0 for imaginary roots
C	(with arguments: upper and lower boundaries and g)
C********************************************************************
C	Integral1 from -inf til -3.7 for large g 

	DOUBLE PRECISION FUNCTION Integral1a(g,xl,xu)

	implicit none
	double precision g,xl,xu,Li_2

	Integral1a=dlog(dabs(xl))/g**2-dlog(dabs(g+xl))/g**2
     C		-dlog(dabs(xu))/g**2+dlog(dabs(g+xu))/g**2
     C		-2.d0*(dlog(dabs(xl))*dlog(dabs(xl/g+1.d0))
     C		+Li_2(-(xl/g)))+2.d0*(dlog(dabs(xu))
     C		*dlog(dabs(xu/g+1.d0))
     C		+Li_2(-(xu/g)))+1.d0/(g*xl)-1.d0/(g*xu)

	return
	end

C********************************************************************
C 	Integral2 from -3.7 til -1.53 around (-12/5) for large g

	DOUBLE PRECISION FUNCTION Integral2a(g,xl,xu)

	implicit none
	double precision g,xl,xu

	Integral2a=(149375.d0*xl**6)/(9788768652.d0*g**2)-
     C		(29875.d0*xl**5)/(1631461442.d0*g)-(2188800.d0*xl**5)
     C		/(815730721.d0*g**2)+(2736000.d0*xl**4)
     C		/(815730721.d0*g)-(166346975.d0*xl**4)
     C		/(3262922884.d0*g**2)+(166346975.d0*xl**3)
     C		/(2447192163.d0*g)-(394813440.d0*xl**3)
     C		/(815730721.d0*g**2)+(592220160.d0*xl**2)
     C		/(815730721.d0*g)-(149375.d0*xu**6)
     C		/(9788768652.d0*g**2)+(29875.d0*xu**5)
     C		/(1631461442.d0*g)+(2188800.d0*xu**5)
     C		/(815730721.d0*g**2)-(2736000.d0*xu**4)
     C		/(815730721.d0*g)+(166346975.d0*xu**4)
     C		/(3262922884.d0*g**2)-(166346975.d0*xu**3)
     C		/(2447192163.d0*g)+(394813440.d0*xu**3)
     C		/(815730721.d0*g**2)-(592220160.d0*xu**2)
     C		/(815730721.d0*g)+(953856.d0*dlog((g+xl)
     C		/(g+xu)))/4826809.d0+dlog(169.d0/25.d0)
     C		*dlog((g+xu)/(g+xl))-(1877064336.d0
     C		*dlog((g+xu)/(g+xl)))/815730721.d0

	return
	end

C********************************************************************
C	Integral3 from -1.53 til -0.68 around (-21/20) for large g

	DOUBLE PRECISION FUNCTION Integral3a(g,xl,xu)

	implicit none
	double precision g,xl,xu

	Integral3a=(28156760000.d0*xl**6)/(1500739238883.d0*g**2)
     C		-(11262704000.d0*xl**5)/(500246412961.d0*g)
     C		+(61601702400.d0*xl**5)/(500246412961.d0*g**2)
     C		-(77002128000.d0*xl**4)/(500246412961.d0*g)
     C		+(146528367400.d0*xl**4)/(500246412961.d0*g**2)
     C		-(586113469600.d0*xl**3)/(1500739238883.d0*g)
     C		-(8930567520.d0*xl**3)/(500246412961.d0*g**2)
     C		+(13395851280.d0*xl**2)/(500246412961.d0*g)
     C		-(28156760000.d0*xu**6)/(1500739238883.d0*g**2)
     C		+(11262704000.d0*xu**5)/(500246412961.d0*g)
     C		-(61601702400.d0*xu**5)/(500246412961.d0*g**2)
     C		+(77002128000.d0*xu**4)/(500246412961.d0*g)
     C		-(146528367400.d0*xu**4)/(500246412961.d0*g**2)
     C		+(586113469600.d0*xu**3)/(1500739238883.d0*g)
     C		+(8930567520.d0*xu**3)/(500246412961.d0*g**2)
     C		-(13395851280.d0*xu**2)/(500246412961.d0*g)
     C		-(98407386.d0*dlog((g+xl)/(g+xu)))/594823321.d0
     C		+dlog(841.d0/400.d0)*dlog((g+xu)/(g+xl))
     C		-(937946162727.d0*dlog((g+xu)/(g+xl)))
     C		/1000492825922.d0

	return
	end

C********************************************************************
C	Integral4 from -0.68 til -0.23 around (-2/5) for large g

	DOUBLE PRECISION FUNCTION Integral4a(g,xl,xu)

	implicit none
	double precision g,xl,xu

	Integral4a=-((25625.d0*xl**6)/(8487372.d0*g**2))
     C		+(5125.d0*xl**5)/(1414562.d0*g)+(193600.d0*xl**5)
     C		/(2121843.d0*g**2)-(242000.d0*xl**4)/(2121843.d0*g)
     C		+(841025.d0*xl**4)/(2829124.d0*g**2)-(841025.d0*xl**3)
     C		/(2121843.d0*g)+(8960.d0*xl**3)/(707281.d0*g**2)
     C		-(13440.d0*xl**2)/(707281.d0*g)+(25625.d0*xu**6)
     C		/(8487372.d0*g**2)-(5125.d0*xu**5)/(1414562.d0*g)
     C		-(193600.d0*xu**5)/(2121843.d0*g**2)+(242000.d0*xu**4)
     C		/(2121843.d0*g)-(841025.d0*xu**4)/(2829124.d0*g**2)
     C		+(841025.d0*xu**3)/(2121843.d0*g)-(8960.d0*xu**3)
     C		/(707281.d0*g**2)+(13440.d0*xu**2)/(707281.d0*g)
     C		-(2272.d0*dlog((g+xl)/(g+xu)))/73167.d0
     C		+dlog(29.d0/25.d0)*dlog((g+xu)/(g+xl))
     C		-(124796.d0*dlog((g+xu)/(g+xl)))/707281.d0

	return
	end

C********************************************************************
C	Integral5 from -0.23 til +0.23 around (0) for large g

	DOUBLE PRECISION FUNCTION Integral5a(g,xl,xu)

	implicit none
	double precision g,xl,xu

	Integral5a=-(xl**6/(12.d0*g**2))+xl**5/(10.d0*g)+xl**4
     C	/(4.d0*g**2)-xl**3/(3.d0*g)+xu**6/(12.d0*g**2)-xu**5
     C	/(10.d0*g)-xu**4/(4.d0*g**2)+xu**3/(3.d0*g)

	return
	end

C********************************************************************
C	Integral6 from 0.23 til 0.68 around (2/5) for large g

	DOUBLE PRECISION FUNCTION Integral6a(g,xl,xu)

	implicit none
	double precision g,xl,xu

	Integral6a=-((25625.d0*xl**6)/(8487372.d0*g**2))
     C		+(5125.d0*xl**5)/(1414562.d0*g)-(193600.d0*xl**5)
     C		/(2121843.d0*g**2)+(242000.d0*xl**4)/(2121843.d0*g)
     C		+(841025.d0*xl**4)/(2829124.d0*g**2)
     C		-(841025.d0*xl**3)
     C		/(2121843.d0*g)-(8960.d0*xl**3)/(707281.d0*g**2)
     C		+(13440.d0*xl**2)/(707281.d0*g)+(25625.d0*xu**6)
     C		/(8487372.d0*g**2)-(5125.d0*xu**5)/(1414562.d0*g)
     C		+(193600.d0*xu**5)/(2121843.d0*g**2)
     C		-(242000.d0*xu**4)
     C		/(2121843.d0*g)-(841025.d0*xu**4)/(2829124.d0*g**2)
     C		+(841025.d0*xu**3)/(2121843.d0*g)+(8960.d0*xu**3)
     C		/(707281.d0*g**2)-(13440.d0*xu**2)/(707281.d0*g)
     C		-(2272.d0*dlog((g+xl)/(g+xu)))/73167.d0
     C		+dlog(29.d0/25.d0)*dlog((g+xu)/(g+xl))
     C		-(124796.d0*dlog((g+xu)/(g+xl)))/707281.d0

	return
	end

C********************************************************************
C	Integral7 from 0.68 til 1.53 around (21/20) for large g

	DOUBLE PRECISION FUNCTION Integral7a(g,xl,xu)

	implicit none
	double precision g,xl,xu
	
	Integral7a=(28156760000.d0*xl**6)/(1500739238883.d0*g**2)
     C		-(11262704000.d0*xl**5)/(500246412961.d0*g)
     C		-(61601702400.d0*xl**5)/(500246412961.d0*g**2)
     C		+(77002128000.d0*xl**4)/(500246412961.d0*g)
     C		+(146528367400.d0*xl**4)/(500246412961.d0*g**2)
     C		-(586113469600.d0*xl**3)/(1500739238883.d0*g)
     C		+(8930567520.d0*xl**3)/(500246412961.d0*g**2)
     C		-(13395851280.d0*xl**2)/(500246412961.d0*g)
     C		-(28156760000.d0*xu**6)/(1500739238883.d0*g**2)
     C		+(11262704000.d0*xu**5)/(500246412961.d0*g)
     C		+(61601702400.d0*xu**5)/(500246412961.d0*g**2)
     C		-(77002128000.d0*xu**4)/(500246412961.d0*g)
     C		-(146528367400.d0*xu**4)/(500246412961.d0*g**2)
     C		+(586113469600.d0*xu**3)/(1500739238883.d0*g)
     C		-(8930567520.d0*xu**3)/(500246412961.d0*g**2)
     C		+(13395851280.d0*xu**2)/(500246412961.d0*g)
     C		-(98407386.d0*dlog((g+xl)/(g+xu)))/594823321.d0
     C	 	+dlog(841.d0/400.d0)*dlog((g+xu)/(g+xl))
     C		-(937946162727.d0*dlog((g+xu)/(g+xl)))
     C		/1000492825922.d0

	return
	end

C********************************************************************
C	Integral8 from 1.53 til 3.7 around (12/5) for large g

	DOUBLE PRECISION FUNCTION Integral8a(g,xl,xu)

	implicit none
	double precision g,xl,xu

	Integral8a=(149375.d0*xl**6)/(9788768652.d0*g**2)
     C		-(29875.d0*xl**5)/(1631461442.d0*g)
     C	 	+(2188800.d0*xl**5)
     C		/(815730721.d0*g**2)-(2736000.d0*xl**4)
     C		/(815730721.d0*g)-(166346975.d0*xl**4)
     C		/(3262922884.d0*g**2)+(166346975.d0*xl**3)
     C		/(2447192163.d0*g)+(394813440.d0*xl**3)
     C		/(815730721.d0*g**2)-(592220160.d0*xl**2)
     C		/(815730721.d0*g)-(149375.d0*xu**6)
     C		/(9788768652.d0*g**2)+(29875.d0*xu**5)
     C		/(1631461442.d0*g)-(2188800.d0*xu**5)
     C		/(815730721.d0*g**2)+(2736000.d0*xu**4)
     C		/(815730721.d0*g)+(166346975.d0*xu**4)
     C		/(3262922884.d0*g**2)-(166346975.d0*xu**3)
     C		/(2447192163.d0*g)-(394813440.d0*xu**3)
     C		/(815730721.d0*g**2)+(592220160.d0*xu**2)
     C 	  	/(815730721.d0*g)+(953856.d0*dlog((g+xl)/(g+xu)))
     C		/4826809.d0+dlog(169.d0/25.d0)*dlog((g+xu)/(g+xl))
     C		-(1877064336.d0*dlog((g+xu)/(g+xl)))/815730721.d0

	return
	end

C********************************************************************

        DOUBLE PRECISION function fhilf(x)
	
	implicit none
        INTEGER i
	DOUBLE PRECISION b(12),x,z,cCc,sum
        z=-dlog(1.d0-x)
        b(1)=-.5d0
        b(2)=1.d0/6.d0
        b(3)=0.d0
        b(4)=-1.d0/30.d0
        b(5)=0.d0
        b(6)=1.d0/42.d0
        b(7)=0.d0
        b(8)=-1.d0/30.d0
        b(9)=0.d0
        b(10)=5.d0/66.d0
        b(11)=0.d0
        b(12)=-691.d0/2730.d0
        cCc=z
        sum=z   
        do i=1,12
        cCc=cCc*z/(i+1.d0)
        sum=sum+b(i)*cCc
        enddo
        fhilf=sum
        end

C*********************************************************************
C	from: arXive: 0810.0751
C	only valid: if a==b!=0

	double precision function DD0(a,b,c,d)

	implicit none
	double precision a,b,c,d,aux1,aux2

	aux1=(c**2-d**2)/d**2
	aux2=(a**2-c**2)/d**2
	
	If(a.ne.b)then
	DD0=0.d0 !dlog(0.d0)
	else
	DD0=1.d0/d**4*((-(1.d0+aux1)*dlog(1.d0+aux1))/(aux1*aux2**2)
     C		+(-aux2*(aux1+aux2)+((aux1+aux2)*(1.d0+aux1+aux2)+aux2)
     C		*dlog(1.d0+aux1+aux2))/(aux2**2*(aux1+aux2)**2))
	endif

	return
	end
C*********************************************************************
C	from: arXive: 0810.0751
C	only valid: if a==b

	double precision function D27(a,b,c,d)

	implicit none
	double precision a,b,c,d,aux1,aux2

	aux1=(c**2-d**2)/d**2
	aux2=(a**2-c**2)/d**2
	
	If(a.ne.b)then
	D27=0.d0 !dlog(0.d0)
	elseif(b.eq.0.d0)then
	D27=-1.d0/(2.d0*d**2)*((1.d0+aux1)**2*dlog(1.d0+aux1))/(2.d0
     C	    *aux1*aux2**2)
	else
	D27=-1.d0/(2.d0*d**2)*(((1.d0+aux1)**2*dlog(1.d0+aux1))/(2.d0
     C	    *aux1*aux2**2)-((1.d0+aux1+aux2)*(-aux2*(aux1+aux2)+((aux1
     C	    +aux2)*(1.d0+aux1)+aux2)*dlog(1.d0+aux1+aux2)))/(2.d0
     C	    *aux2**2*(aux1+aux2)**2))
	endif

	return
	end

C*********************************************************************
C*********************************************************************
C	->derivatives of some two-point-functions (here indefinit 
C		expressions can occur in L/beta it is needed a case 
C		discrimination) 
C*********************************************************************
C*********************************************************************
C		->B22prime
	DOUBLE PRECISION function B22pr(x,a,b)

	implicit none
	DOUBLE PRECISION x,a,b,B5pr
	
	B22pr=-1.d0/4.d0*B5pr(x,a,b)

	return
	end
C*********************************************************************
C		->B21prime

	DOUBLE PRECISION function B21pr(x,a,b)

	implicit none
	DOUBLE PRECISION x,a,b,aux
	DOUBLE PRECISION sigma,delta,beta,L,Lpr,F0pr,F3pr,FApr

	if(x.eq.0.d0)then
	
	if(a**2-b**2.le.a**2*5.d-6)then
	aux=1.d0/(20.d0*a**2)
	else
		if(a.eq.0.d0)then
		aux=1.d0/(12.d0*b**2)
		elseif(b.eq.0.d0)then
		aux=1.d0/(4.d0*a**2)
		else
	   aux=-(-3.d0*a**8-10.d0*b**2*a**6
     C		-24.d0*b**2*dlog(dabs(b/a))*a**6+18.d0*b**4*a**4
     C		-6.d0*b**6*a**2+b**8)/(12.d0*(a**2-b**2)**5)

		endif
	endif

	else
C		definition of sigma and delta	
	sigma=(a**2+b**2)/x**2
	delta=(a**2-b**2)/x**2

	If(dabs(a**2-b**2).le.a**2*7.d-6)Then
	delta=0.d0
	ENDIF
C		various cases for beta
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.(x**2)
     C	         .lt.((dabs(a)-dabs(b))**2)) then
	beta=dsqrt(1.d0-2.d0*sigma+delta**2)
	else
	beta=dsqrt(2.d0*sigma-delta**2-1.d0)
	endif

C		various cases for L 
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.(x**2)
     C	 	 .lt.((dabs(a)-dabs(b))**2)) then
	L=1.d0/2.d0*dlog(dabs((1.d0+beta-sigma)/(1.d0-beta-sigma)))
	Lpr=1.d0/2.d0*dlog(dabs((1.d0+beta-sigma)/(1.d0-beta-sigma)))
	else
	L=datan((1.d0-delta)/dabs(beta))
     C	+datan((1.d0+delta)/dabs(beta))
	Lpr=-datan((1.d0-delta)/dabs(beta))
     C	-datan((1.d0+delta)/dabs(beta))
	endif
C********************************************************************
	
	if(beta.eq.0.d0)then
	F0pr=1.d0/x**2*(1.d0+delta*dlog(dabs(b/a))-(delta**2-sigma)
     C		*2.d0/(1.d0-delta**2))
	F3pr=1.d0/x**2*(1.d0/6.d0-delta**2+sigma/2.d0+delta
     C		*(sigma-delta**2)*dlog(dabs(b/a))
     C 		+((sigma**2+delta**2)/2.d0+delta**2*(delta**2
     C		-2.d0*sigma))*2.d0/(1.d0-delta**2))
	FApr=1.d0/x**2*((sigma-2.d0*delta**2)*dlog(dabs(b/a))
     C		-2.d0*delta
     C		+delta*(1.d0-3.d0*sigma+2.d0*delta**2)*2.d0
     C		/(1.d0-delta**2))


	else
	F0pr=1.d0/x**2*(1.d0+delta*dlog(dabs(b/a))-(delta**2-sigma)
     C		*Lpr/beta)
	F3pr=1.d0/x**2*(1.d0/6.d0-delta**2+sigma/2.d0+delta
     C		*(sigma-delta**2)*dlog(dabs(b/a))
     C		+((sigma**2+delta**2)
     C		/2.d0+delta**2*(delta**2-2.d0*sigma))*Lpr/beta)
	FApr=1.d0/x**2*((sigma-2.d0*delta**2)*dlog(dabs(b/a))
     C		-2.d0*delta+delta*(1.d0-3.d0*sigma+2.d0*delta**2)
     C		*Lpr/beta)

	endif
	aux=-((F0pr-FApr)/2.d0-F3pr)
	endif
	
	B21pr=aux

	return
	end
C*********************************************************************
C		->B5prime
	DOUBLE PRECISION function B5pr(x,a,b)

	implicit none
	DOUBLE PRECISION x,a,b,aux,div
	DOUBLE PRECISION sigma,delta,beta,L,Lpr,F0,F3,F0pr,F3pr,FApr
	div=0.d0
	if(x.eq.0.d0)then
	
	if(a**2-b**2.le.a**2*5.d-6)then
	aux=-2.d0*dlog(dabs(a))/3.d0
	else
		if(a.eq.0.d0)then
		aux=1.d0/18.d0*(-12.d0*dlog(dabs(b))+5.d0)
		elseif(b.eq.0.d0)then
		aux=1.d0/18.d0*(-12.d0*dlog(dabs(a))+5.d0)
		else
	   aux=-(1.d0/(18.d0*(a**2-b**2)**3))*(-5.d0*a**6+27.d0*b**2
     C		*a**4+12.d0*(a**2-3.d0*b**2)*dlog(dabs(a))*a**4
     C		-27.d0*b**4*a**2+5.d0*b**6-12.d0*b**4*(b**2-3.d0*a**2)
     C		*dlog(dabs(b)))

		endif
	endif

	else
C		definition of sigma and delta	
	sigma=(a**2+b**2)/x**2
	delta=(a**2-b**2)/x**2

	If(dabs(a**2-b**2).le.a**2*7.d-6)Then
	delta=0.d0
	ENDIF
C		various cases for beta
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.(x**2)
     C	         .lt.((dabs(a)-dabs(b))**2)) then
	beta=dsqrt(1.d0-2.d0*sigma+delta**2)
	else
	beta=dsqrt(2.d0*sigma-delta**2-1.d0)
	endif

C		various cases for L 
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.(x**2)
     C	 	 .lt.((dabs(a)-dabs(b))**2)) then
	L=1.d0/2.d0*dlog(dabs((1.d0+beta-sigma)/(1.d0-beta-sigma)))
	Lpr=1.d0/2.d0*dlog(dabs((1.d0+beta-sigma)/(1.d0-beta-sigma)))
	else
	L=datan((1.d0-delta)/dabs(beta))
     C	+datan((1.d0+delta)/dabs(beta))
	Lpr=-datan((1.d0-delta)/dabs(beta))
     C	-datan((1.d0+delta)/dabs(beta))
	endif

	F0=dlog(dabs(a*b))-delta*dlog(dabs(b/a))-2.d0+beta*L
	F3=1.d0/6.d0*dlog(dabs(a*b))-(3.d0*sigma-2.d0*delta**2)
     C	*delta*dlog(dabs(b/a))/6.d0-5.d0/18.d0-(sigma-delta**2)
     C	/3.d0+(1.d0+sigma-2.d0*delta**2)*beta*L/6.d0
	
	if(beta.eq.0.d0)then
	F0pr=1.d0/x**2*(1.d0+delta*dlog(dabs(b/a))-(delta**2-sigma)
     C		*2.d0/(1.d0-delta**2))
	F3pr=1.d0/x**2*(1.d0/6.d0-delta**2+sigma/2.d0+delta
     C		*(sigma-delta**2)*dlog(dabs(b/a))+((sigma**2
     C		+delta**2)/2.d0+delta**2*(delta**2-2.d0*sigma))
     C 		*2.d0/(1.d0-delta**2))
	FApr=1.d0/x**2*((sigma-2.d0*delta**2)*dlog(dabs(b/a))-2.d0
     C		*delta+delta*(1.d0-3.d0*sigma+2.d0*delta**2)*2.d0
     C		/(1.d0-delta**2))


	else
	F0pr=1.d0/x**2*(1.d0+delta*dlog(dabs(b/a))-(delta**2-sigma)
     C		*Lpr/beta)
	F3pr=1.d0/x**2*(1.d0/6.d0-delta**2+sigma/2.d0+delta
     C		*(sigma-delta**2)*dlog(dabs(b/a))
     C		+((sigma**2+delta**2)
     C		/2.d0+delta**2*(delta**2-2.d0*sigma))*Lpr/beta)
	FApr=1.d0/x**2*((sigma-2.d0*delta**2)*dlog(dabs(b/a))
     C		-2.d0*delta+delta*(1.d0-3.d0*sigma+2.d0*delta**2)
     C		*Lpr/beta)

	endif
	aux=-(F0-4.d0*F3+x**2*(F0pr-4.d0*F3pr)+(a**2-b**2)*FApr)
	endif
	
	B5pr=aux+1.d0/3.d0*div

	return
	end

C*********************************************************************
C		B0prime
	DOUBLE PRECISION function B0pr(x,a,b)

	implicit none
	DOUBLE PRECISION x,a,b,aux,sigma,delta,beta,L,F0pr

	if(x.eq.0.d0)then
		if(a**2-b**2.le.a**2*5.d-6)then
		aux=1.d0/(6.d0*a**2)
		else
		if(a.eq.0.d0)then
		aux=1.d0/(2.d0*b**2)
		elseif(b.eq.0.d0)then
		aux=1.d0/(2.d0*a**2)
		else
		aux=(a**4+4.d0*b**2*dlog(dabs(b/a))*a**2-b**4)
     C	     	    /(2.d0*(a**2-b**2)**3)
		endif
		endif
	else
C		definition of sigma and delta	
	sigma=(a**2+b**2)/x**2
	delta=(a**2-b**2)/x**2

	If(dabs(a**2-b**2).le.a**2*7.d-6)Then
	delta=0.d0
	ENDIF
C		various cases for beta
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.(x**2)
     C	         .lt.((dabs(a)-dabs(b))**2)) then
	beta=dsqrt(1.d0-2.d0*sigma+delta**2)
	else
	beta=dsqrt(2.d0*sigma-delta**2-1.d0)
	endif

C		various cases for L 
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.(x**2)
     C	   	 .lt.((dabs(a)-dabs(b))**2)) then
	L=1.d0/2.d0*dlog(dabs((1.d0+beta-sigma)/(1.d0-beta-sigma)))
	else
	L=-datan((1.d0-delta)/dabs(beta))
     C	-datan((1.d0+delta)/dabs(beta))	
	endif

	if(beta.eq.0.d0)then
	F0pr=1.d0/x**2*(1.d0+delta*dlog(dabs(b/a))-(delta**2-sigma)
     C	*(2.d0)/(1.d0-delta**2))
	aux=-F0pr
	else
	F0pr=1.d0/x**2*(1.d0+delta*dlog(dabs(b/a))-(delta**2-sigma)
     C		*L/beta)
	aux=-F0pr
	endif
	endif

	B0pr=aux
	
	return
	end
C*********************************************************************
C		->B1prime (finite)
	DOUBLE PRECISION function B1pr(x,a,b)

	implicit none
	DOUBLE PRECISION x,a,b,aux
	DOUBLE PRECISION sigma,delta,beta,L,Lpr,F0pr,FApr

	if(x.eq.0.d0)then
	
	if(a**2-b**2.le.5.d-6*a**2)then
	aux=-1.d0/(12.d0*a**2)
	else
		if(a.eq.0.d0)then
		aux=-1.d0/(6.d0*b**2)
		elseif(b.eq.0.d0)then
		aux=-1.d0/(3.d0*a**2)
		else
	   aux=-(2.d0*a**6+3.d0*b**2*a**4+12.d0*b**2*a**4
     C		*dlog(dabs(b/a))
     C           -6.d0*b**4*a**2+b**6)/(6.d0*(a**2-b**2)**4)

		endif
	endif

	else
C		definition of sigma and delta	
	sigma=(a**2+b**2)/x**2
	delta=(a**2-b**2)/x**2

	If(dabs(a**2-b**2).le.a**2*7.d-6)Then
	delta=0.d0
	ENDIF
C		various cases for beta
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.(x**2)
     C		 .lt.((dabs(a)-dabs(b))**2)) then
	beta=dsqrt(1.d0-2.d0*sigma+delta**2)
	else
	beta=dsqrt(2.d0*sigma-delta**2-1.d0)
	endif

C		various cases for L 
	if((x**2).gt.((dabs(a)+dabs(b))**2).or.(x**2)
     C		 .lt.((dabs(a)-dabs(b))**2)) then
	L=1.d0/2.d0*dlog(dabs((1.d0+beta-sigma)/(1.d0-beta-sigma)))
	Lpr=1.d0/2.d0*dlog(dabs((1.d0+beta-sigma)/(1.d0-beta-sigma)))
	else
	L=datan((1.d0-delta)/dabs(beta))
     C	+datan((1.d0+delta)/dabs(beta))
	Lpr=-datan((1.d0-delta)/dabs(beta))
     C	-datan((1.d0+delta)/dabs(beta))
	endif

	
	if(beta.eq.0.d0)then
	F0pr=1.d0/x**2*(1.d0+delta*dlog(dabs(b/a))-(delta**2-sigma)
     C		*2.d0/(1.d0-delta**2))
	FApr=1.d0/x**2*((sigma-2.d0*delta**2)*dlog(dabs(b/a))-2.d0
     C		*delta+delta*(1.d0-3.d0*sigma+2.d0*delta**2)*2.d0
     C		/(1.d0-delta**2))


	else
	F0pr=1.d0/x**2*(1.d0+delta*dlog(dabs(b/a))-(delta**2-sigma)
     C		*Lpr/beta)
	FApr=1.d0/x**2*((sigma-2.d0*delta**2)*dlog(dabs(b/a))-2.d0
     C		*delta+delta*(1.d0-3.d0*sigma+2.d0*delta**2)*Lpr/beta)

	endif
	aux=1.d0/2.d0*(F0pr-FApr)
	endif
	
	B1pr=aux

	return
	end

C*********************************************************************
C*********************************************************************
C	new Physics Boson selfenergies in Higgs-sector
C	(np = newphysics = minus sm plus nmssm)
C	Higgsmass (Standard Model)=125.5GeV	
C*********************************************************************
C*********************************************************************
C		->W-Boson NP-contribution in Higgs-sector
	
	DOUBLE PRECISION function SigmaWHiggs(k,tanb,MW)
	
	implicit none
	INTEGER i,j
	DOUBLE PRECISION k,A0zdec,B0zdec,B22,aux1,aux2,pi,Gmu
	DOUBLE PRECISION tanb,sinb,cosb,MZ,MW,C2TW,S2TW,g2
	DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2)
	DOUBLE PRECISION CMASS
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0


	pi=4.d0*datan(1.d0) 

C Weinberg angle
	MZ=MZ0   !91.1876d0
	Gmu=GF0  !1.16637d-5
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW
	g2=4.d0*sqrt(2.d0)*Gmu*MW**2

C       Trig. Functions of Betab 
        sinb=tanb/dsqrt(1.d0+tanb**2)
        cosb=sinb/tanb

C		->summation over CP-even Higgs contribution	

	aux1=0.d0
	do i=1,3
	aux1=aux1+(
     C	        -(SCOMP(i,1)**2+SCOMP(i,2)**2)*A0zdec(SMASS(i))
     C		+4.d0*(sinb*SCOMP(i,2)-cosb
     C		*SCOMP(i,1))**2*B22(k,SMASS(i),CMASS)
     C		+4.d0*(cosb*SCOMP(i,2)+sinb*SCOMP(i,1))**2
     C		*B22(k,SMASS(i),MW)
     C		-4.d0*MW**2*(sinb*SCOMP(i,1)
     C		+cosb*SCOMP(i,2))**2*B0zdec(k,SMASS(i),MW)
     C		)

	enddo

C		->summation over CP-odd Higgs contribution

	aux2=0.d0
	do j=1,2
	aux2=aux2+(
     C		-PCOMP(j,1)**2*A0zdec(PMASS(j))
     C		+4.d0*PCOMP(j,1)**2*B22(k,PMASS(j),CMASS)
     C		)

	enddo 



	SigmaWHiggs=g2/(64.d0*pi**2)*(-A0zdec(125.5d0)+4.d0
     C		*B22(k,MW,125.5d0)-4.d0*MW**2*B0zdec(k,MW,125.5d0))
     C		-g2/(64.d0*pi**2)*(aux1+aux2
     C		-2.d0*A0zdec(CMASS)
     C		)


	return
	end

C*********************************************************************
C		->Z-Boson NP-contribution in Higgs-sector

	DOUBLE PRECISION function SigmaZHiggs(k,tanb,MW)

	implicit none
	INTEGER i,j
	DOUBLE PRECISION k,A0zdec,B0zdec,B22,aux1,aux2,aux3,Gmu
	DOUBLE PRECISION tanb,sinb,cosb,MZ,MW,C2TW,S2TW,g2,pi
	DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2)
	DOUBLE PRECISION CMASS
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0


	pi=4.d0*datan(1.d0) 

C Weinberg angle
	MZ=MZ0  !91.1876d0
	Gmu=GF0  !1.16637d-5
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW
	g2=4.d0*sqrt(2.d0)*Gmu*MW**2

C       Trig. Functions of Betab 
        sinb=tanb/dsqrt(1.d0+tanb**2)
        cosb=sinb/tanb

C		->summation  over CP-even Higgs contribution	

	aux1=0.d0
	do i=1,3
	aux1=aux1+(-(SCOMP(i,1)**2+SCOMP(i,2)**2)*A0zdec(SMASS(i))
     C		+4.d0*(sinb*SCOMP(i,1)+cosb
     C		*SCOMP(i,2))**2*B22(k,SMASS(i),MZ)-4.d0*MW**2/C2TW
     C		*(sinb*SCOMP(i,1)+cosb*SCOMP(i,2))**2
     C		*B0zdec(k,SMASS(i),MZ))

	enddo

C		->summation over CP-odd Higgs contribution

	aux2=0.d0
	do j=1,2
	aux2=aux2+(-PCOMP(j,1)**2*A0zdec(PMASS(j)))
	enddo

C		->sum over mixed CP-odd(-even) Higgs contribution
	
	aux3=0.d0
	do j=1,2
		do i=1,3
		aux3=aux3+(4.d0*PCOMP(j,1)**2*(cosb*SCOMP(i,1)
     C		-sinb*SCOMP(i,2))**2*B22(k,SMASS(i),PMASS(j)))
		enddo
	enddo



	SigmaZHiggs=g2/(16.d0*4.d0*pi**2*C2TW)*(-A0zdec(125.5d0)+4.d0
     C		*B22(k,125.5d0,MZ)-4.d0*MW**2/C2TW*B0zdec(k,125.5d0,MZ))
     C 		-g2/(16.d0*4.d0*pi**2*C2TW)*(-2.d0*(S2TW-C2TW)**2
     C		*A0zdec(CMASS)+4.d0*(S2TW-C2TW)**2*B22(k,CMASS,CMASS)+
     C		aux1+aux2+aux3)


	return
	end

C*********************************************************************
C		->Photon NP-contribution in Higgs-sector

	DOUBLE PRECISION function SigmaGamHiggs(k,MW)

	implicit none
	DOUBLE PRECISION k,B5,pi,g2,S2TW,MZ,MW,Gmu
	DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2)
	DOUBLE PRECISION CMASS
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0


	pi=4.d0*datan(1.d0) 

C Weinberg angle
	Gmu=GF0 !1.16637d-5
	MZ=MZ0  !91.1876d0
	S2TW=1.d0-MW**2/MZ**2
	g2=4.d0*sqrt(2.d0)*Gmu*MW**2

	SigmaGamHiggs=g2*S2TW/(16.d0*pi**2)*B5(k,CMASS,CMASS)


	return
	end

C*********************************************************************
C		->Photon-Z-Mixing NP-contribution in Higgs-sector



	DOUBLE PRECISION function SigmaGamZHiggs(k,MW)

	implicit none
	DOUBLE PRECISION k,MW,B5,pi,MZ,S2TW,C2TW,g2,Gmu
	DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2)
	DOUBLE PRECISION CMASS
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0



	pi=4.d0*datan(1.d0) 

C Weinberg angle
	MZ=MZ0  !91.1876d0
	Gmu=GF0 !1.16637d-5
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW
	g2=4.d0*sqrt(2.d0)*Gmu*MW**2

	SigmaGamZHiggs=g2*dsqrt(S2TW/C2TW)/(32.d0*pi**2)*(S2TW-C2TW)
     C		*B5(k,CMASS,CMASS)


	return
	end

C*********************************************************************
C		->Photon Prime NP-contribution in Higgs-sector

	DOUBLE PRECISION function SigmaGamPrHiggs(k,MW)

	implicit none
	DOUBLE PRECISION k,pi,B5pr,MW,MZ,S2TW,g2,Gmu
	DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2)
	DOUBLE PRECISION CMASS
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0


	pi=4.d0*datan(1.d0) 

C Weinberg angle
	MZ=MZ0  !91.1876d0
	Gmu=GF0 !1.16637d-5
	S2TW=1.d0-MW**2/MZ**2
	g2=4.d0*sqrt(2.d0)*Gmu*MW**2

	SigmaGamPrHiggs=g2*S2TW/(16.d0*pi**2)*B5pr(k,CMASS,CMASS)


	return
	end
C*********************************************************************
C		->Z-Boson-Prime NP-contribution in Higgs-sector

	DOUBLE PRECISION function SigmaZPrHiggs(k,tanb,MW)

	implicit none
	INTEGER i,j
	DOUBLE PRECISION aux1,aux3
	DOUBLE PRECISION k,pi,B22pr,B0pr,Gmu
	DOUBLE PRECISION tanb,sinb,cosb,MW,MZ,S2TW,C2TW,g2
	DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2)
	DOUBLE PRECISION CMASS
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0



	pi=4.d0*datan(1.d0) 

C Weinberg angle
	MZ=MZ0  !91.1876d0
	Gmu=GF0 !1.16637d-5
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW
	g2=4.d0*sqrt(2.d0)*Gmu*MW**2

C       Trig. Functions of Betab 
        sinb=tanb/dsqrt(1.d0+tanb**2)
        cosb=sinb/tanb

C		->summation over pure CP-even Higgs contribution	

	aux1=0.d0
	do i=1,3
	aux1=aux1+(4.d0*(sinb*SCOMP(i,1)+cosb
     C		*SCOMP(i,2))**2*B22pr(k,SMASS(i),MZ)-4.d0*MW**2/C2TW
     C		*(sinb*SCOMP(i,1)+cosb*SCOMP(i,2))**2
     C		*B0pr(k,SMASS(i),MZ))
	enddo


C		->sum over mixed CP-odd(-even) Higgs contribution
	
	aux3=0.d0
	do j=1,2
		do i=1,3
		aux3=aux3+(4.d0*PCOMP(j,1)**2*(cosb*SCOMP(i,1)
     C		-sinb*SCOMP(i,2))**2*B22pr(k,SMASS(i),PMASS(j)))
		enddo
	enddo
	


	SigmaZPrHiggs=(g2)/(16.d0*4.d0*pi**2*C2TW)*(4.d0
     C		*B22pr(k,125.5d0,MZ)-4.d0*MW**2/C2TW*B0pr(k,125.5d0,MZ))
     C 		-g2/(16.d0*4.d0*pi**2*C2TW)*(4.d0*(S2TW-C2TW)**2
     C		*B22pr(k,CMASS,CMASS)
     C          +aux1+aux3)



	return
	end


C*********************************************************************
C*********************************************************************
C	renormalized Boson selfenergies (see for relations e.g. [5])
C*********************************************************************
C*********************************************************************
C		-> W-Boson NP-contribution in Higgs-sector
	DOUBLE PRECISION function SigmaWrenHiggs(k,tanb,MW)

	implicit none
	DOUBLE PRECISION k,SigmaWHiggs,dMZ2,dMW2,dzw2,tanb
	DOUBLE PRECISION SigmaZHiggs,SigmaGamPrHiggs,SigmaGamZHiggs
	DOUBLE PRECISION MZ,MW,S2TW,C2TW
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0

C Weinberg angle
	MZ=MZ0  !91.1876d0
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW

	dMZ2=SigmaZHiggs(MZ,tanb,MW)
	dMW2=SigmaWHiggs(MW,tanb,MW)
	dzw2=-SigmaGamPrHiggs(0.d0,MW)-2.d0*dsqrt(C2TW/S2TW)
     C		*SigmaGamZHiggs(0.d0,MW)/MZ**2+C2TW/S2TW
     C		*(dMZ2/MZ**2-dMW2/MW**2)

	SigmaWrenHiggs=SigmaWHiggs(k,tanb,MW)-dMW2+dzw2*(k**2-MW**2)


	return
	end
C*********************************************************************
C		->Z-Boson NP-contribution in Higgs-sector
	DOUBLE PRECISION function SigmaZrenHiggs(k,tanb,MW)

	implicit none
	DOUBLE PRECISION k,SigmaZHiggs,SigmaWHiggs,SigmaGamZHiggs
	DOUBLE PRECISION SigmaGamPrHiggs,MW,MZ,S2TW,C2TW
	DOUBLE PRECISION dMZ2,dMW2,dzz2,tanb
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0

C Weinberg angle
	MZ=MZ0 !91.1876d0
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW

	dMZ2=SigmaZHiggs(MZ,tanb,MW)
	dMW2=SigmaWHiggs(MW,tanb,MW)
	dzz2=-SigmaGamPrHiggs(0.d0,MW)-2.d0*(C2TW-S2TW)
     C   /dsqrt(S2TW*C2TW)*SigmaGamZHiggs(0.d0,MW)/MZ**2
     C		+(C2TW-S2TW)/S2TW*(dMZ2/MZ**2-dMW2/MW**2)

	SigmaZrenHiggs=SigmaZHiggs(k,tanb,MW)-dMZ2+dzz2*(k**2-MZ**2)

	return
	end
C*********************************************************************
C		->Photon NP-contribution in Higgs-sector
	DOUBLE PRECISION function SigmaGamrenHiggs(k,MW)

	implicit none
	DOUBLE PRECISION k,MW,SigmaGamHiggs,SigmaGamPrHiggs

	SigmaGamrenHiggs=SigmaGamHiggs(k,MW)
     C                    -SigmaGamPrHiggs(0.d0,MW)*k**2

	return
	end
C*********************************************************************
C		->Photon-Z-mixing NP-contribution in Higgs-sector
	DOUBLE PRECISION function SigmaGamZrenHiggs(k,tanb,MW)

	implicit none
	DOUBLE PRECISION k,SigmaGamZHiggs,dzgaz2,dzgaz1,tanb
	DOUBLE PRECISION dMZ2,dMW2,dzz1,dzz2,dzga1,dzga2
	DOUBLE PRECISION SigmaGamPrHiggs,SigmaZHiggs,SigmaWHiggs
	DOUBLE PRECISION MZ,MW,S2TW,C2TW
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	
C Weinberg angle
	MZ=MZ0 !91.1876d0
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW

	dMZ2=SigmaZHiggs(MZ,tanb,MW)
	dMW2=SigmaWHiggs(MW,tanb,MW)

	dzga1=-SigmaGamPrHiggs(0.d0,MW)-dsqrt(S2TW/C2TW)
     C	      *SigmaGamZHiggs(0.d0,MW)/MZ**2
	dzga2=-SigmaGamPrHiggs(0.d0,MW)

	dzz1=-SigmaGamPrHiggs(0.d0,MW)-(3.d0*C2TW-2.d0*S2TW)
     C	        /dsqrt(C2TW*S2TW)*SigmaGamZHiggs(0.d0,MW)/MZ**2
     C	        +(C2TW-S2TW)/S2TW*(dMZ2/MZ**2-dMW2/MW**2)


	dzz2=-SigmaGamPrHiggs(0.d0,MW)-2.d0*(C2TW-S2TW)
     C	/dsqrt(S2TW*C2TW)*SigmaGamZHiggs(0.d0,MW)/MZ**2+(C2TW-S2TW)
     C		/S2TW*(dMZ2/MZ**2-dMW2/MW**2)

	dzgaz1=dsqrt(C2TW*S2TW)/(C2TW-S2TW)*(dzz1-dzga1)

	dzgaz2=dsqrt(C2TW*S2TW)/(C2TW-S2TW)*(dzz2-dzga2)


	SigmaGamZrenHiggs=SigmaGamZHiggs(k,MW)-dzgaz2*k**2
     C          	  +(dzgaz1-dzgaz2)*MZ**2

	return
	end
C*********************************************************************
C*********************************************************************
C	SUSY-contribution to gauge bosons selfenergies 
C	(all SUSY contribution besides the one from the Higgs Sector->
C	already included  e.g. SigmaWHiggs...)
C*********************************************************************
C*********************************************************************
C		->W-Selfenergy SUSY-contribution
	
	DOUBLE PRECISION function SigmaWSusy(k,MW)
	
	implicit none
	INTEGER i,j
	DOUBLE PRECISION k,A0zdec,B1zdec,B0zdec,B22,B21
	DOUBLE PRECISION aux,C2TW,pi,Cf,MW,MZ,S2TW,Gmu
	DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
	DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
	DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	DOUBLE PRECISION MSU(2),MSD(2),MSE(2),MST(2),MSB(2),MSL(2)
	DOUBLE PRECISION UT(2,2),UB(2,2),UL(2,2),DELT(2,2)
	DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
	DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C		MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C		CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
        COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0

	pi=4.d0*datan(1.d0) 
	Gmu=GF0  !1.16637d-5

C Weinberg angle
	MZ=MZ0  !91.1876d0
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW

	Cf=3.d0

	
C	    Additionnal contributions from the SUSY sector

C       0) Masses and couplings

C	everything for 1.+ 2. Generation	
        MSU(1)=MUL
        MSU(2)=MUR
        MSD(1)=MDL
        MSD(2)=MDR
        MSE(1)=MLL
        MSE(2)=MLR
	
C	everything for Stop	
        MST(1)=dsqrt(RMST1)
	MST(2)=dsqrt(RMST2)
        UT(1,1)=CST
        UT(1,2)=+dsqrt(1.d0-CST**2)
        UT(2,1)=-dsqrt(1.d0-CST**2)
        UT(2,2)=CST

C	everything for Sbottom
        MSB(1)=dsqrt(RMSB1)
	MSB(2)=dsqrt(RMSB2)
        UB(1,1)=CSB
        UB(1,2)=+dsqrt(1.d0-CSB**2)
        UB(2,1)=-dsqrt(1.d0-CSB**2)
        UB(2,2)=CSB

C	everthing for Stau
        MSL(1)=MSL1
	MSL(2)=MSL2
        UL(1,1)=CSL
        UL(1,2)=+dsqrt(1.d0-CSL**2)
        UL(2,1)=-dsqrt(1.d0-CSL**2)
        UL(2,2)=CSL

C	Kronecker Delta		
        DELT(1,1)=1.d0
        DELT(1,2)=0.d0
        DELT(2,1)=0.d0
        DELT(2,2)=1.d0

C	1) sfermion contributions
      
        aux=0.d0
        do i=1,2
C	3. Generation
         aux=aux
     C		+Cf*UT(i,1)**2*A0zdec(MST(i))
     C		+Cf*UB(i,1)**2*A0zdec(MSB(i))
     C       +UL(i,1)**2*A0zdec(MSL(i))+DELT(i,1)**2*A0zdec(MSNT)
C	1.+2. Generation -> therefore "2*"
         aux=aux+2.d0*(
     C        Cf*DELT(i,1)**2*A0zdec(MSU(i))
     C	     +Cf*DELT(i,1)**2*A0zdec(MSD(i))
     C       +DELT(i,1)**2*A0zdec(MSE(i))+DELT(i,1)**2*A0zdec(MNL))
C	3. Generation
         do j=1,2
          aux=aux
     C		-4.d0*Cf*(UT(i,1)*UB(j,1))**2*B22(k,MST(i),MSB(j))
     C          -4.d0*(UL(i,1)*DELT(j,1))**2*B22(k,MSNT,MSL(i))
         enddo
C	1.+2. Generation -> therefore "2*"
         do j=1,2
          aux=aux
     C       -8.d0*Cf*(DELT(i,1)*DELT(j,1))*B22(k,MSU(i),MSD(j))
     C       -8.d0*(DELT(i,1)*DELT(j,1))**2*B22(k,MNL,MSE(i))
         enddo
        enddo

	SigmaWSusy=4.d0*dsqrt(2.d0)*Gmu*MW**2/(2.d0*16.d0*Pi**2)*aux

	
C	2) chargino/neutralino contributions

	aux=0.d0
	do i=1,2
	do j=1,5
	aux=aux+
     C      (-2.d0*k**2*(B1zdec(k,MCHA(i),MNEU(j))
     C	                   +B21(k,MCHA(i),MNEU(j)))
     C      +MNEU(j)**2*B1zdec(k,MCHA(i),MNEU(j))
     C	    +MCHA(i)**2*B1zdec(k,MNEU(j),MCHA(i)))
     C	    *((-NEU(j,2)*V(i,1)+NEU(j,3)*V(i,2)/dsqrt(2.d0))**2
     C      +(U(i,1)*NEU(j,2)+U(i,2)*NEU(j,4)/dsqrt(2.d0))**2)
     C	    +2.d0*MNEU(j)*MCHA(i)*B0zdec(k,MCHA(i),MNEU(j))
     C	    *(NEU(j,2)*V(i,1)-NEU(j,3)*V(i,2)/dsqrt(2.d0))
     C	    *(U(i,1)*NEU(j,2)+U(i,2)*NEU(j,4)/dsqrt(2.d0))
	enddo
	enddo  	

	If(k.le.1.d-3)THEN
	aux=0.d0
	do i=1,2
	do j=1,5
	aux=aux+
     C      (MNEU(j)**2*B1zdec(0.d0,MCHA(i),MNEU(j))
     C	    +MCHA(i)**2*B1zdec(0.d0,MNEU(j),MCHA(i)))
     C	    *((-NEU(j,2)*V(i,1)+NEU(j,3)*V(i,2)/dsqrt(2.d0))**2
     C      +(U(i,1)*NEU(j,2)+U(i,2)*NEU(j,4)/dsqrt(2.d0))**2)
     C	    +2.d0*MNEU(j)*MCHA(i)*B0zdec(0.d0,MCHA(i),MNEU(j))
     C	    *(NEU(j,2)*V(i,1)-NEU(j,3)*V(i,2)/dsqrt(2.d0))
     C	    *(U(i,1)*NEU(j,2)+U(i,2)*NEU(j,4)/dsqrt(2.d0))
	enddo
	enddo  	
	ENDIF

       SigmaWSusy=SigmaWSusy
     C         +4.d0*dsqrt(2.d0)*Gmu*MW**2*2.d0/(16.d0*Pi**2)*aux

c	print*,'k,SigW',k,
c     c 4.d0*dsqrt(2.d0)*Gmu*MW**2*2.d0/(16.d0*Pi**2)*aux,SigmaWSusy

	return
	end
C*********************************************************************
C		->Z-Selfenergy SUSY-contribution
	
	DOUBLE PRECISION function SigmaZSusy(k,MW)
	
	implicit none
	INTEGER i,j
	DOUBLE PRECISION k,A0zdec,B0zdec,B1zdec,B22,B21
	DOUBLE PRECISION aux,C2TW,pi,Cf,Gmu,MW,MZ,S2TW
	DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
	DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
	DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	DOUBLE PRECISION MSU(2),MSD(2),MSE(2),MST(2),MSB(2),MSL(2)
	DOUBLE PRECISION UT(2,2),UB(2,2),UL(2,2),DELT(2,2)
	DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
	DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C		MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C		CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
        COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0

	pi=4.d0*datan(1.d0) 
	Gmu=GF0 !1.16637d-5

C Weinberg angle
	MZ=MZ0  !91.1876d0
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW

	Cf=3.d0

	
C	    Additionnal contributions from the SUSY sector

C       0) Masses and couplings

C	everything for 1.+ 2. Generation	
        MSU(1)=MUL
        MSU(2)=MUR
        MSD(1)=MDL
        MSD(2)=MDR
        MSE(1)=MLL
        MSE(2)=MLR
	
C	everything for Stop	
        MST(1)=dsqrt(RMST1)
	MST(2)=dsqrt(RMST2)
        UT(1,1)=CST
        UT(1,2)=+dsqrt(1.d0-CST**2)
        UT(2,1)=-dsqrt(1.d0-CST**2)
        UT(2,2)=CST

C	everything for Sbottom
        MSB(1)=dsqrt(RMSB1)
	MSB(2)=dsqrt(RMSB2)
        UB(1,1)=CSB
        UB(1,2)=+dsqrt(1.d0-CSB**2)
        UB(2,1)=-dsqrt(1.d0-CSB**2)
        UB(2,2)=CSB

C	everthing for Stau
        MSL(1)=MSL1
	MSL(2)=MSL2
        UL(1,1)=CSL
        UL(1,2)=+dsqrt(1.d0-CSL**2)
        UL(2,1)=-dsqrt(1.d0-CSL**2)
        UL(2,2)=CSL

C	Kronecker Delta		
        DELT(1,1)=1.d0
        DELT(1,2)=0.d0
        DELT(2,1)=0.d0
        DELT(2,2)=1.d0


C	1) sfermion contributions
      
        aux=0.d0
        do i=1,2
C	3. Generation	
         aux=aux
     C	    +Cf*((1.d0/2.d0-2.d0*S2TW/3.d0)**2*UT(i,1)**2
     C           +4.d0/9.d0*S2TW**2*UT(i,2)**2)*A0zdec(MST(i))
     C      +Cf*((1.d0/2.d0-S2TW/3.d0)**2*UB(i,1)**2
     C           +1.d0/9.d0*S2TW**2*UB(i,2)**2)*A0zdec(MSB(i))
     C      +((1.d0/2.d0-S2TW)**2*UL(i,1)**2
     C           +S2TW**2*UL(i,2)**2)*A0zdec(MSL(i))
C	1.+ 2. Generation         
	 aux=aux+2.d0*(
     C      Cf*((1.d0/2.d0-2.d0*S2TW/3.d0)**2*DELT(i,1)**2
     C           +4.d0/9.d0*S2TW**2*DELT(i,2)**2)*A0zdec(MSU(i))
     C      +Cf*((1.d0/2.d0-S2TW/3.d0)**2*DELT(i,1)**2
     C           +1.d0/9.d0*S2TW**2*DELT(i,2)**2)*A0zdec(MSD(i))
     C      +((1.d0/2.d0-S2TW)**2*DELT(i,1)**2
     C           +S2TW**2*DELT(i,2)**2)*A0zdec(MSE(i)))
C	3. Generation
         do j=1,2
          aux=aux
     C		 -2.d0*Cf*((1.d0/2.d0-2.d0*S2TW/3.d0)*UT(i,1)*UT(j,1)
     C                  -2.d0/3.d0*S2TW*UT(i,2)*UT(j,2))**2
     C             *B22(k,MST(i),MST(j))
     C           -2.d0*Cf*((1.d0/2.d0-S2TW/3.d0)*UB(i,1)*UB(j,1)
     C                  -1.d0/3.d0*S2TW*UB(i,2)*UB(j,2))**2
     C             *B22(k,MSB(i),MSB(j))
     C           -2.d0*((1.d0/2.d0-S2TW)*UL(i,1)*UL(j,1)
     C                               -S2TW*UL(i,2)*UL(j,2))**2
     C             *B22(k,MSL(i),MSL(j))
         enddo
C	1. + 2. Generation
         do j=1,2
          aux=aux
     C      -4.d0*Cf*((1.d0/2.d0-2.d0*S2TW/3.d0)*DELT(i,1)*DELT(j,1)
     C                  -2.d0/3.d0*S2TW*DELT(i,2)*DELT(j,2))**2
     C             *B22(k,MSU(i),MSU(j))
     C      -4.d0*Cf*((1.d0/2.d0-S2TW/3.d0)*DELT(i,1)*DELT(j,1)
     C               -1.d0/3.d0*S2TW*DELT(i,2)*DELT(j,2))**2
     C             *B22(k,MSD(i),MSD(j))
     C      -4.d0*((1.d0/2.d0-S2TW)*DELT(i,1)*DELT(j,1)
     C                         -S2TW*DELT(i,2)*DELT(j,2))**2
     C             *B22(k,MSE(i),MSE(j))
         enddo
        enddo
        aux=aux+1.d0/4.d0*(A0zdec(MSNT)-2.d0*B22(k,MSNT,MSNT))
     C     +1.d0/2.d0*(A0zdec(MNL)-2.d0*B22(k,MNL,MNL))

	SigmaZSusy=4.d0*dsqrt(2.d0)*Gmu*MW**2
     C                *2.d0/((1.d0-S2TW)*16.d0*Pi**2)*aux



C	2) chargino/neutralino contributions

C	->chargino
        aux=0.d0
        do i=1,2
        do j=1,2
	aux=aux+
     C      (-2.d0*k**2*(B1zdec(k,MCHA(j),MCHA(i))
     C       +B21(k,MCHA(j),MCHA(i)))
     C	     +MCHA(i)**2*B1zdec(k,MCHA(j),MCHA(i))
     C       +MCHA(j)**2*B1zdec(k,MCHA(i),MCHA(j)))
     C	     *((V(j,1)*V(i,1)+1.d0/2.d0*V(j,2)*V(i,2)
     C		-S2TW*delt(i,j))**2
     C	     +(U(i,1)*U(j,1)+1.d0/2.d0*U(i,2)*U(j,2)
     C		-S2TW*delt(i,j))**2)	
     C       +2.d0*MCHA(i)*MCHA(j)*B0zdec(k,MCHA(j),MCHA(i))
     C       *(V(j,1)*V(i,1)+1.d0/2.d0*V(j,2)*V(i,2)-S2TW*delt(i,j))
     C	     *(U(i,1)*U(j,1)+1.d0/2.d0*U(i,2)*U(j,2)-S2TW*delt(i,j))	
        enddo
        enddo

C	->Neutralino


C	not sure about 1/2 from Majorana-Loop
        do i=1,5
        do j=1,5
	aux=aux+1.d0/4.d0*(
     C  (-2.d0*k**2*(B1zdec(k,MNEU(j),MNEU(i))+B21(k,MNEU(j),MNEU(i)))
     C   +MNEU(i)**2*B1zdec(k,MNEU(j),MNEU(i))
     C   +MNEU(j)**2*B1zdec(k,MNEU(i),MNEU(j)))
     C  *(NEU(j,4)*NEU(i,4)-NEU(j,3)*NEU(i,3))**2
     C  -MNEU(i)*MNEU(j)*B0zdec(k,MNEU(j),MNEU(i))
     C   *(NEU(j,4)*NEU(i,4)-NEU(j,3)*NEU(i,3))**2)
        enddo
        enddo


        SigmaZSusy=SigmaZSusy
     C		+4.d0*dsqrt(2.d0)*Gmu*MW**2*2.d0/(16.d0*Pi**2
     C		*(1.d0-S2TW))*aux	
c      print*,'k,SigZ',k,4.d0*dsqrt(2.d0)*Gmu*MW**2*2.d0/(16.d0*Pi**2
c     C		*(1.d0-S2TW))*aux,SigmaZSusy

	return
	end
C*********************************************************************
C		->Gamma-Selfenergy SUSY-contribution
	
	DOUBLE PRECISION function SigmaGamSusy(k,MW)
	
	implicit none
	INTEGER i
	DOUBLE PRECISION k,B0zdec,B1zdec,B21,B5,aux,pi,Cf
	DOUBLE PRECISION MW,MZ,Gmu,C2TW,S2TW
	DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
	DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
	DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	DOUBLE PRECISION MSU(2),MSD(2),MSE(2),MST(2),MSB(2),MSL(2)
	DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
	DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C		MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C		CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
        COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0

	pi=4.d0*datan(1.d0) 
	Gmu=GF0 !1.16637d-5

C Weinberg angle
	MZ=MZ0  !91.1876d0
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW

	Cf=3.d0

	
C	    Additionnal contributions from the SUSY sector

C       0) Masses and couplings

C	everything for 1.+ 2. Generation	
        MSU(1)=MUL
        MSU(2)=MUR
        MSD(1)=MDL
        MSD(2)=MDR
        MSE(1)=MLL
        MSE(2)=MLR
	
C	everything for Stop	
        MST(1)=dsqrt(RMST1)
	MST(2)=dsqrt(RMST2)

C	everything for Sbottom
        MSB(1)=dsqrt(RMSB1)
	MSB(2)=dsqrt(RMSB2)

C	everthing for Stau
        MSL(1)=MSL1
	MSL(2)=MSL2
 

C	1) sfermion contributions
      
        aux=0.d0
        do i=1,2
C	3. Generation	
         aux=aux
     C	    +Cf*4.d0/9.d0*B5(k,MST(i),MST(i))
     C      +Cf*1.d0/9.d0*B5(k,MSB(i),MSB(i))
     C      +B5(k,MSL(i),MSL(i))
C	1. + 2. Generation
         aux=aux
     C     		+2.d0*(
     C	     Cf*4.d0/9.d0*B5(k,MSU(i),MSU(i))
     C      +Cf*1.d0/9.d0*B5(k,MSD(i),MSD(i))
     C      +B5(k,MSE(i),MSE(i)))
        enddo

	SigmaGamSusy=4.d0*dsqrt(2.d0)*Gmu*MW**2*S2TW/(16.d0*Pi**2)*aux


C	2) chargino/neutralino contributions

	aux=0.d0
        do i=1,2
	aux=aux+(-2.d0*k**2*(B1zdec(k,MCHA(i),MCHA(i))
     C      +B21(k,MCHA(i),MCHA(i)))
     C      +2.d0*MCHA(i)**2*B1zdec(k,MCHA(i),MCHA(i))
     C      +MCHA(i)**2*B0zdec(k,MCHA(i),MCHA(i)))
        enddo

        SigmaGamSusy=SigmaGamSusy
     C        +4.d0*4.d0*dsqrt(2.d0)*Gmu*MW**2*S2TW/(16.d0*Pi**2)*aux
	

	return
	end
C*********************************************************************
C		->Gamma-Z-Selfenergy SUSY-contribution
	
	DOUBLE PRECISION function SigmaGamZSusy(k,MW)
	
	implicit none
	INTEGER i
	DOUBLE PRECISION k,B0zdec,B1zdec,B21,B5,aux,pi,Cf
	DOUBLE PRECISION MZ,MW,C2TW,S2TW,Gmu
	DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
	DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
	DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	DOUBLE PRECISION MSU(2),MSD(2),MSE(2),MST(2),MSB(2),MSL(2)
	DOUBLE PRECISION UT(2,2),UB(2,2),UL(2,2),DELT(2,2)
	DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
	DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C		MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C		CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
        COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0

	pi=4.d0*datan(1.d0) 
	Gmu=GF0 !1.16637d-5

C Weinberg angle
	MZ=MZ0  !91.1876d0
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW

	Cf=3.d0

	
C	    Additionnal contributions from the SUSY sector

C       0) Masses and couplings

C	everything for 1.+ 2. Generation	
        MSU(1)=MUL
        MSU(2)=MUR
        MSD(1)=MDL
        MSD(2)=MDR
        MSE(1)=MLL
        MSE(2)=MLR
	
C	everything for Stop	
        MST(1)=dsqrt(RMST1)
	MST(2)=dsqrt(RMST2)
        UT(1,1)=CST
        UT(1,2)=+dsqrt(1.d0-CST**2)
        UT(2,1)=-dsqrt(1.d0-CST**2)
        UT(2,2)=CST

C	everything for Sbottom
        MSB(1)=dsqrt(RMSB1)
	MSB(2)=dsqrt(RMSB2)
        UB(1,1)=CSB
        UB(1,2)=+dsqrt(1.d0-CSB**2)
        UB(2,1)=-dsqrt(1.d0-CSB**2)
        UB(2,2)=CSB

C	everthing for Stau
        MSL(1)=MSL1
	MSL(2)=MSL2
        UL(1,1)=CSL
        UL(1,2)=+dsqrt(1.d0-CSL**2)
        UL(2,1)=-dsqrt(1.d0-CSL**2)
        UL(2,2)=CSL

C	Kroenecker-Delta		
        DELT(1,1)=1.d0
        DELT(1,2)=0.d0
        DELT(2,1)=0.d0
        DELT(2,2)=1.d0


C	1) sfermion contributions

        aux=0.d0
        do i=1,2
C	3. Generation	
         aux=aux
     C	    -2.d0/3.d0*Cf*((1-4.d0/3.d0*S2TW)/2.d0*UT(i,1)**2
     C                          -2.d0/3.d0*S2TW*UT(i,2)**2)
     C                        *B5(k,MST(i),MST(i))
     C      -1.d0/3.d0*Cf*((1-2.d0/3.d0*S2TW)/2.d0*UB(i,1)**2
     C                          -1.d0/3.d0*S2TW*UB(i,2)**2)
     C                        *B5(k,MSB(i),MSB(i))
     C      -((1-2.d0*S2TW)/2.d0*UL(i,1)**2-S2TW*UL(i,2)**2)
     C                        *B5(k,MSL(i),MSL(i))
C	1. + 2. Generation
         aux=aux
     C      -4.d0/3.d0*Cf*((1-4.d0/3.d0*S2TW)/2.d0*DELT(i,1)**2
     C                         -2.d0/3.d0*S2TW*DELT(i,2)**2)
     C                        *B5(k,MSU(i),MSU(i))
     C      -2.d0/3.d0*Cf*((1-2.d0/3.d0*S2TW)/2.d0*DELT(i,1)**2
     C                         -1.d0/3.d0*S2TW*DELT(i,2)**2)
     C                        *B5(k,MSD(i),MSD(i))
     C      -2.d0*((1-2.d0*S2TW)/2.d0*DELT(i,1)**2
     C               -S2TW*DELT(i,2)**2)*B5(k,MSE(i),MSE(i))
        enddo

	SigmaGamZSusy=4.d0*dsqrt(2.d0)*Gmu*MW**2*dsqrt(S2TW)
     C                    /(dsqrt(1.d0-S2TW)*16.d0*Pi**2)*aux


C	2) chargino/neutralino contributions

        aux=0.d0
        do i=1,2
	aux=aux+(-2.d0*k**2*(B1zdec(k,MCHA(i),MCHA(i))
     C		 +B21(k,MCHA(i),MCHA(i)))
     C           +2.d0*MCHA(i)**2*B1zdec(k,MCHA(i),MCHA(i))
     C	         +MCHA(i)**2*B0zdec(k,MCHA(i),MCHA(i)))
     C           *(2.d0*V(i,1)**2+V(i,2)**2+2.d0*U(i,1)**2
     C             +U(i,2)**2-4.d0*S2TW)
        enddo

        SigmaGamZSusy=SigmaGamZSusy
     C  -4.d0*dsqrt(2.d0)*Gmu*MW**2*dsqrt(S2TW/(1.d0-S2TW))
     C	/(16.d0*Pi**2)*aux

c	print*,'SiggamZ',-4.d0*dsqrt(2.d0)*Gmu*
c     C MW**2*dsqrt(S2TW/(1.d0-S2TW))
c     C	/(16.d0*Pi**2)*aux*(-2.d0)*dsqrt(C2TW/S2TW)/MZ**2,C2TW/S2TW,
c     C SigmaGamZSusy*(-2.d0)*dsqrt(C2TW/S2TW)/MZ**2

	return
	end
C*********************************************************************
C		->Gamma-Prime-Selfenergy SUSY-contribution
	
	DOUBLE PRECISION function SigmaGamPrSusy(k,MW)
	
	implicit none
	INTEGER i
	DOUBLE PRECISION k,B1zdec,B21,B5pr,B0pr,B1pr,B21pr,Cf
	DOUBLE PRECISION aux,C2TW,pi,MW,MZ,S2TW,Gmu
	DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
	DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
	DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	DOUBLE PRECISION MSU(2),MSD(2),MSE(2),MST(2),MSB(2),MSL(2)
	DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
	DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C		MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C		CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
        COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0

	pi=4.d0*datan(1.d0) 
	Gmu=GF0 !1.16637d-5

C Weinberg angle
	MZ=MZ0  !91.1876d0
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW
 
	Cf=3.d0
	
C	    Additionnal contributions from the SUSY sector

C       0) Masses and couplings

C	everything for 1.+ 2. Generation	
        MSU(1)=MUL
        MSU(2)=MUR
        MSD(1)=MDL
        MSD(2)=MDR
        MSE(1)=MLL
        MSE(2)=MLR
	
C	everything for Stop	
        MST(1)=dsqrt(RMST1)
	MST(2)=dsqrt(RMST2)

C	everything for Sbottom
        MSB(1)=dsqrt(RMSB1)
	MSB(2)=dsqrt(RMSB2)

C	everthing for Stau
        MSL(1)=MSL1
	MSL(2)=MSL2
 

C	1) sfermion contributions
      
        aux=0.d0
        do i=1,2
C	3. Generation	
         aux=aux
     C	    +4.d0/9.d0*Cf*B5pr(k,MST(i),MST(i))
     C      +1.d0/9.d0*Cf*B5pr(k,MSB(i),MSB(i))
     C      +B5pr(k,MSL(i),MSL(i))
C	1. + 2. Generation
         aux=aux
     C		+2.d0*(
     C	     4.d0/9.d0*Cf*B5pr(k,MSU(i),MSU(i))
     C      +1.d0/9.d0*Cf*B5pr(k,MSD(i),MSD(i))
     C      +B5pr(k,MSE(i),MSE(i)))
        enddo

	SigmaGamPrSusy=4.d0*dsqrt(2.d0)*Gmu*MW**2*S2TW/(16.d0*Pi**2)
     C 			*aux



C	2) chargino/neutralino contributions

	aux=0.d0
        do i=1,2
	aux=aux+(-2.d0*(B1zdec(k,MCHA(i),MCHA(i))
     C      +B21(k,MCHA(i),MCHA(i)))-2.d0*k**2
     C		*(B1pr(k,MCHA(i),MCHA(i))
     C      +B21pr(k,MCHA(i),MCHA(i)))
     C      +2.d0*MCHA(i)**2*B1pr(k,MCHA(i),MCHA(i))
     C      +MCHA(i)**2*B0pr(k,MCHA(i),MCHA(i)))
        enddo


        SigmaGamPrSusy=SigmaGamPrSusy
     C        +4.d0*4.d0*dsqrt(2.d0)*Gmu*MW**2*S2TW/(16.d0*Pi**2)*aux
	
c	print*,'-SiggamPr',
c     c -4.d0*4.d0*dsqrt(2.d0)*Gmu*MW**2*S2TW/(16.d0*Pi**2)*aux,
c     c -SigmaGamPrSusy
	return
	end
C*********************************************************************
C		->Z-Prime-Selfenergy SUSY-contribution
	
	DOUBLE PRECISION function SigmaZPrSusy(k,MW)
	
	implicit none
	INTEGER i,j
	DOUBLE PRECISION k,B1zdec,B21,B22pr,B1pr,B0pr,B21pr,Cf
	DOUBLE PRECISION aux,C2TW,pi,MW,MZ,S2TW,Gmu
	DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
	DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
	DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	DOUBLE PRECISION MSU(2),MSD(2),MSE(2),MST(2),MSB(2),MSL(2)
	DOUBLE PRECISION UT(2,2),UB(2,2),UL(2,2),DELT(2,2)
	DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
	DOUBLE PRECISION RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
	COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     C		MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     C		CST,CSB,CSL,MSMU1,MSMU2,MSMUNT,CSMU
	COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
        COMMON/RADCOR/RMST1,RMST2,S2T,RMSB1,RMSB2,S2B,XT,XB
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0

	pi=4.d0*datan(1.d0) 
	Gmu=GF0  !1.16637d-5

C Weinberg angle
	MZ=MZ0  !91.1876d0
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW

	Cf=3.d0

	
C	    Additionnal contributions from the SUSY sector

C       0) Masses and couplings

C	everything for 1.+ 2. Generation	
        MSU(1)=MUL
        MSU(2)=MUR
        MSD(1)=MDL
        MSD(2)=MDR
        MSE(1)=MLL
        MSE(2)=MLR
	
C	everything for Stop	
        MST(1)=dsqrt(RMST1)
	MST(2)=dsqrt(RMST2)
        UT(1,1)=CST
        UT(1,2)=+dsqrt(1.d0-CST**2)
        UT(2,1)=-dsqrt(1.d0-CST**2)
        UT(2,2)=CST

C	everything for Sbottom
        MSB(1)=dsqrt(RMSB1)
	MSB(2)=dsqrt(RMSB2)
        UB(1,1)=CSB
        UB(1,2)=+dsqrt(1.d0-CSB**2)
        UB(2,1)=-dsqrt(1.d0-CSB**2)
        UB(2,2)=CSB

C	everthing for Stau
        MSL(1)=MSL1
	MSL(2)=MSL2
        UL(1,1)=CSL
        UL(1,2)=+dsqrt(1.d0-CSL**2)
        UL(2,1)=-dsqrt(1.d0-CSL**2)
        UL(2,2)=CSL

C	"Mixing matrices" for 1. and 2. Generation
C	(as mixing is absent here it is of course the unity matrix)		
        DELT(1,1)=1.d0
        DELT(1,2)=0.d0
        DELT(2,1)=0.d0
        DELT(2,2)=1.d0



C	1) sfermion contributions


        aux=0.d0
        do i=1,2
C	3. Generation
         do j=1,2
          aux=aux
     C	   	 -2.d0*Cf*((1.d0/2.d0-2.d0*S2TW/3.d0)*UT(i,1)
     C		        *UT(j,1)-2.d0/3.d0*S2TW*UT(i,2)*UT(j,2))**2
     C             *B22pr(k,MST(i),MST(j))
     C           -2.d0*Cf*((1.d0/2.d0-S2TW/3.d0)*UB(i,1)*UB(j,1)
     C                  -1.d0/3.d0*S2TW*UB(i,2)*UB(j,2))**2
     C             *B22pr(k,MSB(i),MSB(j))
     C           -2.d0*((1.d0/2.d0-S2TW)*UL(i,1)*UL(j,1)
     C                               -S2TW*UL(i,2)*UL(j,2))**2
     C             *B22pr(k,MSL(i),MSL(j))
         enddo
C	1. + 2. Generation
         do j=1,2
          aux=aux
     C      -4.d0*Cf*((1.d0/2.d0-2.d0*S2TW/3.d0)*DELT(i,1)*DELT(j,1)
     C                  -2.d0/3.d0*S2TW*DELT(i,2)*DELT(j,2))**2
     C             *B22pr(k,MSU(i),MSU(j))
     C      -4.d0*Cf*((1.d0/2.d0-S2TW/3.d0)*DELT(i,1)*DELT(j,1)
     C               -1.d0/3.d0*S2TW*DELT(i,2)*DELT(j,2))**2
     C             *B22pr(k,MSD(i),MSD(j))
     C      -4.d0*((1.d0/2.d0-S2TW)*DELT(i,1)*DELT(j,1)
     C                         -S2TW*DELT(i,2)*DELT(j,2))**2
     C             *B22pr(k,MSE(i),MSE(j))
         enddo
        enddo
        aux=aux-1.d0/2.d0*B22pr(k,MSNT,MSNT)
     C     -B22pr(k,MNL,MNL)

	SigmaZPrSusy=4.d0*dsqrt(2.d0)*Gmu*MW**2
     C                *2.d0/((1.d0-S2TW)*16.d0*Pi**2)*aux


C	->chargino
        aux=0.d0
        do i=1,2
        do j=1,2
	aux=aux+
     C      (-2.d0*(B1zdec(k,MCHA(j),MCHA(i))
     C       +B21(k,MCHA(j),MCHA(i)))-2.d0*k**2
     C		*(B1pr(k,MCHA(j),MCHA(i))
     C       +B21pr(k,MCHA(j),MCHA(i)))
     C	     +MCHA(i)**2*B1pr(k,MCHA(j),MCHA(i))
     C       +MCHA(j)**2*B1pr(k,MCHA(i),MCHA(j)))
     C	     *((V(j,1)*V(i,1)+1.d0/2.d0*V(j,2)*V(i,2)
     C	     -S2TW*delt(i,j))**2+(U(i,1)*U(j,1)+1.d0/2.d0
     C	     *U(i,2)*U(j,2)-S2TW*delt(i,j))**2)	
     C       +MCHA(i)*MCHA(j)*B0pr(k,MCHA(j),MCHA(i))
     C       *(V(j,1)*V(i,1)+1.d0/2.d0*V(j,2)*V(i,2)-S2TW*delt(i,j))
     C	     *(U(i,1)*U(j,1)+1.d0/2.d0*U(i,2)*U(j,2)-S2TW*delt(i,j))	
        enddo
        enddo
	
C	->Neutralino
        do i=1,5
        do j=1,5
	aux=aux+
     C	1.d0/4.d0*(
     C  (-2.d0*(B1zdec(k,MNEU(j),MNEU(i))+B21(k,MNEU(j),MNEU(i)))
     C	-2.d0*k**2*(B1pr(k,MNEU(j),MNEU(i))+B21pr(k,MNEU(j),MNEU(i)))
     C   +MNEU(i)**2*B1pr(k,MNEU(j),MNEU(i))
     C   +MNEU(j)**2*B1pr(k,MNEU(i),MNEU(j)))
     C  *(NEU(j,4)*NEU(i,4)-NEU(j,3)*NEU(i,3))**2
     C  -MNEU(i)*MNEU(j)*B0pr(k,MNEU(j),MNEU(i))
     C   *(NEU(j,4)*NEU(i,4)-NEU(j,3)*NEU(i,3))**2)
        enddo
        enddo

        SigmaZPrSusy=SigmaZPrSusy
     C          +4.d0*dsqrt(2.d0)*Gmu*MW**2*2.d0/(16.d0*Pi**2
     C		*(1.d0-S2TW))*aux	


	return
	end
C*********************************************************************
C*********************************************************************
C	renormalized selfenergies with Susy-contributions (see e.g [5])
C*********************************************************************
C*********************************************************************
C		-> renormalized W-Selfenergy -> Susy-contributions 
	DOUBLE PRECISION function SigmaWrenSusy(k,MW)

	implicit none
	DOUBLE PRECISION k,C2TW,dMZ2,dMW2,dzw2
	DOUBLE PRECISION SigmaWSusy,SigmaZSusy,SigmaGamZSusy
	DOUBLE PRECISION SigmaGamPrSusy,MW,MZ,S2TW
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0

C Weinberg angle
	MZ=MZ0 !91.1876d0
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW

	dMZ2=SigmaZSusy(MZ,MW)
	dMW2=SigmaWSusy(MW,MW)
	dzw2=-SigmaGamPrSusy(0.d0,MW)-2.d0*dsqrt(C2TW/S2TW)
     C		*SigmaGamZSusy(0.d0,MW)/MZ**2+C2TW/S2TW
     C		*(dMZ2/MZ**2-dMW2/MW**2)


	SigmaWrenSusy=SigmaWSusy(k,MW)-dMW2+dzw2*(k**2-MW**2)
c	print*,'kSigR',k

	return
	end
C*********************************************************************
C		-> renormalized Z-Selfenergy -> Susy-contributions
	DOUBLE PRECISION function SigmaZrenSusy(k,MW)

	implicit none
	DOUBLE PRECISION k,C2TW,MW,MZ,S2TW
	DOUBLE PRECISION SigmaZSusy,SigmaWSusy,SigmaGamZSusy
	DOUBLE PRECISION SigmaGamPrSusy
	DOUBLE PRECISION dMZ2,dMW2,dzz2
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0

C Weinberg angle
	MZ=MZ0 !91.1876d0
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW


	dMZ2=SigmaZSusy(MZ,MW)
	dMW2=SigmaWSusy(MW,MW)
	dzz2=-SigmaGamPrSusy(0.d0,MW)-2.d0*(C2TW-S2TW)/dsqrt(S2TW*C2TW)
     C		*SigmaGamZSusy(0.d0,MW)/MZ**2+(C2TW-S2TW)/S2TW
     C		*(dMZ2/MZ**2-dMW2/MW**2)

	SigmaZrenSusy=SigmaZSusy(k,MW)-dMZ2+dzz2*(k**2-MZ**2)

	return
	end
C*********************************************************************
C	      -> renormalized Gamma-Selfenergy -> Susy-contributions

	DOUBLE PRECISION function SigmaGamrenSusy(k,MW)

	implicit none
	DOUBLE PRECISION k,MW,SigmaGamSusy,SigmaGamPrSusy

	SigmaGamrenSusy=SigmaGamSusy(k,MW)-SigmaGamPrSusy(0.d0,MW)*k**2

	return
	end
C*********************************************************************
C	      -> renormalized Gamma-Z-Selfenergy -> Susy-contributions

	DOUBLE PRECISION function SigmaGamZrenSusy(k,MW)

	implicit none
	DOUBLE PRECISION k,C2TW,MW,MZ,S2TW
	DOUBLE PRECISION SigmaGamZSusy,SigmaGamPrSusy,SigmaZSusy
	DOUBLE PRECISION SigmaWSusy
	DOUBLE PRECISION dzgaz2,dzgaz1,dzz1,dzz2,dzga1,dzga2
	DOUBLE PRECISION dMZ2,dMW2
	DOUBLE PRECISION ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	DOUBLE PRECISION MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0
        COMMON/GAUGE/ALSMZ0,ALEMMZ0,GF0,g10,g20,S2TW0
	COMMON/SMSPEC/MS,MC,MBNP,MB,MT0,MTAU,MMUON,MZ0,MW0

C Weinberg angle
	MZ=MZ0 !91.1876d0
	S2TW=1.d0-MW**2/MZ**2
	C2TW=1.d0-S2TW

	dMZ2=SigmaZSusy(MZ,MW)
	dMW2=SigmaWSusy(MW,MW)

	dzga1=-SigmaGamPrSusy(0.d0,MW)-dsqrt(S2TW/C2TW)
     C	      *SigmaGamZSusy(0.d0,MW)/MZ**2
	dzga2=-SigmaGamPrSusy(0.d0,MW)

	dzz1=-SigmaGamPrSusy(0.d0,MW)-(3.d0*C2TW-2.d0*S2TW)
     C	        /dsqrt(C2TW*S2TW)*SigmaGamZSusy(0.d0,MW)/MZ**2
     C	        +(C2TW-S2TW)/S2TW*(dMZ2/MZ**2-dMW2/MW**2)


	dzz2=-SigmaGamPrSusy(0.d0,MW)-2.d0*(C2TW-S2TW)/dsqrt(S2TW*C2TW)
     C		*SigmaGamZSusy(0.d0,MW)/MZ**2+(C2TW-S2TW)/S2TW
     C		*(dMZ2/MZ**2-dMW2/MW**2)

	dzgaz1=dsqrt(C2TW*S2TW)/(C2TW-S2TW)*(dzz1-dzga1)

	dzgaz2=dsqrt(C2TW*S2TW)/(C2TW-S2TW)*(dzz2-dzga2)


	SigmaGamZrenSusy=SigmaGamZSusy(k,MW)-dzgaz2*k**2
     C          	  +(dzgaz1-dzgaz2)*MZ**2

	return
	end
C*********************************************************************
C	      -> F1(x,y)

	DOUBLE PRECISION function Fone(x,y)

	implicit none
	DOUBLE PRECISION x,y,Li_2,aux
	if(dabs(x-y).le.1.d-5) then
	 aux=0.d0
	else
	 aux=x+y-2.d0*x*y/(x-y)*dlog(x/y)*(2.d0+x/y*dlog(x/y))
     c     +(x+y)*x**2/(x-y)**2*dlog(x/y)**2-2.d0*(x-y)
     c     *Li_2(1.d0-x/y)
	endif
	Fone=aux

	return
	end
C*********************************************************************
C*********************************************************************
C	      -> Phi(x,y,z)

	DOUBLE PRECISION function PhiGL(x,y,z)

	implicit none
	DOUBLE PRECISION x,y,z,Li_2,Cl2,aux,aux1,Pi
	Pi=4.d0*datan(1.d0)
	aux1=(z-x-y)**2-4.d0*x*y
	If(aux1.ge.0.d0)THEN
	aux=dlog((z+x-y-dsqrt(aux1))/2.d0/z)
     C *dlog((z+y-x-dsqrt(aux1))/2.d0/z)
     C -dlog(x/z)*dlog(y/z)/2.d0
     C -Li_2((z+x-y-dsqrt(aux1))/2.d0/z)
     C -Li_2((z+y-x-dsqrt(aux1))/2.d0/z)+Pi**2/6
	ELSE
	aux=Cl2(2.d0*dacos((x+y-z)/(2.d0*dsqrt(x*y))))
     C +Cl2(2.d0*dacos((x-y+z)/(2.d0*dsqrt(x*z))))
     C +Cl2(2.d0*dacos((z-x+y)/(2.d0*dsqrt(y*z))))
	ENDIF
	PhiGL=2*z*aux/dsqrt(dabs(aux1))

	return
	end
C*********************************************************************
C*********************************************************************
C	      -> Cl2(x,y,z)

	DOUBLE PRECISION function Cl2(x)

	implicit none
	DOUBLE PRECISION x,aux,Pi
	Pi=4.d0*datan(1.d0)
	If(x.ge.0.d0.and.x.le.2.d0)then
	 aux=(1.d0-dlog(x))*x+1.388888889d-2*x**3
     C +6.94444445d-5*x**5+7.873519797d-7*x**7
	else 
	 If(x.ge.2.d0.and.x.le.2.d0*(Pi-1.d0))then
	  aux=5.33186646836d-1+1.10292953236d0*x
     C -7.514016413236d-1*x**2+1.5599160887d-3*x**3
     C -1.8207027275269d-2*x**4+1.159095371966d-3*x**5
	 else
	  If(x.ge.2.d0*(Pi-1.d0).and.x.le.2.d0*Pi)then
	   aux=-(1.d0-dlog(2.d0*Pi-x))*(2.d0-x)
     C -1.388888889d-2*(2.d0*Pi-x)**3
     C -6.94444445d-5*(2.d0*Pi-x)**5
     C -7.873519797d-7*(2.d0*Pi-x)**7
	  else 
	   aux=0.d0
	  endif
	 endif
	endif
	Cl2=aux
	return
	end
C*********************************************************************
C*********************************************************************
C	      -> B0Yuk

	DOUBLE PRECISION function B0Yuk(x,y)

	implicit none
	DOUBLE PRECISION x,y
	B0Yuk=-dlog(x**2)+y**2/(y**2-x**2)*dlog((x/y)**2)

	return
	end

