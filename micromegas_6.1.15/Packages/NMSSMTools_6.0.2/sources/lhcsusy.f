      SUBROUTINE LHCSUSY(PAR,PROB,FLAG)

*********************************************************************
*   Subroutine to check SUSY constraints at the LHC using SmodelS
*      PROB(88) =/= 0  excluded with R=RMAX via channel CHAN
*********************************************************************

      IMPLICIT NONE

      CHARACTER CHINL*120,CHAN*20
      CHARACTER(200) PAT,IFILE,OFILE,EFILE,TFILE,SFILE,SMODELS_PATH

      INTEGER CFLAG(6),FLAG,I,N1,N2,N3
      PARAMETER(N1=164,N2=8,N3=19)

      DOUBLE PRECISION PAR(*),PROB(*),V1,V2,V3,V4,R,X,Y,X1(N1),Y1(N1)
      DOUBLE PRECISION X2(N2),Y2H(N2),Y2L(N2),X3(N3),Y3(N3)
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      DOUBLE PRECISION brsellneute(5),brselrneute(5),brsellcharnue(2),
     .         brselrcharnue(2),brsnellneut(5),brsnellchar(5),
     .         brsmu1neutmu(5),brsmu2neutmu(5),brsmu1charnumu(2),
     .         brsmu2charnumu(2),brsnmu1neut(5),brsnmu1char(5),
     .         brstau1neut(5),brstau2neut(5),brstau1char(2),
     .         brstau2char(2),brstau1hcsn(2),brstau2hcsn(2),
     .         brstau1wsn(2),brstau2wsn(2),brstau2ztau,brstau2H(3),
     .         brstau2A(2),brsntauneut(5),brsntauchar(2),
     .         brsntau1hcstau(2),brsntau1wstau(2)

      COMMON/SLEPTON_BR_2BD/brsellneute,brselrneute,brsellcharnue,
     .         brselrcharnue,brsnellneut,brsnellchar,brsmu1neutmu,
     .         brsmu2neutmu,brsmu1charnumu,brsmu2charnumu,brsnmu1neut,
     .         brsnmu1char,brstau1neut,brstau2neut,brstau1char,
     .         brstau2char,brstau1hcsn,brstau2hcsn,brstau1wsn,
     .         brstau2wsn,brstau2ztau,brstau2H,brstau2A,brsntauneut,
     .         brsntauchar,brsntau1hcstau,brsntau1wstau
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     . MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     . CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      COMMON/FILES/PAT,IFILE,OFILE,EFILE,TFILE,SFILE
      COMMON/SMODELS/R,CHAN
      COMMON/CFLAG/CFLAG

      DATA X1/0.97980d2,0.10051d3,0.10109d3,0.10303d3,0.10556d3,
     .0.10573d3,0.10808d3,0.11014d3,0.11061d3,0.11313d3,0.11422d3,
     .0.11566d3,0.11792d3,0.11818d3,0.12071d3,0.12123d3,0.12323d3,
     .0.12417d3,0.12576d3,0.12678d3,0.12828d3,0.12909d3,0.13081d3,
     .0.13117d3,0.13305d3,0.13333d3,0.13480d3,0.13586d3,0.13644d3,
     .0.13804d3,0.13838d3,0.13965d3,0.14091d3,0.14131d3,0.14309d3,
     .0.14343d3,0.14503d3,0.14596d3,0.14717d3,0.14848d3,0.14958d3,
     .0.15101d3,0.15227d3,0.15354d3,0.15527d3,0.15606d3,0.15856d3,
     .0.15859d3,0.16111d3,0.16214d3,0.16364d3,0.16588d3,0.16616d3,
     .0.16869d3,0.16971d3,0.17121d3,0.17348d3,0.17374d3,0.17626d3,
     .0.17711d3,0.17879d3,0.18049d3,0.18131d3,0.18360d3,0.18384d3,
     .0.18636d3,0.18643d3,0.18889d3,0.18901d3,0.19136d3,0.19141d3,
     .0.19354d3,0.19394d3,0.19557d3,0.19646d3,0.19748d3,0.19899d3,
     .0.19929d3,0.20105d3,0.20152d3,0.20276d3,0.20404d3,0.20441d3,
     .0.20604d3,0.20657d3,0.20763d3,0.20909d3,0.20917d3,0.21066d3,
     .0.21162d3,0.21204d3,0.21331d3,0.21414d3,0.21439d3,0.21526d3,
     .0.21580d3,0.21597d3,0.21566d3,0.21481d3,0.21414d3,0.21333d3,
     .0.21162d3,0.21116d3,0.20909d3,0.20825d3,0.20657d3,0.20454d3,
     .0.20404d3,0.20152d3,0.20002d3,0.19899d3,0.19646d3,0.19463d3,
     .0.19394d3,0.19141d3,0.18889d3,0.18839d3,0.18636d3,0.18384d3,
     .0.18133d3,0.18131d3,0.17879d3,0.17626d3,0.17374d3,0.17362d3,
     .0.17121d3,0.16869d3,0.16616d3,0.16553d3,0.16364d3,0.16111d3,
     .0.15859d3,0.15739d3,0.15606d3,0.15354d3,0.15101d3,0.14944d3,
     .0.14848d3,0.14596d3,0.14343d3,0.14180d3,0.14091d3,0.13838d3,
     .0.13586d3,0.13444d3,0.13333d3,0.13081d3,0.12828d3,0.12723d3,
     .0.12576d3,0.12323d3,0.12071d3,0.11992d3,0.11818d3,0.11566d3,
     .0.11313d3,0.11209d3,0.11061d3,0.10808d3,0.10556d3,0.10303d3,
     .0.10286d3,0.10051d3,0.97980d2/
      DATA Y1/0.86180d0,0.88450d0,0.88980d0,0.90870d0,0.93330d0,
     .0.93500d0,0.96010d0,0.98240d0,0.98790d0,0.10186d1,0.10320d1,
     .0.10520d1,0.10844d1,0.10887d1,0.11306d1,0.11395d1,0.11784d1,
     .0.11970d1,0.12337d1,0.12578d1,0.12990d1,0.13213d1,0.13766d1,
     .0.13884d1,0.14588d1,0.14706d1,0.15325d1,0.15827d1,0.16103d1,
     .0.16920d1,0.17104d1,0.17775d1,0.18459d1,0.18677d1,0.19625d1,
     .0.19806d1,0.20616d1,0.21077d1,0.21662d1,0.22269d1,0.22756d1,
     .0.23378d1,0.23911d1,0.24429d1,0.25125d1,0.25433d1,0.26394d1,
     .0.26406d1,0.27346d1,0.27733d1,0.28288d1,0.29141d1,0.29248d1,
     .0.30213d1,0.30613d1,0.31218d1,0.32166d1,0.32278d1,0.33404d1,
     .0.33799d1,0.34618d1,0.35506d1,0.35958d1,0.37308d1,0.37463d1,
     .0.39147d1,0.39201d1,0.41087d1,0.41181d1,0.43271d1,0.43321d1,
     .0.45457d1,0.45899d1,0.47764d1,0.48876d1,0.50188d1,0.52288d1,
     .0.52723d1,0.55399d1,0.56156d1,0.58210d1,0.60492d1,0.61150d1,
     .0.64254d1,0.65328d1,0.67515d1,0.70762d1,0.70925d1,0.74525d1,
     .0.77161d1,0.78289d1,0.82262d1,0.85487d1,0.86437d1,0.90803d1,
     .0.95411d1,0.10025d2,0.10532d2,0.11066d2,0.11316d2,0.11628d2,
     .0.12089d2,0.12215d2,0.12653d2,0.12835d2,0.13128d2,0.13486d2,
     .0.13561d2,0.13941d2,0.14168d2,0.14305d2,0.14642d2,0.14887d2,
     .0.14969d2,0.15276d2,0.15581d2,0.15639d2,0.15867d2,0.16151d2,
     .0.16432d2,0.16432d2,0.16707d2,0.16979d2,0.17250d2,0.17266d2,
     .0.17523d2,0.17795d2,0.18072d2,0.18138d2,0.18349d2,0.18629d2,
     .0.18919d2,0.19059d2,0.19213d2,0.19516d2,0.19829d2,0.20026d2,
     .0.20146d2,0.20479d2,0.20816d2,0.21038d2,0.21164d2,0.21528d2,
     .0.21893d2,0.22105d2,0.22274d2,0.22667d2,0.23062d2,0.23227d2,
     .0.23464d2,0.23873d2,0.24277d2,0.24401d2,0.24689d2,0.25090d2,
     .0.25480d2,0.25639d2,0.25870d2,0.26242d2,0.26589d2,0.26915d2,
     .0.26934d2,0.27227d2,0.27511d2/
      DATA X2/90d0,115d0,140d0,165d0,190d0,210d0,215d0,216d0/
      DATA Y2H/89d0,114d0,138.5d0,162d0,186d0,203d0,206d0,206d0/
      DATA Y2L/62d0,90d0,119d0,147d0,174d0,198d0,204d0,206d0/
      DATA X3/90d0,100d0,110d0,120d0,130d0,140d0,150d0,160d0,170d0,
     .        180d0,190d0,197d0,202d0,203d0,203.5d0,204d0,206.5d0,
     .        211d0,216d0/
      DATA Y3/52.5d0,54d0,58d0,63d0,68d0,71d0,75d0,81d0,84.5d0,86d0,
     .        83.5d0,77.5d0,69.5d0,65d0,53d0,49d0,43d0,31d0,0d0/

      PROB(88)=0d0

      IF(CFLAG(6).EQ.-1)THEN
       PROB(88)=-DDIM(1d0,MIN(MGL,MUR,MUL,MDR,MDL,MST1,MSB1)/1d3)
       RETURN
      ENDIF

      R=0d0
      CHAN=''

      CALL SOUTPUT(PAR,PROB)
      CALL getenv('SMODELS_PATH',SMODELS_PATH)
      if(SMODELS_PATH.eq.' ')  SMODELS_PATH='../smodels'
      CALL SYSTEM(trim(SMODELS_PATH)//'/smodelsTools.py xseccomputer '
     .//'-s 8 13 -e 10000 -p -f '//trim(TFILE)//' -v error')

      CALL SYSTEM(trim(SMODELS_PATH)//'/runSModelS.py -f '//trim(TFILE)
     .//' -p '//trim(SMODELS_PATH)//'parameters.ini -o '//trim(PAT)
     .//' -v error')

      OPEN(19,FILE=trim(TFILE)//'.smodels ',STATUS='UNKNOWN')
*   START TO READ NEW LINE INTO CHINL
 21   CHINL=' '
      READ(19,'(A120)',END=29) CHINL
      IF(CHINL(1:9).ne."#Analysis") GOTO 21
      READ(19,'(A120)',END=29) CHINL
      READ(19,'(A120)',END=29) CHINL
      READ(CHINL,*) CHAN,V1,V2,V3,V4,R
      CLOSE(19)
!      WRITE(0,*)CHAN,R

 29   IF(R.GT.1d0)PROB(88)=R-1d0
      IF(R.EQ.0d0)PROB(88)=-1d0

      CALL SYSTEM('rm '//trim(TFILE))
      IF(FLAG.EQ.0)THEN
       CALL SYSTEM('rm '//trim(TFILE)//'.smodels')
      ELSE
       CALL SYSTEM('mv '//trim(TFILE)//'.smodels '//trim(SFILE))
      ENDIF

      IF(CFLAG(6).EQ.1 .OR.R.EQ.0d0)RETURN

c M_sleptonL vs M_neutralino1

      IF(brsellneute(1).GT.0.5d0)THEN

* ...from ATLAS 1911.12606 fig 16b

       X=MLL
       Y=MLL-DABS(MNEU(1))

       I=1
       DOWHILE(Y1(I).LE.Y .AND. I.LT.N1)
        I=I+1
       ENDDO
       IF(I.GE.2 .AND. Y.LE.Y1(N1))PROB(88)=PROB(88)+
     . DDIM((X1(I-1)+(Y-Y1(I-1))/(Y1(I)-Y1(I-1))*(X1(I)-X1(I-1)))/X,1d0)

* ...and ATLAS-CONF-2022-006 FIG.7b

        X=MLL
        Y=DABS(MNEU(1))

*    orange zone

       I=1
       DOWHILE(X2(I).LE.X .AND. I.LT.N2)
        I=I+1
       ENDDO

       IF(I.GE.2 .AND. X.LE.X2(N2)) PROB(88)=PROB(88)+
     .  DDIM(Y/(Y2L(I-1)+(X-X2(I-1))/(X2(I)-X2(I-1))*(Y2L(I)-Y2L(I-1))),
     .  1d0) * DDIM(1d0,
     .  Y/(Y2H(I-1)+(X-X2(I-1))/(X2(I)-X2(I-1))*(Y2H(I)-Y2H(I-1))))

*    blue curve

       I=1
       DOWHILE(X3(I).LE.X .AND. I.LT.N3)
        I=I+1
       ENDDO

       IF(I.GE.2 .AND. X.LE.X3(N3)) PROB(88)=PROB(88)+DDIM(1d0,
     .  Y/(Y3(I-1)+(X-X3(I-1))/(X3(I)-X3(I-1))*(Y3(I)-Y3(I-1))))

      ENDIF


      END


      SUBROUTINE SOUTPUT(PAR,PROB)

*********************************************************************
*   Subroutine writing the results in the the output file for SmodelS
*********************************************************************

      IMPLICIT NONE

      CHARACTER(200) PAT,IFILE,OFILE,EFILE,TFILE,SFILE

      INTEGER I,NBIN,IFAIL,Q2FIX,NMSFLAG,OMGFLAG,MAFLAG,MOFLAG
      INTEGER PFLAG,VFLAG,OUTFLAG,GRFLAG,MWFLAG,NSUSY,NGUT,NMES
      INTEGER CFLAG(6),PDGH1,PDGH2

      PARAMETER (NSUSY=14,NGUT=21,NMES=21)

      DOUBLE PRECISION PAR(*),PROB(*),SIG(5,8)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),AMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION MGL,MCHA(2),U(2,2),V(2,2),MNEU(5),NEU(5,5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5),CL(5)
      DOUBLE PRECISION ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      DOUBLE PRECISION MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      DOUBLE PRECISION VUS,VCB,VUB,TANB,SINB,COSB
      DOUBLE PRECISION MUR,MUL,MDR,MDL,MLR,MLL,MNL
      DOUBLE PRECISION MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT
      DOUBLE PRECISION CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      DOUBLE PRECISION SST,SSB,SSL,Q2,Q2MIN
      DOUBLE PRECISION MGUT,g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      DOUBLE PRECISION MHUS,MHDS,MSS
      DOUBLE PRECISION G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT
      DOUBLE PRECISION HBOTGUT,HTAUGUT,M1GUT,M2GUT,M3GUT,ALGUT,AKGUT
      DOUBLE PRECISION ATGUT,ABGUT,ATAUGUT,AMUGUT,MHUGUT,MHDGUT
      DOUBLE PRECISION MSGUT,MQ3GUT,MU3GUT,MD3GUT,MQGUT,MUGUT
      DOUBLE PRECISION MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      DOUBLE PRECISION XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      DOUBLE PRECISION XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      DOUBLE PRECISION OMG,OMGMIN,OMGMAX
      DOUBLE PRECISION Xf,sigmaV,vcsll,vcsbb,x(100),dNdx(100),EMIN
      DOUBLE PRECISION sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      DOUBLE PRECISION XSMAX,XENON_SI,XENON_SDn,XENON_SDp,PandaX_SI
      DOUBLE PRECISION LUX_SI,LUX_SDn,LUX_SDp,PICO60_SDp
      DOUBLE PRECISION CRESST_SI,DarkSide50_SI
      DOUBLE PRECISION PRINTCHANNELS,omg_
      DOUBLE PRECISION MHUQ,MHDQ,MSX,LQ,KQ,ALQ,AKQ,MUQ,NUQ
      DOUBLE PRECISION ZHU,ZHD,ZS,H1Q,H2Q,TANBQ,QSTSB
      DOUBLE PRECISION BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     . DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     . BRBtaunumax,BRBtaunumin
      DOUBLE PRECISION BRBSll,BRBSllmin,BRBSllMax,
     .             BRBShll,BRBShllmin,BRBShllMax,
     .             BRDG,BRDGmin,BRDGmax,
     .             BRBdMUMU,BRBdMUMUmin,BRBdMUMUmax,
     .             BRBXsnunu,BRBXsnunumin,BRBXsnunumax,
     .             BRBpKpnunu,BRBpKpnunumin,BRBpKpnunuMax,
     .             BRBKsnunu,BRBKsnunumin,BRBKsnunuMax,
     .             RD_taul,RD_taulmin,RD_taulmax,
     .             RDs_taul,RDs_taulmin,RDs_taulmax
      DOUBLE PRECISION BRKp_Pipnunub,BRKp_Pipnunubmin,BRKp_PipnunubMax,
     .       BRKL_Pi0nunub,BRKL_Pi0nunubmin,BRKL_Pi0nunubMax,
     .       DMK,DMKmin,DMKmax,epsK,epsKmin,epsKmax
      DOUBLE PRECISION eps0,epst0,epst1,epst2,epst3,epsts,epstb
      DOUBLE PRECISION epsY32,epsY31,epsY23,epsY13,epscs,epscb
      DOUBLE PRECISION delmagmu,damumin,damumax,amuthmax,amuthmin
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2),toptot
      DOUBLE PRECISION FTSUSY(NSUSY+2),FTGUT(NGUT+2),FTMES(NMES+2)
      DOUBLE PRECISION PX,PA(6),PB(2),PL(7),PK(8)
      DOUBLE PRECISION MHmin,MHmax
      DOUBLE PRECISION M32,CGR,MPL,DELMB,DELML,DEL1
      DOUBLE PRECISION DSMASS(3),DAMASS(2),DCMASS
      DOUBLE PRECISION xsectot,limtrilep
      DOUBLE PRECISION MWNMSSM,dumw,dMW0,DrNMSSM,MWSM,dMWSM,decztt,
     .      deltadecztt,deczee,deltadeczee,BRZTauTau,BRZTTmin,
     .      BRZTTmax,ratio,deltaratio,S2TWeffTau,deltaS2TWeffTau,
     .      S2TWeffSM

      COMMON/EWPO/MWNMSSM,dumw,dMW0,DrNMSSM,MWSM,dMWSM,decztt,
     .      deltadecztt,deczee,deltadeczee,BRZTauTau,BRZTTmin,
     .      BRZTTmax,ratio,deltaratio,S2TWeffTau,deltaS2TWeffTau,
     .      S2TWeffSM
      COMMON/PFLAG/PFLAG
      COMMON/NMSFLAG/NMSFLAG
      COMMON/FLAGS/OMGFLAG,MAFLAG,MOFLAG
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     . BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     . BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     . BRSUSY,WIDTH
      COMMON/BRC/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     . HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     . HCBRSUSY,HCWIDTH
      COMMON/BRSG/BRSG,BRSGmax,BRSGmin,DMd,DMdmin,DMdmax,DMs,
     .      DMsmax,DMsmin,BRBMUMU,BRBMUMUmax,BRBMUMUmin,BRBtaunu,
     .      BRBtaunumax,BRBtaunumin
      COMMON/FLAV2/BRBSll,BRBSllmin,BRBSllMax,
     .      BRBShll,BRBShllmin,BRBShllMax,
     .      BRDG,BRDGmin,BRDGmax,
     .      BRBdMUMU,BRBdMUMUmin,BRBdMUMUmax,
     .      BRBXsnunu,BRBXsnunumin,BRBXsnunumax,
     .      BRBpKpnunu,BRBpKpnunumin,BRBpKpnunuMax,
     .      BRBKsnunu,BRBKsnunumin,BRBKsnunuMax,
     .      RD_taul,RD_taulmin,RD_taulmax,
     .      RDs_taul,RDs_taulmin,RDs_taulmax
      COMMON/FLAV3/BRKp_Pipnunub,BRKp_Pipnunubmin,BRKp_PipnunubMax,
     .             BRKL_Pi0nunub,BRKL_Pi0nunubmin,BRKL_Pi0nunubMax,
     .             DMK,DMKmin,DMKmax,epsK,epsKmin,epsKmax
      COMMON/EPSCOUP/eps0,epst0,epst1,epst2,epst3,epsts,epstb,
     .               epsY32,epsY31,epsY23,epsY13,epscs,epscb
      COMMON/MAGMU/delmagmu,damumin,damumax,amuthmax,amuthmin
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/topwidth/toptot
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/CB/CB,CL
      COMMON/HIGGSPEC/SMASS,SCOMP,AMASS,PCOMP,CMASS
      COMMON/SUSYSPEC/MGL,MCHA,U,V,MNEU,NEU
      COMMON/GAUGE/ALSMZ,ALEMMZ,GF,g1,g2,S2TW
      COMMON/SMSPEC/MS,MC,MB,MBP,MT,MTAU,MMUON,MZ,MW
      COMMON/CKM/VUS,VCB,VUB
      COMMON/SFSPEC/MUR,MUL,MDR,MDL,MLR,MLL,MNL,
     . MST1,MST2,MSB1,MSB2,MSL1,MSL2,MSNT,
     . CST,CSB,CSL,MSMU1,MSMU2,MSNM,CSMU
      COMMON/Q2FIX/Q2MIN,Q2FIX
      COMMON/RENSCALE/Q2
      COMMON/STSBSCALE/QSTSB
      COMMON/MGUT/MGUT
      COMMON/SUSYCOUP/g1s,g2s,g3s,HTOPS,HBOTS,HTAUS
      COMMON/SUSYMH/MHUS,MHDS,MSS
      COMMON/QMHIGGS/MHUQ,MHDQ,MSX
      COMMON/QPAR/LQ,KQ,ALQ,AKQ,MUQ,NUQ
      COMMON/QHIGGS/ZHU,ZHD,ZS,H1Q,H2Q,TANBQ
      COMMON/GUTCOUP/G1GUT,G2GUT,G3GUT,LGUT,KGUT,HTOPGUT,
     . HBOTGUT,HTAUGUT
      COMMON/GUTPAR/M1GUT,M2GUT,M3GUT,ALGUT,AKGUT,ATGUT,ABGUT,
     . ATAUGUT,AMUGUT,MHUGUT,MHDGUT,MSGUT,MQ3GUT,MU3GUT,MD3GUT,
     . MQGUT,MUGUT,MDGUT,ML3GUT,ME3GUT,MLGUT,MEGUT
      COMMON/GUTEXT/XIFGUT,XISGUT,MUPGUT,MSPGUT,M3HGUT
      COMMON/SUSYEXT/XIFSUSY,XISSUSY,MUPSUSY,MSPSUSY,M3HSUSY
      COMMON/MICROMG/OMG,OMGMIN,OMGMAX,Xf,sigmaV,vcsll,vcsbb,
     .      x,dNdx,EMIN,NBIN
      COMMON/MICROMG2/sigmaPiN,sigmaS,csPsi,csNsi,csPsd,csNsd
      COMMON/FINETUN/FTSUSY,FTGUT,FTMES
      COMMON/EFFCOUP/PX,PA,PB,PL,PK
      COMMON/DELMB/DELMB,DELML,DEL1
      COMMON/LHCSIG/SIG
      COMMON/HIGGSFIT/MHmin,MHmax
      COMMON/VFLAG/VFLAG
      COMMON/OUTFLAG/OUTFLAG
      COMMON/M32/M32,CGR,MPL,GRFLAG
      COMMON/DHIGGSSPEC/DSMASS,DAMASS,DCMASS
      COMMON/XSECTRILEP/xsectot,limtrilep
      COMMON/MWFLAG/MWFLAG
      COMMON/CFLAG/CFLAG
      COMMON/FILES/PAT,IFILE,OFILE,EFILE,TFILE,SFILE
      COMMON/PDGHIGGS/PDGH1,PDGH2

      OPEN(18,FILE=TFILE,STATUS= 'UNKNOWN')

      TANB=PAR(3)
      COSB=1d0/DSQRT(1d0+TANB**2)
      SINB=TANB*COSB

      WRITE(18,899) "# Input parameters"
      WRITE(18,899) "BLOCK MODSEL"
      WRITE(18,921) 3,1,"NMSSM particle content"

      WRITE(18,899) "BLOCK SMINPUTS"
      WRITE(18,901) 1,1d0/ALEMMZ,"ALPHA_EM^-1(MZ)"
      WRITE(18,901) 2,GF,"GF"
      WRITE(18,901) 3,ALSMZ,"ALPHA_S(MZ)"
      WRITE(18,901) 4,MZ,"MZ"
      WRITE(18,901) 5,MB,"MB(MB)"
      WRITE(18,901) 6,MT,"MTOP (POLE MASS)"
      WRITE(18,901) 7,MTAU,"MTAU"

      WRITE(18,899) "BLOCK MINPAR"
      IF(Q2FIX.EQ.1)WRITE(18,901) 0,DSQRT(Q2),"REN. SCALE"
      WRITE(18,901) 3,TANB,"TANBETA(MZ)"

      WRITE(18,899) "BLOCK EXTPAR"
      WRITE(18,901) 1,PAR(20), "M1"
      WRITE(18,901) 2,PAR(21), "M2"
      WRITE(18,901) 3,PAR(22), "M3"
      WRITE(18,901) 11,PAR(12), "ATOP"
      WRITE(18,901) 12,PAR(13), "ABOTTOM"
      WRITE(18,901) 13,PAR(14), "ATAU"
      WRITE(18,901) 16,PAR(25), "AMUON"

      WRITE(18,901) 31,DSQRT(PAR(18)),"LEFT SELECTRON"
      WRITE(18,901) 32,DSQRT(PAR(18)),"LEFT SMUON"
      WRITE(18,901) 33,DSQRT(PAR(10)),"LEFT STAU"

      WRITE(18,901) 34,DSQRT(PAR(19)),"RIGHT SELECTRON"
      WRITE(18,901) 35,DSQRT(PAR(19)),"RIGHT SMUON"
      WRITE(18,901) 36,DSQRT(PAR(11)),"RIGHT STAU"

      WRITE(18,901) 41,DSQRT(PAR(15)),"LEFT 1ST GEN. SQUARKS"
      WRITE(18,901) 42,DSQRT(PAR(15)),"LEFT 2ND GEN. SQUARKS"
      WRITE(18,901) 43,DSQRT(PAR(7)),"LEFT 3RD GEN. SQUARKS"

      WRITE(18,901) 44,DSQRT(PAR(16)),"RIGHT U-SQUARKS"
      WRITE(18,901) 45,DSQRT(PAR(16)),"RIGHT C-SQUARKS"
      WRITE(18,901) 46,DSQRT(PAR(8)),"RIGHT T-SQUARKS"

      WRITE(18,901) 47,DSQRT(PAR(17)),"RIGHT D-SQUARKS"
      WRITE(18,901) 48,DSQRT(PAR(17)),"RIGHT S-SQUARKS"
      WRITE(18,901) 49,DSQRT(PAR(9)),"RIGHT B-SQUARKS"

      WRITE(18,901) 61,PAR(1),"LAMBDA"
      WRITE(18,901) 62,PAR(2),"KAPPA"
      IF(MOD(MAFLAG,3).NE.1)THEN
       WRITE(18,901) 63,PAR(5),"ALAMBDA"
      ELSE
       WRITE(18,920) 63,PAR(5),"ALAMBDA"
      ENDIF
      IF(MAFLAG/3.NE.1)THEN
       WRITE(18,901) 64,PAR(6),"AKAPPA"
      ELSE
       WRITE(18,920) 64,PAR(6),"AKAPPA"
      ENDIF
      WRITE(18,901) 65,PAR(4),"MUEFF"
      IF(MOD(MAFLAG,3).NE.2)THEN
       IF(XIFSUSY.NE.0d0)
     .  WRITE(18,901) 66,XIFSUSY,"XIF"
      ELSE
       WRITE(18,920) 66,XIFSUSY,"XIF"
      ENDIF
      IF(MAFLAG/3.NE.2)THEN
       IF(XISSUSY.NE.0d0)
     .  WRITE(18,901) 67,XISSUSY,"XIS"
      ELSE
       WRITE(18,920) 67,XISSUSY,"XIS"
      ENDIF
      IF(MUPSUSY.NE.0d0)
     .  WRITE(18,901) 68,MUPSUSY,"MUP"
      IF(MSPSUSY.NE.0d0)
     .  WRITE(18,901) 69,MSPSUSY,"MSP"
      IF(M3HSUSY.NE.0d0)
     .  WRITE(18,901) 72,M3HSUSY,"M3H"
      IF(MOD(MAFLAG,3).NE.0)THEN
       WRITE(18,901) 124,PAR(23),"MA AT QSTSB"
      ELSE
       WRITE(18,920) 124,PAR(23),"MA AT QSTSB"
      ENDIF
      IF(MAFLAG/3.NE.0)THEN
       WRITE(18,901) 125,PAR(24),"MP AT QSTSB"
      ELSE
       WRITE(18,920) 125,PAR(24),"MP AT QSTSB"
      ENDIF

      WRITE(18,899) "# "
      WRITE(18,899) "BLOCK MASS   # Mass spectrum "
      WRITE(18,899) "#  PDG Code     mass             particle "
      WRITE(18,902) 5,MB,"MB(MB)"
      WRITE(18,902) 6,MT,"MTOP (POLE MASS)"
      WRITE(18,902) 15,MTAU,"MTAU"
      WRITE(18,902) 23,MZ,"MZ"
      IF(MWFLAG.NE.0)THEN
       WRITE(18,902) 24,MWNMSSM,"MW incl. Delta_MW"
      ELSE
       WRITE(18,902) 24,MW,"MW"
      ENDIF
      WRITE(18,902) PDGH1,SMASS(1),"lightest neutral scalar"
      WRITE(18,902) PDGH2,SMASS(2),"second neutral scalar"
      WRITE(18,902) 45,SMASS(3),"third neutral scalar"
      WRITE(18,902) 36,AMASS(1),"lightest pseudoscalar"
      WRITE(18,902) 46,AMASS(2),"second pseudoscalar"
      WRITE(18,902) 37,CMASS,"charged Higgs"
      WRITE(18,902) 1000001,MDL," ~d_L"
      WRITE(18,902) 2000001,MDR," ~d_R"
      WRITE(18,902) 1000002,MUL," ~u_L"
      WRITE(18,902) 2000002,MUR," ~u_R"
      WRITE(18,902) 1000003,MDL," ~s_L"
      WRITE(18,902) 2000003,MDR," ~s_R"
      WRITE(18,902) 1000004,MUL," ~c_L"
      WRITE(18,902) 2000004,MUR," ~c_R"
      WRITE(18,902) 1000005,MSB1," ~b_1"
      WRITE(18,902) 2000005,MSB2," ~b_2"
      WRITE(18,902) 1000006,MST1," ~t_1"
      WRITE(18,902) 2000006,MST2," ~t_2"
      WRITE(18,902) 1000011,MLL," ~e_L"
      WRITE(18,902) 2000011,MLR," ~e_R"
      WRITE(18,902) 1000012,MNL," ~nue_L"
      WRITE(18,902) 1000013,MSMU1," ~mu_L"
      WRITE(18,902) 2000013,MSMU2," ~mu_R"
      WRITE(18,902) 1000014,MNL," ~numu_L"
      WRITE(18,902) 1000015,MSL1," ~tau_1"
      WRITE(18,902) 2000015,MSL2," ~tau_2"
      WRITE(18,902) 1000016,MSNT," ~nutau_L"
      WRITE(18,902) 1000021,MGL," ~g"
      WRITE(18,902) 1000022,MNEU(1),"neutralino(1)"
      WRITE(18,902) 1000023,MNEU(2),"neutralino(2)"
      WRITE(18,902) 1000025,MNEU(3),"neutralino(3)"
      WRITE(18,902) 1000035,MNEU(4),"neutralino(4)"
      WRITE(18,902) 1000045,MNEU(5),"neutralino(5)"
      WRITE(18,902) 1000024,MCHA(1),"chargino(1)"
      WRITE(18,902) 1000037,MCHA(2),"chargino(2)"
      IF(GRFLAG.EQ.1)WRITE(18,902) 1000039,M32,"gravitino"
      WRITE(18,899) "# "

      WRITE(18,907) "BLOCK HMIX Q=",DSQRT(QSTSB),
     .    " # (STOP/SBOTTOM MASSES)"
      WRITE(18,901) 1,MUQ,"MUEFF"
      WRITE(18,901) 2,TANBQ,"TAN(BETA)"
      WRITE(18,901) 3,DSQRT(2d0*(H1Q**2+H2Q**2)),"V(Q)"
      WRITE(18,901) 4,PAR(23)**2,"MA^2"
      WRITE(18,901) 5,PAR(24)**2,"MP^2"

      WRITE(18,899) "# "
      WRITE(18,899) "# 3*3 Higgs mixing"
      WRITE(18,899) "BLOCK NMHMIX"
      WRITE(18,903) 1,1,SCOMP(1,2),"S_(1,1)"
      WRITE(18,903) 1,2,SCOMP(1,1),"S_(1,2)"
      WRITE(18,903) 1,3,SCOMP(1,3),"S_(1,3)"
      WRITE(18,903) 2,1,SCOMP(2,2),"S_(2,1)"
      WRITE(18,903) 2,2,SCOMP(2,1),"S_(2,2)"
      WRITE(18,903) 2,3,SCOMP(2,3),"S_(2,3)"
      WRITE(18,903) 3,1,SCOMP(3,2),"S_(3,1)"
      WRITE(18,903) 3,2,SCOMP(3,1),"S_(3,2)"
      WRITE(18,903) 3,3,SCOMP(3,3),"S_(3,3)"

      WRITE(18,899) "# "
      WRITE(18,899) "# 3*3 Pseudoscalar Higgs mixing"
      WRITE(18,899) "BLOCK NMAMIX"
      WRITE(18,903) 1,1,SINB*PCOMP(1,1),"P_(1,1)"
      WRITE(18,903) 1,2,COSB*PCOMP(1,1),"P_(1,2)"
      WRITE(18,903) 1,3,PCOMP(1,2),"P_(1,3)"
      WRITE(18,903) 2,1,SINB*PCOMP(2,1),"P_(2,1)"
      WRITE(18,903) 2,2,COSB*PCOMP(2,1),"P_(2,2)"
      WRITE(18,903) 2,3,PCOMP(2,2),"P_(2,3)"

      SST=DSQRT(1-CST**2)
      SSB=DSQRT(1-CSB**2)
      SSL=DSQRT(1-CSL**2)

      WRITE(18,899) "# "
      WRITE(18,899) "# 3rd generation sfermion mixing"
      WRITE(18,899) "BLOCK STOPMIX  # Stop mixing matrix"
      WRITE(18,903) 1,1,CST,"Rst_(1,1)"
      WRITE(18,903) 1,2,SST,"Rst_(1,2)"
      WRITE(18,903) 2,1,-SST,"Rst_(2,1)"
      WRITE(18,903) 2,2,CST,"Rst_(2,2)"
      WRITE(18,899) "BLOCK SBOTMIX  # Sbottom mixing matrix"
      WRITE(18,903) 1,1,CSB,"Rsb_(1,1)"
      WRITE(18,903) 1,2,SSB,"Rsb_(1,2)"
      WRITE(18,903) 2,1,-SSB,"Rsb_(2,1)"
      WRITE(18,903) 2,2,CSB,"Rsb_(2,2)"
      WRITE(18,899) "BLOCK STAUMIX  # Stau mixing matrix"
      WRITE(18,903) 1,1,CSL,"Rsl_(1,1)"
      WRITE(18,903) 1,2,SSL,"Rsl_(1,2)"
      WRITE(18,903) 2,1,-SSL,"Rsl_(2,1)"
      WRITE(18,903) 2,2,CSL,"Rsl_(2,2)"

      WRITE(18,899) "# "
      WRITE(18,899) "# Gaugino-Higgsino mixing"
      WRITE(18,899) "BLOCK NMNMIX  # 5*5 Neutralino Mixing Matrix"
      WRITE(18,903) 1,1,NEU(1,1),"N_(1,1)"
      WRITE(18,903) 1,2,NEU(1,2),"N_(1,2)"
      WRITE(18,903) 1,3,NEU(1,4),"N_(1,3)"
      WRITE(18,903) 1,4,NEU(1,3),"N_(1,4)"
      WRITE(18,903) 1,5,NEU(1,5),"N_(1,5)"
      WRITE(18,903) 2,1,NEU(2,1),"N_(2,1)"
      WRITE(18,903) 2,2,NEU(2,2),"N_(2,2)"
      WRITE(18,903) 2,3,NEU(2,4),"N_(2,3)"
      WRITE(18,903) 2,4,NEU(2,3),"N_(2,4)"
      WRITE(18,903) 2,5,NEU(2,5),"N_(2,5)"
      WRITE(18,903) 3,1,NEU(3,1),"N_(3,1)"
      WRITE(18,903) 3,2,NEU(3,2),"N_(3,2)"
      WRITE(18,903) 3,3,NEU(3,4),"N_(3,3)"
      WRITE(18,903) 3,4,NEU(3,3),"N_(3,4)"
      WRITE(18,903) 3,5,NEU(3,5),"N_(3,5)"
      WRITE(18,903) 4,1,NEU(4,1),"N_(4,1)"
      WRITE(18,903) 4,2,NEU(4,2),"N_(4,2)"
      WRITE(18,903) 4,3,NEU(4,4),"N_(4,3)"
      WRITE(18,903) 4,4,NEU(4,3),"N_(4,4)"
      WRITE(18,903) 4,5,NEU(4,5),"N_(4,5)"
      WRITE(18,903) 5,1,NEU(5,1),"N_(5,1)"
      WRITE(18,903) 5,2,NEU(5,2),"N_(5,2)"
      WRITE(18,903) 5,3,NEU(5,4),"N_(5,3)"
      WRITE(18,903) 5,4,NEU(5,3),"N_(5,4)"
      WRITE(18,903) 5,5,NEU(5,5),"N_(5,5)"

      WRITE(18,899) "# "
      WRITE(18,899) "BLOCK UMIX  # Chargino U Mixing Matrix"
      WRITE(18,903) 1,1,U(1,1),"U_(1,1)"
      WRITE(18,903) 1,2,U(1,2),"U_(1,2)"
      WRITE(18,903) 2,1,U(2,1),"U_(2,1)"
      WRITE(18,903) 2,2,U(2,2),"U_(2,2)"

      WRITE(18,899) "# "
      WRITE(18,899) "BLOCK VMIX  # Chargino V Mixing Matrix"
      WRITE(18,903) 1,1,V(1,1),"V_(1,1)"
      WRITE(18,903) 1,2,V(1,2),"V_(1,2)"
      WRITE(18,903) 2,1,V(2,1),"V_(2,1)"
      WRITE(18,903) 2,2,V(2,2),"V_(2,2)"
      WRITE(18,899) "# "

      WRITE(18,899) "# GAUGE AND YUKAWA COUPLINGS AT THE SUSY SCALE"
      WRITE(18,907) "BLOCK GAUGE Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(18,901) 1,DSQRT(G1S),"g1(Q,DR_bar)"
      WRITE(18,901) 2,DSQRT(G2S),"g2(Q,DR_bar)"
      WRITE(18,901) 3,DSQRT(G3S),"g3(Q,DR_bar)"

      WRITE(18,907) "BLOCK YU Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(18,903) 3,3,HTOPS,"HTOP(Q,DR_bar)"
      WRITE(18,907) "BLOCK YD Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(18,903) 3,3,HBOTS,"HBOT(Q,DR_bar)"
      WRITE(18,907) "BLOCK YE Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(18,903) 3,3,HTAUS,"HTAU(Q,DR_bar)"

      WRITE(18,899) "# "
      WRITE(18,899) "# SOFT TRILINEAR COUPLINGS AT THE SUSY SCALE"
      WRITE(18,899) "# (BOTH SLHA1 AND SLHA2 FORMAT)"
      WRITE(18,907) "BLOCK AU Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(18,903) 3,3,PAR(12),"ATOP"
      WRITE(18,907) "BLOCK TU Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(18,903) 3,3,PAR(12),"ATOP"
      WRITE(18,907) "BLOCK AD Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(18,903) 3,3,PAR(13),"ABOT"
      WRITE(18,907) "BLOCK TD Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(18,903) 3,3,PAR(13),"ABOT"
      WRITE(18,907) "BLOCK AE Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(18,903) 2,2,PAR(25),"AMUON"
      WRITE(18,903) 3,3,PAR(14),"ATAU"
      WRITE(18,907) "BLOCK TE Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(18,903) 2,2,PAR(25),"AMUON"
      WRITE(18,903) 3,3,PAR(14),"ATAU"

      WRITE(18,899) "# "
      WRITE(18,899) "# SOFT MASSES AT THE SUSY SCALE"
      WRITE(18,907) "BLOCK MSOFT Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(18,901) 1,PAR(20),"M1"
      WRITE(18,901) 2,PAR(21),"M2"
      WRITE(18,901) 3,PAR(22),"M3"
      WRITE(18,901) 21,MHDS,"M_HD^2"
      WRITE(18,901) 22,MHUS,"M_HU^2"
      WRITE(18,901) 31,PAR(18)/DSQRT(DABS(PAR(18))),"M_eL"
      WRITE(18,901) 32,PAR(18)/DSQRT(DABS(PAR(18))),"M_muL"
      WRITE(18,901) 33,PAR(10)/DSQRT(DABS(PAR(10))),"M_tauL"
      WRITE(18,901) 34,PAR(19)/DSQRT(DABS(PAR(19))),"M_eR"
      WRITE(18,901) 35,PAR(19)/DSQRT(DABS(PAR(19))),"M_muR"
      WRITE(18,901) 36,PAR(11)/DSQRT(DABS(PAR(11))),"M_tauR"
      WRITE(18,901) 41,PAR(15)/DSQRT(DABS(PAR(15))),"M_q1L"
      WRITE(18,901) 42,PAR(15)/DSQRT(DABS(PAR(15))),"M_q2L"
      WRITE(18,901) 43,PAR(7)/DSQRT(DABS(PAR(7))),"M_q3L"
      WRITE(18,901) 44,PAR(16)/DSQRT(DABS(PAR(16))),"M_uR"
      WRITE(18,901) 45,PAR(16)/DSQRT(DABS(PAR(16))),"M_cR"
      WRITE(18,901) 46,PAR(8)/DSQRT(DABS(PAR(8))),"M_tR"
      WRITE(18,901) 47,PAR(17)/DSQRT(DABS(PAR(17))),"M_dR"
      WRITE(18,901) 48,PAR(17)/DSQRT(DABS(PAR(17))),"M_sR"
      WRITE(18,901) 49,PAR(9)/DSQRT(DABS(PAR(9))),"M_bR"

      WRITE(18,899) "# "
      WRITE(18,899) "# NMSSM SPECIFIC PARAMETERS THE SUSY SCALE"
      WRITE(18,907) "BLOCK NMSSMRUN Q=",DSQRT(Q2)," # (SUSY SCALE)"
      WRITE(18,901) 1,PAR(1),"LAMBDA(Q,DR_bar)"
      WRITE(18,901) 2,PAR(2),"KAPPA(Q,DR_bar)"
      WRITE(18,901) 3,PAR(5),"ALAMBDA"
      WRITE(18,901) 4,PAR(6),"AKAPPA"
      WRITE(18,901) 5,PAR(4),"MUEFF"
      WRITE(18,901) 6,XIFSUSY,"XIF"
      WRITE(18,901) 7,XISSUSY,"XIS"
      WRITE(18,901) 8,MUPSUSY,"MUP"
      WRITE(18,901) 9,MSPSUSY,"MSP"
      WRITE(18,901) 10,MSS,"MS^2"
      WRITE(18,920) 12,M3HSUSY,"M3H"
      WRITE(18,899) "# "
      WRITE(18,899) "BLOCK MSQ2  # Soft l.h. squark masses squared"
      WRITE(18,903) 1,1,PAR(15),"M_q1L"
      WRITE(18,903) 2,2,PAR(15),"M_q2L"
      WRITE(18,903) 3,3,PAR(7),"M_q3L"
      WRITE(18,899) "# "
      WRITE(18,899) "BLOCK MSU2  # Soft r.h. up-squark masses squared"
      WRITE(18,903) 1,1,PAR(16),"M_uR"
      WRITE(18,903) 2,2,PAR(16),"M_cR"
      WRITE(18,903) 3,3,PAR(8),"M_tR"
      WRITE(18,899) "# "
      WRITE(18,899) "BLOCK MSD2  # Soft r.h. down-squark masses squared"
      WRITE(18,903) 1,1,PAR(17),"M_dR"
      WRITE(18,903) 2,2,PAR(17),"M_sR"
      WRITE(18,903) 3,3,PAR(9),"M_bR"
      WRITE(18,899) "# "
      WRITE(18,899) "BLOCK MSL2  # Soft l.h. slepton masses squared"
      WRITE(18,903) 1,1,PAR(18),"M_eL"
      WRITE(18,903) 2,2,PAR(18),"M_muL"
      WRITE(18,903) 3,3,PAR(10),"M_tauL"
      WRITE(18,899) "# "
      WRITE(18,899) "BLOCK MSE2  # Soft r.h. slepton masses squared"
      WRITE(18,903) 1,1,PAR(19),"M_eR"
      WRITE(18,903) 2,2,PAR(19),"M_muR"
      WRITE(18,903) 3,3,PAR(11),"M_tauR"
      WRITE(18,899) "# "
      WRITE(18,899) "BLOCK USQMIX  # Elements of 6x6 up-squark matrix"
      WRITE(18,903) 1,1,1.d0,"R_u_11"
      WRITE(18,903) 2,2,1.d0,"R_u_22"
      WRITE(18,903) 3,3,CST,"R_u_33"
      WRITE(18,903) 3,6,SST,"R_u_36"
      WRITE(18,903) 4,4,1.d0,"R_u_44"
      WRITE(18,903) 5,5,1.d0,"R_u_55"
      WRITE(18,903) 6,6,CST,"R_u_66"
      WRITE(18,903) 6,3,-SST,"R_u_63"
      WRITE(18,899) "# "
      WRITE(18,899) "BLOCK DSQMIX  # Elements of 6x6 down-squark matrix"
      WRITE(18,903) 1,1,1.d0,"R_d_11"
      WRITE(18,903) 2,2,1.d0,"R_d_22"
      WRITE(18,903) 3,3,CSB,"R_d_33"
      WRITE(18,903) 3,6,SSB,"R_u_36"
      WRITE(18,903) 4,4,1.d0,"R_d_44"
      WRITE(18,903) 5,5,1.d0,"R_d_55"
      WRITE(18,903) 6,6,CSB,"R_d_66"
      WRITE(18,903) 6,3,-SSB,"R_u36"
      WRITE(18,899) "# "
      WRITE(18,899) "BLOCK SELMIX  # Elements of 6x6 ch. slepton matrix"
      WRITE(18,903) 1,1,1.d0,"R_e_11"
      WRITE(18,903) 2,2,1.d0,"R_e_22"
      WRITE(18,903) 3,3,CSL,"R_e_33"
      WRITE(18,903) 3,6,SSL,"R_e_36"
      WRITE(18,903) 4,4,1.d0,"R_e_44"
      WRITE(18,903) 5,5,1.d0,"R_e_55"
      WRITE(18,903) 6,6,CSL,"R_e_66"
      WRITE(18,903) 6,3,-SSL,"R_e_63"
      WRITE(18,899) "# "

      WRITE(18,899) "#           PDG          Width"
      WRITE(18,904) PDGH1,WIDTH(1),"Lightest neutral Higgs scalar"
      IF(BRJJ(1).GT.0d0)
     .  WRITE(18,905) BRJJ(1),2,21,21,"BR(H_1 -> hadrons)"
      IF(BREE(1).GT.0d0)
     .  WRITE(18,905) BREE(1),2,11,-11,"BR(H_1 -> e- e+)"
      IF(BRMM(1).GT.0d0)
     .  WRITE(18,905) BRMM(1),2,13,-13,"BR(H_1 -> muon muon)"
      IF(BRLL(1).GT.0d0)
     .  WRITE(18,905) BRLL(1),2,15,-15,"BR(H_1 -> tau tau)"
      IF(BRCC(1).GT.0d0)
     .  WRITE(18,905) BRCC(1),2,4,-4,"BR(H_1 -> c cbar)"
      IF(BRBB(1).GT.0d0)
     .  WRITE(18,905) BRBB(1),2,5,-5,"BR(H_1 -> b bbar)"
      IF(BRTT(1).GT.0d0)
     .  WRITE(18,905) BRTT(1),2,6,-6,"BR(H_1 -> t tbar)"
      IF(BRWW(1).GT.0d0)
     .  WRITE(18,905) BRWW(1),2,24,-24,"BR(H_1 -> W+ W-)"
      IF(BRZZ(1).GT.0d0)
     .  WRITE(18,905) BRZZ(1),2,23,23,"BR(H_1 -> Z Z)"
      IF(BRGG(1).GT.0d0)
     .  WRITE(18,905) BRGG(1),2,22,22,"BR(H_1 -> gamma gamma)"
      IF(BRZG(1).GT.0d0)
     .  WRITE(18,905) BRZG(1),2,23,22,"BR(H_1 -> Z gamma)"
      IF(BRHAA(1,1).GT.0d0)
     .  WRITE(18,905) BRHAA(1,1),2,36,36,"BR(H_1 -> A_1 A_1)"
      IF(BRHAA(1,2).GT.0d0)
     .  WRITE(18,905) BRHAA(1,2),2,36,46,"BR(H_1 -> A_1 A_2)"
      IF(BRHAA(1,3).GT.0d0)
     .  WRITE(18,905) BRHAA(1,3),2,46,46,"BR(H_1 -> A_2 A_2)"
      IF(BRHAZ(1,1).GT.0d0)
     .  WRITE(18,905) BRHAZ(1,1),2,23,36,"BR(H_1 -> A_1 Z)"
      IF(BRHAZ(1,2).GT.0d0)
     .  WRITE(18,905) BRHAZ(1,2),2,23,46,"BR(H_1 -> A_2 Z)"
      IF(BRHCHC(1).GT.0d0)
     .  WRITE(18,905) BRHCHC(1),2,37,-37,"BR(H_1 -> H+ H-)"
      IF(BRHCW(1).GT.0d0)
     .  WRITE(18,905) BRHCW(1),2,24,-37,"BR(H_1 -> W+ H-)"
      IF(BRHCW(1).GT.0d0)
     .  WRITE(18,905) BRHCW(1),2,-24,37,"BR(H_1 -> W- H+)"
      IF(BRNEU(1,1,1).GT.0d0)
     .  WRITE(18,905) BRNEU(1,1,1),2,1000022,1000022,
     .    "BR(H_1 -> neu_1 neu_1)"
      IF(BRNEU(1,1,2).GT.0d0)
     .  WRITE(18,905) BRNEU(1,1,2),2,1000022,1000023,
     .    "BR(H_1 -> neu_1 neu_2)"
      IF(BRNEU(1,1,3).GT.0d0)
     .  WRITE(18,905) BRNEU(1,1,3),2,1000022,1000025,
     .    "BR(H_1 -> neu_1 neu_3)"
      IF(BRNEU(1,1,4).GT.0d0)
     .  WRITE(18,905) BRNEU(1,1,4),2,1000022,1000035,
     .    "BR(H_1 -> neu_1 neu_4)"
      IF(BRNEU(1,1,5).GT.0d0)
     .  WRITE(18,905) BRNEU(1,1,5),2,1000022,1000045,
     .    "BR(H_1 -> neu_1 neu_5)"
      IF(BRNEU(1,2,2).GT.0d0)
     .  WRITE(18,905) BRNEU(1,2,2),2,1000023,1000023,
     .    "BR(H_1 -> neu_2 neu_2)"
      IF(BRNEU(1,2,3).GT.0d0)
     .  WRITE(18,905) BRNEU(1,2,3),2,1000023,1000025,
     .    "BR(H_1 -> neu_2 neu_3)"
      IF(BRNEU(1,2,4).GT.0d0)
     .  WRITE(18,905) BRNEU(1,2,4),2,1000023,1000035,
     .    "BR(H_1 -> neu_2 neu_4)"
      IF(BRNEU(1,2,5).GT.0d0)
     .  WRITE(18,905) BRNEU(1,2,5),2,1000023,1000045,
     .    "BR(H_1 -> neu_2 neu_5)"
      IF(BRNEU(1,3,3).GT.0d0)
     .  WRITE(18,905) BRNEU(1,3,3),2,1000025,1000025,
     .    "BR(H_1 -> neu_3 neu_3)"
      IF(BRNEU(1,3,4).GT.0d0)
     .  WRITE(18,905) BRNEU(1,3,4),2,1000025,1000035,
     .    "BR(H_1 -> neu_3 neu_4)"
      IF(BRNEU(1,3,5).GT.0d0)
     .  WRITE(18,905) BRNEU(1,3,5),2,1000025,1000045,
     .    "BR(H_1 -> neu_3 neu_5)"
      IF(BRNEU(1,4,4).GT.0d0)
     .  WRITE(18,905) BRNEU(1,4,4),2,1000035,1000035,
     .    "BR(H_1 -> neu_4 neu_4)"
      IF(BRNEU(1,4,5).GT.0d0)
     .  WRITE(18,905) BRNEU(1,4,5),2,1000035,1000045,
     .    "BR(H_1 -> neu_4 neu_5)"
      IF(BRNEU(1,5,5).GT.0d0)
     .  WRITE(18,905) BRNEU(1,5,5),2,1000045,1000045,
     .    "BR(H_1 -> neu_5 neu_5)"
      IF(BRCHA(1,1).GT.0d0)
     .  WRITE(18,905) BRCHA(1,1),2,1000024,-1000024,
     .    "BR(H_1 -> cha_1 cha_1bar)"
      IF(BRCHA(1,2).GT.0d0)
     .  WRITE(18,905) BRCHA(1,2),2,1000024,-1000037,
     .    "BR(H_1 -> cha_1 cha_2bar)"
      IF(BRCHA(1,2).GT.0d0)
     .  WRITE(18,905) BRCHA(1,2),2,1000037,-1000024,
     .    "BR(H_1 -> cha_2 cha_1bar)"
      IF(BRCHA(1,3).GT.0d0)
     .  WRITE(18,905) BRCHA(1,3),2,1000037,-1000037,
     .    "BR(H_1 -> cha_2 cha_2bar)"
      IF(BRHSQ(1,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,1),2,1000002,-1000002,
     .    "BR(H_1 -> ~u_L ~ubar_L)"
      IF(BRHSQ(1,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,1),2,1000004,-1000004,
     .    "BR(H_1 -> ~c_L ~cbar_L)"
      IF(BRHSQ(1,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,2),2,2000002,-2000002,
     .    "BR(H_1 -> ~u_R ~ubar_R)"
      IF(BRHSQ(1,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,2),2,2000004,-2000004,
     .    "BR(H_1 -> ~c_R ~cbar_R)"
      IF(BRHSQ(1,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,3),2,1000001,-1000001,
     .    "BR(H_1 -> ~d_L ~dbar_L)"
      IF(BRHSQ(1,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,3),2,1000003,-1000003,
     .    "BR(H_1 -> ~s_L ~sbar_L)"
      IF(BRHSQ(1,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,4),2,2000001,-2000001,
     .    "BR(H_1 -> ~d_R ~dbar_R)"
      IF(BRHSQ(1,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,4),2,2000003,-2000003,
     .    "BR(H_1 -> ~s_R ~sbar_R)"
      IF(BRHSQ(1,5).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,5),2,1000006,-1000006,
     .    "BR(H_1 -> ~t_1 ~tbar_1)"
      IF(BRHSQ(1,6).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,6),2,2000006,-2000006,
     .    "BR(H_1 -> ~t_2 ~tbar_2)"
      IF(BRHSQ(1,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,7),2,1000006,-2000006,
     .    "BR(H_1 -> ~t_1 ~tbar_2)"
      IF(BRHSQ(1,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,7),2,2000006,-1000006,
     .    "BR(H_1 -> ~t_2 ~tbar_1)"
      IF(BRHSQ(1,8).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,8),2,1000005,-1000005,
     .    "BR(H_1 -> ~b_1 ~bbar_1)"
      IF(BRHSQ(1,9).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,9),2,2000005,-2000005,
     .    "BR(H_1 -> ~b_2 ~bbar_2)"
      IF(BRHSQ(1,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,10),2,1000005,-2000005,
     .    "BR(H_1 -> ~b_1 ~bbar_2)"
      IF(BRHSQ(1,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(1,10),2,2000005,-1000005,
     .    "BR(H_1 -> ~b_2 ~bbar_1)"
      IF(BRHSL(1,1).GT.0d0)
     .  WRITE(18,905) BRHSL(1,1),2,1000011,-1000011,
     .    "BR(H_1 -> ~e_L ~ebar_L)"
      IF(BRHSL(1,1).GT.0d0)
     .  WRITE(18,905) BRHSL(1,1),2,1000013,-1000013,
     .    "BR(H_1 -> ~mu_L ~mubar_L)"
      IF(BRHSL(1,2).GT.0d0)
     .  WRITE(18,905) BRHSL(1,2),2,2000011,-2000011,
     .    "BR(H_1 -> ~e_R ~ebar_R)"
      IF(BRHSL(1,2).GT.0d0)
     .  WRITE(18,905) BRHSL(1,2),2,2000013,-2000013,
     .    "BR(H_1 -> ~mu_R ~mubarRL)"
      IF(BRHSL(1,3).GT.0d0)
     .  WRITE(18,905) BRHSL(1,3),2,1000012,-1000012,
     .    "BR(H_1 -> ~nu_e_L ~nu_ebar_L)"
      IF(BRHSL(1,3).GT.0d0)
     .  WRITE(18,905) BRHSL(1,3),2,1000014,-1000014,
     .    "BR(H_1 -> ~nu_mu_L ~nu_mubar_L)"
      IF(BRHSL(1,4).GT.0d0)
     .  WRITE(18,905) BRHSL(1,4),2,1000015,-1000015,
     .    "BR(H_1 -> ~tau_1 ~taubar_1)"
      IF(BRHSL(1,5).GT.0d0)
     .  WRITE(18,905) BRHSL(1,5),2,2000015,-2000015,
     .    "BR(H_1 -> ~tau_2 ~taubar_2)"
      IF(BRHSL(1,6).GT.0d0)
     .  WRITE(18,905) BRHSL(1,6),2,1000015,-2000015,
     .    "BR(H_1 -> ~tau_1 ~taubar_2)"
      IF(BRHSL(1,6).GT.0d0)
     .  WRITE(18,905) BRHSL(1,6),2,2000015,-1000015,
     .    "BR(H_1 -> ~tau_2 ~taubar_1)"
      IF(BRHSL(1,7).GT.0d0)
     .  WRITE(18,905) BRHSL(1,7),2,1000016,-1000016,
     .    "BR(H_1 -> ~nu_tau_L ~nu_taubar_L)"

      WRITE(18,904) PDGH2,WIDTH(2),"2nd neutral Higgs scalar"
      IF(BRJJ(2).GT.0d0)
     .  WRITE(18,905) BRJJ(2),2,21,21,"BR(H_2 -> hadrons)"
      IF(BREE(2).GT.0d0)
     .  WRITE(18,905) BREE(2),2,11,-11,"BR(H_2 -> e- e+)"
      IF(BRMM(2).GT.0d0)
     .  WRITE(18,905) BRMM(2),2,13,-13,"BR(H_2 -> muon muon)"
      IF(BRLL(2).GT.0d0)
     .  WRITE(18,905) BRLL(2),2,15,-15,"BR(H_2 -> tau tau)"
      IF(BRCC(2).GT.0d0)
     .  WRITE(18,905) BRCC(2),2,4,-4,"BR(H_2 -> c cbar)"
      IF(BRBB(2).GT.0d0)
     .  WRITE(18,905) BRBB(2),2,5,-5,"BR(H_2 -> b bbar)"
      IF(BRTT(2).GT.0d0)
     .  WRITE(18,905) BRTT(2),2,6,-6,"BR(H_2 -> t tbar)"
      IF(BRWW(2).GT.0d0)
     .  WRITE(18,905) BRWW(2),2,24,-24,"BR(H_2 -> W+ W-)"
      IF(BRZZ(2).GT.0d0)
     .  WRITE(18,905) BRZZ(2),2,23,23,"BR(H_2 -> Z Z)"
      IF(BRGG(2).GT.0d0)
     .  WRITE(18,905) BRGG(2),2,22,22,"BR(H_2 -> gamma gamma)"
      IF(BRZG(2).GT.0d0)
     .  WRITE(18,905) BRZG(2),2,23,22,"BR(H_2 -> Z gamma)"
      IF(BRHHH(1).GT.0d0)
     .  WRITE(18,905) BRHHH(1),2,PDGH1,PDGH1,"BR(H_2 -> H_1 H_1)"
      IF(BRHAA(2,1).GT.0d0)
     .  WRITE(18,905) BRHAA(2,1),2,36,36,"BR(H_2 -> A_1 A_1)"
      IF(BRHAA(2,2).GT.0d0)
     .  WRITE(18,905) BRHAA(2,2),2,36,46,"BR(H_2 -> A_1 A_2)"
      IF(BRHAA(2,3).GT.0d0)
     .  WRITE(18,905) BRHAA(2,3),2,46,46,"BR(H_2 -> A_2 A_2)"
      IF(BRHAZ(2,1).GT.0d0)
     .  WRITE(18,905) BRHAZ(2,1),2,23,36,"BR(H_2 -> A_1 Z)"
      IF(BRHAZ(2,2).GT.0d0)
     .  WRITE(18,905) BRHAZ(2,2),2,23,46,"BR(H_2 -> A_2 Z)"
      IF(BRHCHC(2).GT.0d0)
     .  WRITE(18,905) BRHCHC(2),2,37,-37,"BR(H_2 -> H+ H-)"
      IF(BRHCW(2).GT.0d0)
     .  WRITE(18,905) BRHCW(2),2,24,-37,"BR(H_2 -> W+ H-)"
      IF(BRHCW(2).GT.0d0)
     .  WRITE(18,905) BRHCW(2),2,-24,37,"BR(H_2 -> W- H+)"
      IF(BRNEU(2,1,1).GT.0d0)
     .  WRITE(18,905) BRNEU(2,1,1),2,1000022,1000022,
     .    "BR(H_2 -> neu_1 neu_1)"
      IF(BRNEU(2,1,2).GT.0d0)
     .  WRITE(18,905) BRNEU(2,1,2),2,1000022,1000023,
     .    "BR(H_2 -> neu_1 neu_2)"
      IF(BRNEU(2,1,3).GT.0d0)
     .  WRITE(18,905) BRNEU(2,1,3),2,1000022,1000025,
     .    "BR(H_2 -> neu_1 neu_3)"
      IF(BRNEU(2,1,4).GT.0d0)
     .  WRITE(18,905) BRNEU(2,1,4),2,1000022,1000035,
     .    "BR(H_2 -> neu_1 neu_4)"
      IF(BRNEU(2,1,5).GT.0d0)
     .  WRITE(18,905) BRNEU(2,1,5),2,1000022,1000045,
     .    "BR(H_2 -> neu_1 neu_5)"
      IF(BRNEU(2,2,2).GT.0d0)
     .  WRITE(18,905) BRNEU(2,2,2),2,1000023,1000023,
     .    "BR(H_2 -> neu_2 neu_2)"
      IF(BRNEU(2,2,3).GT.0d0)
     .  WRITE(18,905) BRNEU(2,2,3),2,1000023,1000025,
     .    "BR(H_2 -> neu_2 neu_3)"
      IF(BRNEU(2,2,4).GT.0d0)
     .  WRITE(18,905) BRNEU(2,2,4),2,1000023,1000035,
     .    "BR(H_2 -> neu_2 neu_4)"
      IF(BRNEU(2,2,5).GT.0d0)
     .  WRITE(18,905) BRNEU(2,2,5),2,1000023,1000045,
     .    "BR(H_2 -> neu_2 neu_5)"
      IF(BRNEU(2,3,3).GT.0d0)
     .  WRITE(18,905) BRNEU(2,3,3),2,1000025,1000025,
     .    "BR(H_2 -> neu_3 neu_3)"
      IF(BRNEU(2,3,4).GT.0d0)
     .  WRITE(18,905) BRNEU(2,3,4),2,1000025,1000035,
     .    "BR(H_2 -> neu_3 neu_4)"
      IF(BRNEU(2,3,5).GT.0d0)
     .  WRITE(18,905) BRNEU(2,3,5),2,1000025,1000045,
     .    "BR(H_2 -> neu_3 neu_5)"
      IF(BRNEU(2,4,4).GT.0d0)
     .  WRITE(18,905) BRNEU(2,4,4),2,1000035,1000035,
     .    "BR(H_2 -> neu_4 neu_4)"
      IF(BRNEU(2,4,5).GT.0d0)
     .  WRITE(18,905) BRNEU(2,4,5),2,1000035,1000045,
     .    "BR(H_2 -> neu_4 neu_5)"
      IF(BRNEU(2,5,5).GT.0d0)
     .  WRITE(18,905) BRNEU(2,5,5),2,1000045,1000045,
     .    "BR(H_2 -> neu_5 neu_5)"
      IF(BRCHA(2,1).GT.0d0)
     .  WRITE(18,905) BRCHA(2,1),2,1000024,-1000024,
     .    "BR(H_2 -> cha_1 cha_1bar)"
      IF(BRCHA(2,2).GT.0d0)
     .  WRITE(18,905) BRCHA(2,2),2,1000024,-1000037,
     .    "BR(H_2 -> cha_1 cha_2bar)"
      IF(BRCHA(2,2).GT.0d0)
     .  WRITE(18,905) BRCHA(2,2),2,1000037,-1000024,
     .    "BR(H_2 -> cha_2 cha_1bar)"
      IF(BRCHA(2,3).GT.0d0)
     .  WRITE(18,905) BRCHA(2,3),2,1000037,-1000037,
     .    "BR(H_2 -> cha_2 cha_2bar)"
      IF(BRHSQ(2,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,1),2,1000002,-1000002,
     .    "BR(H_2 -> ~u_L ~ubar_L)"
      IF(BRHSQ(2,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,1),2,1000004,-1000004,
     .    "BR(H_2 -> ~c_L ~cbar_L)"
      IF(BRHSQ(2,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,2),2,2000002,-2000002,
     .    "BR(H_2 -> ~u_R ~ubar_R)"
      IF(BRHSQ(2,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,2),2,2000004,-2000004,
     .    "BR(H_2 -> ~c_R ~cbar_R)"
      IF(BRHSQ(2,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,3),2,1000001,-1000001,
     .    "BR(H_2 -> ~d_L ~dbar_L)"
      IF(BRHSQ(2,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,3),2,1000003,-1000003,
     .    "BR(H_2 -> ~s_L ~sbar_L)"
      IF(BRHSQ(2,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,4),2,2000001,-2000001,
     .    "BR(H_2 -> ~d_R ~dbar_R)"
      IF(BRHSQ(2,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,4),2,2000003,-2000003,
     .    "BR(H_2 -> ~s_R ~sbar_R)"
      IF(BRHSQ(2,5).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,5),2,1000006,-1000006,
     .    "BR(H_2 -> ~t_1 ~tbar_1)"
      IF(BRHSQ(2,6).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,6),2,2000006,-2000006,
     .    "BR(H_2 -> ~t_2 ~tbar_2)"
      IF(BRHSQ(2,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,7),2,1000006,-2000006,
     .    "BR(H_2 -> ~t_1 ~tbar_2)"
      IF(BRHSQ(2,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,7),2,2000006,-1000006,
     .    "BR(H_2 -> ~t_2 ~tbar_1)"
      IF(BRHSQ(2,8).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,8),2,1000005,-1000005,
     .    "BR(H_2 -> ~b_1 ~bbar_1)"
      IF(BRHSQ(2,9).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,9),2,2000005,-2000005,
     .    "BR(H_2 -> ~b_2 ~bbar_2)"
      IF(BRHSQ(2,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,10),2,1000005,-2000005,
     .    "BR(H_2 -> ~b_1 ~bbar_2)"
      IF(BRHSQ(2,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(2,10),2,2000005,-1000005,
     .    "BR(H_2 -> ~b_2 ~bbar_1)"
      IF(BRHSL(2,1).GT.0d0)
     .  WRITE(18,905) BRHSL(2,1),2,1000011,-1000011,
     .    "BR(H_2 -> ~e_L ~ebar_L)"
      IF(BRHSL(2,1).GT.0d0)
     .  WRITE(18,905) BRHSL(2,1),2,1000013,-1000013,
     .    "BR(H_2 -> ~mu_L ~mubar_L)"
      IF(BRHSL(2,2).GT.0d0)
     .  WRITE(18,905) BRHSL(2,2),2,2000011,-2000011,
     .    "BR(H_2 -> ~e_R ~ebar_R)"
      IF(BRHSL(2,2).GT.0d0)
     .  WRITE(18,905) BRHSL(2,2),2,2000013,-2000013,
     .    "BR(H_2 -> ~mu_R ~mubarRL)"
      IF(BRHSL(2,3).GT.0d0)
     .  WRITE(18,905) BRHSL(2,3),2,1000012,-1000012,
     .    "BR(H_2 -> ~nu_e_L ~nu_ebar_L)"
      IF(BRHSL(2,3).GT.0d0)
     .  WRITE(18,905) BRHSL(2,3),2,1000014,-1000014,
     .    "BR(H_2 -> ~nu_mu_L ~nu_mubar_L)"
      IF(BRHSL(2,4).GT.0d0)
     .  WRITE(18,905) BRHSL(2,4),2,1000015,-1000015,
     .    "BR(H_2 -> ~tau_1 ~taubar_1)"
      IF(BRHSL(2,5).GT.0d0)
     .  WRITE(18,905) BRHSL(2,5),2,2000015,-2000015,
     .    "BR(H_2 -> ~tau_2 ~taubar_2)"
      IF(BRHSL(2,6).GT.0d0)
     .  WRITE(18,905) BRHSL(2,6),2,1000015,-2000015,
     .    "BR(H_2 -> ~tau_1 ~taubar_2)"
      IF(BRHSL(2,6).GT.0d0)
     .  WRITE(18,905) BRHSL(2,6),2,2000015,-1000015,
     .    "BR(H_2 -> ~tau_2 ~taubar_1)"
      IF(BRHSL(2,7).GT.0d0)
     .  WRITE(18,905) BRHSL(2,7),2,1000016,-1000016,
     .    "BR(H_2 -> ~nu_tau_L ~nu_taubar_L)"

      WRITE(18,904) 45,WIDTH(3),"3rd neutral Higgs scalar"
      IF(BRJJ(3).GT.0d0)
     .  WRITE(18,905) BRJJ(3),2,21,21,"BR(H_3 -> hadrons)"
      IF(BREE(3).GT.0d0)
     .  WRITE(18,905) BREE(3),2,11,-11,"BR(H_3 -> e- e+)"
      IF(BRMM(3).GT.0d0)
     .  WRITE(18,905) BRMM(3),2,13,-13,"BR(H_3 -> muon muon)"
      IF(BRLL(3).GT.0d0)
     .  WRITE(18,905) BRLL(3),2,15,-15,"BR(H_3 -> tau tau)"
      IF(BRCC(3).GT.0d0)
     .  WRITE(18,905) BRCC(3),2,4,-4,"BR(H_3 -> c cbar)"
      IF(BRBB(3).GT.0d0)
     .  WRITE(18,905) BRBB(3),2,5,-5,"BR(H_3 -> b bbar)"
      IF(BRTT(3).GT.0d0)
     .  WRITE(18,905) BRTT(3),2,6,-6,"BR(H_3 -> t tbar)"
      IF(BRWW(3).GT.0d0)
     .  WRITE(18,905) BRWW(3),2,24,-24,"BR(H_3 -> W+ W-)"
      IF(BRZZ(3).GT.0d0)
     .  WRITE(18,905) BRZZ(3),2,23,23,"BR(H_3 -> Z Z)"
      IF(BRGG(3).GT.0d0)
     .  WRITE(18,905) BRGG(3),2,22,22,"BR(H_3 -> gamma gamma)"
      IF(BRZG(3).GT.0d0)
     .  WRITE(18,905) BRZG(3),2,23,22,"BR(H_3 -> Z gamma)"
      IF(BRHHH(2).GT.0d0)
     .  WRITE(18,905) BRHHH(2),2,PDGH1,PDGH1,"BR(H_3 -> H_1 H_1)"
      IF(BRHHH(3).GT.0d0)
     .  WRITE(18,905) BRHHH(3),2,PDGH1,PDGH2,"BR(H_3 -> H_1 H_2)"
      IF(BRHHH(4).GT.0d0)
     .  WRITE(18,905) BRHHH(4),2,PDGH2,PDGH2,"BR(H_3 -> H_2 H_2)"
      IF(BRHAA(3,1).GT.0d0)
     .  WRITE(18,905) BRHAA(3,1),2,36,36,"BR(H_3 -> A_1 A_1)"
      IF(BRHAA(3,2).GT.0d0)
     .  WRITE(18,905) BRHAA(3,2),2,36,46,"BR(H_3 -> A_1 A_2)"
      IF(BRHAA(3,3).GT.0d0)
     .  WRITE(18,905) BRHAA(3,3),2,46,46,"BR(H_3 -> A_2 A_2)"
      IF(BRHAZ(3,1).GT.0d0)
     .  WRITE(18,905) BRHAZ(3,1),2,23,36,"BR(H_3 -> A_1 Z)"
      IF(BRHAZ(3,2).GT.0d0)
     .  WRITE(18,905) BRHAZ(3,2),2,23,46,"BR(H_3 -> A_2 Z)"
      IF(BRHCHC(3).GT.0d0)
     .  WRITE(18,905) BRHCHC(3),2,37,-37,"BR(H_3 -> H+ H-)"
      IF(BRHCW(3).GT.0d0)
     .  WRITE(18,905) BRHCW(3),2,24,-37,"BR(H_3 -> W+ H-)"
      IF(BRHCW(3).GT.0d0)
     .  WRITE(18,905) BRHCW(3),2,-24,37,"BR(H_3 -> W- H+)"
      IF(BRNEU(3,1,1).GT.0d0)
     .  WRITE(18,905) BRNEU(3,1,1),2,1000022,1000022,
     .    "BR(H_3 -> neu_1 neu_1)"
      IF(BRNEU(3,1,2).GT.0d0)
     .  WRITE(18,905) BRNEU(3,1,2),2,1000022,1000023,
     .    "BR(H_3 -> neu_1 neu_2)"
      IF(BRNEU(3,1,3).GT.0d0)
     .  WRITE(18,905) BRNEU(3,1,3),2,1000022,1000025,
     .    "BR(H_3 -> neu_1 neu_3)"
      IF(BRNEU(3,1,4).GT.0d0)
     .  WRITE(18,905) BRNEU(3,1,4),2,1000022,1000035,
     .    "BR(H_3 -> neu_1 neu_4)"
      IF(BRNEU(3,1,5).GT.0d0)
     .  WRITE(18,905) BRNEU(3,1,5),2,1000022,1000045,
     .    "BR(H_3 -> neu_1 neu_5)"
      IF(BRNEU(3,2,2).GT.0d0)
     .  WRITE(18,905) BRNEU(3,2,2),2,1000023,1000023,
     .    "BR(H_3 -> neu_2 neu_2)"
      IF(BRNEU(3,2,3).GT.0d0)
     .  WRITE(18,905) BRNEU(3,2,3),2,1000023,1000025,
     .    "BR(H_3 -> neu_2 neu_3)"
      IF(BRNEU(3,2,4).GT.0d0)
     .  WRITE(18,905) BRNEU(3,2,4),2,1000023,1000035,
     .    "BR(H_3 -> neu_2 neu_4)"
      IF(BRNEU(3,2,5).GT.0d0)
     .  WRITE(18,905) BRNEU(3,2,5),2,1000023,1000045,
     .    "BR(H_3 -> neu_2 neu_5)"
      IF(BRNEU(3,3,3).GT.0d0)
     .  WRITE(18,905) BRNEU(3,3,3),2,1000025,1000025,
     .    "BR(H_3 -> neu_3 neu_3)"
      IF(BRNEU(3,3,4).GT.0d0)
     .  WRITE(18,905) BRNEU(3,3,4),2,1000025,1000035,
     .    "BR(H_3 -> neu_3 neu_4)"
      IF(BRNEU(3,3,5).GT.0d0)
     .  WRITE(18,905) BRNEU(3,3,5),2,1000025,1000045,
     .    "BR(H_3 -> neu_3 neu_5)"
      IF(BRNEU(3,4,4).GT.0d0)
     .  WRITE(18,905) BRNEU(3,4,4),2,1000035,1000035,
     .    "BR(H_3 -> neu_4 neu_4)"
      IF(BRNEU(3,4,5).GT.0d0)
     .  WRITE(18,905) BRNEU(3,4,5),2,1000035,1000045,
     .    "BR(H_3 -> neu_4 neu_5)"
      IF(BRNEU(3,5,5).GT.0d0)
     .  WRITE(18,905) BRNEU(3,5,5),2,1000045,1000045,
     .    "BR(H_3 -> neu_5 neu_5)"
      IF(BRCHA(3,1).GT.0d0)
     .  WRITE(18,905) BRCHA(3,1),2,1000024,-1000024,
     .    "BR(H_3 -> cha_1 cha_1bar)"
      IF(BRCHA(3,2).GT.0d0)
     .  WRITE(18,905) BRCHA(3,2),2,1000024,-1000037,
     .    "BR(H_3 -> cha_1 cha_2bar)"
      IF(BRCHA(3,2).GT.0d0)
     .  WRITE(18,905) BRCHA(3,2),2,1000037,-1000024,
     .    "BR(H_3 -> cha_2 cha_1bar)"
      IF(BRCHA(3,3).GT.0d0)
     .  WRITE(18,905) BRCHA(3,3),2,1000037,-1000037,
     .    "BR(H_3 -> cha_2 cha_2bar)"
      IF(BRHSQ(3,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,1),2,1000002,-1000002,
     .    "BR(H_3 -> ~u_L ~ubar_L)"
      IF(BRHSQ(3,1).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,1),2,1000004,-1000004,
     .    "BR(H_3 -> ~c_L ~cbar_L)"
      IF(BRHSQ(3,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,2),2,2000002,-2000002,
     .    "BR(H_3 -> ~u_R ~ubar_R)"
      IF(BRHSQ(3,2).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,2),2,2000004,-2000004,
     .    "BR(H_3 -> ~c_R ~cbar_R)"
      IF(BRHSQ(3,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,3),2,1000001,-1000001,
     .    "BR(H_3 -> ~d_L ~dbar_L)"
      IF(BRHSQ(3,3).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,3),2,1000003,-1000003,
     .    "BR(H_3 -> ~s_L ~sbar_L)"
      IF(BRHSQ(3,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,4),2,2000001,-2000001,
     .    "BR(H_3 -> ~d_R ~dbar_R)"
      IF(BRHSQ(3,4).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,4),2,2000003,-2000003,
     .    "BR(H_3 -> ~s_R ~sbar_R)"
      IF(BRHSQ(3,5).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,5),2,1000006,-1000006,
     .    "BR(H_3 -> ~t_1 ~tbar_1)"
      IF(BRHSQ(3,6).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,6),2,2000006,-2000006,
     .    "BR(H_3 -> ~t_2 ~tbar_2)"
      IF(BRHSQ(3,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,7),2,1000006,-2000006,
     .    "BR(H_3 -> ~t_1 ~tbar_2)"
      IF(BRHSQ(3,7).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,7),2,2000006,-1000006,
     .    "BR(H_3 -> ~t_2 ~tbar_1)"
      IF(BRHSQ(3,8).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,8),2,1000005,-1000005,
     .    "BR(H_3 -> ~b_1 ~bbar_1)"
      IF(BRHSQ(3,9).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,9),2,2000005,-2000005,
     .    "BR(H_3 -> ~b_2 ~bbar_2)"
      IF(BRHSQ(3,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,10),2,1000005,-2000005,
     .    "BR(H_3 -> ~b_1 ~bbar_2)"
      IF(BRHSQ(3,10).GT.0d0)
     .  WRITE(18,905) BRHSQ(3,10),2,2000005,-1000005,
     .    "BR(H_3 -> ~b_2 ~bbar_1)"
      IF(BRHSL(3,1).GT.0d0)
     .  WRITE(18,905) BRHSL(3,1),2,1000011,-1000011,
     .    "BR(H_3 -> ~e_L ~ebar_L)"
      IF(BRHSL(3,1).GT.0d0)
     .  WRITE(18,905) BRHSL(3,1),2,1000013,-1000013,
     .    "BR(H_3 -> ~mu_L ~mubar_L)"
      IF(BRHSL(3,2).GT.0d0)
     .  WRITE(18,905) BRHSL(3,2),2,2000011,-2000011,
     .    "BR(H_3 -> ~e_R ~ebar_R)"
      IF(BRHSL(3,2).GT.0d0)
     .  WRITE(18,905) BRHSL(3,2),2,2000013,-2000013,
     .    "BR(H_3 -> ~mu_R ~mubarRL)"
      IF(BRHSL(3,3).GT.0d0)
     .  WRITE(18,905) BRHSL(3,3),2,1000012,-1000012,
     .    "BR(H_3 -> ~nu_e_L ~nu_ebar_L)"
      IF(BRHSL(3,3).GT.0d0)
     .  WRITE(18,905) BRHSL(3,3),2,1000014,-1000014,
     .    "BR(H_3 -> ~nu_mu_L ~nu_mubar_L)"
      IF(BRHSL(3,4).GT.0d0)
     .  WRITE(18,905) BRHSL(3,4),2,1000015,-1000015,
     .    "BR(H_3 -> ~tau_1 ~taubar_1)"
      IF(BRHSL(3,5).GT.0d0)
     .  WRITE(18,905) BRHSL(3,5),2,2000015,-2000015,
     .    "BR(H_3 -> ~tau_2 ~taubar_2)"
      IF(BRHSL(3,6).GT.0d0)
     .  WRITE(18,905) BRHSL(3,6),2,1000015,-2000015,
     .    "BR(H_3 -> ~tau_1 ~taubar_2)"
      IF(BRHSL(3,6).GT.0d0)
     .  WRITE(18,905) BRHSL(3,6),2,2000015,-1000015,
     .    "BR(H_3 -> ~tau_2 ~taubar_1)"
      IF(BRHSL(3,7).GT.0d0)
     .  WRITE(18,905) BRHSL(3,7),2,1000016,-1000016,
     .    "BR(H_3 -> ~nu_tau_L ~nu_taubar_L)"

      WRITE(18,904) 36,WIDTH(4),"Lightest pseudoscalar"
      IF(BRJJ(4).GT.0d0)
     .  WRITE(18,905) BRJJ(4),2,21,21,"BR(A_1 -> hadrons)"
      IF(BREE(4).GT.0d0)
     .  WRITE(18,905) BREE(4),2,11,-11,"BR(A_1 -> e- e+)"
      IF(BRMM(4).GT.0d0)
     .  WRITE(18,905) BRMM(4),2,13,-13,"BR(A_1 -> muon muon)"
      IF(BRLL(4).GT.0d0)
     .  WRITE(18,905) BRLL(4),2,15,-15,"BR(A_1 -> tau tau)"
      IF(BRCC(4).GT.0d0)
     .  WRITE(18,905) BRCC(4),2,4,-4,"BR(A_1 -> c cbar)"
      IF(BRBB(4).GT.0d0)
     .  WRITE(18,905) BRBB(4),2,5,-5,"BR(A_1 -> b bbar)"
      IF(BRTT(4).GT.0d0)
     .  WRITE(18,905) BRTT(4),2,6,-6,"BR(A_1 -> t tbar)"
      IF(BRGG(4).GT.0d0)
     .  WRITE(18,905) BRGG(4),2,22,22,"BR(A_1 -> gamma gamma)"
      IF(BRZG(4).GT.0d0)
     .  WRITE(18,905) BRZG(4),2,23,22,"BR(A_1 -> Z gamma)"
      IF(BRAHZ(1,1).GT.0d0)
     .  WRITE(18,905) BRAHZ(1,1),2,23,PDGH1,"BR(A_1 -> Z H_1)"
      IF(BRAHZ(1,2).GT.0d0)
     .  WRITE(18,905) BRAHZ(1,2),2,23,PDGH2,"BR(A_1 -> Z H_2)"
      IF(BRAHZ(1,3).GT.0d0)
     .  WRITE(18,905) BRAHZ(1,3),2,23,45,"BR(A_1 -> Z H_3)"
      IF(BRHCW(4).GT.0d0)
     .  WRITE(18,905) BRHCW(4),2,24,-37,"BR(A_1 -> W+ H-)"
      IF(BRHCW(4).GT.0d0)
     .  WRITE(18,905) BRHCW(4),2,-24,37,"BR(A_1 -> W- H+)"
      IF(BRNEU(4,1,1).GT.0d0)
     .  WRITE(18,905) BRNEU(4,1,1),2,1000022,1000022,
     .    "BR(A_1 -> neu_1 neu_1)"
      IF(BRNEU(4,1,2).GT.0d0)
     .  WRITE(18,905) BRNEU(4,1,2),2,1000022,1000023,
     .    "BR(A_1 -> neu_1 neu_2)"
      IF(BRNEU(4,1,3).GT.0d0)
     .  WRITE(18,905) BRNEU(4,1,3),2,1000022,1000025,
     .    "BR(A_1 -> neu_1 neu_3)"
      IF(BRNEU(4,1,4).GT.0d0)
     .  WRITE(18,905) BRNEU(4,1,4),2,1000022,1000035,
     .    "BR(A_1 -> neu_1 neu_4)"
      IF(BRNEU(4,1,5).GT.0d0)
     .  WRITE(18,905) BRNEU(4,1,5),2,1000022,1000045,
     .    "BR(A_1 -> neu_1 neu_5)"
      IF(BRNEU(4,2,2).GT.0d0)
     .  WRITE(18,905) BRNEU(4,2,2),2,1000023,1000023,
     .    "BR(A_1 -> neu_2 neu_2)"
      IF(BRNEU(4,2,3).GT.0d0)
     .  WRITE(18,905) BRNEU(4,2,3),2,1000023,1000025,
     .    "BR(A_1 -> neu_2 neu_3)"
      IF(BRNEU(4,2,4).GT.0d0)
     .  WRITE(18,905) BRNEU(4,2,4),2,1000023,1000035,
     .    "BR(A_1 -> neu_2 neu_4)"
      IF(BRNEU(4,2,5).GT.0d0)
     .  WRITE(18,905) BRNEU(4,2,5),2,1000023,1000045,
     .    "BR(A_1 -> neu_2 neu_5)"
      IF(BRNEU(4,3,3).GT.0d0)
     .  WRITE(18,905) BRNEU(4,3,3),2,1000025,1000025,
     .    "BR(A_1 -> neu_3 neu_3)"
      IF(BRNEU(4,3,4).GT.0d0)
     .  WRITE(18,905) BRNEU(4,3,4),2,1000025,1000035,
     .    "BR(A_1 -> neu_3 neu_4)"
      IF(BRNEU(4,3,5).GT.0d0)
     .  WRITE(18,905) BRNEU(4,3,5),2,1000025,1000045,
     .    "BR(A_1 -> neu_3 neu_5)"
      IF(BRNEU(4,4,4).GT.0d0)
     .  WRITE(18,905) BRNEU(4,4,4),2,1000035,1000035,
     .    "BR(A_1 -> neu_4 neu_4)"
      IF(BRNEU(4,4,5).GT.0d0)
     .  WRITE(18,905) BRNEU(4,4,5),2,1000035,1000045,
     .    "BR(A_1 -> neu_4 neu_5)"
      IF(BRNEU(4,5,5).GT.0d0)
     .  WRITE(18,905) BRNEU(4,5,5),2,1000045,1000045,
     .    "BR(A_1 -> neu_5 neu_5)"
      IF(BRCHA(4,1).GT.0d0)
     .  WRITE(18,905) BRCHA(4,1),2,1000024,-1000024,
     .    "BR(A_1 -> cha_1 cha_1bar)"
      IF(BRCHA(4,2).GT.0d0)
     .  WRITE(18,905) BRCHA(4,2),2,1000024,-1000037,
     .    "BR(A_1 -> cha_1 cha_2bar)"
      IF(BRCHA(4,2).GT.0d0)
     .  WRITE(18,905) BRCHA(4,2),2,1000037,-1000024,
     .    "BR(A_1 -> cha_2 cha_1bar)"
      IF(BRCHA(4,3).GT.0d0)
     .  WRITE(18,905) BRCHA(4,3),2,1000037,-1000037,
     .    "BR(A_1 -> cha_2 cha_2bar)"
      IF(BRASQ(1,1).GT.0d0)
     .  WRITE(18,905) BRASQ(1,1),2,1000006,-2000006,
     .    "BR(A_1 -> ~t_1 ~tbar_2)"
      IF(BRASQ(1,1).GT.0d0)
     .  WRITE(18,905) BRASQ(1,1),2,2000006,-1000006,
     .    "BR(A_1 -> ~t_2 ~tbar_1)"
      IF(BRASQ(1,2).GT.0d0)
     .  WRITE(18,905) BRASQ(1,2),2,1000005,-2000005,
     .    "BR(A_1 -> ~b_1 ~bbar_2)"
      IF(BRASQ(1,2).GT.0d0)
     .  WRITE(18,905) BRASQ(1,2),2,2000005,-1000005,
     .    "BR(A_1 -> ~b_2 ~bbar_1)"
      IF(BRASL(1).GT.0d0)
     .  WRITE(18,905) BRASL(1),2,1000015,-2000015,
     .    "BR(A_1 -> ~tau_1 ~taubar_2)"
      IF(BRASL(1).GT.0d0)
     .  WRITE(18,905) BRASL(1),2,2000015,-1000015,
     .    "BR(A_1 -> ~tau_2 ~taubar_1)"

      WRITE(18,904) 46,WIDTH(5),"2nd pseudoscalar"
      IF(BRJJ(5).GT.0d0)
     .  WRITE(18,905) BRJJ(5),2,21,21,"BR(A_2 -> hadrons)"
      IF(BREE(5).GT.0d0)
     .  WRITE(18,905) BREE(5),2,11,-11,"BR(A_2 -> e- e+)"
      IF(BRMM(5).GT.0d0)
     .  WRITE(18,905) BRMM(5),2,13,-13,"BR(A_2 -> muon muon)"
      IF(BRLL(5).GT.0d0)
     .  WRITE(18,905) BRLL(5),2,15,-15,"BR(A_2 -> tau tau)"
      IF(BRCC(5).GT.0d0)
     .  WRITE(18,905) BRCC(5),2,4,-4,"BR(A_2 -> c cbar)"
      IF(BRBB(5).GT.0d0)
     .  WRITE(18,905) BRBB(5),2,5,-5,"BR(A_2 -> b bbar)"
      IF(BRTT(5).GT.0d0)
     .  WRITE(18,905) BRTT(5),2,6,-6,"BR(A_2 -> t tbar)"
      IF(BRGG(5).GT.0d0)
     .  WRITE(18,905) BRGG(5),2,22,22,"BR(A_2 -> gamma gamma)"
      IF(BRZG(5).GT.0d0)
     .  WRITE(18,905) BRZG(5),2,23,22,"BR(A_2 -> Z gamma)"
      IF(BRAHA(1).GT.0d0)
     .  WRITE(18,905) BRAHA(1),2,36,PDGH1,"BR(A_2 -> A_1 H_1)"
      IF(BRAHA(2).GT.0d0)
     .  WRITE(18,905) BRAHA(2),2,36,PDGH2,"BR(A_2 -> A_1 H_2)"
      IF(BRAHA(3).GT.0d0)
     .  WRITE(18,905) BRAHA(3),2,36,45,"BR(A_2 -> A_1 H_3)"
      IF(BRAHZ(2,1).GT.0d0)
     .  WRITE(18,905) BRAHZ(2,1),2,23,PDGH1,"BR(A_2 -> Z H_1)"
      IF(BRAHZ(2,2).GT.0d0)
     .  WRITE(18,905) BRAHZ(2,2),2,23,PDGH2,"BR(A_2 -> Z H_2)"
      IF(BRAHZ(2,3).GT.0d0)
     .  WRITE(18,905) BRAHZ(2,3),2,23,45,"BR(A_2 -> Z H_3)"
      IF(BRHCW(5).GT.0d0)
     .  WRITE(18,905) BRHCW(5),2,24,-37,"BR(A_2 -> W+ H-)"
      IF(BRHCW(5).GT.0d0)
     .  WRITE(18,905) BRHCW(5),2,-24,37,"BR(A_2 -> W- H+)"
      IF(BRNEU(5,1,1).GT.0d0)
     .  WRITE(18,905) BRNEU(5,1,1),2,1000022,1000022,
     .    "BR(A_2 -> neu_1 neu_1)"
      IF(BRNEU(5,1,2).GT.0d0)
     .  WRITE(18,905) BRNEU(5,1,2),2,1000022,1000023,
     .    "BR(A_2 -> neu_1 neu_2)"
      IF(BRNEU(5,1,3).GT.0d0)
     .  WRITE(18,905) BRNEU(5,1,3),2,1000022,1000025,
     .    "BR(A_2 -> neu_1 neu_3)"
      IF(BRNEU(5,1,4).GT.0d0)
     .  WRITE(18,905) BRNEU(5,1,4),2,1000022,1000035,
     .    "BR(A_2 -> neu_1 neu_4)"
      IF(BRNEU(5,1,5).GT.0d0)
     .  WRITE(18,905) BRNEU(5,1,5),2,1000022,1000045,
     .    "BR(A_2 -> neu_1 neu_5)"
      IF(BRNEU(5,2,2).GT.0d0)
     .  WRITE(18,905) BRNEU(5,2,2),2,1000023,1000023,
     .    "BR(A_2 -> neu_2 neu_2)"
      IF(BRNEU(5,2,3).GT.0d0)
     .  WRITE(18,905) BRNEU(5,2,3),2,1000023,1000025,
     .    "BR(A_2 -> neu_2 neu_3)"
      IF(BRNEU(5,2,4).GT.0d0)
     .  WRITE(18,905) BRNEU(5,2,4),2,1000023,1000035,
     .    "BR(A_2 -> neu_2 neu_4)"
      IF(BRNEU(5,2,5).GT.0d0)
     .  WRITE(18,905) BRNEU(5,2,5),2,1000023,1000045,
     .    "BR(A_2 -> neu_2 neu_5)"
      IF(BRNEU(5,3,3).GT.0d0)
     .  WRITE(18,905) BRNEU(5,3,3),2,1000025,1000025,
     .    "BR(A_2 -> neu_3 neu_3)"
      IF(BRNEU(5,3,4).GT.0d0)
     .  WRITE(18,905) BRNEU(5,3,4),2,1000025,1000035,
     .    "BR(A_2 -> neu_3 neu_4)"
      IF(BRNEU(5,3,5).GT.0d0)
     .  WRITE(18,905) BRNEU(5,3,5),2,1000025,1000045,
     .    "BR(A_2 -> neu_3 neu_5)"
      IF(BRNEU(5,4,4).GT.0d0)
     .  WRITE(18,905) BRNEU(5,4,4),2,1000035,1000035,
     .    "BR(A_2 -> neu_4 neu_4)"
      IF(BRNEU(5,4,5).GT.0d0)
     .  WRITE(18,905) BRNEU(5,4,5),2,1000035,1000045,
     .    "BR(A_2 -> neu_4 neu_5)"
      IF(BRNEU(5,5,5).GT.0d0)
     .  WRITE(18,905) BRNEU(5,5,5),2,1000045,1000045,
     .    "BR(A_2 -> neu_5 neu_5)"
      IF(BRCHA(5,1).GT.0d0)
     .  WRITE(18,905) BRCHA(5,1),2,1000024,-1000024,
     .    "BR(A_2 -> cha_1 cha_1bar)"
      IF(BRCHA(5,2).GT.0d0)
     .  WRITE(18,905) BRCHA(5,2),2,1000024,-1000037,
     .    "BR(A_2 -> cha_1 cha_2bar)"
      IF(BRCHA(5,2).GT.0d0)
     .  WRITE(18,905) BRCHA(5,2),2,1000037,-1000024,
     .    "BR(A_2 -> cha_2 cha_1bar)"
      IF(BRCHA(5,3).GT.0d0)
     .  WRITE(18,905) BRCHA(5,3),2,1000037,-1000037,
     .    "BR(A_2 -> cha_2 cha_2bar)"
      IF(BRASQ(2,1).GT.0d0)
     .  WRITE(18,905) BRASQ(2,1),2,1000006,-2000006,
     .    "BR(A_2 -> ~t_1 ~tbar_2)"
      IF(BRASQ(2,1).GT.0d0)
     .  WRITE(18,905) BRASQ(2,1),2,2000006,-1000006,
     .    "BR(A_2 -> ~t_2 ~tbar_1)"
      IF(BRASQ(2,2).GT.0d0)
     .  WRITE(18,905) BRASQ(2,2),2,1000005,-2000005,
     .    "BR(A_2 -> ~b_1 ~bbar_2)"
      IF(BRASQ(2,2).GT.0d0)
     .  WRITE(18,905) BRASQ(2,2),2,2000005,-1000005,
     .    "BR(A_2 -> ~b_2 ~bbar_1)"
      IF(BRASL(2).GT.0d0)
     .  WRITE(18,905) BRASL(2),2,1000015,-2000015,
     .    "BR(A_2 -> ~tau_1 ~taubar_2)"
      IF(BRASL(2).GT.0d0)
     .  WRITE(18,905) BRASL(2),2,2000015,-1000015,
     .    "BR(A_2 -> ~tau_2 ~taubar_1)"

      WRITE(18,904) 37,HCWIDTH,"Charged Higgs"
      IF(HCBRM.GT.0d0)
     .  WRITE(18,905) HCBRM,2,-13,14,"BR(H+ -> muon nu_muon)"
      IF(HCBRL.GT.0d0)
     .  WRITE(18,905) HCBRL,2,-15,16,"BR(H+ -> tau nu_tau)"
      IF(HCBRSU.GT.0d0)
     .  WRITE(18,905) HCBRSU,2,2,-3,"BR(H+ -> u sbar)"
      IF(HCBRSC.GT.0d0)
     .  WRITE(18,905) HCBRSC,2,4,-3,"BR(H+ -> c sbar)"
      IF(HCBRBU.GT.0d0)
     .  WRITE(18,905) HCBRBU,2,2,-5,"BR(H+ -> u bbar)"
      IF(HCBRBC.GT.0d0)
     .  WRITE(18,905) HCBRBC,2,4,-5,"BR(H+ -> c bbar)"
      IF(HCBRBT.GT.0d0)
     .  WRITE(18,905) HCBRBT,2,6,-5,"BR(H+ -> t bbar)"
      IF(HCBRWH(1).GT.0d0)
     .  WRITE(18,905) HCBRWH(1),2,24,PDGH1,"BR(H+ -> W+ H_1)"
      IF(HCBRWH(2).GT.0d0)
     .  WRITE(18,905) HCBRWH(2),2,24,PDGH2,"BR(H+ -> W+ H_2)"
      IF(HCBRWH(4).GT.0d0)
     .  WRITE(18,905) HCBRWH(4),2,24,36,"BR(H+ -> W+ A_1)"
      IF(HCBRNC(1,1).GT.0d0)
     .  WRITE(18,905) HCBRNC(1,1),2,1000024,1000022,
     .    "BR(H+ -> cha_1 neu_1)"
      IF(HCBRNC(2,1).GT.0d0)
     .  WRITE(18,905) HCBRNC(2,1),2,1000024,1000023,
     .    "BR(H+ -> cha_1 neu_2)"
      IF(HCBRNC(3,1).GT.0d0)
     .  WRITE(18,905) HCBRNC(3,1),2,1000024,1000025,
     .    "BR(H+ -> cha_1 neu_3)"
      IF(HCBRNC(4,1).GT.0d0)
     .  WRITE(18,905) HCBRNC(4,1),2,1000024,1000035,
     .    "BR(H+ -> cha_1 neu_4)"
      IF(HCBRNC(5,1).GT.0d0)
     .  WRITE(18,905) HCBRNC(5,1),2,1000024,1000045,
     .    "BR(H+ -> cha_1 neu_5)"
      IF(HCBRNC(1,2).GT.0d0)
     .  WRITE(18,905) HCBRNC(1,2),2,1000037,1000022,
     .    "BR(H+ -> cha_2 neu_1)"
      IF(HCBRNC(2,2).GT.0d0)
     .  WRITE(18,905) HCBRNC(2,2),2,1000037,1000023,
     .    "BR(H+ -> cha_2 neu_2)"
      IF(HCBRNC(3,2).GT.0d0)
     .  WRITE(18,905) HCBRNC(3,2),2,1000037,1000025,
     .    "BR(H+ -> cha_2 neu_3)"
      IF(HCBRNC(4,2).GT.0d0)
     .  WRITE(18,905) HCBRNC(4,2),2,1000037,1000035,
     .    "BR(H+ -> cha_2 neu_4)"
      IF(HCBRNC(5,2).GT.0d0)
     .  WRITE(18,905) HCBRNC(5,2),2,1000037,1000045,
     .    "BR(H+ -> cha_2 neu_5)"
      IF(HCBRSQ(1).GT.0d0)
     .  WRITE(18,905) HCBRSQ(1),2,1000002,-1000001,
     .    "BR(H+ -> ~u_L ~dbar_L)"
      IF(HCBRSQ(1).GT.0d0)
     .  WRITE(18,905) HCBRSQ(1),2,1000004,-1000003,
     .    "BR(H+ -> ~c_L ~sbar_L)"
      IF(HCBRSQ(2).GT.0d0)
     .  WRITE(18,905) HCBRSQ(2),2,1000006,-1000005,
     .    "BR(H+ -> ~t_1 ~bbar_1)"
      IF(HCBRSQ(3).GT.0d0)
     .  WRITE(18,905) HCBRSQ(3),2,1000006,-2000005,
     .    "BR(H+ -> ~t_1 ~bbar_2)"
      IF(HCBRSQ(4).GT.0d0)
     .  WRITE(18,905) HCBRSQ(4),2,2000006,-1000005,
     .    "BR(H+ -> ~t_2 ~bbar_1)"
      IF(HCBRSQ(5).GT.0d0)
     .  WRITE(18,905) HCBRSQ(5),2,2000006,-2000005,
     .    "BR(H+ -> ~t_2 ~bbar_2)"
      IF(HCBRSL(1).GT.0d0)
     .  WRITE(18,905) HCBRSL(1),2,1000012,-1000011,
     .    "BR(H+ -> ~nu_e_L ~ebar_L)"
      IF(HCBRSL(1).GT.0d0)
     .  WRITE(18,905) HCBRSL(1),2,1000014,-1000013,
     .    "BR(H+ -> ~nu_mu_L ~mubar_L)"
      IF(HCBRSL(2).GT.0d0)
     .  WRITE(18,905) HCBRSL(2),2,1000016,-1000015,
     .    "BR(H+ -> ~nu_tau_L ~taubar_1)"
      IF(HCBRSL(3).GT.0d0)
     .  WRITE(18,905) HCBRSL(3),2,1000016,-2000015,
     .    "BR(H+ -> ~nu_tau_L ~taubar_2)"

      WRITE(18,904) 6,toptot,'Top Quark'
      IF(brtopbw.ne.0.D0)
     .  WRITE(18,905) brtopbw,2,5,24,'BR(t ->  b    W+)'
      IF(brtopbh.ne.0.D0)
     .  WRITE(18,905) brtopbh,2,5,37,'BR(t ->  b    H+)'
      IF(brtopneutrstop(1,1).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(1,1),2,1000006,1000022,
     . 'BR(t -> ~t_1 ~chi_10)'
      IF(brtopneutrstop(2,1).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(2,1),2,1000006,1000023,
     . 'BR(t -> ~t_1 ~chi_20)'
      IF(brtopneutrstop(3,1).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(3,1),2,1000006,1000025,
     . 'BR(t -> ~t_1 ~chi_30)'
      IF(brtopneutrstop(4,1).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(4,1),2,1000006,1000025,
     . 'BR(t -> ~t_1 ~chi_40)'
      IF(brtopneutrstop(5,1).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(5,1),2,1000006,1000025,
     . 'BR(t -> ~t_1 ~chi_50)'
      IF(brtopneutrstop(1,2).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(1,2),2,2000006,1000022,
     . 'BR(t -> ~t_2 ~chi_10)'
      IF(brtopneutrstop(2,2).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(2,2),2,2000006,1000023,
     . 'BR(t -> ~t_2 ~chi_20)'
      IF(brtopneutrstop(3,2).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(3,2),2,2000006,1000025,
     .'BR(t -> ~t_2 ~chi_30)'
      IF(brtopneutrstop(4,2).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(4,2),2,2000006,1000025,
     . 'BR(t -> ~t_2 ~chi_40)'
      IF(brtopneutrstop(5,2).ne.0.D0)
     .  WRITE(18,905) brtopneutrstop(5,2),2,2000006,1000025,
     . 'BR(t -> ~t_2 ~chi_50)'

      CALL NS_OUTPUT(18)

      CLOSE(18)

 899  FORMAT(A)
 901  FORMAT(1X,I5,3X,1P,E16.8,0P,3X,'#',1X,A)
 902  FORMAT(1X,I9,3X,1P,E16.8,0P,3X,'#',1X,A)
 903  FORMAT(1X,I2,1X,I2,3X,1P,E16.8,0P,3X,'#',1X,A)
 904  FORMAT('DECAY',1X,I9,3X,1P,E16.8,0P,3X,'#',1X,A)
 905  FORMAT(3X,1P,E16.8,0P,3X,I2,3X,I9,1X,I9,1X,2X,'#',1X,A)
 906  FORMAT('#',1X,A,3X,E16.8)
 907  FORMAT(A,1P,E16.8,A)
 920  FORMAT('#',0P,I5,3X,1P,E16.8,0P,3X,'#',1X,A)
 921  FORMAT(1X,I2,1X,I2,3X,'#',1X,A)

      END
