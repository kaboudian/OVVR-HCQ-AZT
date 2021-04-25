c Fortran Implementation of the O'Hara-Rudy dynamic (ORd) model for the
c undiseased human ventricular action potential and calcium transient
c
c The ORd model is described in the article "Simulation of the Undiseased
c Human Cardiac Ventricular Action Potential: Model Formulation and
c Experimental Validation"
c by Thomas O'Hara, Laszlo Virag, Andras Varro, and Yoram Rudy
      program ohara0dMid
      implicit none
      integer nx,ii
      parameter(nx=1)
c variables
      real*8 v(nx),Nai(nx),Nass(nx),Ki(nx),Kss(nx),
     &     Cai(nx),Cass(nx),Cansr(nx),Cajsr(nx),xm(nx),
     &     xhf(nx),xhs(nx),xj(nx),hCaMK_slow(nx),jCaMK(nx),
     &     xmL(nx),xhL(nx),xhLCaMK(nx),xa(nx),i_fast(nx),
     &     i_slow(nx),a_CaMK(nx),iCaMK_fast(nx),iCaMK_slow(nx),xd(nx),
     &     xff(nx),xfs(nx),xfcaf(nx),xfcas(nx),xjca(nx),
     &     xnca(nx),xffp(nx),xfcafp(nx),xrf(nx),xrs(nx),
     &     xs1(nx),xs2(nx),xk1(nx),Jrel_NP(nx),Jrel_CaMK(nx),
     &     CaMK_trap(nx),JNakNa(nx),JnakK(nx),vt(nx)
c fixed parameters and constants
      real*8 nao,cao,ko,R,T,F,FoRT,RToF,L,rad,vcell,Ageo,Acap,vmyo,vnsr,
     &       vjsr,vss,KmCaMK,aCaMK,bCaMK,CaMKo,KmCaM,PKNa,
     &       Ahf,Ahs,GNa,thL,tjca,Kmn,k2n,zca,kna1,kna2,kna3,kasymm,
     &       wna,wca,wnaca,kcaon,kcaoff,qna,qca,KmCaAct,zna,
     &       k1p,k1m,k2p,k2m,k3m,k4m,KNai0,Knao0,delta,
     &       KKi,Kko,MgADP,MgATP,Kmgatp,H,eP,Khp,Knap,Kxkur,zk,
     &       PNab,PCab,GpCa,amp,duration,bt,kmcmdn,
     &       trpnmax,kmtrpn,BSRmax,KmBSR,BSLmax,KmBSL,csqnmax,kmcsqn
      parameter(nao=140.0,cao=1.8,ko=5.4,
     &     R=8314.0,T=310.0,F=96485.0,
     &     L=0.01,rad=0.0011,
     &     vcell=1000*3.14*rad*rad*L,
     &     Ageo=2*3.14*rad*rad+2*3.14*rad*L,
     &     Acap=2*Ageo,vmyo=0.68*vcell,vnsr=0.0552*vcell,
     &     vjsr=0.0048*vcell,vss=0.02*vcell,
     &     KmCaMK=0.15,aCaMK=0.05,bCaMK=0.00068,CaMKo=0.05,
     &     KmCaM=0.0015,
     &     PKNa=0.01833,
     &     Ahf=0.99,Ahs=1.0-Ahf,
     &     GNa=75,thL=200.0,tjca=75.0,
     &     Kmn=0.002,k2n=1000.0,zca=2.0,
     &     kna1=15.0,kna2=5.0,kna3=88.12,kasymm=12.5,
     &     wna=6.0e4,wca=6.0e4,wnaca=5.0e3,kcaon=1.5e6,
     &     kcaoff=5.0e3,qna=0.5224,qca=0.1670,KmCaAct=150.0e-6,
     &     zna=1.0,
     &     k1p=949.5,k1m=182.4,k2p=687.2,k2m=39.4,
     &     k3m=79300.0,k4m=40.0,KNai0=9.073,Knao0=27.78,
     &     delta=-0.1550,
     &     KKi=0.5,Kko=0.3582,MgADP=0.05,MgATP=9.8,
     &     Kmgatp=1.698e-7,H=1.0e-7,eP=4.2,Khp=1.698e-7,
     &     Knap=224.0,Kxkur=292.0,zk=1.0,
     &     PNab=3.75e-10,PCab=2.5e-8,GpCa=0.0005,
     &     amp=-35.0,duration=2.0,
     &     bt=4.7,kmcmdn=0.00238,trpnmax=0.07,
     &     kmtrpn=0.0005,BSRmax=0.047,KmBSR=0.00087,BSLmax=1.124,
     &     KmBSL=0.0087,csqnmax=10.0,kmcsqn=0.8)
c currents
      real*8 INa,INaL,Ito,ICaL,ICaNa,ICaK,IKr,IKs,IK1,
     &     INaCa_i,INaCa_ss,INaK,IKb,INab,ICab,IpCa,Istim,
     &     JdiffNa,JdiffK,Jdiff,Jrel,Jleak,Jup,Jtr
c other values
      real*8 CaMKb,CaMKa,dCaMK_trap,ENa,EK,EKs,
     &     xh,hp,fINap,tauhLCaMKt,GNaL,fINaLp,tauat,delta_epi,xi,ip,Gto,
     &     fItop,Aff,Afs,xf,fca,fp,fcap,km2n,anca,dnca,PhiCaL,
     &     PhiCaNa,PhiCaK,PCa,PCap,PCaNa,PCaK,PCaNap,PCaKp,fICaLp,
     &     xr,GKr,KsCa,GKs,GK1,
     &     h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,k1,k2,k3p,k3pp,
     &     k3,k4p,k4pp,k4,k5,k6,k7,k8,x1,x2,x3,x4,E1,E2,E3,E4,
     &     allo,JncxNa,JncxCa,Gncx,P,a1,b1,a2,b2,a3,b3,a4,b4,!ak1,
     &     Pnak,GKb,a_rel,Jrel_inf,tau_rel,dJrel_NP,btp,
     &     a_relp,Jrel_infp,tau_relp,dJrel_CaMK,fJrel_CaMK,Jupnp,Jupp,
     &     fJupp,cmdnmax,dNai,dNass,dKi,dKss,BCai,dCai,BCass,
     &     dCass,dCansr,BCajsr,dCajsr

      real*8 tauit_fast,tauit_slow,tauaCaMKt,tauiCaMK_fast,
     &       tauiCaMK_slow,taudt,tauft_fast,tauft_slow,taufCat_fast,
     &       taufCat_slow,taujCat,taufCaMKt_fast,taufCa_CaMKt_fast,
     &       tauXrt_fast,tauXrt_slow,tauXs1t,tauXs2t,tauXk1t


      integer ntime,nend,i,icelltype
      real*8 cyclelength,dt,endtime

      character*40 voltfile,output
      character*40 APDDI
      integer idig1,idig2,idig3,idig4,istp,intcl

c threshold used to measure apd from di
      real*8 vr,vru

c APD90 threshold for 0d, 30s, CL=1000ms
c Maximum of the Action Potential value = 39.99
c Resting Potential value = -88.00
c %90 of the action potential will be (-88.00*0.9) + (39.99*0.1)
c APD_90 = -79.200 + 3.999 = -75.21
      parameter(vr=-75.21,vru=vr)

      integer ncl,ncycles,icycle
      parameter(ncl=50)
      real*4 cls(ncl)
      integer nups1,ndowns1,i1,maxbeats,icl
      parameter(maxbeats=1000)
      real*8 ups1(maxbeats),downs1(maxbeats)
      real*8 time,apd,di,cl,prevapd,prevdi,prevcl
!      real*8 aaa,bbb
!      parameter(aaa=1.0)

c tables
      integer nvt,iv1
      real*8 zindexv,vlo,vhi,dvt,vv
      parameter(vlo=-100.0,vhi=100.0,nvt=20000)
      real*8 vffrt(0:nvt),vfrt(0:nvt),Knao(0:nvt),KNai(0:nvt),
     &       mss(0:nvt),hss(0:nvt),jss(0:nvt),hssp(0:nvt),tm(0:nvt),
     &       thf(0:nvt),ths(0:nvt),tj(0:nvt),thsp(0:nvt),tjp(0:nvt),
     &       exptaumt(0:nvt),exptauht_fast(0:nvt),
     &       exptauht_slow(0:nvt),exptaujt(0:nvt),
     &       exptauhCaMK_slow(0:nvt),exptaujCaMK(0:nvt),
     &       exptaumLt(0:nvt),mLss(0:nvt),tmL(0:nvt),hLinft(0:nvt),
     &       exptauhLt(0:nvt),hLCaMKinft(0:nvt),exptauhLCaMK(0:nvt),
     &       ainft(0:nvt),exptauat(0:nvt),iinft(0:nvt),
     &       exptauit_fast(0:nvt),exptauit_slow(0:nvt),Ai_fast(0:nvt),
     &       Ai_slow(0:nvt),aCaMKinft(0:nvt),exptauaCaMKt(0:nvt),
     &       iCaMKinft(0:nvt),DeltaCaMK_develop(0:nvt),
     &       DeltaCaMK_recover(0:nvt),exptauiCaMK_fast(0:nvt),
     &       exptauiCaMK_slow(0:nvt),AiCaMK_fast(0:nvt),
     &       AiCaMK_slow(0:nvt),dinft(0:nvt),exptaudt(0:nvt),
     &       finft(0:nvt),exptauft_fast(0:nvt),exptauft_slow(0:nvt),
     &       fCainft(0:nvt),exptaufCat_fast(0:nvt),
     &       exptaufCat_slow(0:nvt),AfCa_fast(0:nvt),AfCa_slow(0:nvt),
     &       jCainft(0:nvt),exptaujCat(0:nvt),fCaMKinft(0:nvt),
     &       exptaufCaMKt_fast(0:nvt),fCa_CaMKinft(0:nvt),
     &       exptaufCa_CaMKt_fast(0:nvt),AfCa_CaMK_fast(0:nvt),
     &       AfCa_CaMK_slow(0:nvt),Xrinft(0:nvt),exptauXrt_fast(0:nvt),
     &       exptauXrt_slow(0:nvt),AXrt_fast(0:nvt),AXrt_slow(0:nvt),
     &       RKr(0:nvt),Xs1inft(0:nvt),exptauXs1t(0:nvt),Xs2inft(0:nvt),
     &       exptauXs2t(0:nvt),Xk1inft(0:nvt),exptauXk1t(0:nvt),
     &       RK1(0:nvt),hCa(0:nvt),hNa(0:nvt),X_Kb(0:nvt)


!      cls(1)=4000.
!      cls(2)=3750.
!      cls(3)=3500.
!      cls(4)=3250.
!      cls(5)=3000.
!      cls(6)=2800.
!      cls(7)=2600.
!      cls(8)=2400.
!      cls(9)=2200.
!      cls(10)=2000.
!      cls(11)=1800.
!      cls(12)=1600.
!      cls(13)=1400.
!      cls(14)=1200.
!      cls(15)=1000.
!      cls(16)=950.
!      cls(17)=900.
!      cls(18)=850.
!      cls(19)=800.
!      cls(20)=750.



      cls(1)=1000.
      cls(2)=950.
      cls(3)=900.
      cls(4)=850.
      cls(5)=800.
      cls(6)=750.
      cls(7)=700.
      cls(8)=650.
      cls(9)=600.
      cls(10)=575.
      cls(11)=550.
      cls(12)=525.
      cls(13)=500.
      cls(14)=480.
      cls(15)=460.
      cls(16)=440.
      cls(17)=420.
      cls(18)=400.
      cls(19)=390.
      cls(20)=380.
      cls(21)=370.
      cls(22)=360.
      cls(23)=350.
      cls(24)=340.
      cls(25)=330.
      cls(26)=320.
      cls(27)=310.
      cls(28)=300.
      cls(29)=290.
      cls(30)=280.
      cls(31)=270.
      cls(32)=260.
      cls(33)=250.
      cls(34)=245.
      cls(35)=240.
      cls(36)=235.
      cls(37)=230.
      cls(38)=225.
      cls(39)=220.
      cls(40)=215.
      cls(41)=210.
      cls(42)=205.
      cls(43)=200.
      cls(44)=195.
      cls(45)=190.
      cls(46)=185.
      cls(47)=180.
      cls(48)=175.
      cls(49)=170.
      cls(50)=165.

c initial conditions for state variables
      do i=1,nx
         v(i)=-87
         Nai(i)=7
         Nass(i)=Nai(i)
         Ki(i)=145
         Kss(i)=Ki(i)
         Cai(i)=1.0e-4
         Cass(i)=Cai(i)
         Cansr(i)=1.2
         Cajsr(i)=Cansr(i)
         xm(i)=0
         xhf(i)=1
         xhs(i)=1
         xj(i)=1
         hCaMK_slow(i)=1
         jCaMK(i)=1
         xmL(i)=0
         xhL(i)=1
         xhLCaMK(i)=1
         xa(i)=0
         i_fast(i)=1
         i_slow(i)=1
         a_CaMK(i)=0
         iCaMK_fast(i)=1
         iCaMK_slow(i)=1
         xd(i)=0
         xff(i)=1
         xfs(i)=1
         xfcaf(i)=1
         xfcas(i)=1
         xjca(i)=1
         xnca(i)=0
         xffp(i)=1
         xfcafp(i)=1
         xrf(i)=0
         xrs(i)=0
         xs1(i)=0
         xs2(i)=0
         xk1(i)=1
         Jrel_NP(i)=0
         Jrel_CaMK(i)=0
         CaMK_trap(i)=0
         JnakNa(i)=0
         JnakK(i)=0
      enddo

c numerical parameters
c      CL=1000 !pacing cycle length in ms
      dt=0.001
      endtime=30000.0
c      nend=nint(endtime/dt)
      icl=1
      cyclelength=cls(icl)
      icycle=nint(cyclelength/dt)
c     do whole cycles, don't stop in the middle of one
      ncycles=nint(endtime/cyclelength)
      if(ncycles.lt.endtime/cyclelength) then
        ncycles=ncycles+1
      endif
      nend=icycle*ncycles
      write(6,*) 'write n cycles of length x', ncycles,cyclelength

      FoRT = F/(R*T)
      RToF = 1/FoRT
      
      icelltype=2 ! endo = 0, epi = 1, M = 2

         GNaL=0.0075
         if (icelltype.eq.1) then
            GNaL=GNaL*0.6
         endif

         Gto=0.02
         if (icelltype.eq.1) then
            Gto=Gto*4.0
         elseif (icelltype.eq.2) then
            Gto=Gto*4.0
         endif

         Aff=0.6
         Afs=1.0-Aff

         PCa=0.0001
         if(icelltype.eq.1) then
            PCa=PCa*1.2
         elseif (icelltype.eq.2) then
            PCa=PCa*2.5
         endif
         PCap=1.1*PCa
         PCaNa=0.00125*PCa
         PCaK=3.574e-4*PCa
         PCaNap=0.00125*PCap
         PCaKp=3.574e-4*PCap


         GKr=0.046
         if(icelltype.eq.1) then
            GKr=GKr*1.3
         elseif (icelltype.eq.2) then
            GKr=GKr*0.8
         endif

         GKs=0.0034
         if(icelltype.eq.1) then
            GKs=GKs*1.4
         endif

         GK1=0.1908
         if(icelltype.eq.1) then
            GK1=GK1*1.2
         elseif(icelltype.eq.2) then
            GK1=GK1*1.3
         endif

         Gncx=0.0008
         if(icelltype.eq.1) then
            Gncx=Gncx*1.1
         elseif(icelltype.eq.2) then
            Gncx=Gncx*1.4
         endif

         k3p=1899.0
         k4p=639.0

         Pnak=30.
         if(icelltype.eq.1) then
            Pnak=Pnak*0.9
         elseif(icelltype.eq.2) then
            Pnak=Pnak*0.7
         endif

         GKb=0.003
         if(icelltype.eq.1) then
            GKb=GKb*0.6
         endif

c calcium buffer constants
         cmdnmax=0.05
         if(icelltype.eq.1) then
            cmdnmax=cmdnmax*1.3
         endif

c tables setup
      dvt=(vhi-vlo)/float(nvt)

      do ii=0,nvt

         vv=vlo+float(ii)*dvt

         if(abs(vv).lt.1e-10) then
            vv=0.01
            vffrt(ii)=vv*F*FoRT
            vfrt(ii)=vv*FoRT
         else
            vffrt(ii)=vv*F*FoRT
            vfrt(ii)=vv*FoRT
         endif

            mss(ii)=1.0/(1.0+exp((-(vv+39.57))/9.871))
c Equation (4)
            tm(ii)=1.0/(6.765*exp((vv+11.64)/34.77)
     &        +8.552*exp(-(vv+77.42)/5.955))
            exptaumt(ii)=exp(-dt/tm(ii))
!c Equation (5)
            hss(ii)=1.0/(1+exp((vv+82.90)/6.086))
!c Equation (7)
            thf(ii)=1.0/(1.432e-5*exp(-(vv+1.196)/6.285)
     &        +6.149*exp((vv+0.5096)/20.27))
            exptauht_fast(ii)=exp(-dt/thf(ii))
!c Equation (8)
            ths(ii)=1.0/(0.009794*exp(-(vv+17.95)/28.05)
     &        +0.3343*exp((vv+5.730)/56.66))
            exptauht_slow(ii)=exp(-dt/ths(ii))
!c Equation (9)
            jss(ii)=hss(ii)
!c Equation (13)
            tj(ii)=2.038+1.0/(0.02136*exp(-(vv+100.6)/8.281)
     &        +0.3052*exp((vv+0.9941)/38.45))
            exptaujt(ii)=exp(-dt/tj(ii))
!c Equation (14)
             hssp(ii)=1.0/(1+exp((vv+89.1)/6.086))
!            hCaMKinft(ii)=1.0/(1+exp((vv+89.1)/6.086)) !hssp-114
!c Equation (16)
             thsp(ii)=3.0*ths(ii)
!            tauhCaMK_slow=3.0*tauht_slow !thsp-115
            exptauhCaMK_slow(ii)=exp(-dt/thsp(ii))
!c Equation (17)
!            jCaMKinft(ii)=1.0/(1+exp((vv+82.90)/6.086))
!c Equation (21)
             tjp(ii)=1.46*tj(ii)
!            taujCaMK=1.46*taujt !tjp-118
            exptaujCaMK(ii)=exp(-dt/tjp(ii))
!c Equation (22)
            mLss(ii)=1.0/(1.0+exp((-(vv+42.85))/5.264))
!c Equation (26)
            tmL(ii)=1.0/(6.765*exp((vv+11.64)/34.77)
     &           +8.552*exp(-(vv+77.42)/5.955))
            exptaumLt(ii)=exp(-dt/tmL(ii))
!c Equation (27)
            hLinft(ii)=1.0/(1.0+exp((vv+87.61)/7.488))
!c Equation (29)
c            thL=200.0
            exptauhLt(ii)=exp(-dt/thL)
!c Equation (c2)
            hLCaMKinft(ii)=1.0/(1.0+exp((vv+93.81)/7.488)) !hLssp-131
!c Equation (31)
            tauhLCaMKt=3.0*thL !thLp-132
            exptauhLCaMK(ii)=exp(-dt/tauhLCaMKt)
!c Equation (32)
            ainft(ii)=1.0/(1.0+exp((-(vv-14.34))/14.82))
!c Equation (39)
            tauat=1.0515/(1.0/(1.2089*(1.0+exp(-(vv-18.4099)/29.3814)))
     &           +3.5/(1.0+exp((vv+100.0)/29.3814)))
            exptauat(ii)=exp(-dt/tauat)
!c Equation (40)
            iinft(ii)=1.0/(1.0+exp((vv+43.94)/5.711))
!c Equation (42)

            if (icelltype.eq.1) then
               delta_epi=1.0-(0.95/(1.0+exp((vv+70.0)/5.0)))
            else
                delta_epi=1.0
            endif

            tauit_fast=4.562+1/(0.3933*exp((-(vv+100.0))/100.0)
     &                +0.08004*exp((vv+50.0)/16.59))
            tauit_fast=tauit_fast*delta_epi
            exptauit_fast(ii)=exp(-dt/tauit_fast) !tiF-151
!c Equation (43)

            tauit_slow=23.62+1/(0.001416*exp((-(vv+96.52))/59.05)
     &                +1.780e-8*exp((vv+114.1)/8.079)) !tiS-152
            tauit_slow=tauit_slow*delta_epi
            exptauit_slow(ii)=exp(-dt/tauit_slow)
!c Equation (44)
            Ai_fast(ii)=1.0/(1.0+exp((vv-213.6)/151.2))
!c Equation (45)
            Ai_slow(ii)=1.0-Ai_fast(ii)
!c Equation (46)
            aCaMKinft(ii)=1.0/(1.0+exp((-(vv-24.34))/14.82)) !assp-160
!c Equation (50)
         tauaCaMKt=1.0515/(1.0/(1.2089*(1.0+exp(-(vv-18.4099)/29.3814)))
     &            +3.5/(1.0+exp((vv+100.0)/29.3814)))
            exptauaCaMKt(ii)=exp(-dt/tauaCaMKt)
!c Equation (51)
            iCaMKinft(ii)=1.0/(1.0+exp((vv+43.94)/5.711))
!c Equation (53)
            DeltaCaMK_develop(ii)=1.354+1.0e-4/(exp((vv-167.4)/15.89)
     &                            +exp(-(vv-12.23)/0.2154))
!c Equation (54)
            DeltaCaMK_recover(ii)=1.0-0.5/(1.0+exp((vv+70.0)/20.0))
!c Equation (55)
            tauiCaMK_fast=tauit_fast*DeltaCaMK_develop(ii)
     &                        *DeltaCaMK_recover(ii)
            exptauiCaMK_fast(ii)=exp(-dt/tauiCaMK_fast)
!c Equation (56)
            tauiCaMK_slow=tauit_slow*DeltaCaMK_develop(ii)
     &                        *DeltaCaMK_recover(ii)
            exptauiCaMK_slow(ii)=exp(-dt/tauiCaMK_slow)
!c Equation (57)
            AiCaMK_fast(ii)=Ai_fast(ii)
!c Equation (58)
            AiCaMK_slow(ii)=Ai_slow(ii)
!c Equation (59)
            dinft(ii)=1.0/(1.0+exp((-(vv+3.940))/4.230))
!c Equation (65)
            taudt=0.6+1.0/(exp(-0.05*(vv+6.0))+exp(0.09*(vv+14.0)))
            exptaudt(ii)=exp(-dt/taudt)
!c Equation (66)
            finft(ii)=1.0/(1.0+exp((vv+19.58)/3.696))
!c Equation (68)
            tauft_fast=7.0+1.0/(0.0045*exp(-(vv+20.0)/10.0)
     &                +0.0045*exp((vv+20.0)/10.0))
            exptauft_fast(ii)=exp(-dt/tauft_fast)
!c Equation (69)
            tauft_slow=1000.0+1.0/(0.000035*exp(-(vv+5.0)/4.0)
     &                +0.000035*exp((vv+5.0)/6.0))
            exptauft_slow(ii)=exp(-dt/tauft_slow)
!c Equation (70)
            fCainft(ii)=finft(ii)
!c Equation (74)
            taufCat_fast=7.0+1.0/(0.04*exp(-(vv-4.0)/7.0)
     &                  +0.04*exp((vv-4.0)/7.0))
            exptaufCat_fast(ii)=exp(-dt/taufCat_fast)
!c Equation (75)
            taufCat_slow=100.0+1.0/(0.00012*exp(-vv/3.0)
     &                +0.00012*exp(vv/7.0))
            exptaufCat_slow(ii)=exp(-dt/taufCat_slow)
!c Equation (76)
            AfCa_fast(ii)=0.3+0.6/(1.0+exp((vv-10.0)/10.0))
!c Equation (77)
            AfCa_slow(ii)=1.0-AfCa_fast(ii)
!c Equation (78)
            jCainft(ii)=finft(ii)
!c Equation (82)
            taujCat=75.0
            exptaujCat(ii)=exp(-dt/taujCat)
!c Equation (82-83)
            fCaMKinft(ii)=finft(ii)
!c Equation (84)
            taufCaMKt_fast=2.5*tauft_fast
            exptaufCaMKt_fast(ii)=exp(-dt/taufCaMKt_fast)
!c Equation (85)
            fCa_CaMKinft(ii)=finft(ii)
!c Equation (89)
            taufCa_CaMKt_fast=2.5*taufCat_fast
            exptaufCa_CaMKt_fast(ii)=exp(-dt/taufCa_CaMKt_fast)
!c Equation (90)
            AfCa_CaMK_fast(ii)=AfCa_fast(ii)
!c Equation (91)
            AfCa_CaMK_slow(ii)=1.0-AfCa_CaMK_fast(ii)
!c Equation (92)
            Xrinft(ii)=1.0/(1.0+exp((-(vv+8.337))/6.789))
!c Equation (111)
            tauXrt_fast=12.98+1.0/(0.3652*exp((vv-31.66)/3.869)
     &                 +4.123e-5*exp((-(vv-47.78))/20.38))
            exptauXrt_fast(ii)=exp(-dt/tauXrt_fast)
!c Equation (112)
            tauXrt_slow=1.865+1.0/(0.06629*exp((vv-34.70)/7.355)
     &                 +1.128e-5*exp((-(vv-29.74))/25.94))
            exptauXrt_slow(ii)=exp(-dt/tauXrt_slow)
!c Equation (113)
            AXrt_fast(ii) = 1.0/(1.0+exp((vv+54.81)/38.21))
!c Equation (114)
            AXrt_slow(ii) = 1.0 - AXrt_fast(ii)
!c Equation (115)
            RKr(ii) = 1.0/(1.0+exp((vv+55.0)/75.0))*1.0
     &              / (1.0+exp((vv-10.0)/30.0))
!c Equation (119)
            Xs1inft(ii)=1.0/(1.0+exp((-(vv+11.60))/8.932))
!c Equation (121)
            tauXs1t = 817.3+1.0/(2.326e-4*exp((vv+48.28)/17.80)
     &              + 0.001292*exp((-(vv+210.0))/230.0))
            exptauXs1t(ii)=exp(-dt/tauXs1t)
!c Equation (122)
            Xs2inft(ii)=Xs1inft(ii)
!c Equation (124)
            tauXs2t = 1.0/(0.01*exp((vv-50.0)/20.0)
     &              + 0.0193*exp((-(vv+66.54))/31.0))
            exptauXs2t(ii)=exp(-dt/tauXs2t)
!c Equation (125)
            Xk1inft(ii) = 1.0/(1.0+exp(-(vv+2.5538*ko+144.59)
     &                  / (1.5692*ko+3.8115)))
!c Equation (128)
            tauXk1t = 122.2/(exp((-(vv+127.2))/20.36)
     &              + exp((vv+236.8)/69.33))
            exptauXk1t(ii)=exp(-dt/tauXk1t)
!c Equation (129)
            RK1(ii) = 1.0/(1.0+exp((vv+105.8-2.6*ko)/9.493))
!c Equation (131)
            hCa(ii) = exp(qca*vv*FoRT)
!c Equation (133)
            hNa(ii) = exp(qna*vv*FoRT)
!c Equation (134)
            KNai(ii)=KNai0*exp(delta*vv*FoRT/3)
!c Equation (173)
            Knao(ii)=Knao0*exp((1.0-delta)*vv*FoRT/3)
!c Equation (174)
            X_Kb(ii)=1.0/(1.0+exp(-(vv-14.48)/18.34))
!c Equation (197)
      enddo

      voltfile='Voltage.xxxx'
      APDDI='FirstSite.xxxx'
      output='fort.xxxx'

c ========================= TIME LOOP ==========================

      do icl=1,ncl

         intcl=nint(cls(icl))
         write(6,*) 'intcl= ' ,intcl
         idig1=intcl/1000
         istp=mod(intcl,1000)
         idig2=istp/100
         istp=mod(istp,100)
         idig3=istp/10
         idig4=mod(istp,10)

         voltfile(9:9)=char(ichar('0')+idig1)
         voltfile(10:10)=char(ichar('0')+idig2)
         voltfile(11:11)=char(ichar('0')+idig3)
         voltfile(12:12)=char(ichar('0')+idig4)

         APDDI(11:11)=char(ichar('0')+idig1)
         APDDI(12:12)=char(ichar('0')+idig2)
         APDDI(13:13)=char(ichar('0')+idig3)
         APDDI(14:14)=char(ichar('0')+idig4)

         output(6:6)=char(ichar('0')+idig1)
         output(7:7)=char(ichar('0')+idig2)
         output(8:8)=char(ichar('0')+idig3)
         output(9:9)=char(ichar('0')+idig4)


         cyclelength=cls(icl)
         icycle=nint(cyclelength/dt)
         ncycles=nint(endtime/cyclelength)
         if(ncycles.lt.endtime/cyclelength) then
           ncycles=ncycles+1
         endif
         nend=icycle*ncycles
         write(6,*) 'write n cycles of length x', ncycles,cyclelength

         nups1=0
         ndowns1=0
         ntime=0
         time=0.0

c time loop
c      do ntime=1,nend
 707     continue

         ntime=ntime+1
         time=ntime*dt

         if(mod(ntime,25000).eq.1) write(6,*) time


         do i=1,nx

c compute indices into tables
            zindexv=(v(i)-vlo)/dvt
            iv1=nint(zindexv)
            if(iv1.lt.0.or.iv1.gt.nvt) then
               write(6,*) 'Voltage index outside voltage table: ',
     &              ntime,ntime*dt,v(i),iv1
               stop
            endif


         if (isnan(v(i))) then
            write(87,*) cyclelength,time,v(i),Ki(i),IKs,KsCa,Cai(i),
     &      BCai,IpCa,ICab,INaCa_i,Jup,Jdiff,Nai(i),Cass(i),Cansr(i)
            stop
         endif
!         if(mod(ntime,1500).eq.1)
!         write(1115,*) cyclelength,time,v(i),
!     &      Ki(i),Cai(i),Cass(i),Jup,Jdiff,Nai(i)

!         write(1001,*) cyclelength,time,v(i),Ki(i),IKs,KsCa,Cai(i),
!     &   BCai,IpCa,ICab,INaCa_i,Jup,Jdiff,Nai(i),Cass(i),Cansr(i)

c update CaMK
         CaMKb=CaMKo*(1.0-CaMK_trap(i))/(1.0+KmCaM/Cass(i))
         CaMKa=CaMKb+CaMK_trap(i)
         dCaMK_trap=aCaMK*CaMKb*(CaMKb+CaMK_trap(i))-bCaMK*CaMK_trap(i)
         CaMK_trap(i)=CaMK_trap(i)+dt*dCaMK_trap

c reversal potentials
         ENa=RToF*log(nao/Nai(i))
         EK=RToF*log(ko/Ki(i))
         EKs=RToF*log((ko+PKNa*nao)/(Ki(i)+PKNa*Nai(i)))

c calculate INa

         xm(i) = mss(iv1)-(mss(iv1)-xm(i))*exptaumt(iv1)

         xhf(i) = hss(iv1)-(hss(iv1)-xhf(i))*exptauht_fast(iv1)
         xhs(i) = hss(iv1)-(hss(iv1)-xhs(i))*exptauht_slow(iv1)
         xh=Ahf*xhf(i)+Ahs*xhs(i)

         xj(i) = jss(iv1)-(jss(iv1)-xj(i))*exptaujt(iv1)

         hCaMK_slow(i) = hssp(iv1)-(hssp(iv1)-hCaMK_slow(i))
     &                  * exptauhCaMK_slow(iv1)
         hp=Ahf*xhf(i)+Ahs*hCaMK_slow(i)

         jCaMK(i) = jss(iv1)-(jss(iv1)-jCaMK(i))*exptaujCaMK(iv1)
         fINap=(1.0/(1.0+KmCaMK/CaMKa))
         INa=GNa*(v(i)-ENa)*xm(i)*xm(i)*xm(i)*((1.0-fINap)
     &        *xh*xj(i)+fINap*hp*jCaMK(i))

c calculate INaL

         xmL(i) = mLss(iv1)-(mLss(iv1)-xmL(i))*exptaumLt(iv1)

         xhL(i) = hLinft(iv1)-(hLinft(iv1)-xhL(i))*exptauhLt(iv1)

         xhLCaMK(i) = hLCaMKinft(iv1)-(hLCaMKinft(iv1)-xhLCaMK(i))
     &              * exptauhLCaMK(iv1)

         fINaLp=(1.0/(1.0+KmCaMK/CaMKa))

         INaL=GNaL*(v(i)-ENa)*xmL(i)*((1.0-fINaLp)*xhL(i)
     &        +fINaLp*xhLCaMK(i))

c calculate Ito

          xa(i) = ainft(iv1)-(ainft(iv1)-xa(i))*exptauat(iv1)

       i_fast(i) = iinft(iv1)-(iinft(iv1)-i_fast(i))*exptauit_fast(iv1)
       i_slow(i) = iinft(iv1)-(iinft(iv1)-i_slow(i))*exptauit_slow(iv1)

         xi = Ai_fast(iv1)*i_fast(i)+Ai_slow(iv1)*i_slow(i)

       a_CaMK(i) = aCaMKinft(iv1)-(aCaMKinft(iv1)-a_CaMK(i))
     &           * exptauaCaMKt(iv1)



       iCaMK_fast(i) = iCaMKinft(iv1)-(iCaMKinft(iv1)-iCaMK_fast(i))
     &               * exptauiCaMK_fast(iv1)
       iCaMK_slow(i) = iCaMKinft(iv1)-(iCaMKinft(iv1)-iCaMK_slow(i))
     &               * exptauiCaMK_slow(iv1)

      ip = AiCaMK_fast(iv1)*iCaMK_fast(i)+AiCaMK_slow(iv1)*iCaMK_slow(i)

         fItop=(1.0/(1.0+KmCaMK/CaMKa))
         Ito=Gto*(v(i)-EK)*((1.0-fItop)*xa(i)*xi+fItop*a_CaMK(i)*ip)

c calculate ICaL, ICaNa, ICaK

          xd(i) = dinft(iv1)-(dinft(iv1)-xd(i))*exptaudt(iv1)

          xff(i) = finft(iv1)-(finft(iv1)-xff(i))*exptauft_fast(iv1)
          xfs(i) = finft(iv1)-(finft(iv1)-xfs(i))*exptauft_slow(iv1)

         xf=Aff*xff(i)+Afs*xfs(i)

          xfcaf(i) = fCainft(iv1)-(fCainft(iv1)-xfcaf(i))
     &             * exptaufCat_fast(iv1)

          xfcas(i) = fCainft(iv1)-(fCainft(iv1)-xfcas(i))
     &             * exptaufCat_slow(iv1)

         fca=AfCa_fast(iv1)*xfcaf(i)+AfCa_slow(iv1)*xfcas(i)

          xjca(i) = jCainft(iv1)-(jCainft(iv1)-xjca(i))
     &             * exptaujCat(iv1)

          xffp(i) = fCaMKinft(iv1)-(fCaMKinft(iv1)-xffp(i))
     &             * exptaufCaMKt_fast(iv1)

         fp=Aff*xffp(i)+Afs*xfs(i)

          xfcafp(i) = fCa_CaMKinft(iv1)-(fCa_CaMKinft(iv1)-xfcafp(i))
     &             * exptaufCa_CaMKt_fast(iv1)

        fcap=AfCa_CaMK_fast(iv1)*xfcafp(i)+AfCa_CaMK_slow(iv1)*xfcas(i)

         km2n=xjca(i)*1.0
         anca=1.0/(k2n/km2n+(1.0+Kmn/Cass(i))**4)
         dnca=anca*k2n-xnca(i)*km2n
         xnca(i)=xnca(i)+dt*dnca
         PhiCaL=4.0*vffrt(iv1)*(Cass(i)*exp(2.0*vfrt(iv1))-0.341*cao)/
     &        (exp(2.0*vfrt(iv1))-1.0)
         PhiCaNa=1.0*vffrt(iv1)*(0.75*Nass(i)*exp(1.0*vfrt(iv1))
     &        -0.75*nao)/(exp(1.0*vfrt(iv1))-1.0)
         PhiCaK=1.0*vffrt(iv1)*(0.75*Kss(i)*exp(1.0*vfrt(iv1))-0.75*ko)/
     &        (exp(1.0*vfrt(iv1))-1.0)


         fICaLp=(1.0/(1.0+KmCaMK/CaMKa))

         ICaL=(1.0-fICaLp)*PCa*PhiCaL*xd(i)*(xf*(1.0-xnca(i))
     &        +xjca(i)*fca*xnca(i))
     &        +fICaLp*PCap*PhiCaL*xd(i)*(fp*(1.0-xnca(i))
     &        +xjca(i)*fcap*xnca(i))
         ICaNa=(1.0-fICaLp)*PCaNa*PhiCaNa*xd(i)*(xf*(1.0-xnca(i))
     &        +xjca(i)*fca*xnca(i))
     &        +fICaLp*PCaNap*PhiCaNa*xd(i)*(fp*(1.0-xnca(i))
     &        +xjca(i)*fcap*xnca(i))
         ICaK=(1.0-fICaLp)*PCaK*PhiCaK*xd(i)*(xf*(1.0-xnca(i))
     &        +xjca(i)*fca*xnca(i))
     &        +fICaLp*PCaKp*PhiCaK*xd(i)*(fp*(1.0-xnca(i))
     &        +xjca(i)*fcap*xnca(i))

c calculate IKr

         xrf(i) = Xrinft(iv1)-(Xrinft(iv1)-xrf(i))*exptauXrt_fast(iv1)
         xrs(i) = Xrinft(iv1)-(Xrinft(iv1)-xrs(i))*exptauXrt_slow(iv1)

         xr=AXrt_fast(iv1)*xrf(i)+AXrt_slow(iv1)*xrs(i)

         IKr=GKr*sqrt(ko/5.4)*xr*RKr(iv1)*(v(i)-EK)

c calculate IKs

         xs1(i) = Xs1inft(iv1)-(Xs1inft(iv1)-xs1(i))*exptauXs1t(iv1)
         xs2(i) = Xs2inft(iv1)-(Xs2inft(iv1)-xs2(i))*exptauXs2t(iv1)

         KsCa=1.0+0.6/(1.0+(3.8e-5/Cai(i))**1.4)

         IKs=GKs*KsCa*xs1(i)*xs2(i)*(v(i)-EKs)

c calculate IK1

         xk1(i) = Xk1inft(iv1)-(Xk1inft(iv1)-xk1(i))*exptauXk1t(iv1)

         IK1=GK1*sqrt(ko)*RK1(iv1)*xk1(i)*(v(i)-EK)

c calculate INaCa_i

         h1=1+Nai(i)/kna3*(1+hna(iv1))
         h2=(Nai(i)*hna(iv1))/(kna3*h1)
         h3=1.0/h1
         h4=1.0+Nai(i)/kna1*(1+Nai(i)/kna2)
         h5=Nai(i)*Nai(i)/(h4*kna1*kna2)
         h6=1.0/h4
         h7=1.0+nao/kna3*(1.0+1.0/hna(iv1))
         h8=nao/(kna3*hna(iv1)*h7)
         h9=1.0/h7
         h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2)
         h11=nao*nao/(h10*kna1*kna2)
         h12=1.0/h10
         k1=h12*cao*kcaon
         k2=kcaoff
         k3p=h9*wca
         k3pp=h8*wnaca
         k3=k3p+k3pp
         k4p=h3*wca/hca(iv1)
         k4pp=h2*wnaca
         k4=k4p+k4pp
         k5=kcaoff
         k6=h6*Cai(i)*kcaon
         k7=h5*h2*wna
         k8=h8*h11*wna
         x1=k2*k4*(k7+k6)+k5*k7*(k2+k3)
         x2=k1*k7*(k4+k5)+k4*k6*(k1+k8)
         x3=k1*k3*(k7+k6)+k8*k6*(k2+k3)
         x4=k2*k8*(k4+k5)+k3*k5*(k1+k8)
         E1=x1/(x1+x2+x3+x4)
         E2=x2/(x1+x2+x3+x4)
         E3=x3/(x1+x2+x3+x4)
         E4=x4/(x1+x2+x3+x4)
         allo=1.0/(1.0+(KmCaAct/Cai(i))**2)
         JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp
         JncxCa=E2*k2-E1*k1

         INaCa_i=0.8*Gncx*allo*(zna*JncxNa+zca*JncxCa)

c calculate INaCa_ss

         h1=1+Nass(i)/kna3*(1+hna(iv1))
         h2=(Nass(i)*hna(iv1))/(kna3*h1)
         h3=1.0/h1
         h4=1.0+Nass(i)/kna1*(1+Nass(i)/kna2)
         h5=Nass(i)*Nass(i)/(h4*kna1*kna2)
         h6=1.0/h4
         h7=1.0+nao/kna3*(1.0+1.0/hna(iv1))
         h8=nao/(kna3*hna(iv1)*h7)
         h9=1.0/h7
         h10=kasymm+1.0+nao/kna1*(1+nao/kna2)
         h11=nao*nao/(h10*kna1*kna2)
         h12=1.0/h10
         k1=h12*cao*kcaon
         k2=kcaoff
         k3p=h9*wca
         k3pp=h8*wnaca
         k3=k3p+k3pp
         k4p=h3*wca/hca(iv1)
         k4pp=h2*wnaca
         k4=k4p+k4pp
         k5=kcaoff
         k6=h6*Cass(i)*kcaon
         k7=h5*h2*wna
         k8=h8*h11*wna
         x1=k2*k4*(k7+k6)+k5*k7*(k2+k3)
         x2=k1*k7*(k4+k5)+k4*k6*(k1+k8)
         x3=k1*k3*(k7+k6)+k8*k6*(k2+k3)
         x4=k2*k8*(k4+k5)+k3*k5*(k1+k8)
         E1=x1/(x1+x2+x3+x4)
         E2=x2/(x1+x2+x3+x4)
         E3=x3/(x1+x2+x3+x4)
         E4=x4/(x1+x2+x3+x4)
         allo=1.0/(1.0+(KmCaAct/Cass(i))**2)
         JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp
         JncxCa=E2*k2-E1*k1

         INaCa_ss=0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa)

c calculate INaK

         P=eP/(1.0+H/Khp+Nai(i)/Knap+Ki(i)/Kxkur)
         a1=(k1p*(Nai(i)/KNai(iv1))**3)/
     &        ((1.0+Nai(i)/KNai(iv1))**3+(1.0+Ki(i)/KKi)**2-1.0)
         b1=k1m*MgADP
         a2=k2p
         b2=(k2m*(nao/Knao(iv1))**3)/((1.0+nao/Knao(iv1))**3
     &     +(1.0+ko/Kko)**2-1.0)
         a3=(k3p*(ko/Kko)**2)/((1.0+nao/Knao(iv1))**3
     &     +(1.0+ko/Kko)**2-1.0)
         b3=(k3m*P*H)/(1.0+MgATP/Kmgatp)
         a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp)
         b4=(k4m*(Ki(i)/KKi)**2)/((1.0+Nai(i)/KNai(iv1))**3
     &        +(1.0+Ki(i)/KKi)**2-1.0)
         x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2
         x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4
         x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1
         x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1
         E1=x1/(x1+x2+x3+x4)
         E2=x2/(x1+x2+x3+x4)
         E3=x3/(x1+x2+x3+x4)
         E4=x4/(x1+x2+x3+x4)
         JnakNa(i)=3.0*(E1*a3-E2*b3)
         JnakK(i)=2.0*(E4*b1-E3*a1)

         INaK=Pnak*(zna*JnakNa(i)+zk*JnakK(i))

c calculate IKb
c         xkb=1.0/(1.0+exp(-(v(i)-14.48)/18.34))

         IKb=GKb*X_Kb(iv1)*(v(i)-EK)

c calculate INab
        INab=PNab*vffrt(iv1)*(Nai(i)*exp(vfrt(iv1))-nao)
     &      /(exp(vfrt(iv1))-1.0)

c calculate ICab
         ICab=PCab*4.0*vffrt(iv1)*(Cai(i)*exp(2.0*vfrt(iv1))-0.341*cao)/
     &        (exp(2.0*vfrt(iv1))-1.0)

c calculate IpCa
         IpCa=GpCa*Cai(i)/(0.0005+Cai(i))

c calculate the stimulus current, Istim
         if (mod(ntime*dt,cyclelength).le.duration) then
            Istim=amp
         else
            Istim=0.0
         endif



c update the membrane voltage
         vt(i)=v(i)-dt*(INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1
     &        +INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+Istim)

c calculate diffusion fluxes
         JdiffNa=(Nass(i)-Nai(i))/2.0
         JdiffK=(Kss(i)-Ki(i))/2.0
         Jdiff=(Cass(i)-Cai(i))/0.2
c             dJdiff(i)=(Cass(i)-Cai(i))/0.2
c             Jdiff(i)=Jdiff(i)+dt*dJdiff(i)

c calculate ryanodione receptor calcium induced calcium release from the jsr
         a_rel=0.5*bt
         Jrel_inf=a_rel*(-ICaL)/(1.0+(1.5/Cajsr(i))**8)
         if(icelltype.eq.2) then
            Jrel_inf=Jrel_inf*1.7
         endif
         tau_rel=bt/(1.0+0.0123/Cajsr(i))

         if (tau_rel.lt.0.001) then
            tau_rel=0.001
         endif

         dJrel_NP=(Jrel_inf-Jrel_NP(i))/tau_rel
         Jrel_NP(i)=Jrel_NP(i)+dt*dJrel_NP
         btp=1.25*bt
         a_relp=0.5*btp
         Jrel_infp=a_relp*(-ICaL)/(1.0+(1.5/Cajsr(i))**8)
         if(icelltype.eq.2) then
            Jrel_infp=Jrel_infp*1.7
         endif
         tau_relp=btp/(1.0+0.0123/Cajsr(i))

         if(tau_relp.lt.0.001) then
            tau_relp=0.001
         endif

         dJrel_CaMK=(Jrel_infp-Jrel_CaMK(i))/tau_relp
         Jrel_CaMK(i)=Jrel_CaMK(i)+dt*dJrel_CaMK
         fJrel_CaMK=(1.0/(1.0+KmCaMK/CaMKa))
         Jrel=(1.0-fJrel_CaMK)*Jrel_NP(i)+fJrel_CaMK*Jrel_CaMK(i)

c calculate serca pump, ca uptake flux
         Jupnp=0.004375*Cai(i)/(Cai(i)+0.00092)
         Jupp=2.75*0.004375*Cai(i)/(Cai(i)+0.00092-0.00017)
         if(icelltype.eq.1) then
            Jupnp=Jupnp*1.3
            Jupp=Jupp*1.3
         endif
         fJupp=(1.0/(1.0+KmCaMK/CaMKa))
         Jleak=0.0039375*Cansr(i)/15.0
         Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak

c calculate tranlocation flux
         Jtr=(Cansr(i)-Cajsr(i))/100.0


c update intracellular concentrations, using buffers for Cai, Cass, Cajsr
      dNai=-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)
     &     +JdiffNa*vss/vmyo 
      Nai(i)=Nai(i)+dt*dNai

      dNass=-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa
      Nass(i)=Nass(i)+dt*dNass

      dKi=-(Ito+IKr+IKs+IK1+IKb+Istim-2.0*INaK)*Acap/(F*vmyo)
     &     +JdiffK*vss/vmyo
      Ki(i)=Ki(i)+dt*dKi

      dKss=-(ICaK)*Acap/(F*vss)-JdiffK
      Kss(i)=Kss(i)+dt*dKss

      BCai=1.0/(1.0+cmdnmax*kmcmdn/(kmcmdn+Cai(i))**2
     &     +trpnmax*kmtrpn/(kmtrpn+Cai(i))**2)

c      Cai(i)=Cai(i)+dt*(BCai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)
c     &     -Jup*vnsr/vmyo+Jdiff*vss/vmyo))
      dCai=BCai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)
     &     -((1.0-fJupp)*Jupnp+fJupp*Jupp-(0.0039375*Cansr(i)/15.0))
     &     *vnsr/vmyo+((Cass(i)-Cai(i))/0.2)*vss/vmyo)
      Cai(i)=Cai(i)+dt*dCai


cc            Csqn_b(i)=(Csqn_b(i)+dt*kon_csqn*Ca_sr(i)*Bmax_Csqn)/
cc     &           (1.0+dt*kon_csqn*Ca_sr(i)+dt*koff_csqn)
cc            dCsqn_b=(kon_csqn*Ca_sr(i)
cc     &           *(Bmax_Csqn-Csqn_b(i))-koff_csqn*Csqn_b(i))
ccc            Csqn_b(i)=Csqn_b(i)+dt*dCsqn_b


      BCass=1.0/(1.0+BSRmax*KmBSR/(KmBSR+Cass(i))**2
     &     +BSLmax*KmBSL/(KmBSL+Cass(i))**2)
      dCass=BCass*(-(ICaL-2.0*INaCa_ss)*Acap/
     &     (2.0*F*vss)+Jrel*vjsr/vss-Jdiff)
      Cass(i)=Cass(i)+dt*dCass

      dCansr=Jup-Jtr*vjsr/vnsr
      Cansr(i)=Cansr(i)+dt*dCansr

      BCajsr=1.0/(1.0+csqnmax*kmcsqn/(kmcsqn+Cajsr(i))**2)
      dCajsr=BCajsr*(Jtr-Jrel)
      Cajsr(i)=Cajsr(i)+dt*dCajsr

c      if(mod(ntime,20).eq.1) then
c         write(80,*) ntime*dt,v(i)
c         write(81,*) ntime*dt,ICaL,ICaNa,ICaK
c         write(82,*) ntime*dt,IKr,IKs,IK1
c         write(83,*) ntime*dt,INa,INaL,Ito
c         write(84,*) ntime*dt,INaCa_i,INaCa_ss,INaK
c         write(85,*) ntime*dt,IKb,INab,ICab,IpCa
c         write(86,*) ntime*dt,JdiffNa,JdiffK,Jdiff
c         write(87,*) ntime*dt,Jrel,Jup,Jleak,Jtr
c         write(88,*) ntime*dt,Cai(i),Cass(i)
c         write(89,*) ntime*dt,xj(i),hCaMK_slow(i),jCaMK(i)
c      endif

      enddo

       open(15,file=voltfile,status='unknown')
       write(15,*) time,v!,ICaL,ICaNa,INaCa_ss,IKr,IKs,i_fast,i_slow,Ito,
c     &             Cai,Cass,Jrel,Jup,Jleak,Jtr



c check for APD threshold crossings
c APD90 threshold for 0d, 30s, CL=1000ms
c Maximum of the Action Potential value = 40.35
c Resting Potential value = -87.94
c %90 of the action potential will be (-87.94*0.9) + (40.35*0.1)
c APD_90 = -79.146 + 4.035 = -75.111
      i1=1
      if(v(i1).lt.vru.and.vt(i1).ge.vr.and.nups1.eq.ndowns1) then
        nups1=nups1+1
        ups1(nups1)=ntime*dt+dt*(vr-v(i1))/(vt(i1)-v(i1))
      elseif(v(i1).gt.vr.and.vt(i1).le.vr.and.nups1.eq.ndowns1+1) then
            ndowns1=ndowns1+1
            downs1(ndowns1)=ntime*dt+dt*(vr-v(i1))/(vt(i1)-v(i1))
      endif


      do i=1,nx
         v(i)=vt(i)
      enddo

      if(ntime.lt.nend) goto 707
c     end of a time loop for a given cycle length

c finst APD's, DI's, and CV's
c site 1 (i1)
      apd=downs1(ndowns1)-ups1(ndowns1)
      prevapd=downs1(ndowns1-1)-ups1(ndowns1-1)
      di=ups1(ndowns1)-downs1(ndowns1-1)
      prevdi=ups1(ndowns1-1)-downs1(ndowns1-2)
      cl=ups1(ndowns1)-ups1(ndowns1-1)
      prevcl=ups1(ndowns1-1)-ups1(ndowns1-2)

c      enddo

      open(25,file=APDDI,status='unknown')
      write(25,*) time,di,apd,cl,prevdi,prevapd,prevcl
      close(25)

c write out restart information (also useful for S1-S2 protocol)
c      open(8,file=outfile,form='unformatted',status='unknown')

      open(35,file=output,status='unknown')
      write(35,*) v,Nai,Nass,Ki,Kss,Cai,Cass,Cansr,Cajsr,xm,xhf,
     &            xhs,xj,hCaMK_slow,jCaMK,xmL,xhL,xhLCaMK,xa,
     &            i_fast,i_slow,a_CaMK,iCaMK_fast,iCaMK_slow,xd,
     &            xff,xfs,xfcaf,xfcas,xjca,xnca,xffp,xfcafp,xrf,
     &            xrs,xs1,xs2,xk1,Jrel_NP,Jrel_CaMK,CaMK_trap,
     &            JnakNa,JnakK
      close(35)

      enddo
      close(15)
      end
c elshrif:Tusscher mme4362$ cd debug
c elshrif:debug mme4362$ ls
c elshrif:debug mme4362$ ./tusscher
c Properties...> Binary Parser ...> Mach -O 64 Parser
