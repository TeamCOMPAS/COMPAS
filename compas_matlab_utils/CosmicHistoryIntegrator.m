function [SFR, Zlist, Mtlist, etalist, FormationRateByRedshiftByZ, FormationRateByRedshiftByMtByEta, ...
    MergerRateByRedshiftByZ, MergerRateByRedshiftByMtByEta, zlistdetection, Rdetections, DetectableMergerRate]=...
    CosmicHistoryIntegrator(filename, zlistformation, zmaxdetection, Msimulated, makeplots)
% Integrator for the binary black hole merger rate over cosmic history
% COMPAS (Compact Object Mergers: Population Astrophysics and Statistics) 
% software package
%
% USAGE: 
% [SFR, Zlist, Mtlist, etalist, FormationRateByRedshiftByZ, FormationRateByRedshiftByMtByEta, ...
%    MergerRateByRedshiftByZ, MergerRateByRedshiftByMtByEta, zlistdetection, Rdetections, DetectableMergerRate]]=...
%    CosmicHistoryIntegrator(filename, zlistformation, zmaxdetection, Msimulated [,makeplots])
%
% INPUTS:
%   filename: name of population synthesis input file 
%           should be in COMPAS output h5 format
%   zlistformation: vector of redshifts at which the formation rate is
%   computed
%   zmaxdetection:  maximum redshift to which the detection rate is computed
%   Msimulated: total star forming mass represented by the simulation (for
%   normalisation)
%   makeplots:  if set to 1, generates a set of useful plots (default = 0)
%
% OUTPUTS: 
%   SFR is a vector of size length(zlistformation) containing the star formation rate 
% (solar masses per Mpc^3 of comoving volume per year of source time)
%   Zlist is a vector of metallicities, taken from the COMPAS run input file
%   Mtlist is a list of total mass bins
%   etalist is a list of symmetric mass ratio bins
%   FormationRateByRedshiftByZ is a matrix of size length(zformationlist) X length(Zlist) 
% which contains a formation rate of merging double compact objects in the given redshift 
% and metallicity bin, in units of formed DCOs per Mpc^3 of comoving volume per
% year of source time
%   FormationRateByRedshiftByMtByEta is a matrix of size length(zformationlist) 
% X length(Mtlist) X length(etalist) which contains a formation rate of merging double compact objects 
% in the given redshift, total mass and eta bin, in units of formed DCOs per Mpc^3 
% of comoving volume per year of source time
%   MergerRateByRedshiftByZ is a matrix of size length(zformationlist) X length(Zlist) 
% which contains a merger rate of double compact objects in the given redshift 
% and metallicity bin, in units of mergers per Mpc^3 of comoving volume per
% year of source time
%   MergerRateByRedshiftByMtByEta is a matrix of size length(zformationlist) 
% X length(Mtlist) X length(etalist) which contains a merger rate of double compact objects 
% in the given redshift, total mass and eta bin, in units of mergers per Mpc^3 
% of comoving volume per year of source time
%   zlistdetection is a vector of redshifts at which detection rates are
% computed (a subset of zlistformation going up to zmaxdetection)
%   Rdetection is a matrix of size length(zlistdetection) X length(Mtlist) X
% length(etalist) containing the detection rate per year of observer time
% from a given redshift bin and total mass and symmetric mass ratio pixel
%   DetectableMergerRate is a matrix of the same size as Rdetection but 
% containing the intrinsic rate of detectable mergers per Mpc^3 of comoving
% volume per year of source time

%
% EXAMPLE:
% zlist=0:0.01:10;
% [SFR, Zlist, Mtlist, etalist, FormationRateByRedshiftByZ, FormationRateByRedshiftByMtByEta, ...
%    MergerRateByRedshiftByZ, MergerRateByRedshiftByMtByEta, zlistdetection, Rdetections, DetectableMergerRate]=...
% CosmicHistoryIntegrator('~/Work/COMPASresults/runs/Zdistalpha1-031803.h5', zlist, 1.5, 90e6, 1);
% figure(10), semilogy(zlist, sum(MergerRateByRedshiftByZ,2)*1e9,'LineWidth',3), set(gca,'FontSize',20),
% xlabel('Redshift z'), ylabel('Formation rate of merging DCO per Gpc^3 per yr')
% 


%define constants
global Mpcm;
global Mpc;
global yr;
Mpcm=1*10^6 * 3.0856775807e16;  %Mpc in meters
c=299792458;		%speed of light, m/s
Mpc=Mpcm/c;         %Mpc in seconds
yr=3.15569e7;       %year in seconds

if (nargin<4)
    error('Not enough input arguments.');
end;
if (nargin<5), makeplots=0; end;

%cosmology calculator
[tL,Dl,dVc]=Cosmology(zlistformation); 
%load COMPAS data
[M1,M2,Z,Tdelay,maxNS]=DataRead(filename); 
%metallicity-specific SFR
[SFR,Zlist,Zweight]=Metallicity(zlistformation,min(Z),max(Z)); 

%Consider the contribution of every simulated binary to the merger rate 
%in every redshift bin by considering when it would have to be formed to 
%merge at that redshift and normalizing by the relevant 
%metallicity-specific star formation rate
dz=zlistformation(2)-zlistformation(1);
etalist=0.01:0.01:0.25;
Mtlist=1:1:ceil(max(M1+M2));
FormationRateByRedshiftByZ=zeros(length(zlistformation),length(Zlist));
FormationRateByRedshiftByMtByEta=zeros(length(zlistformation),length(Mtlist),length(etalist));
MergerRateByRedshiftByZ=zeros(length(zlistformation),length(Zlist));
MergerRateByRedshiftByMtByEta=zeros(length(zlistformation),length(Mtlist),length(etalist));
x=zeros(size(M1));
for(i=1:length(M1)),
    Zcounter=find(Zlist>=Z(i),1);
    eta=M1(i)*M2(i)/(M1(i)+M2(i))^2;
    etaindex=ceil(eta*100);
    Mtindex=ceil(M1(i)+M2(i));
    FormationRateByRedshiftByZ(:,Zcounter)=transpose(SFR).*Zweight(:,Zcounter)/Msimulated;
    FormationRateByRedshiftByMtByEta(:,Mtindex,etaindex)=transpose(SFR).*Zweight(:,Zcounter)/Msimulated;
    tLform=tL+Tdelay(i);    %lookback time of when binary would have to form in order to merge at lookback time tL
    firsttooearlyindex=find((tLform)>max(tL),1);
    if(isempty(firsttooearlyindex)), firsttooearlyindex=length(tL)+1; end;
    zForm=interp1(tL,zlistformation,tLform(1:firsttooearlyindex-1));
    zFormindex=ceil((zForm-zlistformation(1))./dz)+1;
    if(~isempty(zFormindex))
        x(i)=SFR(zFormindex(1))*Zweight(zFormindex(1),Zcounter)/Msimulated;
        MergerRateByRedshiftByZ(1:firsttooearlyindex-1,Zcounter)=...
            MergerRateByRedshiftByZ(1:firsttooearlyindex-1,Zcounter)+...
            transpose(SFR(zFormindex)).*Zweight(zFormindex,Zcounter)/Msimulated;
        MergerRateByRedshiftByMtByEta(1:firsttooearlyindex-1,Mtindex,etaindex) =...
            MergerRateByRedshiftByMtByEta(1:firsttooearlyindex-1,Mtindex,etaindex) + ...
            transpose(SFR(zFormindex)).*Zweight(zFormindex,Zcounter)/Msimulated;
    end;
end;

zlistdetection=zlistformation(1:find(zlistformation<=zmaxdetection,1,"last"));
fin=load('~/Work/Rai/LIGOfuture_data/freqVector.txt');
%noise=load('~/Work/Rai/LIGOfuture_data/dataNomaLIGO.txt');
noise=load('~/Work/Rai/LIGOfuture_data/dataEarly_low.txt');
[Rdetections,DetectableMergerRate]=...
    DetectionRate(zlistformation,Mtlist,etalist,MergerRateByRedshiftByMtByEta,zlistdetection,fin,noise,Dl,dVc);

if(makeplots==1),   %make a set of default plots
    MakePlots(M1,M2,Z,Tdelay,zlistformation,Zlist,SFR,Zweight,...
        MergerRateByRedshiftByZ, Rdetections, DetectableMergerRate, zlistdetection, Mtlist, etalist, 1);
end;

end %end of CosmicHistoryIntegrator


%Load the data stored in COMPAS .h5 output format from a file
%Select only double compact object mergers of interest, and return the
%component masses, metallicities, and star formation to merger delay times
function [M1,M2,Z,Tdelay, maxNS]=DataRead(file)
    if(exist(file, 'file')~=2), 
        error('Input file does not exist');
    end;    
    type1=h5read(file,'/BSE_Double_Compact_Objects/Stellar_Type(1)');
    type2=h5read(file,'/BSE_Double_Compact_Objects/Stellar_Type(2)');
    mass1=h5read(file,'/BSE_Double_Compact_Objects/Mass(1)');
    mass2=h5read(file,'/BSE_Double_Compact_Objects/Mass(2)');
    seedDCO=h5read(file,'/BSE_Double_Compact_Objects/SEED');
    merges=h5read(file,'/BSE_Double_Compact_Objects/Merges_Hubble_Time');
    a=h5read(file,'/BSE_Double_Compact_Objects/SemiMajorAxis@DCO');
    e=h5read(file,'/BSE_Double_Compact_Objects/Eccentricity@DCO');
    Ttotal=(h5read(file,'/BSE_Double_Compact_Objects/Time')+h5read(file,'/BSE_Double_Compact_Objects/Coalescence_Time'))*1e6; %to years
    %mergingBBH=(type1==14) & (type2==14) & merges;
    %BBH=(type1==14) & (type2==14);
    %mergingBNS=(type1==13) & (type2==13) & merges;
    %BNS=(type1==13) & (type2==13);
    %mergingNSBH=(((type1==13) & (type2==14)) | ((type1==14) & (type2==13))) & merges;
    %NSBH=(((type1==13) & (type2==14)) | ((type1==14) & (type2==13)));
    %mergingDCO=mergingBNS | mergingNSBH | mergingBBH;
    %BNScount=sum(mergingBNS); NSBHcount=sum(mergingNSBH); BBHcount=sum(mergingBBH);
    maxNS=max(max(mass1(type1==13)), max(mass2(type2==13)));
    chirpmass=mass1.^0.6.*mass2.^0.6./(mass1+mass2).^0.2;
    q=mass2./mass1;
    seedCE=h5read(file,'/BSE_Common_Envelopes/SEED');
    [isCE,CEIndex]=ismember(seedDCO,seedCE);
    optCE=h5read(file,'/BSE_Common_Envelopes/Optimistic_CE');
    RLOFCE=h5read(file,'/BSE_Common_Envelopes/Immediate_RLOF>CE');
    OKCE=zeros(size(seedDCO)); OKCE(CEIndex==0)=1; OKCE(CEIndex>0)=(~optCE(CEIndex(CEIndex>0))) & (~RLOFCE(CEIndex(CEIndex>0)));
    %BNSCE=sum(mergingBNS & isCE & OKCE); NSBHCE=sum(mergingNSBH & isCE & OKCE); BBHCE=sum(mergingBBH & isCE & OKCE);
    mergingDCO=merges & OKCE;
    Zsys=h5read(file,'/BSE_System_Parameters/Metallicity@ZAMS(1)');
    seedsys=h5read(file,'/BSE_System_Parameters/SEED');
    [blah,sysIndex]=ismember(seedDCO,seedsys);
    Zdco=Zsys(sysIndex);
    M1=mass1(mergingDCO); M2=mass2(mergingDCO); Z=Zdco(mergingDCO); Tdelay=Ttotal(mergingDCO);
end %end of DataRead

%Compute the lookback time (yr), lumionosity distance (Mpc), and comoving
%volume (Mpc^3) for an array of redshifts
function [tL, Dl, dVc]=Cosmology(zvec)
    global Mpcm
    global Mpc
    global yr
    %zmax=10; dz=0.001; zvec=0:dz:zmax;
    Nz=length(zvec); dz=zvec(2)-zvec(1);
    
    %Planck cosmology
    OmegaM=0.236+0.046;  %2012arXiv1212.5226H
    OmegaL=0.718;
    Ho=100*0.697*1000/Mpcm; %in sec; H=69.7
    Dh=1/Ho;
    E=sqrt(OmegaM.*(1+zvec).^3+OmegaL);	%Hogg, astro-ph/9905116, Eq. 14
    Dc=Dh*dz*cumsum(1./E); %Hogg, Eq. 15
    Dm=Dc;	%Hogg, Eq. 16, k=0;
    Dl=(1+zvec).*Dm/Mpc;  %Hogg, Eq. 20
    %see also Eq. (1.5.46) in Weinberg, "Cosmology", 2008
    dVc=4*pi*Dh^3*(OmegaM*(1+zvec).^3+OmegaL).^(-0.5).*(Dc/Dh).^2*dz/Mpc^3;
    Vc=cumsum(dVc);
    dtL=(1/Ho)*dz./(1+zvec)./E/yr;  %lookback time, (30) of Hogg
    tL=cumsum(dtL);
end %end of Cosmology


%Compute the weight of each star-forming metallicity as a function of redshift
function [SFR,Zvec,Zweight]=Metallicity(zvec,minZ,maxZ)
    %M_/odot per Mpc^3 per year -- Neijssel+ 2019 preferred model 
    %would be SFR=0.015*(1+zvec).^2.7./(1+((1+zvec)/2.9).^5.6) in Madau & Dickinson, 2014, (15)
    SFR=0.01*(1+zvec).^2.77./(1+((1+zvec)/2.9).^4.7); 
    if(maxZ>minZ),
        Zmean=0.035.*10.^(-0.23*zvec);
        Zmu=log(Zmean)-0.39^2/2;
        dlogZ=0.1;
        logZvec=-12:dlogZ:0;  %natural log
        dPdlogZ=1/0.39/sqrt(2*pi)*exp(-(logZvec'-Zmu).^2/2/0.39^2); %size length(logZvec) x length(zvec)
        dPdlogZ=dPdlogZ./(sum(dPdlogZ,1)*dlogZ);    %normalise
        minlogZindex=find(exp(logZvec)>=minZ,1, 'first');
        maxlogZindex=find(exp(logZvec)>=maxZ,1, 'first');
        Zrange=logZvec(maxlogZindex)-logZvec(minlogZindex);   %ugly correction for not including tails
        PdrawZ=1/Zrange;
        Zvec=exp(logZvec(minlogZindex:maxlogZindex));
        Zweight=zeros(length(zvec),length(Zvec));
        dPdlogZ(minlogZindex,:)=dPdlogZ(minlogZindex,:)+sum(dPdlogZ(1:minlogZindex-1,:),1)*dlogZ/dlogZ;
        dPdlogZ(maxlogZindex,:)=dPdlogZ(maxlogZindex,:)+sum(dPdlogZ(maxlogZindex+1:end,:),1)*dlogZ/dlogZ;
        dPdlogZ(1:minlogZindex-1,:)=0; dPdlogZ(maxlogZindex+1:size(dPdlogZ,1),:)=0;
        dPdlogZ=dPdlogZ./(sum(dPdlogZ,1)*dlogZ);    %normalise
        for(i=1:length(Zvec))
            index=find(exp(logZvec)>=Zvec(i), 1, 'first');
            Zweight(:,i)=dPdlogZ(index,:)./PdrawZ; %weight is desired over sampled probabiluty
        end;
    else    %relevant for single-metallicity runs -- just give all binaries the same unit weight
        Zvec=minZ;
        Zweight=ones(length(zvec),1);
    end;
end %end of Metallicity


%Compute detection rates per year of observer time and per year of source time
%per Mpc^3 of comoving volume as a function of total mass and eta
function [Rdetections, DetectableMergerRate]=...
    DetectionRate(zlistformation,Mtlist,etalist,MergerRateByRedshiftByMtByEta,zlistdetection,freqfile,noisefile,Dl,dVc)

    fin=load('~/Work/Rai/LIGOfuture_data/freqVector.txt');
    noise=load('~/Work/Rai/LIGOfuture_data/dataMid_low.txt');

    flow=10;
    df=1;
    f=flow:df:500; %BBH focussed
    Sf=interp1(fin, noise.^2, f);

    Ntheta=1e6;
    psi=rand(1,Ntheta)*pi;
    phi=rand(1,Ntheta)*2*pi;
    costh=rand(1,Ntheta);
    cosiota=rand(1,Ntheta);
    sinth=sqrt(1-costh.^2);
    siniota=sqrt(1-cosiota.^2);
    Fplus=1/2*(1+costh.^2).*cos(2*phi).*cos(2*psi)-costh.*sin(2*phi).*sin(2*psi);
    Fcross=1/2*(1+costh.^2).*cos(2*phi).*sin(2*psi)+costh.*sin(2*phi).*cos(2*psi);
    Theta=1/2*sqrt(Fplus.^2.*(1+cosiota.^2).^2+4*Fcross.^2.*cosiota.^2);
    Thetas=sort(Theta);


    %save time by not doing calculations beyond maximum redshifted total
    %mass corresponding to detection redshift threshold
    Mtzlistdetection=1:1:ceil(max(Mtlist)*(1+max(zlistdetection)));
    SNRat1Mpc=zeros(length(Mtzlistdetection),length(etalist));
    for(i=1:length(Mtzlistdetection)),
        for(j=1:length(etalist)),
            [h,Am,psi]=IMRSAWaveform(f, Mtzlistdetection(i), etalist(j), 0, 0, 0, 1, flow);
            integral=sum(4*Am.^2./Sf*df);
            SNRat1Mpc(i,j)=sqrt(integral);
        end;
    end;

    SNR=zeros(length(zlistdetection),length(Mtlist),length(etalist));

    for(i=1:length(zlistdetection)),
        for(j=1:length(Mtlist)),
            SNR(i,j,:)=SNRat1Mpc(ceil(j*(zlistdetection(i)+1)),:)./Dl(i);
        end;
    end;
    %for(i=1:length(zlistdetection)), SNR(i,:,:)=SNRat1Mpc./Dl(i); end;

    SNR8pre=1:0.1:1000;
    theta=1./SNR8pre;
    pdetect=1-interp1([0,Thetas,1],[(0:Ntheta)/Ntheta,1],theta);
    pdetect(1)=0;   %set of measure zero to exceed threshold, but enforce just in case

    Rdetections=zeros(length(zlistdetection),length(Mtlist),length(etalist));          %Detections per unit observer time
    DetectableMergerRate=zeros(length(zlistdetection),length(Mtlist),length(etalist)); %Detections per unit source time per unit Vc
    SNR8=SNR/8;
    pdetection=zeros(size(Rdetections));
    pdetection=pdetect(max(min(floor(SNR8*10),length(pdetect)),1));

    DetectableMergerRate=MergerRateByRedshiftByMtByEta(1:length(zlistdetection),:,:).*pdetection;
    Rdetections=DetectableMergerRate.*transpose(dVc(1:length(zlistdetection)))./(1+zlistdetection');

end %end of DetectionRate

%Make a set of default plots
function MakePlots(M1,M2,Z,Tdelay,zlistformation,Zlist,SFR,Zweight,...
        MergerRateByRedshiftByZ, Rdetections, DetectableMergerRate, zlistdetection, Mtlist, etazlist, fignumber)

    zvec=zlistformation;

    figure(fignumber), clf(fignumber); %,colormap jet;
    plot(zvec, sum(MergerRateByRedshiftByZ,2)*1e9, 'LineWidth', 3),  hold on;
    plot(zvec, sum(MergerRateByRedshiftByZ(:,Zlist<=0.001),2)*1e9, 'LineWidth', 1);
    plot(zvec, sum(MergerRateByRedshiftByZ(:,Zlist>0.001 & Zlist<0.01),2)*1e9, 'LineWidth', 1);
    plot(zvec, sum(MergerRateByRedshiftByZ(:,Zlist>=0.01),2)*1e9, 'LineWidth', 1); hold off;
    set(gca, 'FontSize', 20); %for labels
    xlabel('z'),
    ylabel('DCO merger rate per Gpc^3 per yr')
    legend('Total rate', 'From Z<=0.001', 'From 0.001<Z<0.01', 'From Z>=0.01'),
    disp(['Total DCO merger rate at z=0: ', ...
        num2str(1e9*sum(MergerRateByRedshiftByZ(1,:))),...
        ' per Gpc^3 per year']);

    figure(2);
    colormap jet;
    scatter(log10(M1),log10(M2),20,log(Z)/log(10),'filled');
    set(gca, 'FontSize', 20); %for labels
    H=colorbar; H.Label.String='log_{10} metallicity'; 
    xlabel('log_{10}(M_1/M_o)'), ylabel('log_{10}(M_2/M_o)');
    
    figure(3);
    colormap jet;
    scatter(M1+M2,log10(Tdelay/1e6),20,log10(Z),'filled');
    set(gca, 'FontSize', 20); %for labels
    H=colorbar; H.Label.String='log_{10} metallicity'; 
    xlabel('Total DCO mass [M_o]'), ylabel('log_{10}(Tdelay/Myr)');

    
    figure(fignumber+3), clf(fignumber+3);
    plot(zvec, SFR*1e9, 'LineWidth', 3); hold on;
    plot(zvec, SFR'.*sum(Zweight(:,Zlist<=0.001),2)./sum(Zweight,2)*1e9, 'LineWidth', 1);
    plot(zvec, SFR'.*sum(Zweight(:,Zlist>0.001&Zlist<0.01),2)./sum(Zweight,2)*1e9, 'LineWidth', 1);
    plot(zvec, SFR'.*sum(Zweight(:,Zlist>=0.01),2)./sum(Zweight,2)*1e9, 'LineWidth', 1); hold off;
    set(gca, 'FontSize', 20); %for labels
    xlabel('z'), ylabel('Star-formation rate, M_o per Gpc^3 per yr');
    legend('Total rate', 'From Z<=0.001', 'From 0.001<Z<0.01', 'From Z>=0.01'),


    figure(fignumber+4), clf(fignumber+4);
    RdetectionsByzMt=sum(Rdetections,3); %sum across eta
    semilogy(zlistdetection, cumsum(sum(RdetectionsByzMt,2)), 'LineWidth', 3),  hold on;
    semilogy(zlistdetection, cumsum(sum(RdetectionsByzMt(:,Mtlist<=5),2)), 'LineWidth', 1);
    semilogy(zlistdetection, cumsum(sum(RdetectionsByzMt(:,Mtlist>5 & Mtlist<20),2)), 'LineWidth', 1);
    semilogy(zlistdetection, cumsum(sum(RdetectionsByzMt(:,Mtlist>=20),2)), 'LineWidth', 1); hold off;
    legend('Total rate', 'From M_t<=5 M_o', 'From 5<M_t/M_o<20', 'From M_t>=20 M_o'),
    set(gca, 'FontSize', 20); %for labels
    xlabel('z'),
    ylabel('cumulative detection rate per observer yr')


end %end of MakePlots



