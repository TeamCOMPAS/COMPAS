function [Zlist, MergerRateByRedshiftByZ, SFR, Zweight]=...
    CosmicHistoryIntegrator(filename, zlistformation, zlistdetection, Msimulated, makeplots)
% Integrator for the binary black hole merger rate over cosmic history
% COMPAS (Compact Object Mergers: Population Astrophysics and Statistics) 
% software package
%
% USAGE: 
% [Zlist, MergerRateByRedshiftByZ]=...
%    CosmicHistoryIntegrator(filename, zlistformation, zlistmerger, [,makeplots, filename2, name1, name2])
%
% INPUTS:
%   filename: name of population synthesis input file 
%           should be in COMPAS output h5 format
%   zlistformation: vector of redshifts at which the formation rate is
%   computed
%   zlistdetection:  vector of redshifts at which the detection rate is computed
%   Msimulated: total star forming mass represented by the simulation (for
%   normalisation)
%   makeplots:  if set to 1, generates a set of useful plots (default = 0)
%
% Zlist is a vector of metallicities, taken from input file
% MergerRateByRedshiftByZ is a matrix of size length(Zlist) X length(zlist)
% which contains a merger rate of binary black holes in the given redshift 
% and metallicity bin, in units of mergers per Gpc^3 of comoving volume per
% year of source time
%
% EXAMPLE:
% CosmicHistoryIntegrator('~/Work/COMPASresults/runs/Zdistalpha1-031803.h5', ...
%   zlist, zlist, 90e6, 1, '~/Work/COMPASresults/runs/Zdist2stage-031803.h5', 'Default', '2 Stage')
% zlist=0:0.1:5;
% filename=['~/Work/COMPASresults/popsynth/runs/',...
%    '20170628-Coen-Pessimistic/AllMergers.dat'];
% [Zlist,MergerRateByRedshiftByZ]=CosmicHistoryIntegrator(filename,zlist);
% figure(1),colormap jet;
% plot(zlist, MergerRateByRedshiftByZ, 'LineWidth', 2), 
% legend(num2str(Zlist)),
% set(gca, 'FontSize', 20); %for labels
% xlabel('z'),
% ylabel('Pessimistic non-RLOF BBH merger rate, per Gpc^3 per yr')
% sum(MergerRateByRedshiftByZ(:,1)) %Total BBH merger rate at z=0
% [Zlist,MergerRateByRedshiftByZLanger]=CosmicHistoryIntegrator(filename,zlist,1);
% [Zlist,MergerRateByRedshiftByZLambert]=CosmicHistoryIntegrator(filename,zlist,2);
% plot(zlist, sum(MergerRateByRedshiftByZLambert,2), 'r', ...
% zlist, sum(MergerRateByRedshiftByZLanger,2), 'b', 'LineWidth', 2),
% set(gca, 'FontSize', 20), xlabel('z'), legend('Langer & Norman','Lambert')
% ylabel('Pessimistic non-RLOF BBH merger rate, per Gpc^3 per yr'),



%define constants
global Mpcm;
global Mpc;
global yr;
Mpcm=1*10^6 * 3.0856775807e16;  %Mpc in meters
c=299792458;		%speed of light, m/s
Mpc=Mpcm/c;         %Gpc in seconds
yr=3.15569e7;       %year in seconds

if (nargin<4)
    error('Not enough input arguments.');
end;
if (nargin<5), makeplots=0; end;

%cosmology calculator
[tL]=Cosmology(zlistformation); 
%load COMPAS data
[M1,M2,Z,Tdelay]=DataRead(filename); 
%metallicity-specific SFR
[SFR,Zlist,Zweight]=Metallicity(zlistformation,min(Z),max(Z)); 


%Consider the contribution of every simulated binary to the merger rate 
%in every redshift bin by considering when it would have to be formed to 
%merge at that redshift and normalizing by the relevant 
%metallicity-specific star formation rate
dz=zlistformation(2)-zlistformation(1);
tLmerge=tL(floor(zlistformation/dz)+1);
MergerRateByRedshiftByZ=zeros(length(zlistformation),length(Zlist));
for(i=1:length(M1)),
    Zcounter=find(Zlist>=Z(i),1);
    zformindex=find(tL>=Tdelay(i),1);
    for(k=1:length(zlistformation)),
        zformindex=find(tL>=(Tdelay(i)+tLmerge(k)),1);
        MergerRateByRedshiftByZ(k,Zcounter)=...
            MergerRateByRedshiftByZ(k,Zcounter)+...
            SFR(zformindex)*Zweight(zformindex,Zcounter)/Msimulated;   
    end;
end;

if(makeplots==1),   %make a set of default plots
    MakePlots(M1,M2,Z,Tdelay,zlistformation,Zlist,SFR,Zweight,...
        MergerRateByRedshiftByZ, 1);
end;

end %end of CosmicHistoryIntegrator


%Load the data stored in COMPAS .h5 output format from a file
%Select only double compact object mergers of interest, and return the
%component masses, metallicities, and star formation to merger delay times
function [M1,M2,Z,Tdelay]=DataRead(file)
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

%Compute the star formation rate and lookback time (in years) 
%for an array of redshifts
function [tL]=Cosmology(zvec)
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
    Dl=(1+zvec).*Dm;  %Hogg, Eq. 20
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
        dPdlogZ=1/0.39/sqrt(2*pi)*exp(-(logZvec'-Zmu).^2/2/0.39^2);
        dPdlogZ=dPdlogZ./(sum(dPdlogZ,1)*dlogZ);    %normalise
        minlogZindex=find(exp(logZvec)>=minZ,1, 'first');
        maxlogZindex=find(exp(logZvec)>=maxZ,1, 'first');
        Zrange=logZvec(maxlogZindex)-logZvec(minlogZindex);   %ugly correction for not including tails
        PdrawZ=1/Zrange;
        Zvec=exp(logZvec(minlogZindex:maxlogZindex));
        dPdlogZ(minlogZindex,:)=dPdlogZ(minlogZindex,:)+sum(dPdlogZ(1:minlogZindex,:),1)*dlogZ/(sum(Zvec==min(Zvec))/length(Zvec))*PdrawZ;
        dPdlogZ(maxlogZindex,:)=dPdlogZ(maxlogZindex,:)+sum(dPdlogZ(maxlogZindex:end,:),1)*dlogZ/(sum(Zvec==max(Zvec))/length(Zvec))*PdrawZ;
        dPdlogZ(1:minlogZindex,:)=0; dPdlogZ(maxlogZindex:size(dPdlogZ,1),:)=0;
        dPdlogZ=dPdlogZ./(sum(dPdlogZ,1)*dlogZ);    %normalise
        for(i=1:length(Zvec))
            index=find(exp(logZvec)>=Zvec(i), 1, 'first');
            Zweight(:,i)=dPdlogZ(index,:)*dlogZ;
        end;
    else    %relevant for single-metallicity runs -- just give all binaries the same unit weight
        Zvec=minZ;
        Zweight=ones(length(zvec),1);
    end;
end %end of Metallicity


%Make a set of default plots
function MakePlots(M1,M2,Z,Tdelay,zvec,Zlist,SFR,Zweight,...
        MergerRateByRedshiftByZ, fignumber)

    figure(fignumber), clf(fignumber); %,colormap jet;
    plot(zvec, sum(MergerRateByRedshiftByZ,2)*1e9, 'LineWidth', 3),  hold on;
    plot(zvec, sum(MergerRateByRedshiftByZ(:,Zlist<=0.001),2)*1e9, 'LineWidth', 1);
    plot(zvec, sum(MergerRateByRedshiftByZ(:,Zlist>0.001 & Zlist<0.01),2)*1e9, 'LineWidth', 1);
    plot(zvec, sum(MergerRateByRedshiftByZ(:,Zlist>=0.01),2)*1e9, 'LineWidth', 1); hold off;
    legend('Total rate', 'From Z<=0.001', 'From 0.001<Z<0.01', 'From Z>=0.01'),
    set(gca, 'FontSize', 20); %for labels
    xlabel('z'),
    ylabel('DCO merger rate per Gpc^3 per yr')
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
    plot(zvec, SFR'.*sum(Zweight(:,Zlist<=0.001),2)*1e9, 'LineWidth', 1);
    plot(zvec, SFR'.*sum(Zweight(:,Zlist>0.001&Zlist<0.01),2)*1e9, 'LineWidth', 1);
    plot(zvec, SFR'.*sum(Zweight(:,Zlist>=0.01),2)*1e9, 'LineWidth', 1); hold off;
    legend('Total rate', 'From Z<=0.001', 'From 0.001<Z<0.01', 'From Z>=0.01'),
    set(gca, 'FontSize', 20); %for labels
    xlabel('z'), ylabel('Star-formation rate, M_o per Gpc^3 per yr');
    
    figure(5), colormap jet;
    plot(zvec, Zweight, 'LineWidth', 3)
    set(gca, 'FontSize', 20); %for labels
    legend(num2str(Zlist))
    xlabel('z'), ylabel('Z-specific SFR weight');

end %end of MakePlots



