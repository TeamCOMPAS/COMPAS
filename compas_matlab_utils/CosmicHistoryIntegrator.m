function [Zlist, MergerRateByRedshiftByZ, SFRfractionZ]=...
    CosmicHistoryIntegrator(filename,zlist,metallicitychoice, makeplots)
% Integrator for the binary black hole merger rate over cosmic history
% COMPAS (Compact Object Mergers: Population Astrophysics and Statistics) 
% software package
%
% USAGE: 
% [Zlist, MergerRateByRedshiftByZ]=...
%    CosmicHistoryIntegrator(filename, zlistformation, zlistmerger, zlistdetection, Msimulated, [,makeplots])
%
% INPUTS:
%   filename: name of population synthesis input file 
%           should be in COMPAS output h5 format
%   zlistformation: vector of redshifts at which the formation rate is
%   computed
%   zlistmerger:  vector of redshifts at which the merger rate is computed
%   zlistdetection: vector of redshifts at which the detection rate is
%   computed
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

if (nargin<2)
    error('Not enough input arguments.');
end;
if (nargin<4), makeplots=0; end;

%cosmology calculator
[zvec,tL]=Cosmology();           
%load COMPAS data
[M1,M2,Z,Tdelay]=DataRead(filename, max(tL));   
%metallicity-specific SFR
[Zlist,SFR,SFRfractionZ]=Metallicity(Z,zvec,metallicitychoice);       


%Consider the contribution of every simulated binary to the merger rate 
%in every redshift bin by considering when it would have to be formed to 
%merge at that redshift and normalizing by the relevant 
%metallicity-specific star formation rate
dz=zvec(2)-zvec(1);
tLmerge=tL(floor(zlist/(zvec(2)-zvec(1)))+1);
MergerRateByRedshiftByZ=zeros(length(zlist),length(Zlist));
for(i=1:length(M1)),
    Zcounter=find(Zlist==Z(i));
    zformindex=find(tL>=Tdelay(i),1);
    for(k=1:length(zlist)),
        zformindex=find(tL>=(Tdelay(i)+tLmerge(k)),1);
        MergerRateByRedshiftByZ(k,Zcounter)=...
            MergerRateByRedshiftByZ(k,Zcounter)+...
            SFR(zformindex)*SFRfractionZ(zformindex,Zcounter)/Msimulated;   
    end;
end;

if(makeplots==1),   %make a set of default plots
    MakePlots(M1,M2,Z,Tdelay,zvec,SFR,SFRfractionZ,...
        zlist,Zlist,MergerRateByRedshiftByZ);
end;

end %end of CosmicHistoryIntegrator


%Load the data stored in COMPAS output format from a file
%Select only binary black hole mergers of interest, and return the
%component masses, metallicities, and star formation to merger delay times
function [M1,M2,Z,Tdelay]=DataRead(filename, tHubble)
    global BH
    if(exist(filename, 'file')~=2), 
        error('Input file does not exist');
    end;
    data = importdata(filename, '\t', 2);
    
    colnames=data.textdata(2,:);
    Z1index=find(strcmp(colnames,'Metallicity1'));
    M1index=find(strcmp(colnames,'M1'));
    M2index=find(strcmp(colnames,'M2'));
    tcindex=find(strcmp(colnames,'tc'));
    tformindex=find(strcmp(colnames,'tform'));
    type1index=find(strcmp(colnames,'stellarType1'));
    type2index=find(strcmp(colnames,'stellarType2'));
    nonRLOF=data.data(:,find(strcmp(colnames,'RLOFSecondaryAfterCEE')))==0;
    Pessimistic=data.data(:,find(strcmp(colnames,'optimisticCEFlag')))==0;
    tc=data.data(:,tcindex);
    tform=data.data(:,tformindex);
    Ttotal=(tc+tform)*1e6;  %convert Megayears to years
    %select only binaries where both members are black holes at the last
    %step, the total delay time is less than the age of the Universe, there
    %is no Roche-lobe overflow immediately after the common-envelope phase, 
    %and the Pessimistic common-envelope conditions are met
    select=Ttotal<tHubble & nonRLOF & Pessimistic & ...
        data.data(:,type1index)==14 & data.data(:,type2index)==14;

    M1=data.data(select,M1index);
    M2=data.data(select,M2index);
    Z=data.data(select,Z1index);
    Tdelay=Ttotal(select);
end %end of DataRead

%Compute the star formation rate and lookback time (in years) 
%for an array of redshifts
function [zvec,tL]=Cosmology()
    global Mpcm
    global Mpc
    global yr
    zmax=10; dz=0.001; zvec=0:dz:zmax;
    Nz=length(zvec);
    
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


%Compute the fraction of star formation that happens in a given metallicity
%bin as a function of redshift
function [Zlist,SFR,SFRfractionZ]=Metallicity(zvec,Z)
    %M_/odot per Mpc^3 per year -- Neijssel+ 2019 preferred model 
    %would be SFR=0.015*(1+zvec).^2.7./(1+((1+zvec)/2.9).^5.6) in Madau & Dickinson, 2014, (15)
    SFR=0.01*(1+zvec).^2.77./(1+((1+zvec)/2.9).^4.7); 
    Zmean=0.035.*10.^(-0.23*zvec);
    Zmu=log(Zmean)-0.39^2/2;
    dlogZ=0.01;
    logZvec=-12:dlogZ:0;  %natural log
    dPdlogZ=1/0.39/sqrt(2*pi)*exp(-(logZvec'-Zmu).^2/2/0.39^2);
    dPdlogZ=dPdlogZ./(sum(dPdlogZ,1)*dlogZ);    %normalise
    Zrange=log(max(Z))-log(min(Z));   %ugly correction for not including tails
    PdrawZ=1/Zrange;
    minlogZindex=find(exp(logZvec)>=min(Z),1, 'first');
    maxlogZindex=find(exp(logZvec)<=max(Z),1, 'last');
    dPdlogZ(minlogZindex,:)=dPdlogZ(minlogZindex,:)+sum(dPdlogZ(1:minlogZindex,:),1)*dlogZ/(sum(Z==min(Z))/length(Z))*PdrawZ;
    dPdlogZ(maxlogZindex,:)=dPdlogZ(maxlogZindex,:)+sum(dPdlogZ(maxlogZindex:end,:),1)*dlogZ/(sum(Z==max(Z))/length(Z))*PdrawZ;
    dPdlogZ=dPdlogZ./(sum(dPdlogZ,1)*dlogZ);    %normalise
end %end of Metallicity


%Make a set of default plots
function MakePlots(M1,M2,Z,Tdelay,zvec,SFR,SFRfractionZ,...
    zlist,Zlist,MergerRateByRedshiftByZ)
    
    figure(1),colormap jet;
    plot(zlist, MergerRateByRedshiftByZ, 'LineWidth', 2), 
    legend(num2str(Zlist)),
    set(gca, 'FontSize', 20); %for labels
    xlabel('z'),
    ylabel('Pessimistic non-RLOF BBH merger rate, per Gpc^3 per yr')
    disp(['Total BBH merger rate at z=0: ', ...
        num2str(sum(MergerRateByRedshiftByZ(:,1))),...
        ' per Gpc^3 per year']);


    figure(2), colormap jet;
    scatter(M1,M2,20,log(Z)/log(10),'filled');
    set(gca, 'FontSize', 20); %for labels
    H=colorbar; H.Label.String='log_{10} metallicity'; 
    xlabel('Primary BH mass [M_o]'), ylabel('Secondary BH mass [M_o]');
    
    figure(3), colormap jet;
    scatter(M1+M2,log(Tdelay/1e6)/log(10),20,log(Z)/log(10),'filled');
    set(gca, 'FontSize', 20); %for labels
    H=colorbar; H.Label.String='log_{10} metallicity'; 
    xlabel('Total BBH mass [M_o]'), ylabel('log(Tdelay/Myr)');
    
    figure(4), colormap jet;
    plot(zvec, SFR)
    set(gca, 'FontSize', 20); %for labels
    xlabel('z'), ylabel('Star-formation rate, M_o per Gpc^3 per yr');

    
    figure(5), colormap jet;
    plot(zvec, SFRfractionZ)
    set(gca, 'FontSize', 20); %for labels
    legend(num2str(Zlist))
    xlabel('z'), ylabel('Z-specific SFR fraction');

end %end of MakePlots



