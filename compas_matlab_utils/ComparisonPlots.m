function ComparisonPlots(filename1, name1, filename2, name2)
% Carries out some basic analysis and makes plots comparing two COMPAS runs
%
% USAGE: 
% ComparisonPlots(filename1, name1, filename2, name2)
%
% INPUTS:
%   filename1: name of population synthesis input file 1 in COMPAS h5 format
%   name1: name of data set 1 for plot legends
%   filename2: name of population synthesis input file 2 in COMPAS h5 format
%   name2: name of data set 2 for plot legends
%   
% the last two arguments are optional; as single output can be plotted 
% by running ComparisonPlots(filename1, name1)
%
% example: 
%       ComparisonPlots('~/Work/COMPASresults/runs/Zsolaralpha1-031803.h5', 'Default', ...
%       '~/Work/COMPASresults/runs/Zsolar2stage-031803.h5', '2 stage')
%
% Warning: --switch-log must be used for all runs to be analysed
% It is recommended, but not required, to use the same random seed for the 
% runs being compared in order to compare individual binary evolution


    global Msunkg G AU Rsun
    Msunkg=1.98892e30;	%Msun in kg
    G=6.67428e-11;		%G in m^3/kg/s^2
    AU=149597871e3;     %m
    Rsun = 695500000;   %m

    %Plot DCO mass distribution, BNS P-e distribution, chirp mass vs period
    %at DCO formation, BH mass vs secondary core mass for 2->1 CEs leading
    %to merging BBH formation
    figure(1); clf(1); figure(2); clf(2); figure(3); clf(3); figure(4); clf(4); figure(5); clf(5); 
    [BNS,NSBH,BBH,BNSCE,NSBHCE,BBHCE,CEBBH1]=DCOplot(filename1, name1, 1, 'r', 40);
    fprintf('\nDCOs:\t\t#Merging DNS\t#Merging NSBH\t#Merging BBH\t%% BNS via CE\t%% NSBH via CE\t%% BBH via CE\t%% BBH via CE with BH primary\n');
    fprintf('%s:\t%d\t\t%d\t\t%d\t\t%.0f\t\t%.0f\t\t%.0f\t\t%.0f\n', ...
        name1, BNS, NSBH, BBH, BNSCE/BNS*100, NSBHCE/NSBH*100, BBHCE/BBH*100, CEBBH1/BBH*100);
    if(nargin==4),
        [BNS,NSBH,BBH,BNSCE,NSBHCE,BBHCE,CEBBH1]=DCOplot(filename2, name2, 1, 'b', 20);
        fprintf('%s:\t%d\t\t%d\t\t%d\t\t%.0f\t\t%.0f\t\t%.0f\t\t%.0f\n', ...
            name2, BNS, NSBH, BBH, BNSCE/BNS*100, NSBHCE/NSBH*100, BBHCE/BBH*100, CEBBH1/BBH*100);
    end;
    figure(1), hold off; figure(2), hold off; figure(3), hold off; figure(4), hold off;
    figure(1), set(gca,'FontSize',20), xlabel('$M_1$ (M$_\odot$)', 'Interpreter', 'latex'), 
        ylabel('$M_2$ (M$_\odot$)', 'Interpreter', 'latex'), title('Merging DCO masses');  legend;
    figure(2), set(gca,'FontSize',20), xlabel('$\log_{10}$ (Orbital period/hr)', 'Interpreter', 'latex');
    ylabel('Eccentricity'), title('DNS at formation'); legend;
    figure(3), set(gca,'FontSize',20), xlabel('Chirp Mass (M$_\odot$)', 'Interpreter', 'latex'), 
        ylabel('$\log_{10} (P_\mathrm{orb}/\mathrm{d})$', 'Interpreter', 'latex'), title('Merging BBH at formation'); legend;
    figure(4), set(gca,'FontSize',20); xlabel('$M_1$ (M$_\odot$)', 'Interpreter', 'latex'), 
        ylabel('$M_{\textrm{core},2}$ (M$_\odot$)', 'Interpreter', 'latex'), title('CE from 2->1 en route to merging BBH'), legend;
    figure(5), set(gca,'FontSize',20); xlabel('$q \equiv M_2/M_1$', 'Interpreter', 'latex'), 
        ylabel('CDF'), title('CDF of merging BBH'), legend;        
    fignumber=6;

    %Plot BH HMXBs
    figure(fignumber); clf(fignumber);
    HMXBplot(filename1, name1, fignumber, 'r', 40);
    if(nargin==4),
        HMXBplot(filename2, name2, fignumber, 'b', 20);
    end;
    figure(fignumber); hold off;  set(gca,'FontSize',20), legend; title('HMXB masses');
    xlabel('BH mass (M$_\odot$)', 'Interpreter', 'latex'), ylabel('Companion mass (M$_\odot$)', 'Interpreter', 'latex');
    fignumber=fignumber+1;

    %Plot BeXRBs
    figure(fignumber); clf(fignumber);
    BeXRBplot(filename1, name1, fignumber, 'r', 40);
    if(nargin==4),
        BeXRBplot(filename2, name2, fignumber, 'b', 20);
    end;
    figure(fignumber); hold off;  set(gca,'FontSize',20), legend;
    xlabel('Companion mass (M$_\odot$)', 'Interpreter', 'latex');
    ylabel('Formation time, Myr'), title('BeXRBs just after SN');
    fignumber=fignumber+1;

    %Plot LMXBs/IMXBs
    figure(fignumber); clf(fignumber);
    [LMXBcount, NSLMXBcount]=LMXBplot(filename1, name1, fignumber, 'r', 40);
    fprintf('\nLMXBs:\t\t#LMXB\t\t#NS LMXB\n');
    fprintf('%s:\t%d\t\t%d\n', name1, LMXBcount, NSLMXBcount);
    if(nargin==4),
        [LMXBcount, NSLMXBcount]=LMXBplot(filename2, name2, fignumber, 'b', 20);
        fprintf('%s:\t%d\t\t%d\n', name2, LMXBcount, NSLMXBcount);
    end;
    figure(fignumber), hold off; set(gca,'FontSize', 20); legend; title('LMXB on first MT onto CO');
    xlabel('Compact object mass (M$_\odot$)', 'Interpreter', 'latex'), ylabel('Companion mass (M$_\odot$)', 'Interpreter', 'latex'); 
    fignumber=fignumber+1;

    %Plot DWDs (just as a sanity check)
    figure(fignumber); clf(fignumber);
    DWDplot(filename1, name1, fignumber, 'r', 'm', 10);
    if(nargin==4),
        DWDplot(filename2, name2, fignumber, 'b', 'g', 5);
    end;
    figure(fignumber); hold off;  axis([-3 5 -3 5]); set(gca,'FontSize',20); title('Double White Dwarfs'); legend;
    xlabel('$log_{10}(M*a) [M_\odot * R_\odot]$ @ ZAMS', 'Interpreter', 'latex'); 
    ylabel('$log_{10}(M*a) [M_\odot * R_\odot]$ @ end', 'Interpreter', 'latex');


    [binariescount1, SNcount1, BHcompletecount1, SNbothcount1, SNonecount1, ...
            unboundcount1, unbound1count1, unbound2count1, ...
            USSNcount1, ECSNcount1, AICcount1, PISNcount1, PPISNcount1, ...
            MTSNcount1, strippedSNcount1, strippednonECSNcount1, ...
            strippedwindscount1, strippedRLOFcount1, strippedCEcount1, CESNcount1]=SNstats(filename1);
    if(nargin==4),
        [binariescount2, SNcount2, BHcompletecount2, SNbothcount2, SNonecount2, ...
            unboundcount2, unbound1count2, unbound2count2, ...
            USSNcount2, ECSNcount2, AICcount2, PISNcount2, PPISNcount2, ...
            MTSNcount2, strippedSNcount2, strippednonECSNcount2, ...
            strippedwindscount2, strippedRLOFcount2, strippedCEcount2, CESNcount2]=SNstats(filename2);
    end;
    if(nargin==2),
    fprintf('\nSupernova statistics:\t\t\t%s\n\n', name1);
    fprintf('Total number of binaries simulated:\t%d\n', binariescount1);
    fprintf('Total number of supernovae:\t\t%d\n', SNcount1);
    fprintf('Of these, complete collapse to BH:\t%d\n', BHcompletecount1);
    fprintf('Number of binaries with two SNe:\t%d\n', SNbothcount1);
    fprintf('Number of binaries with one SN:\t\t%d\n', SNonecount1);
    fprintf('Number of binaries unbound by SN:\t%d\n', unboundcount1);
    fprintf('Unbound by first SN:\t\t\t%d\n', unbound1count1);
    fprintf('Unbound by second SN:\t\t\t%d\n\n', unbound2count1);
    fprintf('Total number of USSN:\t\t\t%d\n', USSNcount1);
    fprintf('Total number of ECSN:\t\t\t%d\n', ECSNcount1);
    fprintf('Total number of AIC:\t\t\t%d\n', AICcount1);
    fprintf('Total number of PISN:\t\t\t%d\n', PISNcount1);
    fprintf('Total number of PPISN:\t\t\t%d\n\n', PPISNcount1);
    fprintf('#SNe from mass-transfering progenitors:\t%d\n', MTSNcount1);
    fprintf('#Stripped-envelope SNe:\t\t\t%d\n', strippedSNcount1);
    fprintf('Of these, not AIC or ECSN:\t\t%d\n', strippednonECSNcount1);
    fprintf('Of these, stripped by winds, no RLOF:\t%d\n', strippedwindscount1);
    fprintf('Of these, # stripped by RLOF:\t\t%d\n', strippedRLOFcount1);
    fprintf('Of these, previous CE:\t\t\t%d\n', strippedCEcount1);
    fprintf('Or double-core CE simultaneous with SN:\t%d\n', CESNcount1);
    elseif(nargin==4),
    fprintf('\nSupernova statistics:\t\t\t%s\t\t%s\n\n', name1, name2);
    fprintf('Total number of binaries simulated:\t%d\t\t%d\n', binariescount1, binariescount2);
    fprintf('Total number of supernovae:\t\t%d\t\t%d\n', SNcount1, SNcount2);
    fprintf('Of these, complete collapse to BH:\t%d\t\t%d\n', BHcompletecount1, BHcompletecount2);
    fprintf('Number of binaries with two SNe:\t%d\t\t%d\n', SNbothcount1, SNbothcount2);
    fprintf('Number of binaries with one SN:\t\t%d\t\t%d\n', SNonecount1, SNonecount2);
    fprintf('Number of binaries unbound by SN:\t%d\t\t%d\n', unboundcount1, unboundcount2);
    fprintf('Unbound by first SN:\t\t\t%d\t\t%d\n', unbound1count1, unbound1count2);
    fprintf('Unbound by second SN:\t\t\t%d\t\t%d\n\n', unbound2count1, unbound2count2);
    fprintf('Total number of USSN:\t\t\t%d\t\t%d\n', USSNcount1, USSNcount2);
    fprintf('Total number of ECSN:\t\t\t%d\t\t%d\n', ECSNcount1, ECSNcount2);
    fprintf('Total number of AIC:\t\t\t%d\t\t%d\n', AICcount1, AICcount2);
    fprintf('Total number of PISN:\t\t\t%d\t\t%d\n', PISNcount1, PISNcount2);
    fprintf('Total number of PPISN:\t\t\t%d\t\t%d\n\n', PPISNcount1, PPISNcount2);
    fprintf('#SNe from mass-transfering progenitors:\t%d\t\t%d\n', MTSNcount1, MTSNcount2);
    fprintf('#Stripped-envelope SNe:\t\t\t%d\t\t%d\n', strippedSNcount1, strippedSNcount2);
    fprintf('Of these, not AIC or ECSN:\t\t%d\t\t%d\n', strippednonECSNcount1, strippednonECSNcount2);
    fprintf('Of these, stripped by winds, no RLOF:\t%d\t\t%d\n', strippedwindscount1, strippedwindscount2);
    fprintf('Of these, # stripped by RLOF:\t\t%d\t\t%d\n', strippedRLOFcount1, strippedRLOFcount2);
    fprintf('Of these, previous CE:\t\t\t%d\t\t%d\n', strippedCEcount1, strippedCEcount2);
    fprintf('Or double-core CE simultaneous with SN:\t%d\t\t%d\n', CESNcount1, CESNcount2);        
    end;

end %end of ComparisonPlots


%Plot double compact objects; returns DCO counts,  and counts of CEs leading
%to DCOs
function [BNScount, NSBHcount, BBHcount, BNSCE, NSBHCE, BBHCE, CEBBH1count] = ...
        DCOplot(file, name, fignumber, colour, point)
    global Msunkg G AU
    type1=h5read(file,'/BSE_Double_Compact_Objects/Stellar_Type(1)');
    type2=h5read(file,'/BSE_Double_Compact_Objects/Stellar_Type(2)');
    mass1=h5read(file,'/BSE_Double_Compact_Objects/Mass(1)');
    mass2=h5read(file,'/BSE_Double_Compact_Objects/Mass(2)');
    seedDCO=h5read(file,'/BSE_Double_Compact_Objects/SEED');
    merges=h5read(file,'/BSE_Double_Compact_Objects/Merges_Hubble_Time');
    a=h5read(file,'/BSE_Double_Compact_Objects/SemiMajorAxis@DCO');
    e=h5read(file,'/BSE_Double_Compact_Objects/Eccentricity@DCO');
    mergingBBH=(type1==14) & (type2==14) & merges;
    BBH=(type1==14) & (type2==14);
    mergingBNS=(type1==13) & (type2==13) & merges;
    BNS=(type1==13) & (type2==13);
    mergingNSBH=(((type1==13) & (type2==14)) | ((type1==14) & (type2==13))) & merges;
    NSBH=(((type1==13) & (type2==14)) | ((type1==14) & (type2==13)));
    mergingDCO=mergingBNS | mergingNSBH | mergingBBH;
    BNScount=sum(mergingBNS); NSBHcount=sum(mergingNSBH); BBHcount=sum(mergingBBH);
    chirpmass=mass1.^0.6.*mass2.^0.6./(mass1+mass2).^0.2;
    q=mass2./mass1;
    seedCE=h5read(file,'/BSE_Common_Envelopes/SEED');
    [isCE,CEIndex]=ismember(seedDCO,seedCE);
    optCE=h5read(file,'/BSE_Common_Envelopes/Optimistic_CE');
    RLOFCE=h5read(file,'/BSE_Common_Envelopes/Immediate_RLOF>CE');
    OKCE=zeros(size(mergingDCO)); OKCE(CEIndex==0)=1; OKCE(CEIndex>0)=(~optCE(CEIndex(CEIndex>0))) & (~RLOFCE(CEIndex(CEIndex>0)));
    BNSCE=sum(mergingBNS & isCE & OKCE); NSBHCE=sum(mergingNSBH & isCE & OKCE); BBHCE=sum(mergingBBH & isCE & OKCE);
    type1CE = h5read(file,'/BSE_Common_Envelopes/Stellar_Type(1)<CE');
    type2CE = h5read(file,'/BSE_Common_Envelopes/Stellar_Type(2)<CE');
    mass1CE = h5read(file,'/BSE_Common_Envelopes/Mass(1)<CE');
    mass2coreCE = h5read(file,'/BSE_Common_Envelopes/Mass(2)<CE')-h5read(file,'/BSE_Common_Envelopes/Mass_Env(2)');
    [CEtowardmergingBBH,BBHindex]=ismember(seedCE,seedDCO(mergingBBH));  %pick CEs leading to merging BBHs
    %fprintf('Validation: %d valid CE episodes en route to merging BBH, including %d unique CEs, which should match %d merging BBHs produced via CE\n', ...
    %    sum(CEtowardmergingBBH & ~optCE & ~RLOFCE), length(unique(seedCE(CEtowardmergingBBH & ~optCE & ~RLOFCE))), sum(mergingBBH & isCE & OKCE))
    CEBBH1 = CEtowardmergingBBH & ~optCE & ~RLOFCE & type1CE == 14;  CEBBH1count=sum(CEBBH1);
    %Merging DCO mass distribution
    figure(fignumber), hold on;
    scatter(mass1(mergingDCO & isCE & OKCE), mass2(mergingDCO & isCE & OKCE), point, ...
        'filled', colour, 'DisplayName', ['CE, ', name]); 
    scatter(mass1(mergingDCO & ~isCE), mass2(mergingDCO & ~isCE), point, colour,  'DisplayName', ['Stable, ', name]);
    %P-e of Galactic DNS at formation
    P=2*pi*(a*AU).^(3/2)./sqrt(G*Msunkg*(mass1+mass2))/3600; %orbital period at DCO, in hours
    figure(fignumber+1), hold on; 
    scatter(log10(P(BNS & isCE & OKCE)), e(BNS & isCE & OKCE), point, 'filled', colour, ...
         'DisplayName', ['CE, ', name]); hold on;
    scatter(log10(P(BNS & ~isCE)), e(BNS & ~isCE), point, colour,  'DisplayName', ['Stable, ', name]);
    %Chirp mass vs period at BBH formation
    figure(fignumber+2), hold on;
    scatter(chirpmass(mergingDCO & isCE & OKCE), log10(P(mergingDCO & isCE & OKCE)/24), point, ...
        'filled', colour, 'DisplayName', ['CE, ', name]); 
    scatter(chirpmass(mergingDCO & ~isCE), log10(P(mergingDCO & ~isCE)/24), point, colour,  'DisplayName', ['Stable, ', name]);
    %Masses at the time of CE
    figure(fignumber+3), hold on;
    scatter(mass1CE(CEBBH1), mass2coreCE(CEBBH1), point, 'filled', colour, 'DisplayName', name);
    %CSF of q distribution
    figure(fignumber+4), hold on;
    p=cdfplot(q(mergingBBH & isCE & OKCE)), p.Color=colour, p.LineStyle='-', p.LineWidth=point/10, p.DisplayName=['CE, ', name];
    p=cdfplot(q(mergingBBH & ~isCE)), p.Color=colour, p.LineStyle=':', p.LineWidth=point/10, p.DisplayName=['Stable, ', name];
end %end of DCOplot

%Plot BH HMXBs
%BH HMXBs, according to Hirai & Mandel, consist of a BH and a MS O star
%that is more than 80% Roche lobe filling.
function HMXBplot(file, name, fignumberHMXB, colour, point)
    %We can look for binaries that experience Roche lobe overflow from an O
    %star onto a BH, since these must have been HMXBs just previously;
    isinRLOF1RLOF=h5read(file,'/BSE_RLOF/RLOF(1)>MT');
    isinRLOF2RLOF=h5read(file,'/BSE_RLOF/RLOF(2)>MT');
    wasinRLOF1RLOF=h5read(file,'/BSE_RLOF/RLOF(1)<MT');
    wasinRLOF2RLOF=h5read(file,'/BSE_RLOF/RLOF(2)<MT');
    M1RLOF=h5read(file,'/BSE_RLOF/Mass(1)<MT');
    M2RLOF=h5read(file,'/BSE_RLOF/Mass(2)<MT');
    star1RLOF=h5read(file,'/BSE_RLOF/Stellar_Type(1)<MT');
    star2RLOF=h5read(file,'/BSE_RLOF/Stellar_Type(2)<MT');
    %(don't expect any where second-born-star is BH, but check just in case)
    relevantbinaryFrom1=(isinRLOF1RLOF & ~wasinRLOF1RLOF & ~wasinRLOF2RLOF & star1RLOF==1 & M1RLOF>15 & star2RLOF==14);
    %if this is zero (no such systems), can ignore this possibility, focus on the other one
    if(sum(relevantbinaryFrom1)>0) fprintf('Number of HMXBs with second-born BH is non-zero (%d), please check\n', sum(relevantbinaryFrom1)); end;
    relevantbinaryFrom2=(isinRLOF2RLOF & ~wasinRLOF1RLOF & ~wasinRLOF2RLOF & star2RLOF==1 & M2RLOF>15 & star1RLOF==14);
    MBH=[M1RLOF(relevantbinaryFrom2); M2RLOF(relevantbinaryFrom1)]; 
    MO=[M2RLOF(relevantbinaryFrom2); M1RLOF(relevantbinaryFrom1)];
    %this misses binaries in which the companion overflows its RL shortly after 
    %evolving onto the HG; 
    %Now check the binaries at the time of the switch onto the HG, are they 
    %at least 80% RL filling at that point, if so, include them
    M1Switch=h5read(file,'/BSE_Switch_Log/Mass(1)');
    M2Switch=h5read(file,'/BSE_Switch_Log/Mass(2)');
    whichSwitch=h5read(file,'/BSE_Switch_Log/Star_Switching');
    star1Switch=h5read(file,'/BSE_Switch_Log/Stellar_Type(1)');
    star2Switch=h5read(file,'/BSE_Switch_Log/Stellar_Type(2)');
    fromSwitch=h5read(file,'/BSE_Switch_Log/Switching_From');
    toSwitch=h5read(file,'/BSE_Switch_Log/Switching_To');
    roche1Switch=h5read(file,'/BSE_Switch_Log/RocheLobe(1)');
    roche2Switch=h5read(file,'/BSE_Switch_Log/RocheLobe(2)');
    radius1Switch=h5read(file,'/BSE_Switch_Log/Radius(1)');
    radius2Switch=h5read(file,'/BSE_Switch_Log/Radius(2)');
    isMergerSwitch=h5read(file,'/BSE_Switch_Log/Is_Merger');
    relevantBinaryFrom1 = whichSwitch==1 & fromSwitch==1 & toSwitch==2 & star2Switch==14 & radius1Switch>0.8*roche1Switch & M1Switch>15 & ~isMergerSwitch;
    relevantBinaryFrom2 = whichSwitch==2 & fromSwitch==1 & toSwitch==2 & star1Switch==14 & radius2Switch>0.8*roche2Switch & M2Switch>15 & ~isMergerSwitch;
    MBH=[MBH; M2Switch(relevantBinaryFrom1); M1Switch(relevantBinaryFrom2)];
    MO=[MO; M1Switch(relevantBinaryFrom2); M2Switch(relevantBinaryFrom1)];
    figure(fignumberHMXB), scatter(MBH, MO, point, 'filled', colour, 'DisplayName', name); hold on;
end %end of HMXBplot



%Plot BeXRBs
%BeXRBs at time of SN consist of an NS and a MS B (or O) star that previously experienced
%mass accretion to spin it up. No other separation threshold is applied.
function BeXRBplot(file, name, fignumberBeXRB, colour, point)
    starSNSN=h5read(file,'/BSE_Supernovae/Stellar_Type(SN)');
    starCPSN=h5read(file,'/BSE_Supernovae/Stellar_Type(CP)');
    MCPSN=h5read(file,'/BSE_Supernovae/Mass(CP)');
    unboundSN=h5read(file,'/BSE_Supernovae/Unbound');
    timeSN=h5read(file,'/BSE_Supernovae/Time');
    seedSN=h5read(file,'/BSE_Supernovae/SEED');
    seedRLOF=h5read(file,'/BSE_RLOF/SEED');
    timeRLOF=h5read(file,'/BSE_RLOF/Time>MT');
    seedCE=h5read(file,'/BSE_Common_Envelopes/SEED');
    timeCE=h5read(file,'/BSE_Common_Envelopes/Time');
    [experiencedRLOF,indexRLOF]=ismember(seedSN,seedRLOF);
    precededbyRLOF=false(size(indexRLOF)); 
    for(i=1:length(indexRLOF)), if(experiencedRLOF(i)), precededbyRLOF(i)=(timeRLOF(indexRLOF(i))<timeSN(i)); end; end;
    [isCE,indexCE]=ismember(seedSN,seedCE); 
    for(i=1:length(isCE)), if(isCE(i)), isCE(i)=(timeCE(indexCE(i)) < timeSN(i)); end; end;
    relevantbinary=starSNSN==13 & starCPSN==1 & MCPSN>3 & ~unboundSN & experiencedRLOF & precededbyRLOF; %haven't checked direction of RLOF
    figure(fignumberBeXRB), scatter(MCPSN(relevantbinary & isCE), timeSN(relevantbinary & isCE), ...
        point, 'filled', colour, 'DisplayName', ['CE, ', name]); hold on;
    scatter(MCPSN(relevantbinary & ~isCE), timeSN(relevantbinary & ~isCE), ...
        point, colour, 'DisplayName', ['Stable, ', name]);
end %end of BeXRBplot



%Plot LMXBs and IMXBs
%Find all LMXBs/IMXBs at the time of the first mass transfer episode onto the
%compact object.  They consist of an NS or a BH and a companion with a mass
%below 5 Msun.
function [LMXBcount, NSLMXBcount] = LMXBplot(file, name, fignumberLMXB, colour, point)
    seedRLOF=h5read(file,'/BSE_RLOF/SEED');
    timeRLOF=h5read(file,'/BSE_RLOF/Time>MT');
    M1RLOF=h5read(file,'/BSE_RLOF/Mass(1)<MT');
    M2RLOF=h5read(file,'/BSE_RLOF/Mass(2)<MT');
    star1RLOF=h5read(file,'/BSE_RLOF/Stellar_Type(1)<MT');
    star2RLOF=h5read(file,'/BSE_RLOF/Stellar_Type(2)<MT');
    RLOFisCE=h5read(file, '/BSE_RLOF/CEE>MT');
    seedCE=h5read(file,'/BSE_Common_Envelopes/SEED');
    timeCE=h5read(file,'/BSE_Common_Envelopes/Time');
    relevantbinary=(star1RLOF>=13 & star1RLOF<=14 & star2RLOF==1 & M2RLOF<=5 & ~RLOFisCE) ...
        | (star2RLOF>=13 & star2RLOF<=14 & star1RLOF==1 & M1RLOF<=5 & ~RLOFisCE);
    relevantseedRLOF=cast(seedRLOF, 'int64');  relevantseedRLOF(~relevantbinary)=-1;
    uniqueseeds=unique(seedRLOF(relevantbinary));
    [blah,indexlist]=ismember(cast(uniqueseeds,'int64'), relevantseedRLOF);
    M1relevant=M1RLOF(indexlist); M2relevant=M2RLOF(indexlist);
    MCO=M1relevant; MCO(star2RLOF(indexlist)>=13)=M2relevant(star2RLOF(indexlist)>=13);
    Mcomp=M1relevant; Mcomp(star2RLOF(indexlist)==1)=M2relevant(star2RLOF(indexlist)==1);
    LMXBcount=length(MCO); NSLMXBcount=sum(MCO<2); %number of LMXBs, NS LXMBs
    [isCE,CEIndex]=ismember(uniqueseeds,seedCE);  %doesn't check if CE happened before or after
    precededbyCE=false(size(indexlist)); 
    for(i=1:length(indexlist)), if(isCE(i)), precededbyCE(i)=(timeCE(CEIndex(i)) < timeRLOF(indexlist(i))); end; end;
    figure(fignumberLMXB), hold on;
    scatter(MCO(precededbyCE), Mcomp(precededbyCE), point, 'filled', colour, 'DisplayName', ['CE, ', name]);
    scatter(MCO(~precededbyCE), Mcomp(~precededbyCE), point, colour, 'DisplayName', ['Stable, ', name]); 
end %end of LMXBplot


%Plot a*Mtot for white dwarfs
%Should be constant for non-mass-transferring WDs in the absence of tides
%and GR (and supernovae, which is why we focus on WDs as a sanity check
function DWDplot(file, name, fignumberDWD, colourNMT, colourMT, point)
    global AU Rsun
    M1ZAMS=h5read(file, '/BSE_System_Parameters/Mass@ZAMS(1)');
    M2ZAMS=h5read(file, '/BSE_System_Parameters/Mass@ZAMS(2)');
    aZAMS=h5read(file, '/BSE_System_Parameters/SemiMajorAxis@ZAMS');
    SeedZAMS=h5read(file, '/BSE_System_Parameters/SEED');
    Type1=h5read(file, '/BSE_Switch_Log/Stellar_Type(1)');
    Type2=h5read(file, '/BSE_Switch_Log/Stellar_Type(2)');
    M1Switch=h5read(file, '/BSE_Switch_Log/Mass(1)');
    M2Switch=h5read(file, '/BSE_Switch_Log/Mass(2)');
    aSwitch=h5read(file, '/BSE_Switch_Log/SemiMajorAxis');
    SeedSwitch=h5read(file, '/BSE_Switch_Log/SEED');
    isMergerSwitch=h5read(file,'/BSE_Switch_Log/Is_Merger');
    WDIndex=find(Type1>=10 & Type1<=12 & Type2>=10 & Type2<=12 & ~isMergerSwitch);
    ind=unique(SeedSwitch(WDIndex));
    [~,ZAMSIndex]=ismember(ind,SeedZAMS);
    MAZAMS=(M1ZAMS(ZAMSIndex)+M2ZAMS(ZAMSIndex)).*aZAMS(ZAMSIndex);
    MAWD=(M1Switch(WDIndex)+M2Switch(WDIndex)).*aSwitch(WDIndex)*Rsun/AU;
    SeedRLOF=h5read(file, '/BSE_RLOF/SEED');
    [hadRLOF,RLOFIndex]=ismember(ind,SeedRLOF);
    figure(fignumberDWD), hold on;
    scatter(log10(MAZAMS(~hadRLOF)), log10(MAWD(~hadRLOF)), point, 'filled', colourNMT, 'DisplayName', ['DWDs without mass transfer, ', name]);
    scatter(log10(MAZAMS(hadRLOF)), log10(MAWD(hadRLOF)), point, colourMT, 'DisplayName', ['DWDs after mass transfer, ', name]);
end %end of DWDplot



%SN varieties -- just count
function [binariescount, SNcount, BHcompletecount, SNbothcount, SNonecount, ...
            unboundcount, unbound1count, unbound2count, ...
            USSNcount, ECSNcount, AICcount, PISNcount, PPISNcount, ...
            MTSNcount, strippedSNcount, strippednonECSNcount, ...
            strippedwindscount, strippedRLOFcount, strippedCEcount, CESNcount] = SNstats(file)
    starSNSN=h5read(file,'/BSE_Supernovae/Stellar_Type(SN)');
    starCPSN=h5read(file,'/BSE_Supernovae/Stellar_Type(CP)');
    MSNSN=h5read(file,'/BSE_Supernovae/Mass(SN)');
    unboundSN=h5read(file,'/BSE_Supernovae/Unbound');
    MpreSN=h5read(file,'/BSE_Supernovae/Mass_Total@CO(SN)');
    IsStrippedSN=h5read(file,'/BSE_Supernovae/Is_Hydrogen_Poor(SN)');
    HadRLOFSN=h5read(file,'/BSE_Supernovae/Experienced_RLOF(SN)');
    SNtypeSN=h5read(file,'/BSE_Supernovae/SN_Type(SN)');
    timeSN=h5read(file,'/BSE_Supernovae/Time');
    seedSN=h5read(file,'/BSE_Supernovae/SEED');
    prevtypeSN=h5read(file,'/BSE_Supernovae/Stellar_Type_Prev(SN)');
    seedRLOF=h5read(file,'/BSE_RLOF/SEED');
    timeRLOF=h5read(file,'/BSE_RLOF/Time>MT');
    seedCE=h5read(file,'/BSE_Common_Envelopes/SEED');
    timeCE=h5read(file,'/BSE_Common_Envelopes/Time');
    seedAll=h5read(file, '/BSE_System_Parameters/SEED');
    [isCE,indexCE]=ismember(seedSN,seedCE); 
    [isCE,lastindexCE]=ismember(seedSN,seedCE,'legacy'); %last index of CE seed matching given SN seed (for binaries with 2+ CE events)
    simultaneouswithCE = false(size(indexCE)); 
    simultaneouswithCE(isCE) = (timeCE(indexCE(isCE)) == timeSN(isCE)) | (timeCE(lastindexCE(isCE)) == timeSN(isCE)); 
    %Note: could (very rarely) miss a coincidence if >2 CE events and only intermediate CE matches SN in time
    precedingCE = false(size(indexCE)); 
    precedingCE(isCE) = (timeCE(indexCE(isCE)) < timeSN(isCE)); 
    [experiencedRLOF,indexRLOF]=ismember(seedSN,seedRLOF);
    [experiencedRLOF,lastindexRLOF]=ismember(seedSN,seedRLOF, 'legacy');
    simultaneouswithRLOF = false(size(indexRLOF)); 
    simultaneouswithRLOF(experiencedRLOF) = (timeRLOF(indexRLOF(experiencedRLOF)) == timeSN(experiencedRLOF) ...
                                            | timeRLOF(lastindexRLOF(experiencedRLOF)) == timeSN(experiencedRLOF)); 
    %Note: could miss a coincidence if >2 RLOF events and only intermediate RLOF matches SN in time
    binariescount=length(seedAll);
    SNcount=length(seedSN);
    BHcompletecount=sum(starSNSN==14 & MpreSN==MSNSN);
    SNbothcount=length(seedSN)-length(unique(seedSN));
    SNonecount=2*length(unique(seedSN))-length(seedSN);
    firstSN=(starCPSN~=13 & starCPSN~=14);
    secondSN=(starCPSN==13 | starCPSN==14);
    unboundcount=length(unique(seedSN(unboundSN==1)));
    unbound1count=sum(unboundSN(firstSN));
    unbound2count=length(unique(seedSN(unboundSN==1))) - sum(unboundSN(firstSN));
    USSNcount=sum(SNtypeSN==16);
    ECSNcount=sum(SNtypeSN==2);
    AICcount=sum(SNtypeSN==32);
    PISNcount=sum(SNtypeSN==4);
    PPISNcount=sum(SNtypeSN==8);
    MTSNcount=sum(HadRLOFSN);
    strippedSNcount=sum(IsStrippedSN);
    strippednonECSNcount=sum(IsStrippedSN & SNtypeSN~=2 & SNtypeSN~=32);
    strippedwindscount=sum(IsStrippedSN & ~HadRLOFSN & SNtypeSN~=2 & SNtypeSN~=32);
    strippedRLOFcount=sum(IsStrippedSN & HadRLOFSN & prevtypeSN>=7 & prevtypeSN<=9 & SNtypeSN~=2);
    strippedCEcount=sum(IsStrippedSN & HadRLOFSN & precedingCE & SNtypeSN~=2 & SNtypeSN~=32);
    CESNcount=sum(IsStrippedSN & HadRLOFSN & simultaneouswithCE & SNtypeSN~=2 & SNtypeSN~=32);
end %end of SNstats



