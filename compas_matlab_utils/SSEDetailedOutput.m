function SSEDetailedOutput(filename, nfiles)
% Carries out some basic analysis and makes plots for a single SSE run
%
% USAGE: 
% SSEDetailedOutput(filename [, nfiles])
%
% INPUTS:
%   filename: name of detailed output file in COMPAS h5 format
%   nfiles: optional if multiple files are being compared; in that case,
%   the files are all assumed to start with filename, followed by numbers
%   0...nfiles-1, followed by .h5
%
% examples: 
%       SSEDetailedOutput('~/Work/COMPAS/src/COMPAS_Output/Detailed_Output/SSE_Detailed_Output_0.h5')
%       SSEDetailedOutput('~/Work/COMPAS/src/COMPAS_Output/Detailed_Output/SSE_Detailed_Output_', 5)
%       


if(nargin==1) %parse one file

file=filename;
time=h5read(file,'/Time');
Z=h5read(file,'/Metallicity@ZAMS');
mass=h5read(file,'/Mass');
massCO=h5read(file,'/Mass_CO_Core');
massHe=h5read(file,'/Mass_He_Core');
luminosity=h5read(file,'/Luminosity');
type=h5read(file,'/Stellar_Type');
radius=h5read(file,'/Radius');
temperature=h5read(file,'/Teff');
mdot=h5read(file,'/Mdot');
ml=h5read(file,'/Dominant_Mass_Loss_Rate');
%envMass=h5read(file,'/Mass_Env');
%convEnvMass=h5read(file,'/Mass_Convective_Env');
%lambda=h5read(file,'/Lambda_Convective');
%binding=h5read(file,'/BE_ConvectiveEnvelope');
rec=h5read(file,'/Record_Type');


figure(1),     %HR diagram
scatter(temperature, luminosity, 20, type, 'filled');
set(gca,'FontSize',20); xlabel('Teff, T_o'); ylabel('Luminosity, L_o'); set(gca,'XDir','reverse')
H=colorbar; H.Label.String='Stellar type';


figure(2),     %Radial expansion by phase
scatter(time(type==1),radius(type==1),20,'filled', 'b'); hold on;
scatter(time(type==2),radius(type==2),20,'filled', 'g');
scatter(time(type==3),radius(type==3),20,'filled', 'y');
scatter(time(type==4),radius(type==4),20,'filled', 'r');
scatter(time(type==5),radius(type==5),20,'filled', 'm');
scatter(time(type==6),radius(type==6),20,'filled', 'k');
hold off;
set(gca,'FontSize',20); xlabel('Time, Myr'); ylabel('Radius, Rsun'), legend('MS','HG','FGB','CHeB','EAGB','TPAGB')

%for(i=2:length(radius)), %remove points where radius is less than the radius at previous time 
%    if(~isempty(find(radius(1:i-1)>=radius(i)))),
%        envMass(i)=0; convEnvMass(i)=0; radius(i)=0;
%    end;
%end;
%goodindices=type>1 & type<=6 & radius>0;
%figure(2); 
%subplot(2,1,1), scatter(radius(goodindices), envMass(goodindices)-convEnvMass(goodindices), 20, 'b', 'filled'); hold on;
%scatter(radius(goodindices), convEnvMass(goodindices), 20, 'm', 'filled'); hold on;
%scatter(radius(goodindices), mass(goodindices)-envMass(goodindices), 20, 'y', 'filled'); 
%scatter(radius(goodindices), mass(goodindices), 20, 'k', 'filled'); 
%hold off;
%set(gca,'FontSize',20); xlabel('Radius, Rsun'); ylabel('Mass, Msun'), legend('Radiative intershell','Convective Envelope','Core', 'Total')
%subplot(2,1,2); scatter(radius(goodindices), binding(goodindices), 20, 'm', 'filled')
%set(gca,'FontSize',20); xlabel('Radius, Rsun'); ylabel('Binding E conv env, erg')
%%subplot(2,1,2); scatter(radius(goodindices), lambda(goodindices), 20, 'm', 'filled')
%%set(gca,'FontSize',20); xlabel('Radius, Rsun'); ylabel('\lambda')


%radius=h5read(file,'/BSE_Common_Envelopes/Radius(1)<CE');
%apre=h5read(file,'/BSE_Common_Envelopes/SemiMajorAxis<CE');
%a1=h5read(file,'/BSE_Common_Envelopes/SemiMajorAxisStage1>CE');
%apost=h5read(file,'/BSE_Common_Envelopes/SemiMajorAxis>CE');
%primarydonor=(h5read(file,'/BSE_Common_Envelopes/RLOF(1)')==1);
%figure(3); semilogy(radius(primarydonor), apre(primarydonor), 'y-.', radius(primarydonor), a1(primarydonor), 'y--', radius(primarydonor), apost(primarydonor), 'y-', 'LineWidth', 3); hold off;
%set(gca,'FontSize',20); xlabel('Radius, Rsun'); ylabel('Semimajor axis, Rsun'), legend('Pre CE','After Stage1','After Stage 2')

end; % end of if(nargin==1)

if(nargin==2),

for (k=0:nfiles-1),
file=[filename,int2str(k),'.h5'];    
time=h5read(file,'/Time');
Z=h5read(file,'/Metallicity@ZAMS');
mass=h5read(file,'/Mass');
massCO=h5read(file,'/Mass_CO_Core');
massHe=h5read(file,'/Mass_He_Core');
luminosity=h5read(file,'/Luminosity');
type=h5read(file,'/Stellar_Type');
radius=h5read(file,'/Radius');
binding=h5read(file,'/BE_Nanjing');
rec=h5read(file,'/Record_Type');
for(i=2:length(radius)), %remove points where radius is less than one of previously achieved radii
    if(~isempty(find(radius(1:i-1)>=radius(i)))),
        binding(i)=0; radius(i)=0;
    end;
end;
figure(k+1), 
scatter(radius(type==1),log10(binding(type==1)),20,'filled'); hold on;
scatter(radius(type==2),log10(binding(type==2)),20,'filled');
scatter(radius(type==3),log10(binding(type==3)),20,'filled');
scatter(radius(type==4),log10(binding(type==4)),20,'filled');
scatter(radius(type==5),log10(binding(type==5)),20,'filled');
scatter(radius(type==6),log10(binding(type==6)),20,'filled');
hold off;
set(gca,'FontSize',20); xlabel('Radius, Rsun'); ylabel('log_{10}(Nanjing binding energy / erg)'), legend('MS','HG','FGB','CHeB','EAGB','TPAGB'), title(['ZAMS Mass in solar masses:', int2str(mass(1))])
end;

end; %end of if(nargin==2)
