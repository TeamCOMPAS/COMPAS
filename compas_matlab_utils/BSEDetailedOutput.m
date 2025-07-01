function BSEDetailedOutput(file)
% Carries out some basic analysis and makes plots for a single BSE run
%
% USAGE: 
% BSEDetailedOutput(file)
%
% INPUTS:
%   file: name of detailed output file in COMPAS h5 format
%
% examples: 
%       BSEDetailedOutput('~/Work/COMPAS/src/COMPAS_Output/Detailed_Output/BSE_Detailed_Output_0.h5')
%       

time=h5read(file,'/Time');
MThistory=h5read(file,'/MT_History');
Z=h5read(file,'/Metallicity@ZAMS(1)');
mass1=h5read(file,'/Mass(1)');
massCO1=h5read(file,'/Mass_CO_Core(1)');
mass2=h5read(file,'/Mass(2)');
massCO2=h5read(file,'/Mass_CO_Core(2)');
massHe1=h5read(file,'/Mass_He_Core(1)');
massHe2=h5read(file,'/Mass_He_Core(2)');
luminosity1=h5read(file,'/Luminosity(1)');
luminosity2=h5read(file,'/Luminosity(2)');
type1=h5read(file,'/Stellar_Type(1)');
type2=h5read(file,'/Stellar_Type(2)');
radius1=h5read(file,'/Radius(1)');
radius2=h5read(file,'/Radius(2)');
RL1=h5read(file,'/RocheLobe(1)');
RL2=h5read(file,'/RocheLobe(2)');
a=h5read(file,'/SemiMajorAxis');
e=h5read(file,'/Eccentricity');
rec=h5read(file,'/Record_Type');
AM1=h5read(file, '/Ang_Momentum(1)');
AM2=h5read(file, '/Ang_Momentum(2)');
AMtot=h5read(file, '/Ang_Momentum_Total');
omega1=h5read(file, '/Omega(1)');
omega2=h5read(file, '/Omega(2)');
omegab1=h5read(file, '/Omega_Break(1)');
omegab2=h5read(file, '/Omega_Break(2)');
MTtimescale=h5read(file, '/MassTransferTimescale');


disp('Time (Myr),   Event,  M1 (M_o),   type1,  M2 (M_o),   type2,  a (R_o),    e');
disp([num2str(time(1)), '     start: Z=', num2str(Z(1)), '      ', num2str(mass1(1)), '      ', num2str(type1(1)),...
    '      ', num2str(mass2(1)), '      ', num2str(type2(1)), '      ', num2str(a(1)), '      ', num2str(e(1))]);
for(i=2:length(time)),
    if(rec(i)==4), %only print records at end of timestep
        prev=find(rec(1:i-1)==4, 1, 'last');   %previous timestep index           
        if(MThistory(i)~=0),% & MThistory(i)~=MThistory(i-1)), 
            switch (MThistory(i)), case 1, MTstring='Stable MT from 1 TO 2'; case 2, MTstring='Stable MT from 2 TO 1';...
                case 3, MTstring='CE from 1 to 2'; case 4, MTstring='CE from 2 to 1'; case 5, MTstring='CE Double Core'; ...
                case 6, MTstring='CE merger';   end;
            switch(MTtimescale(i)), case 1, MTtimescalestring='Nuclear'; case 2, MTtimescalestring='Thermal'; case 3, ...
                    MTtimescalestring='Dynamical'; end;
            disp([num2str(time(i)), '  ', MTstring, ', @', MTtimescalestring, '    ', num2str(mass1(i)), '      ', num2str(type1(i)),...
                '      ', num2str(mass2(i)), '      ', num2str(type2(i)), '      ', num2str(a(i)), '      ', num2str(e(i))]);
        end;
        if(type1(i)~=type1(prev))
            disp([num2str(time(i)), '     star 1 ', num2str(type1(prev)), '->', num2str(type1(i)), '   ', num2str(mass1(i)), '      ', num2str(type1(i)),...
                '      ', num2str(mass2(i)), '      ', num2str(type2(i)), '      ', num2str(a(i)), '      ', num2str(e(i))]);
        end;
        if(type2(i)~=type2(prev)),
            disp([num2str(time(i)), '     star 2 ', num2str(type2(prev)), '->', num2str(type2(i)), '   ', num2str(mass1(i)), '      ', num2str(type1(i)),...
                '      ', num2str(mass2(i)), '      ', num2str(type2(i)), '      ', num2str(a(i)), '      ', num2str(e(i))]);
        end;
    end;
end;
if((type1(i)==13 || type1(i)==14) && (type2(i)==13 || type2(i)==14) ),
    Msunkg=1.98892e30;	c=299792458;		G=6.67428e-11;		Rsun = 695500000; 
    beta=64/5*G^3*mass1(i)*mass2(i)*(mass1(i)+mass2(i))*Msunkg^3/c^5; T0=(a(i)*Rsun)^4/4/beta;
    Tdelay=T0*(1-e(i)^2).^(7/2).*(1+0.27*e(i)^10 + 0.33*e(i)^20 +  0.2*e(i)^1000)/3.15e7/1e6;
    disp([num2str(time(i)+Tdelay), '     GW merger in ', num2str(Tdelay,'%.0f'), ' Myr      ', num2str(mass1(i)), '      ', num2str(type1(i)),...
        '      ', num2str(mass2(i)), '      ', num2str(type2(i))]);
end;

figure(1), plot(time,  mass1, 'LineWidth', 3),  set(gca,'FontSize',20), xlabel('Time, Myr'), ylabel('Total mass 1, M_o')
%figure(2), scatter(time,  radius1, 30, type1, 'filled'),  set(gca,'FontSize',20), xlabel('Time, Myr'), ylabel('Radius 1, R_o')
figure(2), plot(time, radius2, 'b', time, RL2, 'r'); set(gca,'FontSize',20), xlabel('Time, Myr'), ylabel('Radius, R_o'), legend('R_2', 'RL_2')
