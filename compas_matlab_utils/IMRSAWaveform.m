function [h, Aeff, psiEff] = IMRSAWaveform (fVec, totalMass, eta, chi, startPhase, ...
    startTime, dEff, fLower) 
%
% IMRSAWaveform: Generate the parametrised frequency domain BBH coalescence 
% waveforms proposed in the paper P. Ajith et al (2009)
%
% usage: [h, Aeff, psiEff] = IMRSAWaveform (f, totalMass, eta, chi, startPhase, 
%   startTime, dEff, fLower)
% 
%   fVec        : a vector of frequencies at which the waveform is generated
%   totalMass   : total mass of the binary (in solar masses)
%   eta         : symmetric mass ratio 
%   chi         : spin parameter 
%   startPhase  : start phase of the waveform (in radian)
%   startTime   : start time of the waveform (in seconds)
%   dEff        : effective distance (in Mpc)
%   fLower      : low-frequency cutoff (Hz)
%
%   h           : waveform in Fourier domain (complex vector)
%   Aeff        : amplitude in the Fourier domain 
%   psiEff      : phase in the Fourier domain 
%
% P. Ajith, 30.06.2010
%
% $Id: IMRSAWaveform.m 73 2010-07-01 21:39:26Z ajith $

    MSOLAR_TIME = 4.92579497077314e-06;
    PARSEC_SEC = 1.0292712503e8;

    % parameters
    M = totalMass*MSOLAR_TIME;
    piM = pi*M;
    d = dEff*PARSEC_SEC*1e6;
    phi  = startPhase;
    shft = 2*pi*startTime;
    delta = sqrt(1.-4.*eta);
    mergPower = -2/3;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute the phenomenological parameters 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    psi0 = 3./(128*eta);

    psi1 = 0.;

    psi2 = 3715./756. + ...
        -9.2091e+02*eta^1 + 4.9213e+02*eta^1*chi^1 + 1.3503e+02*eta^1*chi^2 + ...
        6.7419e+03*eta^2 + -1.0534e+03*eta^2*chi^1 + ...
        -1.3397e+04*eta^3 ;    
    
    psi3 = -16.*pi + 113.*chi/3. + ...
        1.7022e+04*eta^1 + -9.5659e+03*eta^1*chi^1 + -2.1821e+03*eta^1*chi^2 + ...
        -1.2137e+05*eta^2 + 2.0752e+04*eta^2*chi^1 + ...
        2.3859e+05*eta^3 ;    
    
    psi4 = 15293365./508032. - 405.*chi^2/8. + ...
        -1.2544e+05*eta^1 + 7.5066e+04*eta^1*chi^1 + 1.3382e+04*eta^1*chi^2 + ...
        8.7354e+05*eta^2 + -1.6573e+05*eta^2*chi^1 + ...
        -1.6936e+06*eta^3 ;    
    
    psi5 = 0.;

    psi6 = -8.8977e+05*eta^1 + 6.3102e+05*eta^1*chi^1 + 5.0676e+04*eta^1*chi^2 + ...
        5.9808e+06*eta^2 + -1.4148e+06*eta^2*chi^1 + ...
        -1.1280e+07*eta^3 ;    
        
    psi7 = 8.6960e+05*eta^1 + -6.7098e+05*eta^1*chi^1 + -3.0082e+04*eta^1*chi^2 + ...
        -5.8379e+06*eta^2 + 1.5145e+06*eta^2*chi^1 + ...
        1.0891e+07*eta^3 ;    
    
    psi8 = -3.6600e+05*eta^1 + 3.0670e+05*eta^1*chi^1 + 6.3176e+02*eta^1*chi^2 + ...
        2.4265e+06*eta^2 + -7.2180e+05*eta^2*chi^1 + ...
        -4.5524e+06*eta^3 ; 
        
    f1 =  1. - 4.4547*(1.-chi).^0.217 + 3.521*(1.-chi).^0.26 + ...
        6.4365e-01*eta^1 + 8.2696e-01*eta^1*chi^1 + -2.7063e-01*eta^1*chi^2 + ...
        -5.8218e-02*eta^2 + -3.9346e+00*eta^2*chi^1 + ...
        -7.0916e+00*eta^3 ;    
    
    f2 = (1. - 0.63*(1.-chi)^0.3)/2. + ...
        1.4690e-01*eta^1 + -1.2281e-01*eta^1*chi^1 + -2.6091e-02*eta^1*chi^2 + ...
        -2.4900e-02*eta^2 + 1.7013e-01*eta^2*chi^1 + ...
        2.3252e+00*eta^3 ;    
    
    sigma = (1. - 0.63*(1.-chi)^0.3)*(1.-chi)^0.45/4. + ...
        -4.0979e-01*eta^1 + -3.5226e-02*eta^1*chi^1 + 1.0082e-01*eta^1*chi^2 + ...
        1.8286e+00*eta^2 + -2.0169e-02*eta^2*chi^1 + ...
        -2.8698e+00*eta^3 ;    
    
    f3 = 3.2361e-01 + 4.8935e-02*chi^1 + 1.3463e-02*chi^2 + ...
        -1.3313e-01*eta^1 + -8.1719e-02*eta^1*chi^1 + 1.4512e-01*eta^1*chi^2 + ...
        -2.7140e-01*eta^2 + 1.2788e-01*eta^2*chi^1 + ...
        4.9220e+00*eta^3 ;    
    
    f3 = f3/piM;
	f1 = f1/piM;
	f2 = f2/piM;
    sigma = sigma/piM;

    %fprintf('\nPhenomenological parameters:\n');
    %fprintf('psi0 = %5.4e \t psi1 = %5.4e \t psi2 = %5.4e \t psi3 = %5.4e \t psi4 = %5.4e\n',...
    %    psi0, psi1, psi2, psi3, psi4);
    %fprintf('psi5 = %5.4e  \t psi6 = %5.4e \t psi7 = %5.4e \t psi8 = %5.4e\n',...
    %    psi5, psi6, psi7, psi8);
    %fprintf('f1 = %5.4e \t f2 = %5.4e \t f3 = %5.4e  \t sigma = %5.4e\n',...
    %    f1, f2, f3, sigma);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now generate the waveforms 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % PN corrections to the inspiral amplitude
    alpha2 = -323./224. + 451.*eta/168.;
    alpha3 = (27./8. - 11.*eta/6.)*chi;

    % correction to the merger power law
    epsilon_1 =  1.4547*chi - 1.8897; 
    epsilon_2 = -1.8153*chi + 1.6557;

    % normalisation constants for the merger and ring down amplitude 
    vMerg = power(pi*M*f1, 1./3.);
    vRing = power(pi*M*f2, 1./3.);
    w1 = 1. + alpha2*power(vMerg,2.) + alpha3*power(vMerg,3.);
    w1 = w1/(1. + epsilon_1*vMerg + epsilon_2*vMerg*vMerg);
    w2 = w1*(pi*sigma/2.)*power(f2/f1, mergPower)*(1. ...
       + epsilon_1*vRing + epsilon_2*vRing*vRing);

    % amplitude scale 
    C = M^(5/6)*(f1)^(-7/6)*sqrt(5*eta/24)/(d*pi^(2/3));

    Aeff = zeros(size(fVec));
    psiEff = zeros(size(fVec));

    inspIdx = find(fVec >= fLower & fVec < f1);
    mergIdx = find(fVec >= f1 & fVec < f2);
    ringIdx = find(fVec >= f2 & fVec <= f3);
    outOfBand = find(fVec < fLower | fVec > f3);
    
    v = power(pi*M*fVec, 1./3.);

    % effective amplitude - inspiral, merger and ring down 
    fNorm = fVec/f1;
    pnCorr = 1 + alpha2*power(v(inspIdx),2.) + alpha3*power(v(inspIdx),3.);

    Aeff(inspIdx) = power(fNorm(inspIdx), -7./6.).*pnCorr;
    Aeff(mergIdx) = w1*power(fNorm(mergIdx), mergPower).*(1. ...
                        + epsilon_1*v(mergIdx) + epsilon_2*power(v(mergIdx),2));
    
    Lorentzian = sigma./(power(fVec(ringIdx)-f2, 2.) + sigma.*sigma/4.);
    Aeff(ringIdx) = w2*Lorentzian/(2*pi);

    Aeff = C*Aeff;
    psiEff = shft*fVec + phi + psi0*v.^-5.*(1 + psi2*v.^2 + psi3*v.^3 + ...
                psi4*v.^4 + psi5*v.^5 + psi6*v.^6 + psi7*v.^7 + psi8*v.^8);

    h = Aeff.*exp(1i*(-psiEff)); 

