%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that generates a Landau distribution of charge %
% and computes the  collected charge                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WorkTransportTotal = Total "Work-Transport" matrix
% x, y       = Axes
% Step       = Unit step of the lattice on which the field is computed [um]
% subStep    = Unit step of the lattice on which the field is interpolated [um]
% Bulk       = Bulk thickness [um]
% Pitch      = Strip pitch [um]
% Radius     = Unit step of the movements [um]
% NParticles = Total number of particles to be simulated

function [] = ComputeSpectra(WorkTransportTotal,x,y,...
    NParticles,Pitch,Step,subStep,Bulk,Radius,particle)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(particle,'beta') == true
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Charge spectrum for Beta %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ehLength = 60; % Electron-holes pairs per unit length [electrons]
    depth    = Bulk; % Source penetration depth [um]
    eMax     = ehLength*Bulk; % Maximum released charge [electrons]
elseif strcmp(particle,'alpha') == true
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Charge spectrum for Alpha %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    energyPair = 13; % Energy to create an electron-hole pair [eV]
    depth      = 10; % Source penetration depth [um]
    eMax       = 4520000 / energyPair; % Maximum released charge [electrons]
    mean       = eMax; % Alpha spectrum mean [electrons]
    sigma      = 17000 / energyPair; % Alpha spectrum sigma [electrons]
elseif strcmp(particle,'gamma') == true
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Charge spectrum for Gamma %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Particle not yet implemented\n');
    return;
else
    fprintf('Unknown particle: %s\n',particle);
    return;
end

eNoise      = 1000; % Electronic noise [electrons]
nBins       = 100;  % Spectrum's number of bins
EnergyScale = 0:eMax/nBins:eMax; % Spectrum energy axis [electrons]


%%%%%%%%%%%%%%%%%%%
% Start algorithm %
%%%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m calculating the spectra of %d %s particles @@@\n',NParticles,particle);
HistoCharge = zeros(NParticles,1);

for i = 1:NParticles

    if strcmp(particle,'beta') == true
    % Particles enter randomly between -Pitch/2 -- +Pitch/2
    % Particles exit  randomly between -Pitch/2 -- +Pitch/2
%        enter = Pitch/2 * (2*rand(1,1) - 1); % x-coordinate entering particle
%        exit  = Pitch/2 * (2*rand(1,1) - 1); % x-coordinate exiting  particle
        enter = 0;
        exit  = 0;
        mean  = ehLength*sqrt(depth^2 + (exit-enter)^2); % Landau MPV [electrons]
        sigma = mean / 10; % 8 = scale factor between MPV and sigma of Landau [electrons]
        ChargeDensity = LandauRND(mean,sigma);
    elseif strcmp(particle,'alpha') == true
    % Particles enter randomly between -Pitch/2 -- +Pitch/2
    % Particles exit  randomly at Xin + Depth*sin(theta)
    % where theta is chosen randomly between -pi/4 and pi/4
        enter = Pitch/2 * (2*rand(1,1) - 1);
        exit  = enter + depth*tan(pi/2*rand(1,1) - pi/4);
        ChargeDensity = normrnd(mean,sigma);
    elseif strcmp(particle,'gamma') == true
        fprintf('Particle not yet implemented\n');
        return;
    end
    
    Noise         = normrnd(0,eNoise);
    ChargeDensity = (ChargeDensity + Noise) / sqrt(depth^2 + (exit-enter)^2);
    [Charge] = ComputeSignal(WorkTransportTotal,x,y,depth,...
        enter,exit,Step,subStep,Bulk,Radius,ChargeDensity);
    
    HistoCharge(i) = Charge;
end


%%%%%%%%%
% Plots %
%%%%%%%%%
figure (8);
hold off
histo = hist(HistoCharge,EnergyScale);
plot(EnergyScale,histo,'-o','LineWidth',1,...
    'MarkerEdgeColor','black','MarkerFaceColor','blue','MarkerSize',7);
title('Collected charge histogram');
xlabel('Energy [e-]');
ylabel('Entries [#]');
grid on

fprintf('CPU time --> %d[min]\n\n',(cputime-TStart)/60);
end