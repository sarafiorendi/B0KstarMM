%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that generates a Landau distribution of charge %
% and computes the  collected charge                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WorkTransportTotal = Total "Work-Transport" matrix
% x, y       = Axes
% Bulk       = Bulk thickness [um]
% Pitch      = Strip pitch [um]
% Radius     = Unit step of the movements and field interpolation [um]
% NParticles = Total number of particles to be simulated
% ItFigIn    = Figure iterator input

function [ItFigOut] = ComputeSpectra(WorkTransportTotal,x,y,...
    NParticles,Pitch,Bulk,Radius,particle,ItFigIn)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
eNoise  = 1000;  % Electronic noise [electrons]
nBins   = 100;   % Spectrum's number of bins
IsGamma = false; % Define whether the particle is a gamma or not

if strcmp(particle,'beta') == true
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Charge spectrum for Beta %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ehLength = 60;            % Electron-holes pairs per unit length [electrons]
    depth    = Bulk;          % Source penetration depth [um]
    eMax     = ehLength*Bulk; % Maximum released charge [electrons]
elseif strcmp(particle,'alpha') == true || strcmp(particle,'gamma') == true
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Charge spectrum for Alpha/Gamma %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    energyPair = 3.2;                  % Energy to create an electron-hole
                                       % pair [eV] [3.2 Silicon, 13 Diamond]
    eMax       = 4520000 / energyPair; % Maximum released charge [electrons]
    mean       = eMax;                 % Particle spectrum mean [electrons]
    sigma      = 17000 / energyPair;   % Particle spectrum sigma [electrons]
    % Particle penetration depth [um]
    if strcmp(particle,'alpha') == true
        depth = 10;
    elseif strcmp(particle,'gamma') == true
        depth   = Bulk * rand(1,1);
        IsGamma = true;
    end
else
    fprintf('Unknown particle: %s\n',particle);
    return;
end

EnergyScale = 0:eMax/nBins:2*eMax; % Spectrum energy axis [electrons]


%%%%%%%%%%%%%%%%%%%
% Start algorithm %
%%%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m calculating the spectra of %d %s particles @@@\n',NParticles,particle);
HistoCharge = zeros(NParticles,1);

for i = 1:NParticles

    if strcmp(particle,'beta') == true
    % Particles enter randomly between -Pitch/2 -- +Pitch/2
    % Particles exit  randomly between -Pitch/2 -- +Pitch/2
% @TMP@
%        enter = Pitch/2 * (2*rand(1,1) - 1); % x-coordinate entering particle
%        exit  = Pitch/2 * (2*rand(1,1) - 1); % x-coordinate exiting  particle
        enter = 0;
        exit  = enter;
        mean  = ehLength*sqrt(depth^2 + (exit-enter)^2); % Landau MPV [electrons]
        sigma = mean / 10; % 8 or 10: scale factor between MPV and sigma of Landau [electrons]
        ChargeDensity = LandauRND(mean,sigma);
    elseif strcmp(particle,'alpha') == true || strcmp(particle,'gamma') == true
    % Particles enter randomly between -Pitch/2 -- +Pitch/2
    % Particles exit  randomly at Xin + Depth*tan(theta)
    % where theta is chosen randomly between -pi/4 and pi/4
        enter = Pitch/2 * (2*rand(1,1) - 1);
        exit  = enter + depth*tan(pi/2*rand(1,1) - pi/4);
        ChargeDensity = normrnd(mean,sigma);
    end
    
    Noise         = normrnd(0,eNoise);
    ChargeDensity = ChargeDensity + Noise;
    if strcmp(particle,'gamma') == false
        ChargeDensity = ChargeDensity / sqrt(depth^2 + (exit-enter)^2);
    else
        ChargeDensity = ChargeDensity / Radius;
    end
    [Charge] = ComputeSignal(WorkTransportTotal,x,y,depth,...
        enter,exit,Bulk,Radius,ChargeDensity,IsGamma);
    
    HistoCharge(i) = Charge;
end


%%%%%%%%%
% Plots %
%%%%%%%%%
figure(ItFigIn);
histo = hist(HistoCharge,EnergyScale);

% @TMP@
%ft = fittype('Landau(x,a,b,c)');
%f  = fit(EnergyScale,histo,ft);

plot(EnergyScale,histo,'-o','LineWidth',1,...
    'MarkerEdgeColor','black','MarkerFaceColor','blue','MarkerSize',7);
title('Collected charge histogram');
xlabel('Energy [e-]');
ylabel('Entries [#]');
grid on


%%%%%%%%%%%%%%%%%%%%%%%
% Save data into file %
%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen('Spectrum.txt','w');
fprintf(fileID,'%14s %14s\n','Energy [e-]','Entries [#]');
fprintf(fileID,'%14.2f %14.2f\n',[EnergyScale; histo]);
fclose(fileID);

ItFigOut = ItFigIn + 1;
fprintf('CPU time --> %d[min]\n\n',(cputime-TStart)/60);
end