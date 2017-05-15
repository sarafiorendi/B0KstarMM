%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO DO                                         %
% - Define Sensor&Air volumes in SolvePoisson3D %
%   (not available in MATLAB 2016)              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Clean up everything
close all;
clear;
clc;

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Program to caculate the signal in particle detectors %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
BiasV = -300; % Sensor backplane voltage [V]

Fluence = 0.12; % Irradiation fluence [10^16 1MeV n.eq./cm^2]
                % 1/tau = c*Fluence/(1 + c*Fluence/t), extracted from fit to data [ns^-1]
ce = 5.36;
te = 0.8295;
ch = 3.361;
th = 107.6;
TauBe = (1 + ce*Fluence/te)/(ce*Fluence); % Life-time on the backplane side [ns]
TauSe = (1 + ce*Fluence/te)/(ce*Fluence); % Life-time on the strip side [ns]
TauBh = (1 + ch*Fluence/th)/(ch*Fluence); % Life-time on the backplane side [ns]
TauSh = (1 + ch*Fluence/th)/(ch*Fluence); % Life-time on the strip side [ns]

Bulk   = 120; % Bulk thickness [um]
PitchX = 100; % Pitch along X [um] (for 2D&3D geometry)
PitchY = 150; % Pitch along Y [um] (for 3D geometry)

qe    = -1.6e-19; % Electron charge [Coulomb]
eps0  = 8.85e-18; % Vacuum permittivity [F/um]
epsR  = 3.9;      % Relative permittivity [3.9 Silicon, 5.7 Diamond]
dN_dPhi = 30;     % dN/dPhi extracted from data [#/(um^3 10^16)]
DeplV = qe*Bulk^2/(2*epsR*eps0)*dN_dPhi*Fluence - 20; % Sensor full depletion voltage [V]
rho   = 2*DeplV*epsR*eps0/(qe*Bulk^2); % Bulk doping concentration [#/um^3]

BField = 0.0; % Magnetic field (orthogonal+outgoing from 2D geometry) [T]

mu_e   = 140; % Electron mobility [um^2/(V*ns)] [140 Silicon, 180 Diamond]
RH_e   = 1;   % Relative Hall electron mobility [1 Silicon, 1 Diamond]
vs_e   = 110; % Saturation velocity of the electrons [um/ns] [110 Silicon, 260 Diamond]
beta_e = 1;   % Exponent for the electric field dependence of the mobility [0.81 Silicon, 0.81 Diamond]

mu_h   = 45;  % Hole mobility in [um^2/(V*ns)] [45 Silicon, 120 Diamond]
RH_h   = 1;   % Relative Hall hole mobility in [1 Silicon, 1 Diamond] 
vs_h   = 95;  % Saturation velocity of the holes [um/ns] [95 Silicon, 160 Diamond]
beta_h = 1;   % Exponent for the electric field dependence of the mobility [0.42 Silicon, 0.42 Diamond]

Step   = 2;       % Unit step of the lattice on which the field is computed [um]
Radius = Step/10; % Unit step of the movements and field interpolation [um]

XQ = 0; % Coordinate for potential query along z [um]
YQ = 0; % Coordinate for potential query along z [um]

NAverage   = 10;     % Generate NAverage "Work-Transport" matrices and average them
NParticles = 10000;  % Total number of particles to be simulated
PType      = 'beta'; % Particle type ['alpha' 'beta' 'gamma']

fprintf('@@@ Derived parameters @@@\n');
fprintf('\t- Electron''s life-time --> %.2f ns, %.2f [ns]\n',TauBe,TauSe);
fprintf('\t- Hole''s life-time --> %.2f ns, %.2f [ns]\n',TauBh,TauSh);
fprintf('\t- Full depletion voltage --> %.1f [V]\n',DeplV);
fprintf('\t- Doping concentration --> %.1E [#/cm^3]\n',rho*1e12);
fprintf('\t- Resistivity --> %.1E [Ohm cm]\n\n',-1/(qe*mu_h*rho)*1e-13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rng default; % Reset random seed
ItFig = 1;   % Figure iterator


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the potentials %
%%%%%%%%%%%%%%%%%%%%%%%%%%
[TotalPot,  ~,       ~, ItFig] = SolvePoissonPDE2D(Bulk,PitchX,BiasV,0,0,epsR,rho*qe/eps0,XQ,ItFig);
[WeightPot, Sq2D, xq2D, ItFig] = SolvePoissonPDE2D(Bulk,PitchX,0,0,1,epsR,0,XQ,ItFig);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare the potential in 2D and 3D %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m comparing the potential for 2D and 3D @@@\n');
[~, Sq3D, xq3D, ItFig] = SolvePoissonPDE3D(Bulk,PitchX,PitchY,0,1,epsR,0,XQ,YQ,ItFig);
Diff2D3D = ((Sq2D ./ xq2D' + Sq2D ./ (Bulk - xq2D')) ./...
    (Sq3D ./ xq3D' + Sq3D ./ (Bulk - xq3D')) - 1) * 100;
fprintf('@@@ Weighted potential difference %.1f%% @@@\n\n',mean(Diff2D3D,'omitnan'));

figure(ItFig);
plot(xq2D,Diff2D3D);
title(sprintf('Weighted potential difference at x = %.2f um y = %.2f um',XQ,YQ));
xlabel('Z [\mum]');
ylabel('Percentage [%]');
grid on;
ItFig = ItFig + 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the velocity field %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[VFieldx_e, VFieldy_e, VFieldx_h, VFieldy_h, x, y, ItFig] =...
    VelocityField(TotalPot,Step,Bulk,BField,PitchX,...
    mu_e,RH_e,vs_e,beta_e,mu_h,RH_h,vs_h,beta_h,ItFig);


%%%%%%%%%%%%%%%%%%%%
% Compute the work %
%%%%%%%%%%%%%%%%%%%%
[WorkTransportTotal, x, y, ItFig] =...
    ManyWorkTransport(WeightPot,VFieldx_e,VFieldy_e,VFieldx_h,VFieldy_h,...
    x,y,Step,Bulk,Radius,TauBe,TauSe,TauBh,TauSh,NAverage,ItFig);


%%%%%%%%%%%%%%%%%
% Fit/Smoothing %
%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m interpolating the work transport matrix with a step of %.2f um @@@\n\n',Radius);
subx = x(1):Radius:x(length(x));
suby = y(1):Radius:y(length(y));
[ssubx,ssuby] = meshgrid(subx,suby);
[x0, y0, z0] = prepareSurfaceData(x,y,WorkTransportTotal);
[WorkTransportTotal_, goodness, output] = fit([x0 y0],z0,'cubicinterp');
subWorkTransportTotal = WorkTransportTotal_(ssubx, ssuby);

figure(ItFig);
colormap jet;
surf(ssubx,ssuby,subWorkTransportTotal,'EdgeColor','none');
title('Interpolated Total <Work-Transport>');
xlabel('X [\mum]');
ylabel('Z [\mum]');
zlabel('Work / q [#charges * V]');
ItFig = ItFig + 1;


%%%%%%%%%%%%%%%%%%%%%%%
% Compute the spectra %
%%%%%%%%%%%%%%%%%%%%%%%
[ItFig] = ComputeSpectra(subWorkTransportTotal,subx,suby,NParticles,...
    PitchX,Bulk,Radius,PType,ItFig);


%%%%%%%%%%%%%%%%%%%%%%%%
% Alert for completion %
%%%%%%%%%%%%%%%%%%%%%%%%
load handel;
sound(y,Fs);
