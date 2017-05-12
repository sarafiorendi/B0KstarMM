%%%%%%%%%
% TO DO %
%%%%%%%%%
% Define Sensor&Air volumes in SolvePoisson3D - not available in MATLAB 2016


% Clean up everything
close all;
clear;
clc;

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Program to caculate the signal in particle detectors %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Variable initialization %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsR = 3.9;                    % Relative dielectric constant [3.9 Silicon, 5.7 Diamond]
rho = -4 * 1.6e-19 / 8.85e-18; % Charge denisty in the bulk [(Coulomb / um^3) / eps0 [F/um]]

Step   = 1;       % Unit step of the lattice on which the field is computed [um]
Radius = Step/10; % Unit step of the movements and field interpolation [um]

BiasV  = -300; % Sensor backplane voltage [V]
Bulk   = 120;  % Bulk thickness [um]
PitchX = 100;  % Pitch along X [um] (for 2D geometry)
PitchY = 150;  % Pitch along Y [um] (for 3D geometry)

XQ = 0; % Coordinate for potential query along z [um]
YQ = 0; % Coordinate for potential query along z [um]

BField = 0.0; % Magnetic field (orthogonal+outgoing from the 2D geometry) [T]

Fluence = 0.12; % Irradiation fluence [10^16 1 MeV eq. n. / cm^2]
                % 1/tau = c*(1-exp(-(Fluence-s)/t)), extracted from fit to data [ns^-1]
ce = 0.6735;
se = -0.01235;
te = 0.1772;
ch = 2.586;
sh = 0.008732;
th = 0.6264;
TauBe = 1/(ce*(1 - exp(-(Fluence - se)/te))); % Life-time on the backplane side [ns]
TauSe = 1/(ce*(1 - exp(-(Fluence - se)/te))); % Life-time on the strip side [ns]
TauBh = 1/(ch*(1 - exp(-(Fluence - sh)/th))); % Life-time on the backplane side [ns]
TauSh = 1/(ch*(1 - exp(-(Fluence - sh)/th))); % Life-time on the strip side [ns]
fprintf('@@@ Electron''s life-time: %.2f ns, %.2f ns @@@\n',TauBe,TauSe);
fprintf('@@@ Hole''s life-time: %.2f ns, %.2f ns @@@\n\n',TauBh,TauSh);

NAverage   = 10;     % Generate NAverage "Work-Transport" matrices and average them
NParticles = 10000;  % Total number of particles to be simulated
PType      = 'beta'; % Particle type ['alpha' 'beta' 'gamma']
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rng default; % Reset random seed
ItFig = 1;   % Figure iterator


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the potentials %
%%%%%%%%%%%%%%%%%%%%%%%%%%
[TotalPot,  ~,       ~, ItFig] = SolvePoissonPDE2D(Bulk,PitchX,BiasV,0,0,epsR,rho,XQ,ItFig);
[WeightPot, Sq2D, xq2D, ItFig] = SolvePoissonPDE2D(Bulk,PitchX,0,0,1,epsR,0,XQ,ItFig);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare the potential in 2D and 3D %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m comparing the potential in 2D and 3D @@@\n');
[~, Sq3D, xq3D, ItFig] = SolvePoissonPDE3D(Bulk,PitchX,PitchY,0,1,epsR,0,XQ,YQ,ItFig);
Diff2D3D = (Sq2D - Sq3D) ./ Sq3D * 100;

figure(ItFig);
plot(xq2D,Diff2D3D);
title(sprintf('Potential difference (2D - 3D) / 3D at x = %.2f um y = %.2f um',XQ,YQ));
xlabel('Z [\mum]');
ylabel('Percentage [%]');
grid on;
ItFig = ItFig + 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the velocity field %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[VFieldx_e, VFieldy_e, VFieldx_h, VFieldy_h, x, y, ItFig] =...
    VelocityField(TotalPot,Step,Bulk,BField,PitchX,ItFig);


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
