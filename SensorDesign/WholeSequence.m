%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to caculate the signal in particle-detection sensors %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step   = Unit step of the lattice on which the field is computed [um]
% Radius = Unit step of the movements and field interpolation [um]

% BiasV = Sensor backplane voltage [V]
% Bulk  = Bulk thickness [um]
% Pitch = Strip pitch [um]

% BField = Magnetic field (orthogonal+outgoing from the 2D geometry) [T]

% TauBe/h = Life-time on the backplane side [ns]
% TauSe/h = Life-time on the strip side [ns]

% NAverage = Generate NAverage "Work-Transport" matrices and average them
% NParticles = Total number of particles to be simulated
% ParticleType = 'alpha', 'beta', 'gamma'


Step   = 5;
Radius = Step/10;

BiasV = -200;
Bulk  = 100;
Pitch = 100;

BField = 0.0;

TauBe = 89;
TauSe = 89;
TauBh = 65;
TauSh = 65;

NAverage     = 10;
NParticles   = 10000;
ParticleType = 'beta';


rng default; % Reset random seed
ItFig = 1;   % Figure iterator

% Run PDE_WeightingField(Pitch,Bulk) to export geometry
[WeightPot, ItFig] = SolvePoissonPDE(gd,sf,ns,Bulk,Pitch,0,0,1,ItFig);
[TotalPot,  ItFig] = SolvePoissonPDE(gd,sf,ns,Bulk,Pitch,BiasV,0,0,ItFig);

[VFieldx_e, VFieldy_e, VFieldx_h, VFieldy_h, x, y, ItFig] =...
    VelocityField(TotalPot,Step,Bulk,BField,Pitch,ItFig);

[WorkTransportTotal, x, y, ItFig] =...
    ManyWorkTransport(WeightPot,VFieldx_e,VFieldy_e,VFieldx_h,VFieldy_h,...
    x,y,Step,Bulk,Radius,TauBe,TauSe,TauBh,TauSh,NAverage,ItFig);


%%%%%%%%%%%%%%%%%
% Fit/Smoothing %
%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m interpolating the work transport matrix with a step of %d um @@@\n\n',Radius);
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
ylabel('Y [\mum]');
zlabel('Work / q [#charges * V]');
ItFig = ItFig + 1;


[ItFig] = ComputeSpectra(subWorkTransportTotal,subx,suby,NParticles,Pitch,...
    Bulk,Radius,ParticleType,ItFig);


%%%%%%%%%
% To do %
%%%%%%%%%
% 1. Definition of multiple dielectrics
% 2. Automatic generation of the geometry