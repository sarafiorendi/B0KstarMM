         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Program to caculate the signal in particle detectors %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Variable initialization %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsR = 3.9;                       % Relative dielectric constant [3.9 Silicon, 5.7 Diamond]
rho  = (-4 * 1.6e-19) / 8.85e-18; % Charge denisty in the bulk [(Coulomb / um^3) / eps0 [F/um]]

Step   = 3;       % Unit step of the lattice on which the field is computed [um]
Radius = Step/10; % Unit step of the movements and field interpolation [um]

BiasV = -200; % Sensor backplane voltage [V]
Bulk  = 100;  % Bulk thickness [um]
Pitch = 100;  % Strip pitch [um]

BField = 0.0; % Magnetic field (orthogonal+outgoing from the 2D geometry) [T]

TauBe = 3.79; % Life-time on the backplane side [ns]
TauSe = 3.79; % Life-time on the strip side [ns]
TauBh = 4.46; % Life-time on the backplane side [ns]
TauSh = 4.46; % Life-time on the strip side [ns]

NAverage   = 100;    % Generate NAverage "Work-Transport" matrices and average them
NParticles = 20000;  % Total number of particles to be simulated
PType      = 'beta'; % Particle type ['alpha' 'beta' 'gamma']
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rng default; % Reset random seed
ItFig = 1;   % Figure iterator


[WeightPot, ItFig] = SolvePoissonPDE(Bulk,Pitch,0,0,1,epsR,0,ItFig);
[TotalPot,  ItFig] = SolvePoissonPDE(Bulk,Pitch,BiasV,0,0,epsR,rho,ItFig);

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


[ItFig] = ComputeSpectra(subWorkTransportTotal,subx,suby,NParticles,...
    Pitch,Bulk,Radius,PType,ItFig);
