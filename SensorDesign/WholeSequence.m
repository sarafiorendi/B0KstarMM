%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that run the whole sequence of steps %
% to caculate the signal in particle sensors    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step     = Unit step of the lattice on which the field is computed [um]
% Bulk     = Bulk thickness [um]
% Pitch    = Strip pitch [um]
% Radius   = Unit step of the movements and field interpolation [um]
% BField   = Magnetic field (orthogonal+outgoing from the 2D geometry) [T]
% TauBe/h  = Life-time on the backplane side [ns]
% TauSe/h  = Life-time on the strip side [ns]
% NAverage = Generate NAverage "Work-Transport" matrices and average them
% NParticles   = Total number of particles to be simulated
% ParticleType = 'alpha', 'beta', 'gamma'

Step         = 5;
Bulk         = 100;
Pitch        = 100;
Radius       = Step/10;
BField       = 0.0;
TauBe        = 89;
TauSe        = 89;
TauBh        = 65;
TauSh        = 65;
NAverage     = 1;
NParticles   = 10000;
ParticleType = 'beta';


rng default; % Reset random seed

[pt,et,tt,ut] = PDE_AllStrips(Pitch,Bulk);
[pw,ew,tw,uw] = PDE_WeightingField(Pitch,Bulk);

[VFieldx_e, VFieldy_e, VFieldx_h, VFieldy_h, x, y] =...
    VelocityField(ut,pt,tt,Step,Bulk,BField,2 * Pitch);

[WorkTransportTotal, x, y] =...
    ManyWorkTransport(uw,pw,tw,VFieldx_e,VFieldy_e,VFieldx_h,VFieldy_h,...
    x,y,Step,Bulk,Radius,TauBe,TauSe,TauBh,TauSh,NAverage);

fprintf('@@@ I''m interpolating the work transport matrix with a step of %d um @@@\n\n',Radius);
[xx,yy] = meshgrid(x,y);
subx = x(1):Radius:x(length(x));
suby = y(1):Radius:y(length(y));
[subxx,subyy] = meshgrid(subx,suby);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolation/Extrapolation: linear, nearest, cubic, spline %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%WorkTransportTotal_ = griddedInterpolant(xx', yy', WorkTransportTotal', 'spline', 'spline');
%subWorkTransportTotal = WorkTransportTotal_(subxx', subyy');
%subWorkTransportTotal = subWorkTransportTotal';
%%%%%%%%%%%%%%%%%
% Fit/Smoothing %
%%%%%%%%%%%%%%%%%
[x0, y0, z0] = prepareSurfaceData(x,y,WorkTransportTotal);
[WorkTransportTotal_, goodness, output] = fit([x0 y0],z0,'cubicinterp');
subWorkTransportTotal = WorkTransportTotal_(subxx, subyy);

figure (7);
colormap jet;
surf(subxx,subyy,subWorkTransportTotal,'EdgeColor','none');
%plot(WorkTransportTotal_,[x0 y0],z0); % To show and compare fit resuts
title('Interpolated Total <Work-Transport>');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Work / q [#charges * V]');

ComputeSpectra(subWorkTransportTotal,subx,suby,NParticles,Pitch,...
    Bulk,Radius,ParticleType);


%%%%%%%%%%%%%%%%%%%%%%%%
% Interesting commands %
%%%%%%%%%%%%%%%%%%%%%%%%
%[ux,uy] = pdegrad(p,t,u);
%ugrad = [ux;uy];
%pdeplot(p,e,t,'flowdata',ugrad);

%pdeplot(p,e,t,'zdata',u,'xygrid','on');

%dlmwrite('WeightingField.txt',uwxy,'\t'); Write a matrix into a txt file
