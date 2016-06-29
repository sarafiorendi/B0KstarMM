%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that calculates the "Work-Transport" matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% uw = weighting potential
% pw (points), tw (triangles) = mesh
% VFieldx_e, VFieldy_e = xy Velocity-Field components for electrons [um/ns]
% VFieldx_h, VFieldy_h = xy Velocity-Field components for holes [um/ns]
% x, y    = Axes
% Step    = Unit step of the lattice on which the field is computed [um]
% Bulk    = Bulk thickness [um]
% Radius  = Unit step of the movements and field interpolation [um]
% TauBe/h = Life-time on the backplane side [ns]
% TauSe/h = Life-time on the strip side [ns]

function [WorkTransportTotal, x, y] = ...
    WorkTransport(uw,pw,tw,VFieldx_e,VFieldy_e,VFieldx_h,VFieldy_h,...
    x,y,Step,Bulk,Radius,TauBe,TauSe,TauBh,TauSh)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
ReSample   = 2;   % Used in order to make nice plots
ContLevel  = 40;  % Contour plot levels
MagnVector = 1.2; % Vector field magnification


%%%%%%%%%%%%%%%%%%%
% Start algorithm %
%%%%%%%%%%%%%%%%%%%
uwxy = tri2grid(pw,tw,uw,x,y);
[Ewx, Ewy] = gradient(uwxy,Step,Step);

fprintf('@@@ I''m calculating the work for electrons @@@\n');
[WTransport_e] =...
    WorkTransportMatrix(VFieldx_e,VFieldy_e,Ewx,Ewy,x,y,...
    TauBe,TauSe,Step,Bulk,Radius,-1);

fprintf('@@@ I''m calculating the work for holes @@@\n');
[WTransport_h] =...
    WorkTransportMatrix(VFieldx_h,VFieldy_h,Ewx,Ewy,x,y,...
    TauBh,TauSh,Step,Bulk,Radius,+1);

WorkTransportTotal = WTransport_e + WTransport_h;


%%%%%%%%%
% Plots %
%%%%%%%%%
Ewx_ReSample = Ewx(1:ReSample:length(y), 1:ReSample:length(x));
Ewy_ReSample = Ewy(1:ReSample:length(y), 1:ReSample:length(x));
xp           = x(1:ReSample:length(x));
yp           = y(1:ReSample:length(y));
[xx,yy]      = meshgrid(x,y);

figure (3);
colormap jet;
subplot(1,2,1);
surf(xx,yy,uwxy,'EdgeColor','none');
title('Weighting potential');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Potential [V]');
subplot(1,2,2);
contour(x,y,uwxy,ContLevel);
hold on
quiver(xp,yp,Ewx_ReSample,Ewy_ReSample,MagnVector);
colormap jet;
hold off
title('Weighting field');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Weighting Field [V / \mum]');

figure (4);
colormap jet;
surf(xx,yy,WorkTransportTotal,'EdgeColor','none');
title('Total Work-Transport');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Work / q [#charges * V]');

figure (5);
colormap jet;
subplot(1,2,1);
surf(xx,yy,WTransport_e,'EdgeColor','none');
title('Electron Work-Transport');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Work / q [#charges * V]');
subplot(1,2,2);
surf(xx,yy,WTransport_h,'EdgeColor','none');
title('Hole Work-Transport');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Work / q [#charges * V]');

fprintf('CPU time --> %d[min]\n',(cputime-TStart)/60);
end
