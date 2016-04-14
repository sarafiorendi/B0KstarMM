%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to plot a filed and its gradient obtained from the PDE tool %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u = potential
% p (points), t (triangles) = mesh
% Step  = Unit step of the lattice on which the field are computed [um]
% Bulk  = Bulk thickness [um]
% Pitch = Strip pitch [um]

function [] = PlotVectorField(u,p,t,Step,Bulk,Pitch)


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
ReSample   = 2;   % Used in order to have nice plots
ContLevel  = 40;  % Contour plot levels
MagnVector = 1.2; % Vector field magnification
x          = -Pitch:Step:Pitch;
y          = 0:Step:Bulk+Bulk/2;


%%%%%%%%%%%%%%%%%%%
% Start algorithm %
%%%%%%%%%%%%%%%%%%%
uxy      = tri2grid(p,t,u,x,y);
[Fx, Fy] = gradient(uxy,Step,Step);

Fx_ReSample = Fx(1:ReSample:length(y), 1:ReSample:length(x));
Fy_ReSample = Fy(1:ReSample:length(y), 1:ReSample:length(x));
xp          = x(1:ReSample:length(x));
yp          = y(1:ReSample:length(y));


%%%%%%%%%
% Plots %
%%%%%%%%%
figure (9);
contour(x,y,uxy,ContLevel);
hold on
quiver(xp,yp,Fx_ReSample,Fy_ReSample,MagnVector);
colormap jet;
hold off

title('Vector field');
xlabel('X');
ylabel('Y');
zlabel('Vector field amplitude');

end
