%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation to compute the potential %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk    = Bulk thickness [um]
% PitchX  = Pixel pitch along X [um]
% PitchY  = Pixel pitch along Y [um]
% BiasB   = Sensor backplane voltage [V] [0 Weighting; -200 All]
% BiasS   = Sensor piel voltage [V]
% BiasW   = Sensor central pixel voltage [V] [1 Weighting; 0 All]
% epsR    = Relative dielectric constant [3.9 Silicon, 5.7 Diamond]
% rho     = Charge denisty in the bulk [(Coulomb / um^3) / eps0 [F/um]]
% ItFigIn = Figure iterator input

function [potential, ItFigOut] = SolvePoissonPDE3D(Bulk,PitchX,PitchY,...
    BiasB,BiasS,BiasW,epsR,rho,ItFigIn)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
MetalThick  =   5; % Metalization thickness [um]
MetalWidthX =  50; % Metalization width along X [um]
MetalWidthY = 100; % Metalization width along Y [um]
SHeight     =   2; % Sensor height [units of bulk thickness]
Step = 5;


%%%%%%%%%%%%%%%%%%%%
% Create PDE model %
%%%%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m solving Poisson equation to calculate the potential @@@\n');
pdem = createpde(1);


%%%%%%%%%%%%%%%%%%%%%%
% Create 3D geometry %
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
% Row 0 %
%%%%%%%%%
[x00,y00,z00] = meshgrid(-MetalWidthX/2:Step:MetalWidthX/2,...
    -MetalWidthY/2:Step:MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

x00 = x00(:);
y00 = y00(:);
z00 = z00(:);

[xp10,yp10,zp10] = meshgrid(1*PitchX-MetalWidthX/2:Step:1*PitchX+MetalWidthX/2,...
    -MetalWidthY/2:Step:MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xp10 = xp10(:);
yp10 = yp10(:);
zp10 = zp10(:);

[xp20,yp20,zp20] = meshgrid(2*PitchX-MetalWidthX/2:Step:2*PitchX+MetalWidthX/2,...
    -MetalWidthY/2:Step:MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xp20 = xp20(:);
yp20 = yp20(:);
zp20 = zp20(:);

[xm10,ym10,zm10] = meshgrid(-1*PitchX-MetalWidthX/2:Step:-1*PitchX+MetalWidthX/2,...
    -MetalWidthY/2:Step:MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xm10 = xm10(:);
ym10 = ym10(:);
zm10 = zm10(:);

[xm20,ym20,zm20] = meshgrid(-2*PitchX-MetalWidthX/2:Step:-2*PitchX+MetalWidthX/2,...
    -MetalWidthY/2:Step:MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xm20 = xm20(:);
ym20 = ym20(:);
zm20 = zm20(:);

%%%%%%%%%%
% Row +1 %
%%%%%%%%%%
[x0p1,y0p1,z0p1] = meshgrid(-MetalWidthX/2:Step:MetalWidthX/2,...
    1*PitchY-MetalWidthY/2:Step:1*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

x0p1 = x0p1(:);
y0p1 = y0p1(:);
z0p1 = z0p1(:);

[xp1p1,yp1p1,zp1p1] = meshgrid(1*PitchX-MetalWidthX/2:Step:1*PitchX+MetalWidthX/2,...
    1*PitchY-MetalWidthY/2:Step:1*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xp1p1 = xp1p1(:);
yp1p1 = yp1p1(:);
zp1p1 = zp1p1(:);

[xp2p1,yp2p1,zp2p1] = meshgrid(2*PitchX-MetalWidthX/2:Step:2*PitchX+MetalWidthX/2,...
    1*PitchY-MetalWidthY/2:Step:1*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xp2p1 = xp2p1(:);
yp2p1 = yp2p1(:);
zp2p1 = zp2p1(:);

[xm1p1,ym1p1,zm1p1] = meshgrid(-1*PitchX-MetalWidthX/2:Step:-1*PitchX+MetalWidthX/2,...
    1*PitchY-MetalWidthY/2:Step:1*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xm1p1 = xm1p1(:);
ym1p1 = ym1p1(:);
zm1p1 = zm1p1(:);

[xm2p1,ym2p1,zm2p1] = meshgrid(-2*PitchX-MetalWidthX/2:Step:-2*PitchX+MetalWidthX/2,...
    1*PitchY-MetalWidthY/2:Step:1*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xm2p1 = xm2p1(:);
ym2p1 = ym2p1(:);
zm2p1 = zm2p1(:);

%%%%%%%%%%
% Row +2 %
%%%%%%%%%%
[x0p2,y0p2,z0p2] = meshgrid(-MetalWidthX/2:Step:MetalWidthX/2,...
    2*PitchY-MetalWidthY/2:Step:2*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

x0p2 = x0p2(:);
y0p2 = y0p2(:);
z0p2 = z0p2(:);

[xp1p2,yp1p2,zp1p2] = meshgrid(1*PitchX-MetalWidthX/2:Step:1*PitchX+MetalWidthX/2,...
    2*PitchY-MetalWidthY/2:Step:2*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xp1p2 = xp1p2(:);
yp1p2 = yp1p2(:);
zp1p2 = zp1p2(:);

[xp2p2,yp2p2,zp2p2] = meshgrid(2*PitchX-MetalWidthX/2:Step:2*PitchX+MetalWidthX/2,...
    2*PitchY-MetalWidthY/2:Step:2*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xp2p2 = xp2p2(:);
yp2p2 = yp2p2(:);
zp2p2 = zp2p2(:);

[xm1p2,ym1p2,zm1p2] = meshgrid(-1*PitchX-MetalWidthX/2:Step:-1*PitchX+MetalWidthX/2,...
    2*PitchY-MetalWidthY/2:Step:2*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xm1p2 = xm1p2(:);
ym1p2 = ym1p2(:);
zm1p2 = zm1p2(:);

[xm2p2,ym2p2,zm2p2] = meshgrid(-2*PitchX-MetalWidthX/2:Step:-2*PitchX+MetalWidthX/2,...
    2*PitchY-MetalWidthY/2:Step:2*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xm2p2 = xm2p2(:);
ym2p2 = ym2p2(:);
zm2p2 = zm2p2(:);

%%%%%%%%%%
% Row -1 %
%%%%%%%%%%
[x0m1,y0m1,z0m1] = meshgrid(-MetalWidthX/2:Step:MetalWidthX/2,...
    -1*PitchY-MetalWidthY/2:Step:-1*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

x0m1 = x0m1(:);
y0m1 = y0m1(:);
z0m1 = z0m1(:);

[xp1m1,yp1m1,zp1m1] = meshgrid(1*PitchX-MetalWidthX/2:Step:1*PitchX+MetalWidthX/2,...
    -1*PitchY-MetalWidthY/2:Step:-1*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xp1m1 = xp1m1(:);
yp1m1 = yp1m1(:);
zp1m1 = zp1m1(:);

[xp2m1,yp2m1,zp2m1] = meshgrid(2*PitchX-MetalWidthX/2:Step:2*PitchX+MetalWidthX/2,...
    -1*PitchY-MetalWidthY/2:Step:-1*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xp2m1 = xp2m1(:);
yp2m1 = yp2m1(:);
zp2m1 = zp2m1(:);

[xm1m1,ym1m1,zm1m1] = meshgrid(-1*PitchX-MetalWidthX/2:Step:-1*PitchX+MetalWidthX/2,...
    -1*PitchY-MetalWidthY/2:Step:-1*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xm1m1 = xm1m1(:);
ym1m1 = ym1m1(:);
zm1m1 = zm1m1(:);

[xm2m1,ym2m1,zm2m1] = meshgrid(-2*PitchX-MetalWidthX/2:Step:-2*PitchX+MetalWidthX/2,...
    -1*PitchY-MetalWidthY/2:Step:-1*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xm2m1 = xm2m1(:);
ym2m1 = ym2m1(:);
zm2m1 = zm2m1(:);

%%%%%%%%%%
% Row -2 %
%%%%%%%%%%
[x0m2,y0m2,z0m2] = meshgrid(-MetalWidthX/2:Step:MetalWidthX/2,...
    -2*PitchY-MetalWidthY/2:Step:-2*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

x0m2 = x0m2(:);
y0m2 = y0m2(:);
z0m2 = z0m2(:);

[xp1m2,yp1m2,zp1m2] = meshgrid(1*PitchX-MetalWidthX/2:Step:1*PitchX+MetalWidthX/2,...
    -2*PitchY-MetalWidthY/2:Step:-2*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xp1m2 = xp1m2(:);
yp1m2 = yp1m2(:);
zp1m2 = zp1m2(:);

[xp2m2,yp2m2,zp2m2] = meshgrid(2*PitchX-MetalWidthX/2:Step:2*PitchX+MetalWidthX/2,...
    -2*PitchY-MetalWidthY/2:Step:-2*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xp2m2 = xp2m2(:);
yp2m2 = yp2m2(:);
zp2m2 = zp2m2(:);

[xm1m2,ym1m2,zm1m2] = meshgrid(-1*PitchX-MetalWidthX/2:Step:-1*PitchX+MetalWidthX/2,...
    -2*PitchY-MetalWidthY/2:Step:-2*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xm1m2 = xm1m2(:);
ym1m2 = ym1m2(:);
zm1m2 = zm1m2(:);

[xm2m2,ym2m2,zm2m2] = meshgrid(-2*PitchX-MetalWidthX/2:Step:-2*PitchX+MetalWidthX/2,...
    -2*PitchY-MetalWidthY/2:Step:-2*PitchY+MetalWidthY/2,...
    Bulk:Step:Bulk+MetalThick);

xm2m2 = xm2m2(:);
ym2m2 = ym2m2(:);
zm2m2 = zm2m2(:);

%%%%%%%%%%%%%%%%
% Whole volume %
%%%%%%%%%%%%%%%%
[x,y,z] = meshgrid(-2*PitchX-MetalWidthX:Step:2*PitchX+MetalWidthX,...
    -2*PitchY-MetalWidthY:Step:2*PitchY+MetalWidthY,...
    0:Step:Bulk*SHeight);

x = x(:);
y = y(:);
z = z(:);


%%%%%%%%%%%%%%%%%%%%%%
% Combine all pixels %
%%%%%%%%%%%%%%%%%%%%%%
xPx = [xm2m2' xm1m2' x0m2' xp1m2' xp2m2'...
       xm2m1' xm1m1' x0m1' xp1m1' xp2m1'...
       xm20'  xm10'  x00'  xp10'  xp20'...
       xm2p1' xm1p1' x0p1' xp1p1' xp2p1'...
       xm2p2' xm1p2' x0p2' xp1p2' xp2p2']';
yPx = [ym2m2' ym1m2' y0m2' yp1m2' yp2m2'...
       ym2m1' ym1m1' y0m1' yp1m1' yp2m1'...
       ym20'  ym10'  y00'  yp10'  yp20'...
       ym2p1' ym1p1' y0p1' yp1p1' yp2p1'...
       ym2p2' ym1p2' y0p2' yp1p2' yp2p2']';
zPx = [zm2m2' zm1m2' z0m2' zp1m2' zp2m2'...
       zm2m1' zm1m1' z0m1' zp1m1' zp2m1'...
       zm20'  zm10'  z00'  zp10'  zp20'...
       zm2p1' zm1p1' z0p1' zp1p1' zp2p1'...
       zm2p2' zm1p2' z0p2' zp1p2' zp2p2']';

shp = alphaShape(xPx,yPx,zPx);


%%%%%%%%%%%%%%%%%%%%%%%
% Create final volume %
%%%%%%%%%%%%%%%%%%%%%%%
in = inShape(shp,x,y,z);
x = x(~in);
y = y(~in);
z = z(~in);
shp = alphaShape(x,y,z);

[elements,nodes] = boundaryFacets(shp);
nodes = nodes';
elements = elements';


%%%%%%%%%%%%%%%%%%%%%%%%%
% Create final geometry %
%%%%%%%%%%%%%%%%%%%%%%%%%
geometryFromMesh(pdem,nodes,elements);


%%%%%%%%%
% Plots %
%%%%%%%%%
h = pdegplot(pdem,'FaceLabels','on');
h(1).FaceAlpha = 0.5;


%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%
% Bottom edge
applyBoundaryCondition(pdem,'face',1,'h',1,'r',BiasB);
% Top edge
applyBoundaryCondition(pdem,'face',2,'h',1,'r',0);
% Left edge
applyBoundaryCondition(pdem,'face',4,'h',1,'r',0);
% Front edge
applyBoundaryCondition(pdem,'face',9,'h',1,'r',0);
% Back edge
applyBoundaryCondition(pdem,'face',15,'h',1,'r',0);
% Right edge
applyBoundaryCondition(pdem,'face',17,'h',1,'r',0);

% Central pixel
applyBoundaryCondition(pdem,'face',14,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'face',32,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'face',77,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'face',113,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'face',157,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'face',169,'h',1,'r',BiasW);


%%%%%%%%%%%%%%%
% Create mesh %
%%%%%%%%%%%%%%%
msh = generateMesh(pdem,'GeometricOrder','quadratic');


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%
specifyCoefficients(pdem,'m',0,'d',0,'c',1,'a',0,'f',0);
potential = solvepde(pdem);


%%%%%%%%%
% Plots %
%%%%%%%%%
figure(ItFigIn);
h = pdegplot(pdem,'FaceLabels','on');
h(1).FaceAlpha = 0.5;
title('Geometry');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Z [\mum]');

ItFigIn = ItFigIn + 1;
figure(ItFigIn);
u = potential.NodalSolution;
pdeplot3D(pdem,'colormapdata',u);
title('Delaunay mesh');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Z [\mum]');

ItFigIn = ItFigIn + 1;
figure(ItFigIn);
colormap jet;
[x,y,z] = meshgrid(-PitchX:Step:PitchX,...
    -PitchY:Step:PitchY,...
    0:Step:Bulk*3/2);
V = interpolateSolution(potential,x,y,z);
V = reshape(V,size(x));
contourslice(x,y,z,V,[],[],0:Step:Bulk);
title('Potential and its gradient');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Z [\mum]');
colorbar;

ItFigOut = ItFigIn + 1;
fprintf('CPU time --> %.2f[min]\n\n',(cputime-TStart)/60);
end