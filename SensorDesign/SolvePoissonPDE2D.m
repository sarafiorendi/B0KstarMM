%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation to compute the potential %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk    = Bulk thickness [um]
% Pitch   = Strip pitch [um]
% BiasB   = Sensor backplane voltage [V] [0 Weighting; -V All]
% BiasS   = Sensor strip voltage [V]
% BiasW   = Sensor central strip voltage [V] [1 Weighting; 0 All]
% epsR    = Relative permittivity
% rho     = Charge denisty in the bulk [(Coulomb/um^3) / eps0 [F/um]]
% XQ      = Coordinate for potential query along y [um]
% ItFigIn = Figure iterator input

function [potential, Sq, yq, ItFigOut] = SolvePoissonPDE2D(Bulk,Pitch,...
    BiasB,BiasS,BiasW,epsR,rho,XQ,ItFigIn)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
ReSampleFine   = 1;   % Used in order to make nice plots [um]
ReSampleCoarse = 10;  % Used in order to make nice plots [um]
ContLevel      = 40;  % Contour plot levels
MagnVector     = 1.5; % Vector field magnification
MeshMax        = 15;  % Maximum mesh edge length [um]

MetalThick   = 5;        % Metalization thickness [um]
MetalWidthHV = Pitch-10; % Metalization width HV strip [um]
MetalWidthSg = Pitch-10; % Metalization width Signal strip [um]
SHeight      = 2;        % Sensor height [units of bulk thickness]
NStrips      = 13;       % Total number of strips


%%%%%%%%%%%%%%%%%%%%
% Create PDE model %
%%%%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m solving Poisson equation in 2D to calculate the potential @@@\n');
pdem = createpde(1);


%%%%%%%%%%%%%%%%%%%%%%
% Create 2D geometry %
%%%%%%%%%%%%%%%%%%%%%%
% Central strip
R1 = [ 3 4 -MetalWidthSg/2 MetalWidthSg/2 MetalWidthSg/2 -MetalWidthSg/2 ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
% Positive strips
R2 = [ 3 4 1*Pitch-MetalWidthHV/2 1*Pitch+MetalWidthHV/2 1*Pitch+MetalWidthHV/2 1*Pitch-MetalWidthHV/2 ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
R3 = [ 3 4 2*Pitch-MetalWidthHV/2 2*Pitch+MetalWidthHV/2 2*Pitch+MetalWidthHV/2 2*Pitch-MetalWidthHV/2 ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
R4 = [ 3 4 3*Pitch-MetalWidthHV/2 3*Pitch+MetalWidthHV/2 3*Pitch+MetalWidthHV/2 3*Pitch-MetalWidthHV/2 ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
R5 = [ 3 4 4*Pitch-MetalWidthHV/2 4*Pitch+MetalWidthHV/2 4*Pitch+MetalWidthHV/2 4*Pitch-MetalWidthHV/2 ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
R6 = [ 3 4 5*Pitch-MetalWidthHV/2 5*Pitch+MetalWidthHV/2 5*Pitch+MetalWidthHV/2 5*Pitch-MetalWidthHV/2 ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
R7 = [ 3 4 6*Pitch-MetalWidthHV/2 6*Pitch+MetalWidthHV/2 6*Pitch+MetalWidthHV/2 6*Pitch-MetalWidthHV/2 ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
% Negative strips
R8 = [ 3 4 -(1*Pitch-MetalWidthHV/2) -(1*Pitch+MetalWidthHV/2) -(1*Pitch+MetalWidthHV/2) -(1*Pitch-MetalWidthHV/2) ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
R9 = [ 3 4 -(2*Pitch-MetalWidthHV/2) -(2*Pitch+MetalWidthHV/2) -(2*Pitch+MetalWidthHV/2) -(2*Pitch-MetalWidthHV/2) ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
R10 = [ 3 4 -(3*Pitch-MetalWidthHV/2) -(3*Pitch+MetalWidthHV/2) -(3*Pitch+MetalWidthHV/2) -(3*Pitch-MetalWidthHV/2) ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
R11 = [ 3 4 -(4*Pitch-MetalWidthHV/2) -(4*Pitch+MetalWidthHV/2) -(4*Pitch+MetalWidthHV/2) -(4*Pitch-MetalWidthHV/2) ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
R12 = [ 3 4 -(5*Pitch-MetalWidthHV/2) -(5*Pitch+MetalWidthHV/2) -(5*Pitch+MetalWidthHV/2) -(5*Pitch-MetalWidthHV/2) ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
R13 = [ 3 4 -(6*Pitch-MetalWidthHV/2) -(6*Pitch+MetalWidthHV/2) -(6*Pitch+MetalWidthHV/2) -(6*Pitch-MetalWidthHV/2) ...
    Bulk+MetalThick Bulk+MetalThick Bulk Bulk ]';
% Sensor volume
R14 = [ 3 4 -(6*Pitch+Pitch/2) (6*Pitch+Pitch/2) (6*Pitch+Pitch/2) -(6*Pitch+Pitch/2) ...
    Bulk Bulk 0 0 ]';
% Whole volume
R15 = [ 3 4 -(6*Pitch+Pitch/2) (6*Pitch+Pitch/2) (6*Pitch+Pitch/2) -(6*Pitch+Pitch/2) ...
    Bulk*SHeight Bulk*SHeight 0 0 ]';

gd = [R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14,R15];
sf = 'R15-(R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12+R13)+R14';
ns = char('R1','R2','R3','R4','R5','R6','R7','R8','R9','R10','R11','R12','R13','R14','R15');
ns = ns';

dl = decsg(gd,sf,ns);
geometryFromEdges(pdem,dl);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply boundary conditions (only on conductors) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Central strip
applyBoundaryCondition(pdem,'edge',1,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'edge',2,'h',1,'r',BiasW);
% Positive strips
applyBoundaryCondition(pdem,'edge',3,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',4,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',5,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',6,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',7,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',8,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',9,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',10,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',11,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',12,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',13,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',14,'h',1,'r',0);
% Negative strips
applyBoundaryCondition(pdem,'edge',15,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',16,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',17,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',18,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',19,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',20,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',21,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',22,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',23,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',24,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',25,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',26,'h',1,'r',0);
% Top edge
applyBoundaryCondition(pdem,'edge',27,'h',1,'r',BiasS);
% Top all strips
applyBoundaryCondition(pdem,'edge',28,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',29,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',30,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',31,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',32,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',33,'h',1,'r',BiasS);

applyBoundaryCondition(pdem,'edge',34,'h',1,'r',BiasW);

applyBoundaryCondition(pdem,'edge',35,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',36,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',37,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',38,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',39,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',40,'h',1,'r',0);
% Bottom all strips
applyBoundaryCondition(pdem,'edge',42,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',44,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',46,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',48,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',50,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',52,'h',1,'r',BiasS);

applyBoundaryCondition(pdem,'edge',54,'h',1,'r',BiasW);

applyBoundaryCondition(pdem,'edge',56,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',58,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',60,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',62,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',64,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',66,'h',1,'r',0);
% Right edge sensor
applyBoundaryCondition(pdem,'edge',68,'q',0,'g',0);
% Right edge air
applyBoundaryCondition(pdem,'edge',69,'q',0,'g',0);
% Bottom edge
applyBoundaryCondition(pdem,'edge',70,'h',1,'r',BiasB);
% Left edge sensor
applyBoundaryCondition(pdem,'edge',71,'q',0,'g',0);
% Left edge air
applyBoundaryCondition(pdem,'edge',72,'q',0,'g',0);


%%%%%%%%%%%%%%%%%
% Generate mesh %
%%%%%%%%%%%%%%%%%
msh = generateMesh(pdem,'Hmax',MeshMax,'Jiggle','mean',...
    'GeometricOrder','quadratic','MesherVersion','R2013a');


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%
specifyCoefficients(pdem,'m',0,'d',0,'c',1,   'a',0,'f',0,  'face',1); % Air
specifyCoefficients(pdem,'m',0,'d',0,'c',epsR,'a',0,'f',rho,'face',2); % Sensor
potential = solvepde(pdem);


%%%%%%%%%
% Plots %
%%%%%%%%%
figure(ItFigIn);
subplot(1,2,1);
pdegplot(dl,'EdgeLabels','on','SubdomainLabels','on');
xlim([-Pitch * NStrips/2,+Pitch * NStrips/2]);
ylim([0,Bulk * SHeight]);
title('Geometry');
xlabel('X [\mum]');
ylabel('Z [\mum]');
subplot(1,2,2);
pdegplot(pdem);
hold on;
pdemesh(pdem);
xlim([-Pitch * NStrips/2,+Pitch * NStrips/2]);
ylim([0,Bulk * SHeight]);
hold off;
title('Delaunay mesh');
xlabel('X [\mum]');
ylabel('Z [\mum]');

ItFigIn = ItFigIn + 1;
figure(ItFigIn);
subplot(1,2,1);
colormap jet;
pdeplot(pdem,'xydata',potential.NodalSolution);
xlim([-Pitch * NStrips/2,+Pitch * NStrips/2]);
ylim([0,Bulk * SHeight]);
title('Potential');
xlabel('X [\mum]');
ylabel('Z [\mum]');

subplot(1,2,2);
colormap jet;

x = -Pitch:ReSampleFine:Pitch;
y = 0:ReSampleFine:Bulk * 3/2;
[FineMeshX,FineMeshY] = meshgrid(x,y);
FineQuery = [FineMeshX(:),FineMeshY(:)]';

x = -Pitch:ReSampleCoarse:Pitch;
y = 0:ReSampleCoarse:Bulk * 3/2;
[CoarseMeshX,CoarseMeshY] = meshgrid(x,y);
CoarseQuery = [CoarseMeshX(:),CoarseMeshY(:)]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recompute solution on a different mesh %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interp = interpolateSolution(potential,FineQuery);
interp = reshape(interp,size(FineMeshX));


%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate the gradient %
%%%%%%%%%%%%%%%%%%%%%%%%%
[gradx,grady] = evaluateGradient(potential,CoarseQuery);


contour(FineMeshX,FineMeshY,interp,ContLevel);
hold on;
quiver(CoarseMeshX(:),CoarseMeshY(:),gradx,grady,MagnVector);
hold off;

title('Potential and its gradient');
xlabel('X [\mum]');
ylabel('Z [\mum]');

ItFigIn = ItFigIn + 1;
figure(ItFigIn);
yq = 0:ReSampleFine:Bulk*3/2;
xq = XQ * ones(1,length(yq));
Sq = interpolateSolution(potential,xq,yq);
plot(yq,Sq);
title(sprintf('Potential along z at x = %.2f um',XQ));
xlabel('Z [\mum]');
ylabel('Potential');
grid on;

ItFigOut = ItFigIn + 1;
fprintf('CPU time --> %.2f[min]\n\n',(cputime-TStart)/60);
end
