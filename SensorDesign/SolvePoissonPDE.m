%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that solves Poisson equation to compute the potential %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk    = Bulk thickness [um]
% Pitch   = Strip pitch [um]
% BiasB   = Sensor backplane voltage [V] [0 Weighting; -200 All]
% BiasS   = Sensor strip voltage [V]
% BiasW   = Sensor central strip voltage [V] [1 Weighting; 0 All]
% epsR    = Relative dielectric constant [3.9 Silicon, 5.7 Diamond]
% rho     = Charge denisty in the bulk [(Coulomb / um^3) / eps0]
% ItFigIn = Figure iterator input

function [potential, ItFigOut] = SolvePoissonPDE(Bulk,Pitch,...
    BiasB,BiasS,BiasW,epsR,rho,ItFigIn)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
ReSampleFine   = 1;   % Used in order to make nice plots [um]
ReSampleCoarse = 10;  % Used in order to make nice plots [um]
ContLevel      = 40;  % Contour plot levels
MagnVector     = 1.2; % Vector field magnification
MeshMax        = 15;  % Maximum mesh edge length [um]

StrThick   = 5;  % Strip metalization thickness [um]
EleHVwidth = 80; % HV strip metalization width [um]
EleSGwidth = 80; % Signal strip metalization width [um]
NupBulk    = 2;  % Number of bulk thicknesses above sensor (included)
NStrips    = 21; % Total number of strips


%%%%%%%%%%%%%%%%%%%%
% Create PDE model %
%%%%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m solving Poisson equation to calculate the potential @@@\n');
pdem = createpde(1);


%%%%%%%%%%%%%%%%%%%%%%
% Create 2D geometry %
%%%%%%%%%%%%%%%%%%%%%%
gd = zeros(2+8,NStrips+2);
% Central strip
gd(:,1) = [ 3 4 -EleSGwidth/2 EleSGwidth/2 EleSGwidth/2 -EleSGwidth/2 ...
    Bulk+StrThick Bulk+StrThick Bulk Bulk ]';

for i = 1:(NStrips-1) / 2
    % Positive strips
    gd(:,i+1) = [ 3 4 i*Pitch-EleHVwidth/2 i*Pitch+EleHVwidth/2 i*Pitch+EleHVwidth/2 i*Pitch-EleHVwidth/2 ...
        Bulk+StrThick Bulk+StrThick Bulk Bulk ]';

    % Negative strips
    gd(:,(NStrips-1) / 2 + 1 + i) = [ 3 4 -(i*Pitch-EleHVwidth/2) -(i*Pitch+EleHVwidth/2) -(i*Pitch+EleHVwidth/2) -(i*Pitch-EleHVwidth/2) ...
        Bulk+StrThick Bulk+StrThick Bulk Bulk ]';
end

% Sensor volume
gd(:,NStrips+1) = [ 3 4 -((NStrips-1)/2*Pitch+Pitch/2) ((NStrips-1)/2*Pitch+Pitch/2) ((NStrips-1)/2*Pitch+Pitch/2) -((NStrips-1)/2*Pitch+Pitch/2) ...
    Bulk Bulk 0 0 ]';
% Whole volume
gd(:,NStrips+2) = [ 3 4 -((NStrips-1)/2*Pitch+Pitch/2) ((NStrips-1)/2*Pitch+Pitch/2) ((NStrips-1)/2*Pitch+Pitch/2) -((NStrips-1)/2*Pitch+Pitch/2) ...
    Bulk*NupBulk Bulk*NupBulk 0 0 ]';

sf = sprintf('R%d-(',NStrips+2);
ns = repmat(' ',[NStrips+2 3]);

for i = 1:NStrips
    sf = strcat(sf,sprintf('R%d',i),'+');
    ns(i,:) = sprintf('R%-2d',i);
end

sf([length(sf)]) = [];
sf = strcat(sf,sprintf(')+R%d',NStrips+1));

ns(NStrips+1,:) = sprintf('R%-2d',NStrips+1);
ns(NStrips+2,:) = sprintf('R%-2d',NStrips+2);
ns = ns';

dl = decsg(gd,sf,ns);
geometryFromEdges(pdem,dl);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply boundary conditions (only on conductors) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left edge air
applyBoundaryCondition(pdem,'edge',112,'q',0,'g',0);
% Left edge sensor
applyBoundaryCondition(pdem,'edge',111,'q',0,'g',0);
% Bottom edge
applyBoundaryCondition(pdem,'edge',110,'h',1,'r',BiasB);
% Right edge air
applyBoundaryCondition(pdem,'edge',109,'q',0,'g',0);
% Right edge sensor
applyBoundaryCondition(pdem,'edge',108,'q',0,'g',0);
% Bottom all strips
applyBoundaryCondition(pdem,'edge',106,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',104,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',102,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',100,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',98,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',96,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',94,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',92,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',90,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',88,'h',1,'r',BiasS);

applyBoundaryCondition(pdem,'edge',86,'h',1,'r',BiasW);

applyBoundaryCondition(pdem,'edge',84,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',82,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',80,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',78,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',76,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',74,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',72,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',70,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',68,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',66,'h',1,'r',0);
% Top all strips
applyBoundaryCondition(pdem,'edge',64,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',63,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',62,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',61,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',60,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',59,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',58,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',57,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',56,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',55,'h',1,'r',BiasS);

applyBoundaryCondition(pdem,'edge',54,'h',1,'r',BiasW);

applyBoundaryCondition(pdem,'edge',53,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',52,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',51,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',50,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',49,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',48,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',47,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',46,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',45,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',44,'h',1,'r',0);
% Top edge
applyBoundaryCondition(pdem,'edge',43,'h',1,'r',BiasS);
% Negative strips
applyBoundaryCondition(pdem,'edge',42,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',41,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',40,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',39,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',38,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',37,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',36,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',35,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',34,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',33,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',32,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',31,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',30,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',29,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',28,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',27,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',26,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',25,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',24,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',23,'h',1,'r',BiasS);
% Positive strips
applyBoundaryCondition(pdem,'edge',22,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',21,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',20,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',19,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',18,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',17,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',16,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',15,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',14,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',13,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',12,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',11,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',10,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',9,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',8,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',7,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',6,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',5,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',4,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',3,'h',1,'r',BiasS);
% Central strip
applyBoundaryCondition(pdem,'edge',2,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'edge',1,'h',1,'r',BiasW);


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
ylim([0,Bulk*NupBulk]);
title('Geometry');
xlabel('X [\mum]');
ylabel('Y [\mum]');
subplot(1,2,2);
pdegplot(pdem);
hold on;
pdemesh(pdem);
xlim([-Pitch * NStrips/2,+Pitch * NStrips/2]);
ylim([0,Bulk*NupBulk]);
title('Delaunay mesh');
xlabel('X [\mum]');
ylabel('Y [\mum]');
hold off;

ItFigIn = ItFigIn + 1;
figure(ItFigIn);
subplot(1,2,1);
colormap jet;
pdeplot(pdem,'xydata',potential.NodalSolution);
xlim([-Pitch * NStrips/2,+Pitch * NStrips/2]);
ylim([0,Bulk  * NupBulk]);
title('Potential');
xlabel('X [\mum]');
ylabel('Y [\mum]');
subplot(1,2,2);

x = -Pitch:ReSampleFine:Pitch;
y = 0:ReSampleFine:Bulk * 3/2;
[FineMeshX,FineMeshY] = meshgrid(x,y);
FineQuery = [FineMeshX(:),FineMeshY(:)]';

x = -Pitch:ReSampleCoarse:Pitch;
y = 0:ReSampleCoarse:Bulk*3/2;
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
ylabel('Y [\mum]');

ItFigOut = ItFigIn + 1;
fprintf('CPU time --> %d[min]\n\n',(cputime-TStart)/60);
end
