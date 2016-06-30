%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that solves Poisson equation %
% to compute the weighting field        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [gd,sf,ns] = Geometry description
% Bulk       = Bulk thickness [um]
% Pitch      = Strip pitch [um]
% itFig      = Figure iterator

function [potential] = SolveWeightingFieldPDE(gd,sf,ns,Bulk,Pitch,itFig)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
ContLevel  = 30; % Contour plot levels
MagnVector = 2;  % Vector field magnification
MeshMax    = 15; % Maximum mesh edge length
Step    = 5; % Unit step of the lattice on which the field is computed [um]

BiasB   = 0; % Sensor backplane voltage [V]
BiasS   = 0; % Sensor strip voltage [V]
BiasW   = 1; % Sensor central strip voltage for weighting field [V]

NupBulk = 3; % Number of bulk thicknesses above sensor (included)
Nstrips = 21;% Total number of strips


%%%%%%%%%%%%%%%%%%%%
% Create PDE model %
%%%%%%%%%%%%%%%%%%%%
fprintf('@@@ I''m solving Poisson equation to calculate the weighting field @@@\n');
pdem = createpde(1);


%%%%%%%%%%%%%%%%%%%%%%
% Create 2D geometry %
%%%%%%%%%%%%%%%%%%%%%%
dl = decsg(gd,sf,ns);
geometryFromEdges(pdem,dl);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply boundary conditions (only on conductors) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
applyBoundaryCondition(pdem,'edge',112,'q',0,'g',0);
applyBoundaryCondition(pdem,'edge',111,'q',0,'g',0);
applyBoundaryCondition(pdem,'edge',110,'h',1,'r',BiasB);
applyBoundaryCondition(pdem,'edge',109,'q',0,'g',0);
applyBoundaryCondition(pdem,'edge',108,'q',0,'g',0);
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
applyBoundaryCondition(pdem,'edge',51,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',50,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',49,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',48,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',47,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',46,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',45,'h',1,'r',BiasS);
applyBoundaryCondition(pdem,'edge',44,'h',1,'r',0);
applyBoundaryCondition(pdem,'edge',43,'h',1,'r',BiasS);
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
applyBoundaryCondition(pdem,'edge',2,'h',1,'r',BiasW);
applyBoundaryCondition(pdem,'edge',1,'h',1,'r',BiasS);


%%%%%%%%%%%%%%%%%
% Generate mesh %
%%%%%%%%%%%%%%%%%
msh = generateMesh(pdem,'Hmax',MeshMax,'Jiggle','minimum');


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Poisson equation %
%%%%%%%%%%%%%%%%%%%%%%%%%%
specifyCoefficients(pdem,'m',0,'d',0,'c',1,'a',0,'f',0);
potential = solvepde(pdem);


%%%%%%%%%
% Plots %
%%%%%%%%%
figure (itFig);
pdegplot(pdem,'EdgeLabels','on');
hold on;
pdemesh(pdem);
xlim([-Pitch*Nstrips/2,+Pitch*Nstrips/2]);
ylim([0,Bulk*NupBulk]);
title('Delaunay mesh');
xlabel('x');
ylabel('y');
hold off;

figure (itFig+1);
subplot(1,2,1);
colormap jet;
pdeplot(pdem,'xydata',potential.NodalSolution);
xlim([-Pitch*Nstrips/2,+Pitch*Nstrips/2]);
ylim([0,Bulk*NupBulk]);
title('Weighting field');
xlabel('x');
ylabel('y');

subplot(1,2,2);
x = -Pitch:Step:Pitch;
y = 0:Step:Bulk*3/2;
[mshx,mshy] = meshgrid(x,y);
querypoints = [mshx(:),mshy(:)]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recompute solution on a different mesh %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interp = interpolateSolution(potential,querypoints);
interp = reshape(interp,size(mshx));


%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate the gradient %
%%%%%%%%%%%%%%%%%%%%%%%%%
[gradx,grady] = evaluateGradient(potential,querypoints);


contour(mshx,mshy,interp,ContLevel);
hold on;
quiver(mshx(:),mshy(:),gradx,grady,MagnVector);
hold off;

title('Weighting field contour and gradient');
xlabel('x');
ylabel('y');

fprintf('CPU time --> %d[min]\n\n',(cputime-TStart)/60);
end
