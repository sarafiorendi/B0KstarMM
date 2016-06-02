function [p,e,t,u] = PDE_ParallelElectrodes(Bulk)

% Parameters of the sensor:
DConst     = 3.9;    % Relative dielectric constant [3.9 Silicon, 5.7 Diamond]
rho        = 0.0;    % Charge denisty in the bulk [#charges / m^3]
rho        = rho*1.6e-19; % Charge denisty in the bulk [Coulomb / m^3]
NupBulk    = 3;      % Number of bulk thicknesses above sensor (included)
StrThick   = 5;      % Strip metalization thickness [um]
BiasB      = '-200'; % Sensor backplane voltage [V]
Width      = 1000;   % Plane width [um]

[pde_fig,ax]=pdeinit;
pdetool('appl_cb',5);
pdetool('snapon','on');
set(ax,'DataAspectRatio',[1 1 1]);
set(ax,'PlotBoxAspectRatio',[1 1 1]);
set(ax,'XLim',[-Width/2 Width/2]);
set(ax,'YLim',[0 Bulk*NupBulk]);
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');
pdetool('gridon','on');

% Geometry description:
pderect([-Width/2 Width/2 Bulk 0],'R1');
pderect([-Width/2 Width/2 Bulk+StrThick Bulk],'R2');
pderect([-Width/2 Width/2 Bulk*NupBulk Bulk+StrThick],'R3');
pderect([-Width/2 Width/2 Bulk*NupBulk 0],'R4');

set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','R1+R3+R4-R2')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(8,'dir',1,'1','0')
pdesetbd(7,'neu',1,'0','0')
pdesetbd(6,'dir',1,'1',BiasB)
pdesetbd(5,'dir',1,'1','0')
pdesetbd(4,'neu',1,'0','0')
pdesetbd(3,'dir',1,'1','0')
pdesetbd(2,'dir',1,'1','0')
pdesetbd(1,'dir',1,'1','0')

% Mesh generation:
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
pdetool('initmesh')
pdetool('refine')
pdetool('refine')

% PDE coefficients:
stringa1 = sprintf('1.0!%0.1f',DConst);
stringa2 = sprintf('%0.3f!0',rho);
stringa3 = stringa1;
stringa4 = stringa2;
for i = 1:length(stringa2)-length(stringa1)
    stringa3 = [stringa3,' '];
end
for i = 1:length(stringa1)-length(stringa2)
    stringa4 = [stringa4,' '];
end
pdeseteq(1,stringa1,'0.0!0.0',stringa2,'1.0!1.0','0:10','0.0','0.0','[0 100]')
setappdata(pde_fig,'currparam',[stringa3;stringa4])

% Solve parameters:
setappdata(pde_fig,'solveparam',...
str2mat('0','2928','10','pdeadworst',...
'0.5','longest','0','1E-4','','fixed','Inf'))

% Plotflags and user data strings:
setappdata(pde_fig,'plotflags',[1 1 1 1 1 1 1 1 0 0 0 1 1 1 0 1 0 1]);
setappdata(pde_fig,'colstring','');
setappdata(pde_fig,'arrowstring','');
setappdata(pde_fig,'deformstring','');
setappdata(pde_fig,'heightstring','');

% Solve PDE:
pdetool('solve')

% Export mesh and solution:
meshstat = getappdata(pde_fig,'meshstat');
n = length(meshstat);
p = getappdata(pde_fig,['p' int2str(n-1)]);
e = getappdata(pde_fig,['e' int2str(n-1)]);
t = getappdata(pde_fig,['t' int2str(n-1)]);
u = get(findobj(get(pde_fig,'Children'),'flat','Tag','PDEPlotMenu'),'UserData');

end
