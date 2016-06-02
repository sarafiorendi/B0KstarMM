function [pt,et,tt,ut] = PDE_AllStrips(Pitch,Bulk)

% Parameters of the sensor:
DConst     = 3.9;   % Relative dielectric constant [3.9 Silicon, 5.7 Diamond]
epsilon    = 8.854187817e-12; % Vacuum diectric constant [Farad / m]
rho        = -1.0;  % Charge denisty in the bulk [#charges / um^2 / epsilon]
rho        = rho*1.6e-7/epsilon; % Charge denisty in the bulk [Coulomb / m^2]
NupBulk    = 5;     % Number of bulk thicknesses above sensor (included)
StrThick   = 5;     % Strip metalization thickness [um]
EleHVwidth = 80;    % HV strip metalization width [um]
EleSGwidth = 80;    % Signal strip metalization width [um]
% The Pitch shouldn't be changed because a change in the Pitch can be
% emulated with a change in the bulk thickness and strip metalization width
BiasB      = '-200';% Sensor Backplane voltage [V]
BiasS      = '0';   % Sensor Strip voltage [V]
BiasW      = '0';   % Sensor central Strip voltage for weighting field [V]
ShowNstrip = 8;     % Number of strips to show in the plot

[pde_fig,ax]=pdeinit;
pdetool('appl_cb',5);
pdetool('snapon','on');
set(ax,'DataAspectRatio',[2 1 1]);
set(ax,'PlotBoxAspectRatio',[2 1 1]);
set(ax,'XLim',[-ShowNstrip/2*Pitch ShowNstrip/2*Pitch]);
set(ax,'YLim',[0 Bulk*NupBulk]);
set(ax,'XTickMode','auto');
set(ax,'YTickMode','auto');
pdetool('gridon','on');

% Geometry description:
pderect([EleSGwidth/2 -EleSGwidth/2 Bulk+StrThick Bulk],'R1');

pderect([Pitch-EleHVwidth/2 Pitch+EleHVwidth/2 Bulk+StrThick Bulk],'R2');
pderect([2*Pitch-EleSGwidth/2 2*Pitch+EleSGwidth/2 Bulk+StrThick Bulk],'R3');
pderect([3*Pitch-EleHVwidth/2 3*Pitch+EleHVwidth/2 Bulk+StrThick Bulk],'R4');
pderect([4*Pitch-EleSGwidth/2 4*Pitch+EleSGwidth/2 Bulk+StrThick Bulk],'R5');
pderect([5*Pitch-EleHVwidth/2 5*Pitch+EleHVwidth/2 Bulk+StrThick Bulk],'R6');
pderect([6*Pitch-EleSGwidth/2 6*Pitch+EleSGwidth/2 Bulk+StrThick Bulk],'R7');
pderect([7*Pitch-EleHVwidth/2 7*Pitch+EleHVwidth/2 Bulk+StrThick Bulk],'R8');
pderect([8*Pitch-EleSGwidth/2 8*Pitch+EleSGwidth/2 Bulk+StrThick Bulk],'R9');
pderect([9*Pitch-EleHVwidth/2 9*Pitch+EleHVwidth/2 Bulk+StrThick Bulk],'R10');
pderect([10*Pitch-EleSGwidth/2 10*Pitch+EleSGwidth/2 Bulk+StrThick Bulk],'R11');

pderect([-(Pitch-EleHVwidth/2) -(Pitch+EleHVwidth/2) Bulk+StrThick Bulk],'R12');
pderect([-(2*Pitch-EleSGwidth/2) -(2*Pitch+EleSGwidth/2) Bulk+StrThick Bulk],'R13');
pderect([-(3*Pitch-EleHVwidth/2) -(3*Pitch+EleHVwidth/2) Bulk+StrThick Bulk],'R14');
pderect([-(4*Pitch-EleSGwidth/2) -(4*Pitch+EleSGwidth/2) Bulk+StrThick Bulk],'R15');
pderect([-(5*Pitch-EleHVwidth/2) -(5*Pitch+EleHVwidth/2) Bulk+StrThick Bulk],'R16');
pderect([-(6*Pitch-EleSGwidth/2) -(6*Pitch+EleSGwidth/2) Bulk+StrThick Bulk],'R17');
pderect([-(7*Pitch-EleHVwidth/2) -(7*Pitch+EleHVwidth/2) Bulk+StrThick Bulk],'R18');
pderect([-(8*Pitch-EleSGwidth/2) -(8*Pitch+EleSGwidth/2) Bulk+StrThick Bulk],'R19');
pderect([-(9*Pitch-EleHVwidth/2) -(9*Pitch+EleHVwidth/2) Bulk+StrThick Bulk],'R20');
pderect([-(10*Pitch-EleSGwidth/2) -(10*Pitch+EleSGwidth/2) Bulk+StrThick Bulk],'R21');

pderect([-(10*Pitch+Pitch/2) (10*Pitch+Pitch/2) Bulk 0],'R22');
pderect([-(10*Pitch+Pitch/2) (10*Pitch+Pitch/2) Bulk*NupBulk 0],'R23');

set(findobj(get(pde_fig,'Children'),'Tag','PDEEval'),'String','R23-(R1+R2+R3+R4+R5+R6+R7+R8+R9+R10+R11+R12+R13+R14+R15+R16+R17+R18+R19+R20+R21)+R22')

% Boundary conditions:
pdetool('changemode',0)
pdesetbd(112,'neu',1,'0','0')
pdesetbd(111,'neu',1,'0','0')
pdesetbd(110,'dir',1,'1',BiasB)
pdesetbd(109,'neu',1,'0','0')
pdesetbd(108,'neu',1,'0','0')
pdesetbd(106,'dir',1,'1','0')
pdesetbd(104,'dir',1,'1',BiasS)
pdesetbd(102,'dir',1,'1','0')
pdesetbd(100,'dir',1,'1',BiasS)
pdesetbd(98,'dir',1,'1','0')
pdesetbd(96,'dir',1,'1',BiasS)
pdesetbd(94,'dir',1,'1','0')
pdesetbd(92,'dir',1,'1',BiasS)
pdesetbd(90,'dir',1,'1','0')
pdesetbd(88,'dir',1,'1',BiasS)
pdesetbd(86,'dir',1,'1',BiasW)
pdesetbd(84,'dir',1,'1',BiasS)
pdesetbd(82,'dir',1,'1','0')
pdesetbd(80,'dir',1,'1',BiasS)
pdesetbd(78,'dir',1,'1','0')
pdesetbd(76,'dir',1,'1',BiasS)
pdesetbd(74,'dir',1,'1','0')
pdesetbd(72,'dir',1,'1',BiasS)
pdesetbd(70,'dir',1,'1','0')
pdesetbd(68,'dir',1,'1',BiasS)
pdesetbd(66,'dir',1,'1','0')
pdesetbd(64,'dir',1,'1','0')
pdesetbd(63,'dir',1,'1',BiasS)
pdesetbd(62,'dir',1,'1','0')
pdesetbd(61,'dir',1,'1',BiasS)
pdesetbd(60,'dir',1,'1','0')
pdesetbd(59,'dir',1,'1',BiasS)
pdesetbd(58,'dir',1,'1','0')
pdesetbd(57,'dir',1,'1',BiasS)
pdesetbd(56,'dir',1,'1','0')
pdesetbd(55,'dir',1,'1',BiasS)
pdesetbd(54,'dir',1,'1',BiasW)
pdesetbd(53,'dir',1,'1',BiasS)
pdesetbd(52,'dir',1,'1','0')
pdesetbd(51,'dir',1,'1',BiasS)
pdesetbd(50,'dir',1,'1','0')
pdesetbd(49,'dir',1,'1',BiasS)
pdesetbd(48,'dir',1,'1','0')
pdesetbd(47,'dir',1,'1',BiasS)
pdesetbd(46,'dir',1,'1','0')
pdesetbd(45,'dir',1,'1',BiasS)
pdesetbd(44,'dir',1,'1','0')
pdesetbd(43,'dir',1,'1','0')
pdesetbd(42,'dir',1,'1','0')
pdesetbd(41,'dir',1,'1','0')
pdesetbd(40,'dir',1,'1',BiasS)
pdesetbd(39,'dir',1,'1',BiasS)
pdesetbd(38,'dir',1,'1','0')
pdesetbd(37,'dir',1,'1','0')
pdesetbd(36,'dir',1,'1',BiasS)
pdesetbd(35,'dir',1,'1',BiasS)
pdesetbd(34,'dir',1,'1','0')
pdesetbd(33,'dir',1,'1','0')
pdesetbd(32,'dir',1,'1',BiasS)
pdesetbd(31,'dir',1,'1',BiasS)
pdesetbd(30,'dir',1,'1','0')
pdesetbd(29,'dir',1,'1','0')
pdesetbd(28,'dir',1,'1',BiasS)
pdesetbd(27,'dir',1,'1',BiasS)
pdesetbd(26,'dir',1,'1','0')
pdesetbd(25,'dir',1,'1','0')
pdesetbd(24,'dir',1,'1',BiasS)
pdesetbd(23,'dir',1,'1',BiasS)
pdesetbd(22,'dir',1,'1','0')
pdesetbd(21,'dir',1,'1','0')
pdesetbd(20,'dir',1,'1',BiasS)
pdesetbd(19,'dir',1,'1',BiasS)
pdesetbd(18,'dir',1,'1','0')
pdesetbd(17,'dir',1,'1','0')
pdesetbd(16,'dir',1,'1',BiasS)
pdesetbd(15,'dir',1,'1',BiasS)
pdesetbd(14,'dir',1,'1','0')
pdesetbd(13,'dir',1,'1','0')
pdesetbd(12,'dir',1,'1',BiasS)
pdesetbd(11,'dir',1,'1',BiasS)
pdesetbd(10,'dir',1,'1','0')
pdesetbd(9,'dir',1,'1','0')
pdesetbd(8,'dir',1,'1',BiasS)
pdesetbd(7,'dir',1,'1',BiasS)
pdesetbd(6,'dir',1,'1','0')
pdesetbd(5,'dir',1,'1','0')
pdesetbd(4,'dir',1,'1',BiasS)
pdesetbd(3,'dir',1,'1',BiasS)
pdesetbd(2,'dir',1,'1',BiasW)
pdesetbd(1,'dir',1,'1',BiasW)

% Mesh generation:
if (Pitch/Bulk > 1.5)
    setappdata(pde_fig,'trisize',100);
end
setappdata(pde_fig,'Hgrad',1.3);
setappdata(pde_fig,'refinemethod','regular');
setappdata(pde_fig,'jiggle',char('on','mean',''));
pdetool('initmesh')
pdetool('refine')

% PDE coefficients:
stringa1 = sprintf('1.0!%0.1f',DConst);
stringa2 = sprintf('0!%0.3f',rho);
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
str2mat('0','82368','10','pdeadworst','0.5','longest','0','1E-4','','fixed','Inf'))

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
n  = length(meshstat);
pt = getappdata(pde_fig,['p' int2str(n-1)]);
et = getappdata(pde_fig,['e' int2str(n-1)]);
tt = getappdata(pde_fig,['t' int2str(n-1)]);
ut = get(findobj(get(pde_fig,'Children'),'flat','Tag','PDEPlotMenu'),'UserData');

end