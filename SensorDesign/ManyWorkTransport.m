%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate several "Work-Transport" matrices and compute the average %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAverage = Generate NAverage WorkTransportTotal and average them

function [WorkTransportTotal, x, y, ItFigOut] = ...
    ManyWorkTransport(potential,VFieldx_e,VFieldy_e,VFieldx_h,VFieldy_h,...
    x,y,Step,Bulk,Radius,TauBe,TauSe,TauBh,TauSh,NAverage,ItFigIn)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%
% Start algorithm %
%%%%%%%%%%%%%%%%%%%
[WorkTransportTotal_, x, y, ItFigOut] =...
    WorkTransport(potential,VFieldx_e,VFieldy_e,VFieldx_h,VFieldy_h,...
    x,y,Step,Bulk,Radius,TauBe,TauSe,TauBh,TauSh,ItFigIn);
WorkTransportTotal = WorkTransportTotal_;
fprintf('Generated 1 "Work-Transport" matrix\n\n');

for i = 1:NAverage-1
    [WorkTransportTotal_, x, y] =...
        WorkTransport(potential,VFieldx_e,VFieldy_e,VFieldx_h,VFieldy_h,...
        x,y,Step,Bulk,Radius,TauBe,TauSe,TauBh,TauSh,ItFigIn);
    WorkTransportTotal = WorkTransportTotal + WorkTransportTotal_;
    fprintf('Generated %d "Work-Transport" matrices\n\n',i+1);
end

WorkTransportTotal = WorkTransportTotal ./ NAverage;
fprintf('@@@ Average over all "Work-Transport" matrices @@@\n');


%%%%%%%%%
% Plots %
%%%%%%%%%
[xx,yy] = meshgrid(x,y);

figure(ItFigOut);
colormap jet;
surf(xx,yy,WorkTransportTotal,'EdgeColor','none');
title('Total <Work-Transport>');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Work / q [#charges * V]');

ItFigOut = ItFigOut + 1;
fprintf('CPU time --> %.2f[min]\n\n',(cputime-TStart)/60);
end