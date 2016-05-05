%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that generates several "Work-Transport" %
% matrices and computes the average                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAverage = Generate NAverage WorkTransportTotal and average them

function [WorkTransportTotal, x, y] = ...
    ManyWorkTransport(uw,pw,tw,VFieldx_e,VFieldy_e,VFieldx_h,VFieldy_h,...
    x,y,Step,Bulk,Radius,TauBe,TauSe,TauBh,TauSh,NAverage)
TStart = cputime; % CPU time at start


%%%%%%%%%%%%%%%%%%%
% Start algorithm %
%%%%%%%%%%%%%%%%%%%
[WorkTransportTotal_, x, y] =...
    WorkTransport(uw,pw,tw,VFieldx_e,VFieldy_e,...
    VFieldx_h,VFieldy_h,x,y,Step,Bulk,Radius,TauBe,TauSe,TauBh,TauSh);
WorkTransportTotal = WorkTransportTotal_;
fprintf('Generated 1 "Work-Transport" matrix\n\n');

for i = 1:NAverage-1
    [WorkTransportTotal_, x, y] =...
        WorkTransport(uw,pw,tw,VFieldx_e,VFieldy_e,...
        VFieldx_h,VFieldy_h,x,y,Step,Bulk,Radius,TauBe,TauSe,TauBh,TauSh);
    WorkTransportTotal = WorkTransportTotal + WorkTransportTotal_;
    fprintf('Generated %d "Work-Transport" matrices\n\n',i+1);
end

WorkTransportTotal = WorkTransportTotal ./ NAverage;
fprintf('@@@ Average over all "Work-Transport" matrices @@@\n');


%%%%%%%%%
% Plots %
%%%%%%%%%
[xx,yy] = meshgrid(x,y);

figure (6);
colormap jet;
surf(xx,yy,WorkTransportTotal,'EdgeColor','none');
title('Total <Work-Transport>');
xlabel('X [\mum]');
ylabel('Y [\mum]');
zlabel('Work / q [#charges * V]');

fprintf('CPU time --> %d[min]\n\n',(cputime-TStart)/60);
end