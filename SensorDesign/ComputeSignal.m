%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that calculates the signal given the "Work-Transport" matrix %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WorkTransportTotal = Total "Work-Transport" matrix
% x, y          = Axes
% Depth         = Penetration depth along y (positive = from back side,
%                 negative = from strip side) [um]
% XEnteranceParticle = X cordinate of the entrance of the particle [um]
% (back-plane side)
% XExitParticle = X cordinate of the exit of the particle [um] (strip side)
% Step          = Unit step of the lattice on which the field is computed [um]
% subStep       = Unit step of the lattice on which the field is interpolated [um]
% Bulk          = Bulk thickness [um]
% Radius        = Unit step of the movements [um]
% ChargeDensity = Charge density [electrons/um]

function [Charge] = ComputeSignal(WorkTransportTotal,x,y,Depth,...
    XEnteranceParticle,XExitParticle,Step,subStep,Bulk,Radius,ChargeDensity)


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps = Radius/10; % Minimum step to avoid infinities [um]
Charge = 0; % Starting value of the signal [electrons]


%%%%%%%%%%%%%%%%%%%
% Start algorithm %
%%%%%%%%%%%%%%%%%%%
mx = find(x >= XEnteranceParticle);
mx = x(mx(1));
if Depth >= 0
    my = y(1);
else
    my = find(y >= (Bulk + Depth));
    my = y(my(1));
end

while my >= 0 && ((Depth >= 0 && my <= Depth) || (Depth < 0 && my <= Bulk))

    x_near = find(abs(x - mx) <= subStep);
    y_near = find(abs(y - my) <= subStep);

    % Calculate the average Work within a radius = subStep
    W      = 0;
    weight = 0;
    for ii = 1:length(x_near)
        for jj = 1:length(y_near)

            dist = sqrt((mx-x(x_near(ii)))^2 + (my-y(y_near(jj)))^2);

            if dist <= subStep && ~isnan(WorkTransportTotal(y_near(jj),x_near(ii)))
                
                if dist < eps
                    dist = eps;
                end
                
                W = W + WorkTransportTotal(y_near(jj),x_near(ii)) / Step*subStep / dist;
                weight = weight + 1 / dist;
            end
        end
    end
    if weight ~= 0
        W = W / weight;
    end
    
    Charge = Charge + W;
%    fprintf('Signal: %d\tx: %d\ty: %d\n',S,mx,my);

    % Calculate the next movement
    dx   = XExitParticle - XEnteranceParticle;
    sinv = dx   / sqrt(dx^2 + Bulk^2);
    cosv = Bulk / sqrt(dx^2 + Bulk^2);
    mx   = mx + Radius*sinv;
    my   = my + Radius*cosv;
end

Charge = ChargeDensity * Charge * Radius / subStep;
%fprintf('Total signal: %d [e-]\n',S);
end
