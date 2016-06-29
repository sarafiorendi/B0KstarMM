%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that calculates the "Work-Transport" matrix for the carriers %
% given the velocity field and the weighting potential                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [WorkMatrix] =...
    WorkTransportMatrix(Fieldx,Fieldy,WeightingEx,WeightingEy,x,y,...
    TauB,TauS,Step,Bulk,Radius,Charge)


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable initialization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps = Radius/10; % Minimum step to avoid infinities [um]


%%%%%%%%%%%%%%%%%%%
% Start algorithm %
%%%%%%%%%%%%%%%%%%%
WorkMatrix = zeros(length(y), length(x));

for i = 1:length(x)
    for j = 1:length(y)
        
        mx     = x(i);
        my     = y(j);
        mx_old = sqrt(-1);
        my_old = sqrt(-1);
        Path   = 0; % Path length
        LocalEffectiveV = 0; % Starting value of the effective potential
        ContinuePath = 1; % Boolean variable result of the Markov process

        while my >= 0 && my <= Bulk && (mx ~= mx_old || my ~= my_old) &&...
                ContinuePath == 1

            x_near = find(abs(x - mx) <= Step);
            y_near = find(abs(y - my) <= Step);        

            % Calculate the average Effective-Field within a radius = Step
            EFx    = 0;
            EFy    = 0;
            weight = 0;
            for ii = 1:length(x_near)
                for jj = 1:length(y_near)
                    
                    dist = sqrt((mx-x(x_near(ii)))^2 + (my-y(y_near(jj)))^2);
                    
                    if dist <= Step &&...
                            ~isnan(WeightingEx(y_near(jj),x_near(ii)))&&...
                            ~isnan(WeightingEy(y_near(jj),x_near(ii)))
                        
                        if dist < eps
                            dist = eps;
                        end
                        
                        EFx = EFx + WeightingEx(y_near(jj),x_near(ii)) / dist;
                        EFy = EFy + WeightingEy(y_near(jj),x_near(ii)) / dist;
                        weight = weight + 1 / dist;
                    end
                end
            end
            if weight ~= 0
                EFx = EFx / weight;
                EFy = EFy / weight;
            end
         
            % Calculate the average Velocity-Field within a radius = Step
            VFx    = 0;
            VFy    = 0;
            weight = 0;
            for ii = 1:length(x_near)
                for jj = 1:length(y_near)
                    
                    dist = sqrt((mx-x(x_near(ii)))^2 + (my-y(y_near(jj)))^2);
                    
                    if dist <= Step &&...
                            ~isnan(Fieldx(y_near(jj),x_near(ii))) &&...
                            ~isnan(Fieldy(y_near(jj),x_near(ii)))
                        
                        if dist < eps
                            dist = eps;
                        end
                        
                        VFx = VFx + Fieldx(y_near(jj),x_near(ii)) / dist;
                        VFy = VFy + Fieldy(y_near(jj),x_near(ii)) / dist;
                        weight = weight + 1 / dist;
                    end
                end
            end
            if weight ~= 0
                VFx = VFx / weight;
                VFy = VFy / weight;
            end
            
            % Calculate the next movement
            if VFy ~= 0 || VFx ~= 0
                sinv = VFx / sqrt(VFx^2 + VFy^2);
                cosv = VFy / sqrt(VFx^2 + VFy^2);
            else
                sinv = 0;
                cosv = 0;                
            end
            mx_old = mx;
            my_old = my;
            mx     = mx + Radius*sinv;
            my     = my + Radius*cosv;

            % Markov chain for the propagation of the charge
            % Assumption: linear interpolation of carriers life-time
            Coin = rand(1,1);
            CCD  = sqrt(VFx^2+VFy^2)*(TauB+y(j)*(TauS-TauB)/Bulk);
            P    = exp(-Radius / CCD);
            if Coin < P
                ContinuePath = 1;
                Path = Path + Radius;
            else
                ContinuePath = 0;
            end
            
            % Calculate the local potential difference
            if VFx ~= 0 || VFy ~= 0
                LocalEffectiveV = LocalEffectiveV +...
                (EFx*VFx + EFy*VFy)/sqrt(VFx^2 + VFy^2) * Radius;
            end
            
%fprintf('New point --> x: %d\ty: %d\tPath: %d\tVolt: %d\n',...
%    mx,my,Path,LocalEffectiveV);
        end
        WorkMatrix(j,i) = Charge * Radius * LocalEffectiveV;
    end
end
end
