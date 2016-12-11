%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate the analytical efficiency to check if it becomes negative %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [isPOS] = EvalEffFunc3D(CF,showPlot)

isPOS  = true;
nSteps = 200;

startX = -1;
endX   = 1;
stepX  = (endX - startX) / nSteps;

startY = -1;
endY   = 1;
stepY  = (endY - startY) / nSteps;

startZ = -pi;
endZ   = pi;
stepZ  = (endZ - startZ) / nSteps;

x = startX:stepX:endX;
y = startY:stepY:endY;
z = startZ:stepZ:endZ;

T   = zeros(length(y),length(x),length(z));
Txy = zeros(length(y),length(x));
Txz = zeros(length(z),length(x));
Tyz = zeros(length(z),length(y));

for i = 1:length(x)
    for j = 1:length(y)
        for k = 1:length(z)
            T(j,i,k) = (CF(1)+CF(2)*x(i)+CF(3)*x(i)^2+CF(4)*x(i)^3) +...
                (CF(5)+CF(6)*x(i)+CF(7)*x(i)^2+CF(8)*x(i)^3)*y(j) +...
                (CF(9)+CF(10)*x(i)+CF(11)*x(i)^2+CF(12)*x(i)^3)*y(j)^2 +...
                (CF(13)+CF(14)*x(i)+CF(15)*x(i)^2+CF(16)*x(i)^3)*y(j)^3 +...
                (CF(17)+CF(18)*x(i)+CF(19)*x(i)^2+CF(20)*x(i)^3)*y(j)^4 +...
                (CF(21)+CF(22)*x(i)+CF(23)*x(i)^2+CF(24)*x(i)^3)*y(j)^5 +...
                (CF(25)+CF(26)*x(i)^2+CF(27)*y(j)+CF(28)*y(j)^2)*z(k)^2;
            if (T(j,i,k) < 0)
                isPOS = false;
                return
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%
% Make projections %
%%%%%%%%%%%%%%%%%%%%
for i = 1:length(x)
    for j = 1:length(y)
        sum = 0;
        for k = 1:length(z)
            sum = sum + T(j,i,k);
        end
        
        Txy(j,i) = sum;
    end
end

for i = 1:length(x)
    for k = 1:length(z)
        sum = 0;
        for j = 1:length(y)
            sum = sum + T(j,i,k);
        end
        
        Txz(k,i) = sum;
    end
end

for j = 1:length(y)
    for k = 1:length(z)
        sum = 0;
        for i = 1:length(x)
            sum = sum + T(j,i,k);
        end
        
        Tyz(k,j) = sum;
    end
end


%%%%%%%%%%%%%%%%
% Plot options %
%%%%%%%%%%%%%%%%
if showPlot == true
   figure;
   surf(x, y, Txy, 'EdgeColor', 'none');
   shading interp;
   colormap jet;
   light;
   lighting phong;
   colorbar;
   title('Efficiency');
   xlabel('cos(\theta_K)');
   ylabel('cos(\theta_l)');

   figure;
   surf(x, z, Txz, 'EdgeColor', 'none');
   shading interp;
   colormap jet;
   light;
   lighting phong;
   colorbar;
   title('Efficiency');
   xlabel('cos(\theta_K)');
   ylabel('\phi');
   
   figure;
   surf(y, z, Tyz, 'EdgeColor', 'none');
   shading interp;
   colormap jet;
   light;
   lighting phong;
   colorbar;
   title('Efficiency');
   xlabel('cos(\theta_l)');
   ylabel('\phi');
end
end