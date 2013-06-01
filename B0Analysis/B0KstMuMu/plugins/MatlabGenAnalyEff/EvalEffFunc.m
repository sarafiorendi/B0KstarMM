%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate the analytical efficiency to check if it becomes negative %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [isPOS] = EvalEffFunc(CF)

isPOS = 1;

nSteps = 200;

startX = -1;
endX = 1;
stepX = (endX - startX) / nSteps;

startY = -1;
endY = 1;
stepY = (endY - startY) / nSteps;

x = startX:stepX:endX;
y = startY:stepY:endY;
Z = zeros(length(y),length(x));

for i = 1:length(x)
    for j = 1:length(y)
        Z(j,i) = (CF(1)+CF(2)*x(i)+CF(3)*x(i)^2+CF(4)*x(i)^3) +...
            (CF(5)+CF(6)*x(i)+CF(7)*x(i)^2+CF(8)*x(i)^3)*y(j)^2 +...
            (CF(9)+CF(10)*x(i)+CF(11)*x(i)^2+CF(12)*x(i)^3)*y(j)^3 +...
            (CF(13)+CF(14)*x(i)+CF(15)*x(i)^2+CF(16)*x(i)^3)*y(j)^4 +...
            (CF(17)+CF(18)*x(i)+CF(19)*x(i)^2+CF(20)*x(i)^3)*y(j)^6;
        
        if (Z(j,i) < 0)
            isPOS = 0;
            return
        end
    end
end


%%%%%%%%%%%%%%%%
% Plot options %
%%%%%%%%%%%%%%%%
surf(x, y, Z, 'EdgeColor', 'none');
shading interp;
colormap jet;
light;
lighting phong;
colorbar;
title('Efficiency');
xlabel('cos(\theta_k)');
ylabel('cos(\theta_l)');

end