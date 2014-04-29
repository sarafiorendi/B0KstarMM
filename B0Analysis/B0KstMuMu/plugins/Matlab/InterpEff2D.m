%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to read and interpolate 2D efficiency functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;


%%%%%%%%%%%%%%%%%%%%
% Global constants %
%%%%%%%%%%%%%%%%%%%%
fileName    = '../../efficiency/H2Deff_MisTag_q2Bin';
NbinsX      = 5;
NbinsY      = 5;
NstepsX     = 120;  % [120]
NstepsY     = 120;  % [120]
effAtBound  = 1e-5; % [1e-5]
effMinValue = 2e-5; % [2e-5]


%%%%%%%%%%%%%%%%%%%%
% Global variables %
%%%%%%%%%%%%%%%%%%%%
q2Indx = 0; % q^2 bin index
% - q^2 bin 0-2:     interp2
% - q^2 bin 3-5 7-8: scatteredInterpolant natural+linear
% - q^2 bin 6:       scatteredInterpolant natural+nearest
useMethodInterp2 = true; % If true 'interp2' else 'scatteredInterpolant'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ==> Interpolation methods <==                          %
% - for 'interp2': linear, nearest, cubic, spline        %
% - for 'scatteredInterpolant': linear, nearest, natural %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%
% Prepare variables %
%%%%%%%%%%%%%%%%%%%%%
if useMethodInterp2 == true
    interpMethod = 'spline'; % [spline]
    addContraints = true;
else
    interpMethod = 'natural'; % [natural]
    addContraints = false;
end

if addContraints == true
    deltaBins = 2;
else
    deltaBins = 0;
end


%%%%%%%%%%%%%%%%%%%
% Prepare vectors %
%%%%%%%%%%%%%%%%%%%
Xvec   = zeros(1,NbinsX);
Yvec   = zeros(1,NbinsY);
Values = zeros(NbinsY+deltaBins,NbinsX+deltaBins);

XvecAtBound   = zeros(1,NbinsX+2);
YvecAtBound   = zeros(1,NbinsY+2);
ValuesAtBound = zeros(NbinsY+2,NbinsX+2);


%%%%%%%%%%%%%%%%%%%%%%%%%
% Read values from file %
%%%%%%%%%%%%%%%%%%%%%%%%%
fileIDin = fopen(strcat(fileName,sprintf('_%d.txt',q2Indx)),'r');
formatData = '%f   %f   %f   %f   %f   %f';
i = 1;
j = 1;
fprintf('@@@ Reading 2D binned efficiency from txt file %s @@@\n',fileName);
while i <= NbinsX
    row = textscan(fileIDin,formatData,1);

    Xvec(i) = row{1,1} + row{1,2}/2;
    Yvec(j) = row{1,3} + row{1,4}/2;

    XvecAtBound(i+1) = row{1,1} + row{1,2}/2;
    YvecAtBound(j+1) = row{1,3} + row{1,4}/2;

    Values(j+deltaBins/2,i+deltaBins/2) = row{1,5};
    
    j = j + 1;    
    if (j > NbinsY)
        j = 1;
        i = i + 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if original efficiency goes negative %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isPOS = true;
for i = 1:NbinsX
    for j = 1:NbinsY
        if Values(j+deltaBins/2,i+deltaBins/2) < 0
            isPOS = false;
            Values(j+deltaBins/2,i+deltaBins/2) = effMinValue;
        end
    end
end
if isPOS == false
   fprintf('@@@ The efficiency before interpolation is negative --> corrected to %e @@@\n',effMinValue);
else
    fprintf('@@@ The efficiency before interpolation is OK @@@\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add contraints at boundaries %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XvecAtBound(1)        = -1;
XvecAtBound(NbinsX+2) =  1;

YvecAtBound(1)        = -1;
YvecAtBound(NbinsY+2) =  1;

XvecAtBound(2)        = 2*XvecAtBound(2) - XvecAtBound(1);
XvecAtBound(NbinsX+1) = 2*XvecAtBound(NbinsX+1) - XvecAtBound(NbinsX+2);

YvecAtBound(2)        = 2*YvecAtBound(2) - YvecAtBound(1);
YvecAtBound(NbinsX+1) = 2*YvecAtBound(NbinsY+1) - YvecAtBound(NbinsY+2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraint low X bound to nearest values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:NbinsX
    ValuesAtBound(1,i+1) = Values(1+deltaBins/2,i+deltaBins/2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraint high X bound to nearest values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:NbinsX
    ValuesAtBound(NbinsY+2,i+1) = Values(NbinsY+deltaBins/2,i+deltaBins/2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraint low Y bound to nearest values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:NbinsY
    ValuesAtBound(j+1,1) = Values(j+deltaBins/2,1+deltaBins/2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraint high Y bound to nearest values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:NbinsY
    ValuesAtBound(j+1,NbinsX+2) = Values(j+deltaBins/2,NbinsX+deltaBins/2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraint corners to nearest values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ValuesAtBound(1,1)               = Values(1+deltaBins/2,1+deltaBins/2);
ValuesAtBound(1,NbinsX+2)        = Values(1+deltaBins/2,NbinsX+deltaBins/2);
ValuesAtBound(NbinsY+2,1)        = Values(NbinsY+deltaBins/2,1+deltaBins/2);
ValuesAtBound(NbinsY+2,NbinsX+2) = Values(NbinsY+deltaBins/2,NbinsX+deltaBins/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply constraints at boundaries %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if addContraints == true
    for i = 1:NbinsX+2
        Values(1,i)        = effAtBound;
        Values(NbinsY+2,i) = effAtBound;
    end
    
    for j = 1:NbinsY+2
        Values(j,1)        = effAtBound;
        Values(j,NbinsX+2) = effAtBound;
    end
    
    actualX = [-1 Xvec 1];
    actualY = [-1 Yvec 1];
else
    actualX = Xvec;
    actualY = Yvec;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot original efficiency %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
surf(actualX, actualY, Values);
colorbar;
title('Efficiency');
xlabel('cos(\theta_K)');
ylabel('cos(\theta_l)');


%%%%%%%%%%%%%%%%%%%%%%%
% Create interpolants %
%%%%%%%%%%%%%%%%%%%%%%%
[newXvec,newYvec] = meshgrid(linspace(XvecAtBound(1),XvecAtBound(NbinsX+2),NstepsX),...
    linspace(YvecAtBound(1),YvecAtBound(NbinsY+2),NstepsY));

newValuesAtBound = interp2(XvecAtBound, YvecAtBound, ValuesAtBound, newXvec, newYvec, 'nearest', 0.0);

if useMethodInterp2 == true
    newValues = interp2(actualX, actualY, Values, newXvec, newYvec, interpMethod, 0.0);
else
    [newXvec_,newYvec_] = meshgrid(actualX, actualY);

    if q2Indx ~= 6
        newValues_ = scatteredInterpolant(newXvec_(:), newYvec_(:), Values(:), interpMethod, 'linear');
    else
        newValues_ = scatteredInterpolant(newXvec_(:), newYvec_(:), Values(:), interpMethod, 'nearest');
    end
    
    newValues = newValues_(newXvec, newYvec);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cast boundaries to zero %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = 1:length(newXvec)
%     for j = 1:length(newYvec)
%         if newXvec(j,i) < Xvec(1+deltaBins/2) || newXvec(j,i) > Xvec(NbinsX+deltaBins/2) ||...
%                 newYvec(j,i) < Yvec(1+deltaBins/2) || newYvec(j,i) > Yvec(NbinsY+deltaBins/2)
%             newValues(j,i) = 0.0;
%         end
%     end
% end


%%%%%%%%%%%%%%%%%%%
% Plot boundaries %
%%%%%%%%%%%%%%%%%%%
figure;
surf(XvecAtBound, YvecAtBound, ValuesAtBound);
colorbar;
title('Efficiency at Boundaries');
xlabel('cos(\theta_K)');
ylabel('cos(\theta_l)');

figure;
surf(newXvec, newYvec, newValuesAtBound);
colorbar;
title('Interpolated Efficiency at Boundaries');
xlabel('cos(\theta_K)');
ylabel('cos(\theta_l)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot new efficiency without boundaries %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
surf(newXvec, newYvec, newValues);
colorbar;
title('Interpolated Efficiency');
xlabel('cos(\theta_K)');
ylabel('cos(\theta_l)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Merge new efficiency and boundaries %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% newValues = newValues + newValuesAtBound;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if new efficiency goes negative %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isPOS = true;
for i = 1:length(newXvec)
    for j = 1:length(newYvec)
        if newValues(j,i) < 0
            isPOS = false;
            newValues(j,i) = effMinValue;
        end
    end
end
if isPOS == false
   fprintf('@@@ The efficiency after interpolation is negative --> corrected to %e @@@\n',effMinValue);
else
    fprintf('@@@ The efficiency after interpolation is OK @@@\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot new efficiency with boundaries %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
surf(newXvec, newYvec, newValues);
colorbar;
title('Interpolated Efficiency with Boundaries');
xlabel('cos(\theta_K)');
ylabel('cos(\theta_l)');


%%%%%%%%%%%%%%%%%%%%%%%
% Save values to file %
%%%%%%%%%%%%%%%%%%%%%%%
fileIDout = fopen(strcat(fileName,sprintf('_interp_%d.txt',q2Indx)),'w');
i = 1;
j = 1;
fprintf('@@@ Saving 2D binned efficiency to txt file @@@\n');
while i <= NstepsX-1
    formatData = sprintf(' ');
    formatData = strcat(formatData,sprintf('%f   %f',newXvec(j,i),newXvec(j,i+1) - newXvec(j,i)));
    formatData = strcat(formatData,sprintf('   %f   %f',newYvec(j,i),newYvec(j+1,i) - newYvec(j,i)));
    formatData = strcat(formatData,sprintf('   %f   %f',newValues(j,i),0));

    fprintf(fileIDout,'%s\n',formatData);
    
    j = j + 1;    
    if (j > NstepsY-1)
        j = 1;
        i = i + 1;
    end
    
    clear newRow;
end


fprintf('@@@ I''ve generated interpolated efficiency @@@\n');
fclose(fileIDout);
fclose(fileIDin);