%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate new efficiency functions from the mutivariate normal %
% distribution of the parameters                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q2Bin, meanVOut, errVOut, CovMOut, meanVOrig] =...
    ReadCovMatrix(fidVal,fidCov,NcoeffThetaL,NcoeffThetaK)

nRowVal = NcoeffThetaL;
nColVal = NcoeffThetaK*2+1;
nRowCov = NcoeffThetaL*NcoeffThetaK;
nColCov = NcoeffThetaL*NcoeffThetaK+1;

% Read central values from file
myStr = '%f';
for i = 1:nColVal-1
    myStr = [myStr ' %f'];
end
tmpMval = fscanf(fidVal,myStr,[nColVal nRowVal]);

% Read covariance matrix from file
myStr = '%f';
for i = 1:nColCov-1
    myStr = [myStr ' %f'];
end
tmpMcov = fscanf(fidCov,myStr,[nColCov nRowCov]);

CovM = tmpMcov';
ValM = tmpMval';

q2Bin = CovM(1,1);

ValM(:, 1) = [];
meanV = ValM(1, 1:2:end);
for i = 2:nRowVal
    meanV = [meanV ValM(i, 1:2:end)];
end
meanVOrig = meanV;
meanV = meanV';

CovM(:, 1) = [];
errV = diag(CovM);

CovMOut = zeros(size(CovM));
meanVOut = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shifting rows and columns to jump empty rows and columns %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rowIndx = 1;
for j = 1:length(CovM)
    colIndx = 1;
    for k = 1:length(CovM)
        if (errV(j) ~= 0) && (errV(k) ~= 0)
            CovMOut(rowIndx, colIndx) = CovM(j,k);
    		colIndx = colIndx + 1;
        end
    end
    
    if (colIndx ~= 1)
        meanVOut(rowIndx) = meanV(j);
        rowIndx = rowIndx +1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove empty last columns and rows %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CovMOut(:, rowIndx:length(CovM)) = [];
CovMOut(rowIndx:length(CovM), :) = [];

errVOut = errV;

end
