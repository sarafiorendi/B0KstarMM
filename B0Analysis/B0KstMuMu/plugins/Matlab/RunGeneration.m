%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the program to generate the new efficiency functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ==> How to run the program <==           %
% (1) Set the global variables             %
% (2) Run on good-tagged events            %
% (3) Run on mis-tagged events             %
% Output files must remain the same since  %
% good and mis-tagged events should go one %
% after the other in the same file         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%
% Global variables %
%%%%%%%%%%%%%%%%%%%%
nFiles   = 100; % Number of files to generate
nBins    =   9; % Number of q^2 bins
startBin =   1; % Start bin [1...nBins]

NcoeffThetaL = 6;
NcoeffThetaK = 4;
NcoeffPhi    = 0; % If [0] then 2D function, else [4] 3D function

showPlot = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set sandom number seed %
%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(1);


for i = 1:nFiles
    %%%%%%%%%%%%%%
    % Parameters %
    %%%%%%%%%%%%%%
    fidINval    = fopen('../../efficiency/ThetaKThetaLPhi_B0ToKstMuMu_B0ToJPsiKst_B0ToPsi2SKst.txt','r');
    fidINcov    = fopen('../../efficiency/ThetaKThetaLPhiFullCovariance_B0ToKstMuMu_B0ToJPsiKst_B0ToPsi2SKst.txt','r');

    fileNameOut = sprintf('../../efficiency/EffRndGenAnalyFilesSign_JPsi_Psi2S/Efficiency_RndGen_%d.txt',i-1);
    fidOUT      = fopen(fileNameOut,'a+');
    fprintf(fidOUT,'%d\n',nBins*NcoeffThetaL+NcoeffPhi);
    fprintf('\nGenerating file n.%d\n',i);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Skip bins if you want to %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:startBin-1
        fprintf('Skipping bin n.%d\n',j);
        [q2Bin,meanV,errV,CovM,meanVOrig] = ReadCovMatrix(fidINval,...
            fidINcov,NcoeffThetaL,NcoeffThetaK,NcoeffPhi);
    end

    
    for j = startBin:nBins
        fprintf('Making efficiency for bin n.%d\n',j);

        [q2Bin,meanV,errV,CovM,meanVOrig] = ReadCovMatrix(fidINval,...
            fidINcov,NcoeffThetaL,NcoeffThetaK,NcoeffPhi);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check if the efficiency is negative %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ((NcoeffPhi == 0 && EvalEffFunc2D(meanVOrig,showPlot) == false) ||...
            (NcoeffPhi ~= 0 && EvalEffFunc3D(meanVOrig,showPlot) == false))            
            fprintf('@@@ The origianl efficiency is negative @@@\n');
        else
            fprintf('@@@ The origianl efficiency is OK @@@\n');
        end

        
        isPOS = false;
        while (isPOS == false)            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Make multvariate vector %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (size(meanV) ~= 0)
                newV = mvnrnd(meanV,CovM);
            else
                newV = [];
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Add zeros to the new mutivariate vector %
            % in order to fit correct size            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            strechNewV = zeros(size(errV));
            indx = 1;
            for k = 1:length(errV)
                if (errV(k) ~= 0)
                    strechNewV(k) = newV(indx);
                    indx = indx + 1;
                else
                    strechNewV(k) = 0;
                end
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Check if the efficiency is negative %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if NcoeffPhi == 0
                isPOS = EvalEffFunc2D(strechNewV',showPlot);
            else
                isPOS = EvalEffFunc3D(strechNewV',showPlot);
            end
                
            if (isPOS == false)
                fprintf('--> The new efficiency is negative\n');
            else
                fprintf('--> The new efficiency is OK\n');
            end
        end

        
        SaveMVNvecIntoFile(fidOUT,q2Bin,strechNewV,errV,...
            NcoeffThetaL,NcoeffThetaK,NcoeffPhi);
    end
  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save dummy data for mis-tagged events %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(fidOUT,'\n\n');
    fprintf(fidOUT,'%d\n',nBins*NcoeffThetaL+NcoeffPhi);        

    strechNewV = zeros(size(strechNewV));
    errV       = zeros(size(errV));
    
    for j = startBin:nBins        
        SaveMVNvecIntoFile(fidOUT,q2Bin,strechNewV,errV,...
            NcoeffThetaL,NcoeffThetaK,NcoeffPhi);
    end
    
    
fclose(fidINval);
fclose(fidINcov);
fclose(fidOUT);
end


fprintf('@@@ I''ve generated %d efficiencies @@@\n',nFiles);
