%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the program to generate the new efficiency functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%
% Global variables %
%%%%%%%%%%%%%%%%%%%%
nFiles = 200; % Number of files to generate
nBins  = 8;   % Number of q^2 bins
startBin = 1;
NcoeffThetaL = 5;
NcoeffThetaK = 4;


for i = 1:nFiles
    %%%%%%%%%%%%%%
    % Parameters %
    %%%%%%%%%%%%%%
    fidINval    = fopen('../../../Efficiency/ThetaKThetaL_B0ToKstMuMu_B0ToJPsiKst_B0ToPsi2SKst.txt','r');
    fidINcov    = fopen('../../../Efficiency/ThetaKThetaLFullCovariance_B0ToKstMuMu_B0ToJPsiKst_B0ToPsi2SKst.txt','r');
    
    fileNameOut = sprintf('../../../Efficiency/EffRndGenAnalyFilesSign_JPsi_Psi2S/Efficiency_RndGen_%d.txt',i-1);
    fidOUT = fopen(fileNameOut,'w+');

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Skip bins if you want to %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:startBin-1
        [q2Bin,meanV,errV,CovM,meanVOrig] = ReadCovMatrix(fidINval,...
            fidINcov,NcoeffThetaL,NcoeffThetaK);
    end

    
    fprintf('Generating file %d \n',i);
    for j = startBin:startBin+nBins-1
        [q2Bin,meanV,errV,CovM,meanVOrig] = ReadCovMatrix(fidINval,...
            fidINcov,NcoeffThetaL,NcoeffThetaK);
        
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check if the efficiency is negative %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (EvalEffFunc(meanVOrig) == 0)
            fprintf('@@@ The origianl efficiency is negative @@@\n');
        else
            fprintf('@@@ The origianl efficiency is OK @@@\n');
        end

        
        isPOS = 0;
        while (isPOS == 0)            
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
            isPOS = EvalEffFunc(strechNewV');
            if (isPOS == 0)
                fprintf('--> The new efficiency is negative\n');
            else
                fprintf('--> The new efficiency is OK\n');
            end
        end
        
 
        SaveMVNvecIntoFile(fidOUT,q2Bin,strechNewV,errV,...
            NcoeffThetaL,NcoeffThetaK);
    end

    
fclose(fidINval);
fclose(fidINcov);
fclose(fidOUT);
end


fprintf('@@@ I''ve gnerated %d efficiencies @@@\n',nFiles);