%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the program to interpolate the new efficiency functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%
% Global variables %
%%%%%%%%%%%%%%%%%%%%
nFiles    = 100; % Number of files
binVector = [0,1,2,3,5,7,8]; % q^2 bin index
fileName1 = '../../efficiency/EffRndGenBinFilesSign/Efficiency_RndGen_';
fileName2 = '_H2Deff_q2Bin';
showPlot  = false;


for i = 1:nFiles
    for j = 1:length(binVector)
        fileName = sprintf('%s%d%s',fileName1,i-1,fileName2);

        InterpEff2D(fileName,binVector(j),showPlot);
    end
end

fprintf('\n@@@ I''ve interpolated %d efficiencies @@@\n',nFiles);
