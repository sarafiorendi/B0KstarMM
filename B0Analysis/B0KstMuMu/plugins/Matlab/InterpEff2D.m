%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to read and interpolate 2D efficiency functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID = fopen('test.txt');
formatSpec = '%f   %f   %f   %f   %f   %f';

while ~feof(fid)
    row = textscan(fileID,formatSpec,1);
    Xvec[i] = row[1];
    Yvec[i] = row[2];
end

fclose(fileID);