%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save new mutivariate normal vector into file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveMVNvecIntoFile(fid,q2Bin,meanV,errV,NcoeffThetaL,NcoeffThetaK)

for i = 1:NcoeffThetaL
    myS = sprintf('%f',q2Bin);
    
    for j = 1:NcoeffThetaK
        if (errV(j+(i-1)*NcoeffThetaK) ~= 0)
            myS = strcat(myS,sprintf('   %e   0',...
                meanV(j+(i-1)*NcoeffThetaK)));
        else
            myS = strcat(myS,'   0   0');
        end
    end
    
    fprintf(fid,'%s\n',myS);
end
