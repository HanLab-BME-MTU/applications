%% 
% Patchie solution for the different data structures names used to hold
% features
function cellData = getCellData(featsInFname)
load(featsInFname); % 1 x nCells

if exist('dLbpWell','var')
    cellData = dLbpWell;
end

end