%% 
% retreives metastatic efficiency from cell type
% Input: cell type
% Output: high / low
function metEffStr = LCH_getMetEff(cellType)
if sum(strcmpi(cellType,{'m481','m214','um12','m514','m405'}))
    metEffStr = 'High';
else 
    if sum(strcmpi(cellType,{'m530','m610','m528','m498'}))
        metEffStr = 'Low';
    else       
        metEffStr = 'Unknown';
    end
end
end


	




