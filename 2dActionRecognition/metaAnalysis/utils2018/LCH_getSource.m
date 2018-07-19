%% 
% retreives source (tumor, cell line, melanocytes) from cell type
% Input: cell type
% Output: source
function sourceStr = LCH_getSource(cellType)
if sum(strcmpi(cellType,{'m481','m214','m530','m610','um12','m514','m405','m528','m498','ut8','m634','m597','m405c3','m405c2'}))
    sourceStr = 'Tumors';
else
    if sum(strcmpi(cellType,{'atcc','m116'}))
        sourceStr = 'Melanocytes';
    else
        if sum(strcmpi(cellType,{'a375','mv3','wm3670','wm1361','wm1366','skmel2'}))
            sourceStr = 'CellLines';
        else
            warning('Failed finding source %s', cellType);
        end
    end
end
end


	




