%%
% data: n x nFeats
% masterLabels, masterInds - partition criteria
% slaveLabels - labels withint each partition (figure)
function [] = LCH_visualize2D(data,masterUniqueLabels,masterInds,slaveLabels,masterStr,baseDname,suffix,always)

outDname = [baseDname filesep masterStr];

if ~exist(outDname,'dir')
    mkdir(outDname);
end

if isempty(masterUniqueLabels)
    LCH_visualize2D_mapped(data,slaveLabels,[outDname filesep masterStr '_' suffix '.jpg']);
else
    for iMasterLabel = 1 : length(masterUniqueLabels)
        curInds = masterInds{iMasterLabel};
        curData = data(curInds,:);
        curMasterLabel = masterUniqueLabels{iMasterLabel};        
        curSlaveLabels = slaveLabels(curInds);
        
        %         curOutDname = [outDname filesep curMasterLabel];
        %         if ~exist(curOutDname,'dir')
        %             mkdir(curOutDname);
        %         end
        
        outFname = [outDname filesep masterStr '_' curMasterLabel '_' suffix '.jpg'];
        
        if ~exist(outFname,'file') || always
            LCH_visualize2D_mapped(curData,curSlaveLabels,outFname);
        end
        close all;
    end
end

end