function centerDataRevert(obj,cfg)
% function centerDataRevert(obj,cfg)
% SYNOPSIS:
% Reverts the data to its initial position.
%
% REQUIRED INPUTS:
% - cfg
% The configuration data corresponding to the data in obj.data
%
% OPTIONAL INPUTS:
%
% NEEDED PROPERTIES:
% - obj.data.points
% - obj.data.rawPoints
% - obj.data.modelBezCP
% - obj.data.simModelBezCP
%
% MODIFIED PROPERTIES:
% - obj.data.points
% - obj.data.rawPoints
% - obj.data.modelBezCP
% - obj.data.simModelBezCP
%
% OUTPUTS:
%
% Pascal Bérard, October 2011

% Fix the current STORM path
pos = strfind(cfg.path,'_data');
path = [getStormPath() strrep(strrep(cfg.path(pos:end),'\',filesep),'/',filesep)];

% Read and crop the data to find the amount the data has been initially shifted
tmpData = Data.read([path cfg.fileName]);
pro = Processor(tmpData);
pro.cropData(cfg.roiPosition,cfg.roiSize);

% Compute the center of the initial point cloud
center = (max(tmpData.points,[],1)+min(tmpData.points,[],1))/2;

% Translate all the data of the real data set back to the initial position
obj.data.points = obj.data.points+repmat(center,size(obj.data.points,1),1);

if ~isempty(obj.data.rawPoints)
    obj.data.rawPoints = obj.data.rawPoints+repmat(center,size(obj.data.rawPoints,1),1); 
end

if ~isempty(obj.data.modelBezCP);
    nModelBezCP = cellfun(@(a) size(a,1),obj.data.modelBezCP);
    dimModelBezCP = size(obj.data.modelBezCP{1},2);
    modelBezCP = cell2mat(obj.data.modelBezCP);
    modelBezCP = modelBezCP+repmat(center,size(modelBezCP,1),1);
    obj.data.modelBezCP = mat2cell(modelBezCP,nModelBezCP,dimModelBezCP);
end

if ~isempty(obj.data.simModelBezCP)
    nSimModelBezCP = cellfun(@(a) size(a,1),obj.data.simModelBezCP);
    dimSimModelBezCP = size(obj.data.simModelBezCP{1},2);
    simModelBezCP = cell2mat(obj.data.simModelBezCP);
    simModelBezCP = simModelBezCP+repmat(center,size(simModelBezCP,1),1);
    obj.data.simModelBezCP = mat2cell(simModelBezCP,nSimModelBezCP,dimSimModelBezCP);
end

disp('Process: Data translated back to initial position!');

end