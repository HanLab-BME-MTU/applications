function [out] = groupArmaCoefficients(data,options)
%GROUPARMACOEFFICIENTS groups models according to similarity of their ARMA coefficients
%
% SYNOPSIS: [out] = groupArmaCoefficients(data,options)
%
% INPUT 
%
% OUTPUT 
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 22-May-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%============================
%% Test Input / Load Data
%============================

out = [];
def_mode = 4;

% tbd: test input, check options, set defaults

if nargin == 0 || isempty(data)
    % load data
    [fname1, pname1] = uigetfile('*.*','Load ARMA data');
    [fname2, pname2] = uigetfile('*.*','Load strain info');
    if fname1 == 0 || fname2 == 0
        disp('--groupArmaCoefficients aborted')
        return
    end
    
    % armaData is a collection of many files, while strainInfo is already
    % the correct structure that contains the optimal model for each 
    armaData = load(fullfile(pname1,fname1));
    strainInfo = load(fullfile(pname2,fname2));
    
    % strainInfo is a structure in itself
    fn = fieldnames(strainInfo);
    strainInfo = strainInfo.(fn{1});
    
    % create data
    nData = length(strainInfo);
    data(1:nData) = struct;
    
    % read the correct model
    for i=1:nData,
        data(i,1)=...
            armaData.(sprintf('fitLen%s',strainInfo(i).name))...
            (strainInfo(i).orderLen(1)+1,strainInfo(i).orderLen(2)+1);
    end
    
    % read strainInfo into data
    [data.name] = deal(strainInfo.name);
    [data.orderLen] = deal(strainInfo.orderLen);
    
end
    
if nargin < 2 || isempty(options)
    options.mode = def_mode;
end

%============================


%============================
%% Group according to WNV
%============================

[links, groupIdx, linkData] = groupData(data, 'groupArma_distance_WNV', options);

% visualize result
figure,
dendrogram(links,0,'labels',strvcat(data.name),'orientation','right');

% cut off into smaller groups

% for right now: two groups
nGroups = 3;
groupLabels = unique(groupIdx(:,end-nGroups+2));

for iGroup = 1:nGroups
collectedData(iGroup).data = ...
    data(groupIdx(:,end-nGroups+2) == groupLabels(iGroup));
end

%===========================


%===========================
%% Group according to ARMA
%===========================

options.mode = 2;

for iGroup = 1:nGroups
    [collectedData(iGroup).links,collectedData(iGroup).groupIdx,...
        collectedData(iGroup).linkData] = groupData(...
        collectedData(iGroup).data,'groupArma_distance_ARMA',options);
    
    figure,
dendrogram(collectedData(iGroup).links,0,...
    'labels',strvcat(collectedData(iGroup).data.name),...
    'orientation','right');
end