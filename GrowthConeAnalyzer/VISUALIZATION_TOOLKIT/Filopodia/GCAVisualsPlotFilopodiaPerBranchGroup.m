function [ output_args ] = GCAVisualsPlotFilopodiaPerBranchGroup(filoInfo,imgSize,varargin)
%GCAVisualsPlotFilopodiaPerBranchGroup: small function for plotting each
%filo branch group together in a similar color 
%% INPUT

% REQUIRED:
% filoInfo : (REQUIRED) :  Rx1 structure
%   for a given frame where r is the number of filo
%   objects detected. Note more than one filo object can be associated to the
%   same branch group. This information is stored in filoInfo.groupCount
%
% imgSize : (REQUIRED) : 1x2 vector
%   providing the image dimensions
%   imgSize(1) is the ny size of the image and imgSize(2) is the nx size of the image

% PARAMS:
% groupLabels: Rx1 vector of the group labels IDs as indicated in filoInfo.groupCount :
%   where r is the number of groups you would like to plot
%   Default: plot all groups as indicated in filoInfo
%
% plotConnectPoint: logical :
%   flag to plot the recorded connectivity point
%   (ie filoInfo.conXYCoords(1,1) filoInfo.conXYCoords(1,2))
%   on a branch group as a scatter point. Mainly for troubleshooting connectivity
%   Default: false
%% Input check
ip = inputParser;

ip.CaseSensitive = false;


ip.addRequired('filoInfo'); 
ip.addRequired('imgSize'); 

ip.addParameter('groupLabels',[]);
ip.addParameter('plotConnectPoint',false);

ip.parse(filoInfo,imgSize,varargin{:});

%% Initiate 
if isempty(ip.Results.groupLabels); 
    groupLabels = unique(vertcat(filoInfo.groupCount),'stable'); 
else 
    groupLabels = ip.Results.groupLabels; 
end 
 
n = length(groupLabels);
c = linspecer(n);
% c = colormap(jet(n));
reOrder = randperm(n);
c = c(reOrder,:);

hold on
% for each group plot the xy coords of all the filo in that group
for iGroup = 1:length(groupLabels)
    groupLabelsAll = vertcat(filoInfo.groupCount);
     clusterC =  filoInfo(groupLabelsAll == groupLabels(iGroup));
    GCAVisualsMakeOverlaysFilopodia(clusterC,imgSize,1,1,c(iGroup,:),0);
      
    if ip.Results.plotConnectPoint 
        conIdx = vertcat(filoInfo(groupLabelsAll==groupLabels(iGroup)).conXYCoords);
        if ~isempty(conIdx)
        scatter(conIdx(:,1),conIdx(:,2),50,'y','filled'); 
        arrayfun(@(x) text(conIdx(x,1),conIdx(x,2),num2str(conIdx(x,3),2)),1:size(conIdx,1)); 
        end 
    end 
      
    %arrayfun(@(x) plot(x.Ext_coordsXY(:,1),x.Ext_coordsXY(:,2),'color', c(iGroup,:),'Linewidth',2),clusterC);
    % clusterInt = clusterC;
    %clusterInt = clusterInt(arrayfun(@(x) ~isnan(x.Int_coordsXY(1)),clusterC));
    %if ~isempty(clusterInt)
    %arrayfun(@(x) plot(x.Int_coordsXY(:,1),x.Int_coordsXY(:,2),'color',c(iGroup,:),'Linewidth',2),clusterInt);
    %end
end

