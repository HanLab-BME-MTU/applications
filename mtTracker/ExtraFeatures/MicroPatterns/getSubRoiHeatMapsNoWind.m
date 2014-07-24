function [ output_args ] = getSubRoiHeatMapsNoWind(subDir, statDir, segType,colNum,expIdx,refIdx,expCond,refCond )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Check Parameters
if nargin<1 || isempty(subDir)
    subDir=uigetdir(pwd,'select the directory where the subRoi Masks can be found');
end

if nargin<2 || isempty(statDir)
    statDir = uigetdir(pwd,'select the directory where the stats can be found');
end

%%
% Collect Sub-Rois Masks
for i = 1:8
    roiSet(:,:,i)  = imread([subDir filesep 'sub_' num2str(i) filesep 'roiMask.tif']);
end  

% Get Example Cell Outline
roiYX = load([subDir filesep 'roiYX.mat']);

[imL imW] = size(roiSet(:,:,1));
%if YPattern == 1
 %   numRegions = 2;
%else 
    numRegions = 3;
%end 


% Set cell with regionTypes 
    if numRegions == 3 
        regionTypes{1,1} = 'NonAdhesion';
        regionTypes{2,1} = 'Adhesion';
        regionTypes{3,1} = 'AdhesionCorner';
    else 
    regionTypes{1,1} = 'Adhesion';
    regionTypes{2,1} = 'NonAdhesion';
    end % end if numRegions 

 
% does the user want to look at a Growth, Fgap, or Bgap parameter
switch lower(segType)
    case {'growth'}
        statsCell2Load = 'stats_GrowthStats';
        fieldName = 'statsCellGS';
    case {'fgap'} 
        statsCell2Load = 'stats_FGapStats';
        fieldName = 'statsCellFG';
    case {'bgap'} 
        statsCell2Load = 'stats_BGapStats';
        fieldName = 'statsCellBG';
end    
    
    
    
    
% Loop through each region type (adhesion, nonadhesion,etc)
for iRegion = 1:numRegions
   statsCell= load([statDir filesep 'SubRegions' filesep char(regionTypes{iRegion,1}) filesep ...
            statsCell2Load]);
    
    diffValue = (abs(statsCell.(fieldName){expIdx,colNum}) - abs(statsCell.(fieldName){refIdx,colNum}))./(abs(statsCell.(fieldName){refIdx,colNum}))*100;
   

if iRegion == 1
idx = 1 ;     
forFigure = diffValue*roiSet(:,:,idx);

else 
    forFigure = forFigure + diffValue*roiSet(:,:,idx);
   
end    

forFigure = forFigure + diffValue*roiSet(:,:,idx+1); 

if iRegion == 3 
    forFigure = forFigure + diffValue* roiSet(:,:,idx+2);
    forFigure = forFigure + diffValue*roiSet(:,:,idx+3);
end 
idx = idx + 2; 
end 

figure1 = figure;
imagesc(forFigure);
hold all 
plot(roiYX.roiYX(:,2),roiYX.roiYX(:,1),'w');


% Create title

if colNum == 4 ;
    param = ['Percent Difference in Mean ', segType, ' Speed'];
elseif colNum == 7; 
    param = ['Percent Difference in Mean ', segType, ' LifeTime'];
elseif colNum == 10
    param = ['Percent Difference in Mean ', segType, ' Displacement'];
else
end    
 
forTitle = [expCond ' vs ' refCond ':' param];
title({forTitle},...
    'FontWeight','bold',...
    'FontSize',16);
caxis([-50,50]);
colorbar;

if isdir([statDir filesep 'ColorMaps']) ~= 1 
    mkdir([statDir filesep 'ColorMaps']); 
end 
forTitle = strrep(forTitle,' ','_');
forTitle = strrep(forTitle,':','_');
saveas(figure1,[statDir filesep 'ColorMaps' filesep forTitle '.fig']);
saveas(figure1,[statDir filesep 'ColorMaps' filesep forTitle '.eps']);
end

