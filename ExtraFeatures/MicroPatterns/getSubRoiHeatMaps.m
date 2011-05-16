function [ forFigure ] = getSubRoiHeatMaps( subDir, statDir, numWindows, windowSize, YPattern, getConfValues,segType,colNum,condIdx,refIdx,refCond,expCond)
%

% subDir = subRois directory where the subRoiMasks can be found: just has
% to be one example directory
% 
% statDir = statistics directory where one can find the stats_cells and
% where we will output the colorMaps
%
% numWindows = the number of DEFINED size windows from the cell periphery
% (windows of size "windowSize"), note the TOTAL number of windows in each 
% region will be numWindows + 1; as the final window will always be greater than 
% the max windownumber
% 
% windowSize = size of the window in microns away from the cell edge.  
%
% YPattern = default H-pattern, to turn on Y-pattern set this variable = 1
% 
% getConfValues = set to 1 if you want confidence levels designated on the
%                heat maps (default  P-Value < 0.05 = ** 
%                                    P-Value < 0.005 = ***
%                                    P-Value < 0.0005 = ****)
% segType = 'growth', 'fgap','bgap'
% colNum = the column number of the parameter you would like to plot (See Respective Stats
% Cells)
% condIdx = the row of the stats cell to which you would like to compare to
% reference
% refIdx = row of the stats cell to be used as reference value
% refCond = Name of reference Condition for Plot Title (input as String)
% expCond = Name of Experimental Condition for Plot Title (input as String)
%% Check Parameters
if nargin<1 || isempty(subDir)
    subDir=uigetdir(pwd,'select the directory where the subRoi Masks can be found');
end

if nargin<1 || isempty(statDir)
    statDir = uigetdir(pwd,'select the directory for output');
end
%% Set up

% Set up parameters base on whether or not examining a y-pattern vs an
% H-pattern

%If y-pattern only 2 types of regions- adhesion and non-adhesion
% if H-pattern also have adhesion corners
if YPattern == 1
    numRegions = 2;
else 
    numRegions = 3;
end 


% Set cell with regionTypes 
    if numRegions == 3 
        regionTypes{1,1} = 'NonAdhesion';
        regionTypes{2,1} = 'Adhesion';
        regionTypes{3,1} = 'AdhesionCorner';
    else 
    regionTypes{1,1} = 'Adhesion';
    regionTypes{2,1} = 'NonAdhesion';
    end % end if numRegions 
    
%Load List of subRoi Masks
%roiSet = load([statDir filesep 'roiSet.mat']);


% Set number of subRois (Different if Y-Pat or H-Pat)
if YPattern == 1 
    numSubRois = (numWindows+1)*6;
else 
    numSubRois = (numWindows+1)*8;
end % end if Y-Pattern 

% does the user want to look at a growth, fgap, or bgap parameter
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
%% Body

% Collect Sub-Rois Masks
for i = 1:numSubRois
    roiSet(:,:,i)  = imread([subDir filesep 'sub_' num2str(i) filesep 'roiMask.tif']);
end  

% Get Example Cell Outline
roiYX = load([subDir filesep 'roiYX.mat']);

%Collect Statistic of interest for each region and calculate difference 
% between conditions

% Loop through each region type (adhesion, nonadhesion,etc)
for iRegion = 1:numRegions
    
    % loop through each set-size window for each region type
    for iWindow = 1:numWindows

        statsCell= load([statDir filesep 'SubRegions' filesep char(regionTypes{iRegion,1}) filesep ...
           num2str(iWindow*windowSize) 'uM' filesep statsCell2Load]);
        
       diffValue  = statsCell.(fieldName){condIdx,colNum} - statsCell.(fieldName){refIdx,colNum};
       
       % Calculate the number of previously extracted subRegions
       if iRegion == 1
           j = iWindow ; 
           
       else 
            if YPattern  == 1 
                j = iWindow + 3^(iRegion-1)*(numWindows+1);  % number of previously extracted sub_regions 
            else 
                j = iWindow + 2^(iRegion-1)*(numWindows+1);
            end % end if YPattern
           
       end % end if Region
       
       % formulate the new mask with the difference value set in the
       % respective sub-region
        if iRegion == 1 && iWindow == 1
           forFigure = diffValue*roiSet(:,:,j);
        else 
            forFigure  = forFigure + diffValue*roiSet(:,:,j);
        end % End if Region
        
        % Set all equivalent regions (b/c we pooled the data)
       forFigure = forFigure + diffValue*roiSet(:,:,j+numWindows+1); % only one other to set if hpattern ad or nonAD
      
             % more regions to set if Corner Adhesion
             if iRegion == 3 
                    forFigure = forFigure + diffValue*roiSet(:,:,j+2*(numWindows+1));
                    forFigure = forFigure + diffValue*roiSet(:,:,j+3*(numWindows+1));
             else
            end % end if iRegion == 3 
      
            % more regions to set if Y-pattern
            if YPattern == 1 % 
                    forFigure = forFigure + diffValue*roiSet(:,:,j+2*(numWindows+1));
                else 
            end % end if Y-Pattern

     % If user wants information regarding confidence in the statistic 
     % obtain p-value information from the discrimination matrix  
     if getConfValues == 1
         
                discrimMat = load([statDir filesep 'SubRegions' filesep char(regionTypes{iRegion,1}) filesep ...
                     num2str(iWindow*windowSize) 'uM' filesep 'discrimMat.mat']);
                discrimMat = discrimMat.discrimMat;
                
                if colNum == 4 
                    discrimMat2Combine{iRegion,iWindow} = discrimMat.gs_cell;
                elseif colNum == 6
                    discrimMat2Combine{iRegion,iWindow} = discrimMat.gl_cell;
                else 
                    discrimMat2Combine{iRegion,iWindow} = discrimMat.gd_cell;
                   
                end % end if colNum
                
                
                
            else 
      end % end getConfValues  
      
      
      
    end % end iWindows     
    
    % get the stats for the final inner regions (rest of inner portion of
    % the cell after all incremental size windows.
    
    statsCellGS = load([statDir filesep 'SubRegions' filesep char(regionTypes{iRegion,1}) filesep ...
        'Greater' num2str(numWindows*windowSize) 'uM' filesep statsCell2Load]);
   
    diffValue = statsCellGS.statsCellGS{condIdx,colNum} - statsCellGS.statsCellGS{refIdx,colNum};
   
    
    if iRegion == 1 
        j = numWindows+1;
        
        else 
            if YPattern  == 1 
                j = (numWindows+1) + 3^(iRegion-1)*(numWindows+1);  % number of previously extracted sub_regions 
            else 
                j = (numWindows + 1) + 2^(iRegion-1)*(numWindows+1);
            end % end if YPattern
          
            
                
            
            
     end % end if Region
   
      
            % if you want to get the confidence values % NOTE right now
            % only specific for growth parameters
            if getConfValues == 1
                discrimMat = load([statDir filesep 'SubRegions' filesep char(regionTypes{iRegion,1}) filesep ...
                    'Greater' num2str(numWindows*windowSize) 'uM' filesep 'discrimMat.mat']);
                discrimMat = discrimMat.discrimMat;
                
                titles{1,1} = regionTypes{iRegion,1}; 
                
                
                if colNum == 4 
                    discrimMat2Combine{iRegion,numWindows+1} = discrimMat.gs_cell;
                    
                elseif colNum == 6
                    discrimMat2Combine{iRegion,numWindows+1} = discrimMat.gl_cell;
                else 
                    discrimMat2Combine{iRegion,numWindows+1} = discrimMat.gd_cell;
                   
                end % end if colNum
                
                
                
            else 
            end % end getConfValues   
        
        
    
    
    forFigure = forFigure + diffValue*roiSet(:,:,j);
    forFigure = forFigure + diffValue*roiSet(:,:,j+numWindows+1);
    
    if iRegion == 3 
        
    forFigure = forFigure + diffValue*roiSet(:,:,j+2*(numWindows+1));
    forFigure = forFigure + diffValue*roiSet(:,:,j+3*(numWindows+1));
    
        else 
    end % end if iRegion == 3  
    
    if YPattern == 1 
        forFigure = forFigure + diffValue*roiSet(:,:,j+2*(numWindows+1));
    else 
    end 
    
end  

% set up matrices for finalDiscrimMat
if getConfValues == 1;
    titles = cell(1,4);
    spacer = cell(1,4);

for iRegion = 1:length(regionTypes)
for iWindow = 1:numWindows+1
    
    titles{1,1} = regionTypes{iRegion,1};
    
    if iRegion ==1 && iWindow ==1 
        discrimMatFinal = [titles; discrimMat2Combine{iRegion,iWindow};spacer];
    elseif iWindow == 1
        discrimMatFinal = [discrimMatFinal; titles; discrimMat2Combine{iRegion,iWindow};spacer];
    else 
        discrimMatFinal = [discrimMatFinal; discrimMat2Combine{iRegion,iWindow};spacer];
    end % if iRegion ==1 
end % for iWindow 
end % for iRegion

save([statDir filesep 'discrimMatFinal.mat'],'discrimMatFinal');

[numGroups dum] = size(discrimMat.gs);
for iRegion = 1:numRegions
    for iWindow = 1: numWindows + 1
        j = (iRegion-1)*(numWindows+1)*(numGroups+2) + (iRegion-1);
        
        %get appropriate pValue
        pValue = discrimMatFinal{3+(iWindow-1)*(numGroups+2)+j,condIdx};
        
        if pValue < 0.0005  
            subRoiConf{iRegion,iWindow} = '****';
        elseif pValue > 0.0005 && pValue < 0.005
            subRoiConf{iRegion,iWindow} = '***';
        elseif  (pValue > 0.005 && pValue < 0.05)
            subRoiConf{iRegion,iWindow} = '**';
        else 
            subRoiConf{iRegion,iWindow} = [];
        end 
    end 
end 



   
end % if get ConfValues  




figure1 = figure;
imagesc(forFigure);
hold all 
plot(roiYX.roiYX(:,2),roiYX.roiYX(:,1),'w');


% Create title
if colNum == 4 ;
    param = 'Difference in Mean Growth Speed (um)';
elseif colNum == 6; 
    param = 'Difference in Mean Growth LifeTime (sec)';
elseif colNum == 8
    param = 'Difference in Mean Growth Displacement (um)';
else
end



forTitle = [expCond ' vs ' refCond ':' param];
title({forTitle},...
    'FontWeight','bold',...
    'FontSize',16);


if getConfValues == 1 && YPattern == 1 
% Create textbox note only good for 3um windows right now
% need to find better way to get these coordinates


% AD1 
annotation(figure1,'textbox',...
    [0.229090575226929 0.667649239817323 0.0548788904219262 0.0505195727720901],...
    'String',{char(subRoiConf{1,1})},...
    'FontWeight','bold',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none', ...
    'Color', [ 1 1 1]);

% AD2 
% Create textbox
annotation(figure1,'textbox',...
    [0.269029463885735 0.629988075301153 0.0580282116646985 0.0406116840689256],...
    'String',{char(subRoiConf{1,2})},...
    'FontWeight','bold',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'Color', [1 1 1]);

% Create textbox NA1
annotation(figure1,'textbox',...
    [0.467432010214982 0.674358689384468 0.0493030060131902 0.0368968269834719],...
    'String',{char(subRoiConf{2,1})},...
    'FontWeight','bold',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none', ...
    'Color',[1 1 1] );

% Create textbox
annotation(figure1,'textbox',...
    [0.468176251149455 0.6573542357497 0.0442748091603053 0.0181818181818183],...
    'String',{char(subRoiConf{2,2})},...
    'FontWeight','bold',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'Color',[1 1 1]);


% Create textbox NA3
annotation(figure1,'textbox',...
    [0.467888788117716 0.59215895356089 0.0573570166727962 0.0328358208955196],...
    'String',{char(subRoiConf{2,3})},...
    'FontWeight','bold',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'Color',[1 1 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.355422733932345 0.576958842089361 0.0593184027132929 0.0279416772333834],...
    'String',{char(subRoiConf{1,3})},...
    'FontWeight','bold',...
    'FontSize',30,...
    'FitBoxToText','off',...
    'EdgeColor','none',...
    'Color',[1 1 1]);

else 
end 

end

