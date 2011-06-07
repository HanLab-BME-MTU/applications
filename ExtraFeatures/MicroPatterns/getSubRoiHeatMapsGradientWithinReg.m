function [ forFigure ] = getSubRoiHeatMapsGradientWithinReg( subDir, statDir,discrimDir, HPattern, numWindows, windowSize,expCond,expIdx,getConfValues)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%CREATE Gradient HeatMaps for Micropattern Data (H and Y patterns) 

% subDir = subRois directory where the subRoiMasks can be found: just has
% to be one example directory: if leave empty the function will ask you to navigate to 
% to this folder
% 
% statDir = statistics directory where one can find the stats_cells and
% where we will output the colorMaps: if leave empty the function will ask
% you to navigate to this folder
%
% numWindows = the number of DEFINED size windows from the cell periphery
% (windows of size "windowSize"), note the TOTAL number of windows in each 
% region will be numWindows + 1; as the final window will always be greater than 
% the max windownumber
% 
% windowSize = size of the window in microns away from the cell edge.  
%
% HPattern = Value 1 = H-pattern, Y-Pattern =  0
%
% expCond = Name of Experimental Condition for Plot Title (input as String)
%
%
% expIdx = the row of the stats cell to which you would like to compare to
% reference: if leave empty this will default to row 3 (the first
%  experimental condition of the groupList- after control)
%
%
%
%
% getConfValues = set to 1 if you want confidence levels designated on the
%                heat maps (default  P-Value < 0.05 = * 
%                                    P-Value < 0.005 = **
%                                    P-Value < 0.0005 = ***)
%% Check AND Set-up Parameters
if nargin<1 || isempty(subDir)
    subDir=uigetdir(pwd,'Select the directory where the subRoi masks are located');
end

if nargin<1 || isempty(statDir)
    statDir = uigetdir(pwd,'Select the directory where the statistics are located; This will also be the directory to store output ');
end

if nargin<1 || isempty(discrimDir)
    discrimDir = uigetdir(pwd,'Select the directory where the discrimination matrix comparing regions is located  ');
end


% default is to use column 2 (usually control) from stat cell
%if isempty(refIdx)
%   refIdx = 2;
%end

% default is to use column 3 from stat cell
if isempty(expIdx)
    expIdx = 3;
end

doAllSeg = questdlg('Would You Like To Analyze All Segment Types', 'SegmentType', 'Yes','No','Cancel','Yes');
switch doAllSeg
    case 'Yes'
        statsCell2Load = cell(3,1);
        statsCell2Load{1} = 'stats_GrowthStats';
        statsCell2Load{2} = 'stats_FgapStats';
        statsCell2Load{3} = 'stats_BgapStats';
        
        fieldName = cell(3,1);
        fieldName{1} = 'statsCellGS';
        fieldName{2} = 'statsCellFG';
        fieldName{3} = 'statsCellBG';
        
        segType = cell(3,1);
        segType{1} = 'Growth';
        segType{2} = 'FGap';
        segType{3} = 'BGap';
        
    case 'No'
        segChoice = questdlg('What Segment Type Would You Like To Analyze?','SegmentType','Growth','FGap','BGap','Growth');
        
        switch segChoice
            case 'Growth'
                statsCell2Load{1} = 'stats_GrowthStats';
                fieldName{1} = 'statsCellGS';
                segType{1} = 'Growth';
            case 'FGap'
                statsCell2Load{1} = 'stats_FgapStats';
                fieldName{1} = 'statsCellFG';
                segType{1} = 'FGap';
            case 'Bgap'
                statsCell2Load{1} = 'stats_BgapStats';
                fieldName{1} = 'statsCellBG';
                segType{1} = 'Bgap';
                
        end % segChoice
end % doAll

meanOrMedChoice = questdlg('Would you like to analyze mean or median values?','Mean or Median','Mean','Median','Mean');
switch meanOrMedChoice
    case 'Mean'
        mean = 1;
    case 'Median'
        mean = 0;
end

doAllParam = questdlg('Would You Like To Analyze All Parameters?', 'Parameter Type', 'Yes', 'No', 'Cancel', 'Yes');
switch doAllParam
    case 'Yes'
        colNum = zeros(3,1);
        if mean == 1
            colNum(1) = 4;
            colNum(2) = 7;
            colNum(3) = 10;
        else
            colNum(1) = 3;
            colNum(2) = 6;
            colNum(3) = 9;
        end
    case 'No'
        paramChoice = questdlg('What Parameter(s) Would You Like to Analyze?', 'Parameter Type','Velocity','Lifetime','Displacement','Velocity');
        switch paramChoice
            case 'Velocity'
                if mean ==1
                    colNum = 4 ;
                    
                else
                    colNum = 3 ;
                end
            case 'Lifetime'
                if mean == 1
                    colNum = 7;
                else
                    colNum = 6;
                end
            case 'Displacement'
                if mean == 1
                    colNum = 10;
                else
                    colNum = 9;
                end
        end % paramChoice
        
        
end % doAllParam


% Set up parameters base on whether or not examining a y-pattern vs an
% H-pattern

%If y-pattern only 2 types of regions- adhesion and non-adhesion
% if H-pattern also have adhesion corners
if HPattern == 1
    numRegions = 3;
else
    numRegions = 2;
end


% Set cell with regionTypes
if numRegions == 3
    regionTypes{1,1} = 'NonAdhesion';
    regionTypesS{1,1} = 'NonAd';
    regionTypes{2,1} = 'Adhesion';
    regionTypesS{2,1} = 'Ad';
    regionTypes{3,1} = 'AdhesionCorner';
    regionTypesS{3,1} = 'AdCorn';
else
    regionTypes{1,1} = 'Adhesion';
    regionTypesS{1,1} = 'Ad';
    regionTypes{2,1} = 'NonAdhesion';
    regionTypesS{2,1} = 'NonAd';
end % end if numRegions



%Load List of subRoi Masks
%roiSet = load([statDir filesep 'roiSet.mat']);


% Set number of subRois (Different if Y-Pat or H-Pat)
if HPattern == 1
    numSubRois = (numWindows+1)*8;
else
    numSubRois = (numWindows+1)*6;
end % end if H-Pattern

%% Begin Body: Generate the Percent Difference Heat Map

% Collect Sub-Rois Masks
for i = 1:numSubRois
    roiSet(:,:,i)  = imread([subDir filesep 'sub_' num2str(i) filesep 'roiMask.tif']);
end

% Get Example Cell Outline
roiYX = load([subDir filesep 'roiYX.mat']);

%Collect Statistic of interest for each region and calculate difference
% between conditions
for iSeg = 1:length(segType)
    for iColNum = 1:length(colNum);
        
        diffValue = zeros(numRegions,numWindows+1);
        
        % Loop through each region type (adhesion, nonadhesion,etc)
        for iRegion = 1:numRegions
            
            % load the stats for all the cortical regions
            statsCellCortex = cell(numWindows);
            for iWindow = 1:numWindows
                
                statsCellCortex{iWindow} = load([statDir filesep 'SubRegions' filesep char(regionTypes{iRegion,1}) filesep ...
                    num2str(iWindow*windowSize) 'uM' filesep char(statsCell2Load(iSeg))]);
            end
            
            % load the stats for the inner most region to compare
            statsCellInner = load([statDir filesep 'SubRegions' filesep char(regionTypes{iRegion,1}) filesep ...
                'Greater' num2str(numWindows*windowSize) 'uM' filesep char(statsCell2Load(iSeg))]);
            
            innerValue = statsCellInner.(char(fieldName(iSeg))){expIdx,colNum(iColNum)};
            
            
            %for each condition get difference between cortex window and inner window
            %value
            
            % for the number of cortex windows
            for iWindow = 1:numWindows
                
                outerValue = statsCellCortex{iWindow}.(char(fieldName(iSeg))){expIdx,colNum(iColNum)};
                
                
                diffValue(iRegion,iWindow) = ((outerValue - innerValue)./innerValue)*100;
            end
            
            % The difference between the inner most region and itself will be zero
            diffValue(iRegion,numWindows+1) = 0;
        end
        
        %% Create HeatMap
        for iRegion = 1:numRegions
            for iWindow = 1:numWindows
                % Calculate the number of previously extracted subRegions
                if iRegion == 1
                    j = iWindow ;
                    
                else
                    if HPattern  == 0
                        j = iWindow + 3^(iRegion-1)*(numWindows+1);  % number of previously extracted sub_regions
                    else
                        j = iWindow + 2^(iRegion-1)*(numWindows+1);
                    end % end if YPattern
                    
                end % end if Region
                
                % formulate the new mask with the difference value set in the
                % respective sub-region
                if iRegion == 1 && iWindow == 1
                    forFigure = diffValue(iRegion,iWindow)*roiSet(:,:,j);
                else
                    forFigure  = forFigure + diffValue(iRegion,iWindow)*roiSet(:,:,j);
                end % End if Region
                
                % Set all equivalent regions (b/c we pooled the data)
                forFigure = forFigure + diffValue(iRegion,iWindow)*roiSet(:,:,j+numWindows+1); % only one other to set if hpattern ad or nonAD
                
                % more regions to set if Corner Adhesion
                if iRegion == 3
                    forFigure = forFigure + diffValue(iRegion,iWindow)*roiSet(:,:,j+2*(numWindows+1));
                    forFigure = forFigure + diffValue(iRegion,iWindow)*roiSet(:,:,j+3*(numWindows+1));
                else
                end % end if iRegion == 3
                
                % more regions to set if Y-pattern
                if HPattern == 0 %
                    forFigure = forFigure + diffValue(iRegion,iWindow)*roiSet(:,:,j+2*(numWindows+1));
                else
                end % end if Y-Pattern
            end
        end
        %% Extra Feature: Confidence Values
        % If user wants information regarding confidence in the statistic
        % obtain p-value information from the discrimination matrix;
        % note this the statistical test is the perm t-test for the
        % MEANS
        if getConfValues == 1
            pValue = zeros(numRegions,numWindows);
            subRoiConf = cell(numRegions,numWindows);
            for iRegion = 1:numRegions
                for iWindow = 1:numWindows
                    discrimMat = load([discrimDir filesep 'CompareWindowsWithinRegionType' filesep char(regionTypesS{iRegion,1}) filesep ...
                        'discrimMat.mat']);
                    discrimMat = discrimMat.discrimMat;
                    
                    % load the appropriate discrimMat
                    switch char(segType(iSeg))
                        case 'Growth'
                            if colNum(iColNum) == 4
                                pValue(iRegion,iWindow) = discrimMat.gs_cell{iWindow+1,numWindows+2};
                                
                                
                            elseif colNum(iColNum) == 7
                                pValue(iRegion,iWindow)= discrimMat.gl_cell{iWindow+1,numWindows+2};
                            elseif colNum(iColNum) == 10
                                pValue(iRegion,iWindow) = discrimMat.gd_cell{iWindow+1,numWindows+2};
                                
                            end
                            
                        case 'FGap'
                            if colNum(iColNum) == 4
                                pValue(iRegion,iWindow) = discrimMat.fs_cell{iWindow+1,numWindows+2};
                            elseif colNum(iColNum) == 7
                                pValue(iRegion,iWindow) = discrimMat.fl_cell{iWindow+1,numWindows+2};
                            elseif colNum(iColNum) == 10
                                pValue(iRegion,iWindow) = discrimMat.fd_cell{iWindow+1,numWindows+2};
                            end
                            
                        case 'BGap'
                            if colNum(iColNum) == 4
                                pValue(iRegion,iWindow) = discrimMat.bs_cell{iWindow+1,numWindows+2};
                            elseif colNum(iColNum) == 7
                                pValue(iRegion,iWindow) = discrimMat.bl_cell{iWindow+1,numWindows+2};
                            elseif colNum(iColNum) == 10
                                pValue(iRegion,iWindow) = discrimMat.bd_cell{iWindow+1,numWindows+2};
                            end
                            
                    end
                    
                    
                    if pValue(iRegion,iWindow) < 0.0005
                         subRoiConf{iRegion,iWindow} = '***';
                    elseif pValue(iRegion,iWindow) > 0.0005 && pValue(iRegion,iWindow) < 0.005
                        subRoiConf{iRegion,iWindow} = '**';
                    elseif  (pValue(iRegion,iWindow) > 0.005 && pValue(iRegion,iWindow) < 0.05)
                        subRoiConf{iRegion,iWindow} = '*';
                    else
                        subRoiConf{iRegion,iWindow} = 'NotSig';
                    end
                end
            end
            
        end
        
        
        
        
        
        figure1 = figure;
        imagesc(forFigure);
        hold all
        plot(roiYX.roiYX(:,2),roiYX.roiYX(:,1),'w');
        
        
        % Create title
        
        if colNum(iColNum) == 4 ;
            param = ['Percent Difference in Mean ', char(segType(iSeg)), ' Speed '];
        elseif colNum(iColNum) == 7;
            param = ['Percent Difference in Mean ', char(segType(iSeg)), ' LifeTime '];
        elseif colNum(iColNum) == 10
            param = ['Percent Difference in Mean ', char(segType(iSeg)), ' Displacement '];
        elseif colNum(iColNum) == 3;
            param = ['Percent Difference in Median ', char(segType(iSeg)), ' Speed '];
        elseif colNum(iColNum) == 6 ;
            param = ['Percent Difference in Median ', char(segType(iSeg)), ' LifeTime '];
        elseif colNum(iColNum) == 9;
            param = ['Percent Difference in Median ', char(segType(iSeg)), ' Displacement '];
            
        end
        
        forTitle = [expCond ' Gradient : ' param];
        title({forTitle},...
            'FontWeight','bold',...
            'FontSize',16);
        caxis([-50,50]); % Here You Can Change Scale of ColorBar Axis
        colorbar;
        
        
        if getConfValues == 1
            
            
            
            
            
            % Create textbox note only good for 3um windows right now
            % need to find better way to get these coordinates
            
            
            % AD1
            for iRegion = 1:length(regionTypes)
                for iWindow = 1:numWindows
                    if  strcmp(char(subRoiConf{iRegion,iWindow}), 'NotSig') ~= 1;
                        
                        annotation(figure1,'textbox',...
                            [0.229090575226929 0.667649239817323 0.0548788904219262 0.0505195727720901],...
                            'String',{[char(regionTypes{iRegion,1}), ' ' num2str(iWindow*windowSize),'uM ' char(subRoiConf{iRegion,iWindow})]},...
                            'FontWeight','bold',...
                            'FontSize',16,...
                            'FitBoxToText','off',...
                            'EdgeColor','none', ...
                            'Color', [ 1 1 1],...
                            'HorizontalAlignment', 'Center');
                    end % if strcmp
                end % iWindow
                if isdir([statDir filesep 'ColorMaps' filesep 'Gradient']) ~= 1
                    mkdir([statDir filesep 'ColorMaps' filesep 'Gradient']);
                end
            end
            
            forTitle = strrep(forTitle,' ','_');
            forTitle = strrep(forTitle,':','_');
            saveas(figure1,[statDir filesep 'ColorMaps' filesep 'Gradient' filesep  forTitle '.fig']);
            saveas(figure1,[statDir filesep 'ColorMaps' filesep 'Gradient' filesep  forTitle '.eps'],'psc2');
            
            
            
        end
        
    end
    
    
    
    
end
end 

