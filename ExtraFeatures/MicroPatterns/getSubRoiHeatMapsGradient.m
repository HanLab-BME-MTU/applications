function [ forFigure] = getSubRoiHeatMapsGradient( subDir, statDir, discrimDir, HPattern, expCond,  expIdx, getConfValues)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Check AND Set-up Parameters
if nargin<1 || isempty(subDir)
    subDir=uigetdir(pwd,'Select the directory where the subRoi masks are located');
end

if nargin<1 || isempty(statDir)
    statDir = uigetdir(pwd,'Select the directory where the statistics are located; This will also be the directory to store output ');
end

% default is to use column 2 (usually control) from stat cell
%if isempty(refIdx)
%   refIdx = 2;
%end
if nargin<1 || isempty(discrimDir)
    discrimDir = uigetdir(pwd,'Select the directory where the discrimination matrix comparing regions is located  ');
end

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

%If y-pattern only 2 types of regions
% if H-pattern also have adhesion corners
if HPattern == 1
    numRegions = 3;
    statsCell = cell(3,1);
    
else
    numRegions = 2;
    statsCell = cell(2,1);
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
if HPattern == 1
    numSubRois = 8;
    nonAd = 1;
    
else
    numSubRois = 6;
    nonAd = 2;
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
        
        
        for iRegion = 1:numRegions
            
            
            statsCell{iRegion} = load([statDir filesep 'SubRegions' filesep char(regionTypes{iRegion,1}) filesep ...
                char(statsCell2Load(iSeg))]);
            
        end
        
        % calculate the percent difference between regions (relative to the
        % non-adhesion region
        for iRegion = 1:numRegions
            
            diffValue  = (statsCell{iRegion}.(char(fieldName(iSeg))){expIdx,colNum(iColNum)} - ...
                statsCell{nonAd}.(char(fieldName(iSeg))){expIdx,colNum(iColNum)})./(statsCell{nonAd}.(char(fieldName(iSeg))){expIdx,colNum(iColNum)})*100;
            
            if iRegion == 1
                idx = 1;
                forFigure = diffValue*roiSet(:,:,idx);
            else
                forFigure = forFigure + diffValue*roiSet(:,:,idx);
            end
            
            % set all equivalent regions
            forFigure = forFigure + diffValue*roiSet(:,:,idx+1);
            
            if iRegion == 3
                forFigure = forFigure + diffValue* roiSet(:,:,idx+2);
                forFigure = forFigure + diffValue*roiSet(:,:,idx+3);
            end
            
            
            
            
            if HPattern == 0
                forFigure = forFigure + diffValue*roiSet(:,:,idx+2);  % 3-same subRois for y-pattern
                idx = idx + 3;
            else
                idx = idx + 2;
            end
            
        end % end numRegions
        
        
        
        if getConfValues == 1
            discrimMat = load([discrimDir filesep 'discrimMat.mat']);
            discrimMat = discrimMat.discrimMat;
            switch char(segType(iSeg))
                case 'Growth'
                    if colNum(iColNum) == 4
                        discrimMat2Use = discrimMat.gs_cell;
                        
                    elseif colNum(iColNum) == 7
                        discrimMat2Use = discrimMat.gl_cell;
                    elseif colNum(iColNum) == 10
                        discrimMat2Use = discrimMat.gd_cell;
                        
                    end % end if colNum
                    
                case 'FGap'
                    if colNum(iColNum) == 4
                        discrimMat2Use = discrimMat.fs_cell;
                    elseif colNum(iColNum) == 7
                        discrimMat2Use = discrimMat.fl_cell;
                    elseif colNum(iColNum) == 10
                        discrimMat2Use = discrimMat.fd_cell;
                    end
                    
                case 'BGap'
                    if colNum(iColNum) == 4
                        discrimMat2Use = discrimMat.bs_cell;
                    elseif colNum(iColNum) == 7
                        discrimMat2Use= discrimMat.bl_cell;
                    elseif colNum(iColNum) == 10
                        discrimMat2Use = discrimMat.bd_cell;
                    end
                    
            end
            
            
            for iRegion = 1:numRegions
                pValue = discrimMat2Use{2,iRegion+1};
                
                if pValue < 0.0005
                    subRoiConf{iRegion} = '***';
                elseif pValue > 0.0005 && pValue < 0.005
                    subRoiConf{iRegion} = '**';
                elseif  (pValue > 0.005 && pValue < 0.05)
                    subRoiConf{iRegion} = '*';
                else
                    subRoiConf{iRegion} = 'NotSig';
                end
            end
        end
        %%  Make Figure
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
            % AD1
            for iRegion = 1:length(regionTypes)
                
                
                if  strcmp(char(subRoiConf{iRegion}), 'NotSig') ~= 1;
                    if HPattern == 1 
                        name = char(regionTypes{iRegion});
                        sig = char(subRoiConf{iRegion}); 
                    else 
                        name = 'Adhesion'; 
                        sig = char(subRoiConf{2}); 
                    end 
                    annotation(figure1,'textbox',...
                        [0.229090575226929 0.667649239817323 0.0548788904219262 0.0505195727720901],...
                        'String',[name, ' ' sig ],...
                        'FontWeight','bold',...
                        'FontSize',16,...
                        'FitBoxToText','off',...
                        'EdgeColor','none', ...
                        'Color', [ 1 1 1],...
                        'HorizontalAlignment', 'Center');
                end % if strcmp
                
            end
        end
        
        if isdir([statDir filesep 'ColorMaps' filesep 'Gradient']) ~= 1
            mkdir([statDir filesep 'ColorMaps' filesep 'Gradient']);
        end
        
        forTitle = strrep(forTitle,' ','_');
        forTitle = strrep(forTitle,':','_');
        
        saveas(figure1,[statDir filesep 'ColorMaps' filesep 'Gradient' filesep forTitle '.fig']);
        saveas(figure1,[statDir filesep 'ColorMaps' filesep 'Gradient' filesep forTitle '.eps'],'psc2');
    end
end
end

