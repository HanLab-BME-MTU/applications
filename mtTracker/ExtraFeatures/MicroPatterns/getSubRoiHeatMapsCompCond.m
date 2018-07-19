function [ forFigure ] = getSubRoiHeatMapsCompCond( subDir, statDir,HPattern, numWindows, windowSize,expCond,refCond,expIdx,refIdx,getConfValues)
%CREATE HeatMaps for Micropattern Data (H and Y patterns) 

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
% refCond = Name of reference Condition for Plot Title (input as String)
%
% expIdx = the row of the stats cell to which you would like to compare to
% reference: if leave empty this will default to row 3 (the first
%  experimental condition of the groupList- after control)
%
% refIdx = row of the stats cell to be used as reference value
%          if left empty default is row 2 (typcially control)
%
%
% getConfValues = set to 1 if you want confidence levels designated on the
%                heat maps (default  P-Value < 0.05 = ** 
%                                    P-Value < 0.005 = ***
%                                    P-Value < 0.0005 = ****)
%% some extra params
colorText = 'k' ; % black = k white = w 

fontSizeTextBox = 12; 
fontSizeTitle = 12; 

pValueCutOff1 = 0.05;
pValueCutOff2 = 0.005;
pValueCutOff3 = 0.0005; 


%% Check AND Set-up Parameters
if nargin<1 || isempty(subDir)
    subDir=uigetdir(pwd,'Select the directory where the subRoi masks are located');
end

if nargin<1 || isempty(statDir)
    statDir = uigetdir(pwd,'Select the directory where the statistics are located; This will also be the directory to store output ');
end

% default is to use row 2 (usually control) from stat cell
if isempty(refIdx)
    refIdx = 2;
end

% default is to use row 3 from stat cell
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
    
    % Loop through each region type (adhesion, nonadhesion,etc)
    for iRegion = 1:numRegions
        if numWindows ~= 0;
        % loop through each set-size window for each region type
        for iWindow = 1:numWindows
            
            statsCell= load([statDir filesep 'SubRegions' filesep char(regionTypes{iRegion,1}) filesep ...
                num2str(iWindow*windowSize) 'uM' filesep char(statsCell2Load(iSeg))]);
            
            % calculate the percent difference between control and exp condition
            diffValue  = (statsCell.(char(fieldName(iSeg))){expIdx,colNum(iColNum)} - ...
                statsCell.(char(fieldName(iSeg))){refIdx,colNum(iColNum)})./(statsCell.(char(fieldName(iSeg))){refIdx,colNum(iColNum)})*100;
            
            % Calculate the number of previously extracted subRegions
            if iRegion == 1
                j = iWindow ;
                
            else
                if HPattern  == 0
                    j = iWindow + 3^(iRegion-1)*(numWindows+1);  % number of previously extracted sub_regions
                else
                    j = iWindow + 2^(iRegion-1)*(numWindows+1);
                end % end if HPattern
                
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
            if HPattern == 0 %
                forFigure = forFigure + diffValue*roiSet(:,:,j+2*(numWindows+1));
            else
            end % end if Y-Pattern
            
%%          Extra Feature: Confidence Values
            % If user wants information regarding confidence in the statistic
            % obtain p-value information from the discrimination matrix;
            % note this the statistical test is the perm t-test for the
            % MEANS
            if getConfValues == 1
                
                discrimMat = load([statDir filesep 'SubRegions' filesep char(regionTypes{iRegion,1}) filesep ...
                    num2str(iWindow*windowSize) 'uM' filesep 'discrimMat.mat']);
                discrimMat = discrimMat.discrimMat;
             
            switch char(segType(iSeg))
                case 'Growth'
                    if colNum(iColNum) == 4
                        discrimMat2Combine{iRegion,iWindow} = discrimMat.gs_cell;
                        
                    elseif colNum(iColNum) == 7
                        discrimMat2Combine{iRegion,iWindow} = discrimMat.gl_cell;
                    elseif colNum(iColNum) == 10
                        discrimMat2Combine{iRegion,iWindow} = discrimMat.gd_cell;
                        
                    end 
                    
                case 'FGap'
                    if colNum(iColNum) == 4
                        discrimMat2Combine{iRegion,iWindow} = discrimMat.fs_cell;
                    elseif colNum(iColNum) == 7
                        discrimMat2Combine{iRegion,iWindow} = discrimMat.fl_cell;
                    elseif colNum(iColNum) == 10
                        discrimMat2Combine{iRegion,iWindow} = discrimMat.fd_cell;
                    end
                    
                case 'BGap'
                    if colNum(iColNum) == 4
                        discrimMat2Combine{iRegion,iWindow} = discrimMat.bs_cell;
                    elseif colNum(iColNum) == 7
                        discrimMat2Combine{iRegion,iWindow} = discrimMat.bl_cell;
                    elseif colNum(iColNum) == 10
                        discrimMat2Combine{iRegion,iWindow} = discrimMat.bd_cell;
                    end
                    
            end
                
                
            else
            end % end getConfValues
            
            
            
        end % end iWindows
     
        
            %% Code for Greater Than X-um windows 
        % Doesn't fit nicely in window loop so do here. 
        % get the stats for the final inner regions (rest of inner portion of
        % the cell after all incremental size windows.
        
        statsCell = load([statDir filesep 'SubRegions' filesep char(regionTypes{iRegion,1}) filesep ...
            'Greater' num2str(numWindows*windowSize) 'uM' filesep char(statsCell2Load(iSeg))]);
        
        % calculate the percent difference between control and exp condition
            diffValue  = (statsCell.(char(fieldName(iSeg))){expIdx,colNum(iColNum)} - ...
                statsCell.(char(fieldName(iSeg))){refIdx,colNum(iColNum)})./(statsCell.(char(fieldName(iSeg))){refIdx,colNum(iColNum)})*100;
    
       if iRegion == 1
            j = numWindows+1;
            
        else
            if HPattern == 0
                j = (numWindows+1) + 3^(iRegion-1)*(numWindows+1);  % number of previously extracted sub_regions
            else
                j = (numWindows + 1) + 2^(iRegion-1)*(numWindows+1);
            end % end if YPattern
       end     
            
        forFigure = forFigure + diffValue*roiSet(:,:,j);
        forFigure = forFigure + diffValue*roiSet(:,:,j+numWindows+1);
        
        if iRegion == 3 % means H-pattern therefore adhesion corners
            
            forFigure = forFigure + diffValue*roiSet(:,:,j+2*(numWindows+1));
            forFigure = forFigure + diffValue*roiSet(:,:,j+3*(numWindows+1));
            
        else
        end % end if iRegion == 3
        
        if HPattern == 0
            forFigure = forFigure + diffValue*roiSet(:,:,j+2*(numWindows+1)); % y-pattern
        else
        end
       
       
       
        %% if do NOT have windows 
        else 
             statsCell= load([statDir filesep 'SubRegions' filesep char(regionTypes{iRegion,1}) filesep ...
                 char(statsCell2Load(iSeg))]);
            
            % calculate the percent difference between control and exp condition
            diffValue  = (statsCell.(char(fieldName(iSeg))){expIdx,colNum(iColNum)} - ...
                statsCell.(char(fieldName(iSeg))){refIdx,colNum(iColNum)})./(statsCell.(char(fieldName(iSeg))){refIdx,colNum(iColNum)})*100;
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
            
            
            
            
            if HPattern == 0 
                forFigure = forFigure + diffValue*roiSet(:,:,idx+2);  % 3-same subRois for y-pattern 
                idx = idx + 3; 
            else 
                idx = idx + 2; 
            end 
                
            
            
        end % if numWindows 
            
        
        
        %% Again Extra Feature confidence values- GREATER THAN X uM or if no windows
        %
        if getConfValues == 1
            if numWindows ~= 0
                discrimMat = load([statDir filesep 'SubRegions' filesep char(regionTypes{iRegion,1}) filesep ...
                    'Greater' num2str(numWindows*windowSize) 'uM' filesep 'discrimMat.mat']);
            elseif numWindows == 0
                discrimMat = load([statDir filesep 'SubRegions' filesep...
                    char(regionTypes{iRegion,1}) filesep 'discrimMat.mat']);
            end
            
            
            discrimMat = discrimMat.discrimMat;
            
            titles{1,1} = regionTypes{iRegion,1};
            
            switch char(segType(iSeg))
                case 'Growth'
                    if colNum(iColNum) == 4
                        discrimMat2Combine{iRegion,numWindows+1} = discrimMat.gs_cell;
                        
                    elseif colNum(iColNum) == 7
                        discrimMat2Combine{iRegion,numWindows+1} = discrimMat.gl_cell;
                    elseif colNum(iColNum) == 10
                        discrimMat2Combine{iRegion,numWindows+1} = discrimMat.gd_cell;
                        
                    end % end if colNum
                    
                case 'FGap'
                    if colNum(iColNum) == 4
                        discrimMat2Combine{iRegion,numWindows+1} = discrimMat.fs_cell;
                    elseif colNum(iColNum) == 7
                        discrimMat2Combine{iRegion,numWindows+1} = discrimMat.fl_cell;
                    elseif colNum(iColNum) == 10
                        discrimMat2Combine{iRegion,numWindows+1} = discrimMat.fd_cell;
                    end
                    
                case 'BGap'
                    if colNum(iColNum) == 4
                        discrimMat2Combine{iRegion,numWindows+1} = discrimMat.bs_cell;
                    elseif colNum(iColNum) == 7
                        discrimMat2Combine{iRegion,numWindows+1} = discrimMat.bl_cell;
                    elseif colNum(iColNum) == 10
                        discrimMat2Combine{iRegion,numWindows+1} = discrimMat.bd_cell;
                    end
                    
            end
        end % end getConfValues
        
            
            
        
        
        
    end % for iRegion
    
    % set up matrices for finalDiscrimMat
    if getConfValues == 1;
        titles = cell(1,length(discrimMat2Combine{1,1}));
        spacer = cell(1,length(discrimMat2Combine{1,1}));
        
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
        
        if exist('discrimMatFinalAll','var') 
        
        discrimMatFinalAll = [discrimMatFinalAll; discrimMatFinal]; 
        
        else 
            discrimMatFinalAll = discrimMatFinal;
      
        end 
        
        [numGroups dum] = size(discrimMat.gs);
        for iRegion = 1:numRegions
            for iWindow = 1: numWindows + 1
                j = (iRegion-1)*(numWindows+1)*(numGroups+2) + (iRegion-1);
                
                %get appropriate pValue
                pValue = discrimMatFinal{3+(iWindow-1)*(numGroups+2)+j,expIdx};
                
                if pValue < 0.0005
                    subRoiConf{iRegion,iWindow} = '***';
                elseif pValue > pValueCutOff3 && pValue < pValueCutOff2
                    subRoiConf{iRegion,iWindow} = '**';
                elseif  (pValue > pValueCutOff2 && pValue < pValueCutOff1)
                    subRoiConf{iRegion,iWindow} = '*';
                else
                    subRoiConf{iRegion,iWindow} = 'NotSig';
                end
            end
        end
        
        
        
        
    end % if get ConfValues
    
    
    
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
    
    forTitle = [expCond ' vs ' refCond ' : ' param];
    title({forTitle},...
        'FontWeight','bold',...
        'FontSize',fontSizeTitle);
    caxis([-50,50]); % Here You Can Change Scale of ColorBar Axis
    colorbar;
    
    
   
    
    
    if getConfValues == 1 
        
       
            
            
            
        % Create textbox note only good for 3um windows right now
        % need to find better way to get these coordinates
        
        
        % AD1
        for iRegion = 1:length(regionTypes)
            if numWindows ~= 0 
            for iWindow = 1:numWindows
                if  strcmp(char(subRoiConf{iRegion,iWindow}), 'NotSig') ~= 1;
                    
                    annotation(figure1,'textbox',...
                        [0.229090575226929 0.667649239817323 0.0548788904219262 0.0505195727720901],...
                        'String',{[char(regionTypes{iRegion,1}), ' ' num2str(iWindow*windowSize),'uM ' char(subRoiConf{iRegion,iWindow})]},...
                        'FontWeight','bold',...
                        'FontSize',fontSizeTextBox,...
                        'FitBoxToText','off',...
                        'EdgeColor','none', ...
                        'Color', colorText,...
                        'HorizontalAlignment', 'Center');
                end % if strcmp
            end % iWindow
            
            % 
            if strcmp(char(subRoiConf{iRegion,numWindows+1}),'NotSig') ~= 1;
                
                annotation(figure1,'textbox',...
                    [0.229090575226929 0.667649239817323  0.0548788904219262 0.0505195727720901],...
                    'String',{[char(regionTypes{iRegion,1}), ' >' num2str(numWindows*windowSize),'uM ' char(subRoiConf{iRegion,numWindows+1})]},...
                    'FontWeight','bold',...
                    'FontSize',fontSizeTextBox,...
                    'FitBoxToText','off',...
                    'EdgeColor','none', ...
                    'Color', colorText, ...
                    'HorizontalAlignment', 'Center');
            end
            else 
                if  strcmp(char(subRoiConf{iRegion,iWindow}), 'NotSig') ~= 1;
                    
                    annotation(figure1,'textbox',...
                        [0.229090575226929 0.667649239817323 0.0548788904219262 0.0505195727720901],...
                        'String',{[char(regionTypes{iRegion,1}), ' ' char(subRoiConf{iRegion})]},...
                        'FontWeight','bold',...
                        'FontSize',fontSizeTextBox,...
                        'FitBoxToText','off',...
                        'EdgeColor','none', ...
                        'Color', colorText,...
                        'HorizontalAlignment', 'Center');
                end % if strcmp
             
            end % if numWindows 
        end % iRegion
    end % getConfValues
    
    
    
    if isdir([statDir filesep 'ColorMaps' filesep 'CompCond']) ~= 1
        mkdir([statDir filesep 'ColorMaps' filesep 'CompCond']);
    end
    
    forTitle = strrep(forTitle,' ','_');
    forTitle = strrep(forTitle,':','_');
    saveas(figure1,[statDir filesep 'ColorMaps' filesep 'CompCond' filesep  forTitle '.fig']);
    saveas(figure1,[statDir filesep 'ColorMaps' filesep 'CompCond' filesep  forTitle '.eps'],'psc2');
    
    
    
    test = 0;
    if test ==1 
        % AD2
        % Create textbox
        annotation(figure1,'textbox',...
            [0.269029463885735 0.629988075301153 0.0580282116646985 0.0406116840689256],...
            'String',{char(subRoiConf{1,2})},...
            'FontWeight','bold',...
            'FontSize',fontSizeTextBox,...
            'FitBoxToText','off',...
            'LineStyle','none',...
            'Color', colorText);
        
        % Create textbox NA1
        annotation(figure1,'textbox',...
            [0.467432010214982 0.674358689384468 0.0493030060131902 0.0368968269834719],...
            'String',{char(subRoiConf{2,1})},...
            'FontWeight','bold',...
            'FontSize',fontSizeTextBox,...
            'FitBoxToText','off',...
            'EdgeColor','none', ...
            'Color',colorText );
        
        % Create textbox
        annotation(figure1,'textbox',...
            [0.468176251149455 0.6573542357497 0.0442748091603053 0.0181818181818183],...
            'String',{char(subRoiConf{2,2})},...
            'FontWeight','bold',...
            'FontSize',fontSizeTextBox,...
            'FitBoxToText','off',...
            'EdgeColor','none',...
            'Color',colorText);
        
        
        % Create textbox NA3
        annotation(figure1,'textbox',...
            [0.467888788117716 0.59215895356089 0.0573570166727962 0.0328358208955196],...
            'String',{char(subRoiConf{2,3})},...
            'FontWeight','bold',...
            'FontSize',fontSizeTextBox,...
            'FitBoxToText','off',...
            'EdgeColor','none',...
            'Color',colorText);
        
        % Create textbox
        annotation(figure1,'textbox',...
            [0.355422733932345 0.576958842089361 0.0593184027132929 0.0279416772333834],...
            'String',{char(subRoiConf{1,3})},...
            'FontWeight','bold',...
            'FontSize',fontSizeTextBox,...
            'FitBoxToText','off',...
            'EdgeColor','none',...
            'Color',colorText);
    end 
   
   
end % for iColNum
end % enf iSeg    
  save([statDir filesep 'discrimMatFinalAll.mat'],'discrimMatFinalAll');
end

