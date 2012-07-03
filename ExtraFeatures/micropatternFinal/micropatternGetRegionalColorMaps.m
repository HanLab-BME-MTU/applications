function [ forFigure ] = micropatternGetRegionalColorMaps(micropatternOutput,varargin)
%CREATE HeatMaps for Micropattern Data H pattern data using the output from
% 
% REQUIRED INPUT ; 
% micropatternOutput 
% micropatternOutput: output of micropatternSubRoiPrepareData and appended
% by micropatternBtwCompare to including the stat directory name
% the output will automatically be saved in the same file
%
% PARAMETERS : 
% 'useDefaultFieldnames': 1 to use (Default 0); list of all parameter
%  values should pull down and you can choose from them. 
% 'refCondNum' : iGroup you want to use as your reference (Default 1);
% 'expCondNum' : iGroup you want to use as your experimental (Default 2);
% 'pValueCutOff1': cut off p-Value for * (Default = 0.05);
% 'pValueCutOff2' : cut off p-value for ** (Default = 0.01)
% 'pValueCutOff3': cut off p-value for ***(Default = 0.001); 
%% parse input 


%%Input check
ip = inputParser;
ip.addRequired('micropatternOutput',@(x)isempty(x) || iscell(x));

% make time requirements optional for now so don't have to change the
% input structure: 
ip.addParamValue('poolSubTracks',0,@isscalar); % didn't implement the pooled stats yet
ip.addParamValue('useDefaultFieldnames',0,@isscalar); % type 1 if you want to use the default fieldnames
ip.addParamValue('refCondNum',1,@isscalar); % you can compare different 
ip.addParamValue('expCondNum',2,@isscalar); 
ip.addParamValue('fontSizeTextBox',12,@isscalar); 
ip.addParamValue('fontSizeTitle',12,@isscalar); 
ip.addParamValue('pValueCutOff1',0.05,@isscalar); 
ip.addParamValue('pValueCutOff2',0.01,@isscalar);
ip.addParamValue('pValueCutOff3',0.001,@isscalar); 
ip.addParamValue('colorText','k',@ischar); 
ip.parse(micropatternOutput,varargin{:});

poolSubTracks = ip.Results.poolSubTracks; 
useDefaultFieldnames = ip.Results.useDefaultFieldnames; 
refCondNum = ip.Results.refCondNum; 
expCondNum = ip.Results.expCondNum; 
fontSizeTextBox = ip.Results.fontSizeTextBox; 
fontSizeTitle = ip.Results.fontSizeTitle; 
pValueCutOff1 = ip.Results.pValueCutOff1;
pValueCutOff2= ip.Results.pValueCutOff2; 
pValueCutOff3 = ip.Results.pValueCutOff3; 
colorText = ip.Results.colorText; 



%% some extra params
% colorText = 'k' ; % black = k white = w
% 
% fontSizeTextBox = 12;
% fontSizeTitle = 12;
% 
% pValueCutOff1 = 0.05;
% pValueCutOff2 = 0.01;
% pValueCutOff3 = 0.001;

%% Check AND Set-up Parameters
if nargin<1 || isempty(micropatternOutput)
    [filename,pathname]= uigetfile(pwd,'Please Select A micropatternOutput.mat File');
    
    temp = load([pathname filesep filename]);
    
    micropatternOutput = temp.micropatternOutput;
    
end



if useDefaultFieldnames ==1 
    
 
    
    
    params{1} = 'growth_lifetime_mean_INSIDE_REGION';  % keep all local
    params{2} = 'growth_speed_mean_INSIDE_REGION'; %
    
    
    params{3} = 'fgap_lifetime_mean';
    params{4} = 'fgap_speed_mean';
    
    params{5} = 'bgap_speed_mean';
    params{6} = 'bgap_lifetime_mean';
    params{7} = 'nucleationDensity'; 
    params{8} = 'nGrowthTermEvents'; 
    params{9} = 'nGrowths'; 
    params{10} = 'CometDensityPerFrame_mean'; 
else 
  x =   fieldnames(micropatternOutput{1}.groupData{1}{1}{1}.stats{1}{1});
  y = listSelectGUI(x); 
  for i = 1: length(y)
      params{i} = x(y(i)); 
  end 
end


%%

% set up region types
regionTypes{1,1} = 'NonAdhesion';
regionTypes{2,1} = 'Adhesion';
regionTypes{3,1} = 'AdhesionCorner';

%Load List of subRoi Masks

% right now only make for subregions however should make generic because it
% would be helpful
numMasksToAnalyze = numel(micropatternOutput);

subtrackType{1} = {'nGrowths'};
subtrackType{2} = {'nFgaps'};
subtrackType{3} = {'nBgaps'};

%% START MASK LOOP%%


for iMask = 1:numMasksToAnalyze
    
    % first check to see if there are subregions
    
    maskCurrent = micropatternOutput{iMask}.maskParams;
    if maskCurrent.subRegions == 1  % currently only make masks for the the subregion parameters should make more generic
        statDirInd = micropatternOutput{iMask}.statDirInd;
        
        roiSet = micropatternOutput{iMask}.maskParams.exampleSubRoiSet; % load set of subRegions to use as outline %%% FIX
        roiYX = micropatternOutput{iMask}.maskParams.exampleRoiYXWholeCell; % load example whole cell mask
        % colllect the relevent information
        numWindows = maskCurrent.numWindows;
        windowSize = maskCurrent.windowSize;
        names = cell(1,numWindows+2  );
        for iExtract = 1:numel(maskCurrent.extractType)
            
            for iParam = 1:length(params)
                paramName = char(params{iParam});
                
                
                names{1,1} = paramName;
                
                pValueOutput = cell(3,numWindows+1); % reset
                
                
                
                diffValueAll = cell(3,numWindows+1);
                
                %Collect Statistic of interest for each region and calculate difference
                % between conditions
                zoneCount =1 ; % reset zone count
                % Loop through each region type (adhesion, nonadhesion,etc)
                for iRegion = 1:3
                    
                    % loop through each set-size window for each region type
                    for iWindow = 1:numWindows+1
                        
                        
                        groupData = micropatternOutput{iMask}.groupData;
                        
                        groupData = groupData{iExtract}{iRegion}{iWindow};
                        
                       
                        if poolSubTracks ~=1
                            
                            distControl = arrayfun(@(x) groupData.stats{refCondNum}{x}.(paramName),1:numel(groupData.stats{refCondNum}),'uniformOutput',0);
                            meanControlPop = nanmean(cell2mat(distControl));
                            values(1,zoneCount) = nanmean(cell2mat(distControl)); % first value of dataset array is meanControlPop; 
                            
                            
                            
                            distExp = arrayfun(@(x) groupData.stats{expCondNum}{x}.(paramName), 1:numel(groupData.stats{expCondNum}),'uniformOutput',0);
                            
                            meanExpPop = nanmean(cell2mat(distExp));
                            values(2,zoneCount) = nanmean(cell2mat(distExp)); % 2 value is mean Experimental population
                            
                            diffValue = (meanExpPop-meanControlPop)/meanControlPop*100;
                            if isnan(diffValue)
                                diffValue = 0; 
                            end 
                            diffValueAll{iRegion,iWindow} = diffValue;
                          
                            values(3,zoneCount) = diffValue; % 3 is the difference value 
                            
                          
                            
                           
                           
                            
                            
                            loadDir = statDirInd{iExtract}{iRegion}{iWindow};
                            s = load([loadDir filesep 'perCell' filesep 'discrimMat_PerCell.mat']);
                            discrimMat = s.pValues.(paramName);
                            
                            pValueOutput{iRegion,iWindow} = cell2mat(discrimMat(expCondNum+1,refCondNum+1)); % for now just use the t-test as I know there is
                            % bug in the values if there is an NaN.
                            
                            values(4,zoneCount) = cell2mat(discrimMat(expCondNum+1,refCondNum+1)); % 4 is the p-value
                            
                            
                        else % pooled
                            % FIX ME!!!
                        end
                        
                        
                        if iParam ==1
                            NControl{iRegion,iWindow} = numel(groupData.stats{refCondNum});
                            NExp{iRegion,iWindow} = numel(groupData.stats{expCondNum});
                            
                            values(5,zoneCount) = numel(groupData.stats{refCondNum}); % 5 is the number of cells control
                            values(6,zoneCount) = numel(groupData.stats{expCondNum}); % 6 is the number of cells exp
                          
                         
                            
                            for iType = 1:3
                                
                                % always get the total number of subtracks collected and the
                                % average per region
                                
                                numSubtrackControlTot{iRegion,iWindow,iType} = groupData.pooledStats{refCondNum}.(char(subtrackType{iType}));
                                numSubtrackExpTot{iRegion,iWindow,iType} = groupData.pooledStats{expCondNum}.(char(subtrackType{iType}));
                                
                               
                                
                                avgNumPerCellControl{iRegion,iWindow,iType} = arrayfun(@(x) groupData.stats{refCondNum}{x}.(char(subtrackType{iType})),1:NControl{iRegion,iWindow},'uniformOutput',0);
                                avgNumPerCellExp{iRegion,iWindow,iType} = arrayfun(@(x) groupData.stats{expCondNum}{x}.(char(subtrackType{iType})),1:NExp{iRegion,iWindow},'uniformOutput',0);
                                
                                
                            end
                            
                             values(7,zoneCount) = groupData.pooledStats{refCondNum}.nGrowths; 
                             values(8,zoneCount) = groupData.pooledStats{expCondNum}.nGrowths; 
                             values(9,zoneCount) = groupData.pooledStats{refCondNum}.nFgaps; 
                             values(10,zoneCount) = groupData.pooledStats{expCondNum}.nFgaps; 
                             values(11,zoneCount) = groupData.pooledStats{refCondNum}.nBgaps; 
                             values(12,zoneCount) = groupData.pooledStats{expCondNum}.nBgaps; 
                             
                             
                            
                            
                            
                            
                        end
                        
                        
                        
                        % Calculate the number of previously extracted subRegions
                        if iRegion == 1
                            j = iWindow ;
                            
                        else
                            
                            j = iWindow + 2^(iRegion-1)*(numWindows+1);
                            
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
                        
                        
                        
                        zoneCount = zoneCount +1; 
                        % make map
                    end
                    
                end
                
                
                
                %%  Make Figure
                figure1 = figure('Visible', 'off');
                imagesc(forFigure);
                hold all
                plot(roiYX(:,2),roiYX(:,1),'w');
                
                
                
                windows = cell(1,numWindows+1);
                
                for iWindow = 1:numWindows
                    names{1,iWindow+1} = [num2str(iWindow*windowSize) 'uM'];
                    windows{1,iWindow} = [num2str(iWindow*windowSize) 'uM'];
                end
                
                names{1,numWindows+2} = ['Greater Than ' num2str(numWindows*windowSize) 'uM From the Cell Edge'];
                windows{1,numWindows+1} = ['Greater Than ' num2str(numWindows*windowSize) 'uM From the Cell Edge'];
                
                statDir = micropatternOutput{iMask}.primaryStatDir ;
                
                
                statDir = strrep(statDir, 'nucleation',maskCurrent.extractType{iExtract});  
                
                
                expCondDir = groupData.names{expCondNum};
                   
                toremove = ['AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
                
                expCondDir = strrep(expCondDir,toremove,'');
                
                expCondDir = strrep(expCondDir,'grp','');
                
                saveDir = [statDir filesep 'ColorMaps_' expCondDir filesep paramName ] ;
                
                if ~isdir(saveDir);
                    mkdir(saveDir)
                end
                
                
                cellSig = cell(3,1);
                hitIdx = cell(3,3);
                
                
                for iRegion = 1:3
                    % Test PValues % each region will have a different set of
                    % rois
                    forTest = cell2mat(pValueOutput(iRegion,:)); % test over all windows
                    hitIdx{iRegion,1}= find(forTest'< pValueCutOff1 & forTest' > pValueCutOff2);
                    hitIdx{iRegion,2} = find(forTest'< pValueCutOff2 & forTest' > pValueCutOff3);
                    hitIdx{iRegion,3} = find(forTest'< pValueCutOff3);
                end
                
                cellSig{1} = '*';
                cellSig{2} = '**';
                cellSig{3} = '***';
                
                
                addOn(1) = 0;
                addOn(2) = 2*(numWindows+1);
                addOn(3) = 4*(numWindows+1);
                % Add Text for P-Value Hits
                if ~isempty(hitIdx)
                    for k = 1:3 % regiontypes
                        for j = 1:3 % sig level
                            
                            hitIdxRegSig=  hitIdx{k,j}; % get window hit idx for each region/significance level
                            % loop through all hits for
                            for i = 1:length(hitIdxRegSig)
                                tempRoiMask = roiSet(:,:,hitIdxRegSig(i)+addOn(k));% load the appropriate mask
                                
                                % where text should go
                                
                                weightedRoi=bwdist(~tempRoiMask);
                                windowText = windows(hitIdxRegSig(i));
                                
                                [r,c]=find(weightedRoi==max(weightedRoi(:)));
                                text(c(1),r(1),[regionTypes{k} windowText cellSig{j}],'color','k','fontsize',9);
                            end
                        end
                    end
                end
                
                
                
                % Create title
                paramTitle = strrep(paramName,'_',' ');
                
                if poolSubTracks ~=1
                    
                end
                
                refCond = groupData.names{refCondNum};
                expCond = groupData.names{expCondNum};
                
                toremove = ['AdhesionCorn_GreaterThan_' num2str(numWindows*windowSize) 'uM'];
                
                refCond = strrep(refCond, toremove, '');
                
                expCond = strrep(expCond,toremove,'');
                
                refCond = strrep(refCond, 'grp','');
                expCond = strrep(expCond,'grp','');
                
                refCond = strrep(refCond,'_','');
                expCond = strrep(expCond,'_','');
                
                
                forTitle = {[expCond ' vs ' refCond  ' : ' ]; 'Percent Difference in Mean of Regional Distributions'; ['Param: ' paramTitle]};
                title(forTitle,...
                    'FontWeight','bold',...
                    'FontSize',fontSizeTitle);
                caxis([-50,50]); % Here You Can Change Scale of ColorBar Axis
                colorbar;
                
                
                % formulate output
                output.pValues = [regionTypes pValueOutput];
                output.pValues = [names ; output.pValues];
                
                output.percentChanges = [regionTypes diffValueAll];
                output.percentChanges = [names ; output.percentChanges];
                
                
                
                saveas(figure1,[saveDir filesep 'heatMap' paramTitle '.eps'],'psc2')
                print(figure1,'-dtiff', '-r300',[saveDir filesep 'heatMap' paramTitle '.tif']); 
                save([saveDir filesep 'colorMapOutput'],'output');
                
                close(figure1)
                
 %%          Make Dataset array for given parameter     
              
                % get region names
count = 1; 
for iRegion = 1:3
    for iWindow = 1:numWindows+1
       
zone{count,1} = [regionTypes{iRegion} names{iWindow+1}];  
count = count +1;
    end 
end

idxSort = [6,9,3,5,8,2,4,7,1]'; 
zone = zone(idxSort); 
values = values(:,idxSort);
           
%  set obs names 
varNames{1,1} = 'Mean Of Control Population Calculated From Average Of Individual Project Values'; 
varNames{2,1} = 'Mean Of Experimental Population Calculated From Average of Individual Project Values';
varNames{3,1} = 'Percent Difference in Mean of Experimental Relative to Control (Average Project Values)'; 
varNames{4,1} = 'pValue t-test of the Means (Distributions = Project Values)'; 
varNames{5,1} = 'Number of Control Projects'; 
varNames{6,1} = 'Number of Experimental Projects'; 
varNames{7,1} = 'Number of Growth Subtracks Sampled: Pooled Control Population'; 
varNames{8,1} = 'Number of Growth Subtracks Sampled: Pooled Experimental Population'; 
varNames{9,1} = 'Number of Fgap Subtracks Sampled: Pooled Control Population'; 
varNames{10,1} = 'Number of Fgap Subtracks Sampled: Pooled Experimental Population'; 
varNames{11,1} = 'Number of Bgap Subtracks Sampled: Pooled Control Population';
varNames{12,1} = 'Number of Bgap Subtracks Sampled: Pooled Experimental Population'; 


if exist([saveDir filesep 'colorMap_DataSet' paramTitle])~=0;
    delete([saveDir filesep 'colorMap_DataSet' paramTitle]);
end 
    
 colorMapDataSetForText  = dataset({values(:,:),zone{:}},'ObsNames',varNames);                
                
 export(colorMapDataSetForText,'file',[saveDir filesep 'colorMap_DataSet_' paramTitle]);  
  
%% Make bar plots 


% First truncate var names
varNames  = cellfun(@(x) strrep(x,' ',''),varNames,'uniformoutput',0);
varNames{1,1} = 'Control'; 
varNames{2,1} = 'Experimental'; 



% Plot
colorMapDataSetForBarPlots = dataset({values(:,:)',varNames{:}},'ObsNames',zone); 

forbar = ([colorMapDataSetForBarPlots.Control colorMapDataSetForBarPlots.Experimental]); 
figure1 = figure('Visible','off'); 
axes1 = axes('Parent',figure1,'XTick',1:length(colorMapDataSetForBarPlots.Control)); 
hold(axes1,'all'); 
forbar = abs(forbar); 
bar1 = bar(forbar,'Parent',axes1);
set(bar1(1),'DisplayName',refCond); 
set(bar1(2),'DisplayName',expCond); 
xlabel({'Region Number'},'FontSize',14); 
paramName = strrep(paramName,'_',' '); 
ylabel({['Mean of Cellular Population for ' paramName]},'FontSize',14); 

    


% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.818146417445483 0.822383720930233 0.175233644859813 0.0921511627906977]);


% Create textbox
annotation(figure1,'textbox',...
    [0.784475336182471 0.973949101023715 0.001 0.001],'String',{'Periphery'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');

% Create textarrow
annotation(figure1,'textarrow',[0.134384735202492 0.903426791277259],...
    [0.986596899224806 0.986046511627907],'TextEdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.188176082906746 0.868546640667566 0.0879396984924623 0.0394574599260173],...
    'String',{'Adhesion'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create line
annotation(figure1,'line',[0.194790150127585 0.297507788161994],...
    [0.86033685085883 0.86046511627907]);

% Create line
annotation(figure1,'line',[0.425845426509496 0.529798525336966],...
    [0.874300088893987 0.87433831330829]);

% Create textbox
annotation(figure1,'textbox',...
    [0.412163308755617 0.897415392997447 0.0879396984924623 0.0394574599260173],...
    'String',{'Adhesion'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.66430532725935 0.884704699911106 0.0879396984924623 0.0394574599260173],...
    'String',{'Adhesion'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create line
annotation(figure1,'line',[0.683497863147514 0.787551464487547],...
    [0.863493447652912 0.863802942104206]);

% Create textbox
annotation(figure1,'textbox',...
    [0.121628391177069 0.991742236113899 0.001 0.001],...
    'String',{'Central Region'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Arial',...
    'FitBoxToText','off');

%

print(figure1,'-dtiff', '-r300',[saveDir filesep 'values.tif']); 
 %know there is a way to save some of this info but can't quite get it
 %right so just close and start over
close(figure1); 
 
 control = colorMapDataSetForBarPlots.Control; 
 exp = colorMapDataSetForBarPlots.Experimental; 
% make second figure
forbar2 = zeros(9,2);
count = 1;
for j = 3:3:9
for i = j-2:j
forbar2(count,1) = (control(i)-control(j))/control(j); 
count = count +1; 
end
end

count = 1; 
for j = 3:3:9
    for i = j-2:j
        forbar2(count,2) = (exp(i)-exp(j))/exp(j); 
        count = count +1; 
    end 
end 

forbar2 = forbar2*100; 

figure2 = figure('Visible','off'); 
axes1 = axes('Parent',figure2,'XTick',1:length(colorMapDataSetForBarPlots.Control)); 
hold(axes1,'all'); 
 
 

 bar2 = bar(forbar2,'Parent',axes1);
 set(bar2(1),'DisplayName',refCond); 
 set(bar2(2),'DisplayName',expCond); 
    ylabel({['Mean of Cellular Population for ' paramName]},'FontSize',14); 

if (max(forbar2(:))<100 || min(forbar2(:))>-100) 
 axis([0 10 -100 100])
end 


xlabel({'Region Number'},'FontSize',14); 

ylabel({['Percent Change ' paramName 'Relative To'] ; 'Non-Adhesion Region '},'FontSize',14); 

    


% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.818146417445483 0.822383720930233 0.175233644859813 0.0921511627906977]);


% Create textbox
annotation(figure1,'textbox',...
    [0.784475336182471 0.973949101023715 0.001 0.001],'String',{'Periphery'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off');

% Create textarrow
annotation(figure1,'textarrow',[0.134384735202492 0.903426791277259],...
    [0.986596899224806 0.986046511627907],'TextEdgeColor','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.188176082906746 0.868546640667566 0.0879396984924623 0.0394574599260173],...
    'String',{'Adhesion'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create line
annotation(figure1,'line',[0.194790150127585 0.297507788161994],...
    [0.86033685085883 0.86046511627907]);

% Create line
annotation(figure1,'line',[0.425845426509496 0.529798525336966],...
    [0.874300088893987 0.87433831330829]);

% Create textbox
annotation(figure1,'textbox',...
    [0.412163308755617 0.897415392997447 0.0879396984924623 0.0394574599260173],...
    'String',{'Adhesion'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create textbox
annotation(figure1,'textbox',...
    [0.66430532725935 0.884704699911106 0.0879396984924623 0.0394574599260173],...
    'String',{'Adhesion'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');

% Create line
annotation(figure1,'line',[0.683497863147514 0.787551464487547],...
    [0.863493447652912 0.863802942104206]);

% Create textbox
annotation(figure1,'textbox',...
    [0.121628391177069 0.991742236113899 0.001 0.001],...
    'String',{'Central Region'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FontName','Arial',...
    'FitBoxToText','off');


 
                
print(figure2,'-dtiff', '-r300',[saveDir filesep 'PercentChangeAdVsNon.tif']); 
                
            close(figure2);     
                
                
                
                
            end %for iParam
       
     
        
        
        
        
        % create the N-output : I also need to think about putting into text
        % like a more clever way to do this
        % but good enough for now.
        names{1,1} = 'nGrowths';
        outputN.TotalSubTrackControlGrowth = [regionTypes numSubtrackControlTot(:,:,1)];
        outputN.TotalSubTrackControlGrowth = [names ;outputN.TotalSubTrackControlGrowth];
        
        
        outputN.TotalSubTrackExpGrowth = [regionTypes numSubtrackExpTot(:,:,1)];
        outputN.TotalSubTrackExpGrowth = [names ;outputN.TotalSubTrackExpGrowth];
        
        names{1,1} = 'nFGaps';
        outputN.TotalSubTrackControlFgap = [regionTypes numSubtrackControlTot(:,:,2)];
        outputN.TotalSubTrackControlFgap = [names ; outputN.TotalSubTrackControlFgap];
        
        outputN.TotalSubTrackExpFgap = [regionTypes numSubtrackExpTot(:,:,2)];
        outputN.TotalSubTrackExpFgap = [names; outputN.TotalSubTrackExpFgap];
        
        names{1,1} = 'nBGaps';
        outputN.TotalSubTrackControlBgap = [regionTypes numSubtrackControlTot(:,:,3)];
        outputN.TotalSubTrackControlBgap = [names;outputN.TotalSubTrackControlBgap];
        
        outputN.TotalSubTrackExpBgap = [regionTypes numSubtrackExpTot(:,:,3)];
        outputN.TotalSubTrackExpBgap = [names;outputN.TotalSubTrackExpBgap];
        
        save([statDir filesep 'Colormaps_' expCondDir filesep 'outputN'],'outputN');
       end 
         end
            end 



end