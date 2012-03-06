function [ forFigure ] = micropatternGetRegionalColorMaps(micropatternOutput,poolSubTracks,fieldnames )
%CREATE HeatMaps for Micropattern Data H pattern data using the output from
% micropatternOutput
%
% micropatternOutput: output of micropatternSubRoiPrepareData and appended
% by micropatternBtwCompare to including the stat directory name
% the output will automatically be saved in the same file
%
% fieldnames: the name of the fields you would like to generate heatmaps
% for: currently if [] will default to
%
% confValues: will grab the conf values (currently use
% eventually have option to pool subtracks and use permTest
%% some extra params
colorText = 'k' ; % black = k white = w

fontSizeTextBox = 12;
fontSizeTitle = 12;

pValueCutOff1 = 0.05;
pValueCutOff2 = 0.01;
pValueCutOff3 = 0.001;


%% Check AND Set-up Parameters
if nargin<1 || isempty(micropatternOutput)
    [filename,pathname]= uigetfile(pwd,'Please Select A micropatternOutput.mat File');
    
    temp = load([pathname filesep filename]);
    
    micropatternOutput = temp.micropatternOutput;
    
end

if isempty(fieldnames) % set up default fieldnames but eventually choose from a user list
    
    
    fieldnames{1} = 'growth_lifetime_mean_INSIDE_REGION';  % keep all local
    fieldnames{2} = 'growth_speed_mean_INSIDE_REGION'; %
    
    
    fieldnames{3} = 'fgap_lifetime_mean';
    fieldnames{4} = 'fgap_speed_mean';
    
    fieldnames{5} = 'bgap_speed_mean';
    fieldnames{6} = 'bgap_lifetime_mean';
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
            
            for iParam = 1:length(fieldnames)
                paramName = fieldnames{iParam};
                
                
                names{1,1} = paramName;
                
                pValueOutput = cell(3,numWindows+1); % reset
                
                
                
                diffValueAll = cell(3,numWindows+1);
                
                %Collect Statistic of interest for each region and calculate difference
                % between conditions
                
                % Loop through each region type (adhesion, nonadhesion,etc)
                for iRegion = 1:3
                    
                    % loop through each set-size window for each region type
                    for iWindow = 1:numWindows+1
                        
                        
                        groupData = micropatternOutput{iMask}.groupData;
                        
                        groupData = groupData{iExtract}{iRegion}{iWindow};
                        
                        if poolSubTracks ~=1
                            
                            distControl = arrayfun(@(x) groupData.stats{1}{x}.(paramName),1:numel(groupData.stats{1}),'uniformOutput',0);
                            meanControlPop = nanmean(cell2mat(distControl));
                            
                            distExp = arrayfun(@(x) groupData.stats{2}{x}.(paramName), 1:numel(groupData.stats{2}),'uniformOutput',0);
                            
                            meanExpPop = nanmean(cell2mat(distExp));
                            
                            diffValue = (meanExpPop-meanControlPop)/meanControlPop*100;
                            if isnan(diffValue)
                                diffValue = 0; 
                            end 
                            diffValueAll{iRegion,iWindow} = diffValue;
                            
                            
                            loadDir = statDirInd{iExtract}{iRegion}{iWindow};
                            s = load([loadDir filesep 'perCell' filesep 'discrimMat_PerCell.mat']);
                            discrimMat = s.pValues.(paramName);
                            
                            pValueOutput{iRegion,iWindow} = cell2mat(discrimMat(3,2)); % for now just use the t-test as I know there is
                            % bug in the values if there is an NaN.
                            
                            
                            
                            
                        else % pooled
                            % FIX ME!!!
                        end
                        
                        
                        if iParam ==1
                            NControl{iRegion,iWindow} = numel(groupData.stats{1});
                            NExp{iRegion,iWindow} = numel(groupData.stats{1});
                            
                            
                            for iType = 1:3
                                
                                % always get the total number of subtracks collected and the
                                % average per region
                                
                                numSubtrackControlTot{iRegion,iWindow,iType} = groupData.pooledStats{1}.(char(subtrackType{iType}));
                                numSubtrackExpTot{iRegion,iWindow,iType} = groupData.pooledStats{2}.(char(subtrackType{iType}));
                                
                                avgNumPerCellControl{iRegion,iWindow,iType} = arrayfun(@(x) groupData.stats{1}{x}.(char(subtrackType{iType})),1:NControl{iRegion,iWindow},'uniformOutput',0);
                                avgNumPerCellExp{iRegion,iWindow,iType} = arrayfun(@(x) groupData.stats{2}{x}.(char(subtrackType{iType})),1:NExp{iRegion,iWindow},'uniformOutput',0);
                                
                                
                            end
                            
                            
                            
                            
                            
                            
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
                        
                        
                        
                        
                        % make map
                    end
                end
                
                %% Output
                
                %%  Make Figure
                figure1 = figure;
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
                
                saveDir = [statDir filesep 'ColorMaps' filesep paramName] ;
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
                
                refCond = groupData.names{1};
                expCond = groupData.names{2};
                
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
                save([saveDir filesep 'colorMapOutput'],'output');
                
                close(figure1)
                
                
                
                
                
            end
        end
        % create the N-output : I also need to think about putting into text
        % likely a more clever way to do this
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
        
        save([statDir filesep 'Colormaps' filesep 'outputN'],'outputN');
        
        
    end
end
end